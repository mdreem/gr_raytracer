#!/usr/bin/env bash
# Offload a CPU render to a RunPod pod, using Backblaze B2 as the durable hub.
#
# Flow (fully hands-off, driven by .env):
#   build       -> cross-build a static Linux binary in Docker (musl)
#   push-data   -> upload the binary + Gaia parquet to B2 (one-time-ish)
#   render      -> create a CPU pod that pulls inputs from B2, renders, pushes
#                  the result + a DONE marker to B2, then TERMINATES ITSELF
#   fetch       -> download the result from B2 once DONE appears
#
# The pod is ephemeral and self-terminating, so an abort on either side is safe:
# the result lands in B2 regardless, and nothing keeps billing. Re-run `fetch`
# any time.
#
# Requires: docker, rclone, curl, jq locally. Secrets come from .env (gitignored):
#   RUNPOD_API_KEY, RCLONE_CONFIG_B2_{TYPE,ACCOUNT,KEY}, B2_BUCKET
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

# --- config you may want to tweak ---------------------------------------------
# Build in a Debian (bookworm) rust container and run on the matching Debian
# base, so the glibc the binary links against is present on the pod. (musl would
# give a fully static binary, but zstd-sys, pulled in by parquet, fails to
# assemble under musl-gcc; matching the glibc version is simpler and robust.)
BINARY="target/amd64/release/gr_raytracer"         # glibc build, matches POD_IMAGE
                                                   # (own target dir so it never mixes
                                                   # with the host's arm64 artifacts)
POD_IMAGE="${POD_IMAGE:-debian:bookworm-slim}"     # same base as rust:latest (bookworm)
POD_DISK_GB="${POD_DISK_GB:-20}"
# RunPod CPU pod sizing. cpuFlavorIds picks the CPU family (see `types`); vCPU
# count is separate. cpu3c = 3rd-gen compute-optimized; 8 vCPU is a good render box.
POD_CPU_FLAVOR="${POD_CPU_FLAVOR:-cpu3c}"          # one of cpu3c/3g/3m, cpu5c/5g/5m
POD_VCPU="${POD_VCPU:-8}"                           # vCPUs allocated to the pod
POD_CLOUD="${POD_CLOUD:-SECURE}"                    # SECURE or COMMUNITY (cheaper)
RUNPOD_REST="https://rest.runpod.io/v1"
PARQUET="data/gaia_mag12.parquet"                  # matches [star_catalog].path in scenes
# ------------------------------------------------------------------------------

load_env() {
  [ -f .env ] || { echo "missing .env (see script header)"; exit 1; }
  set -a; . ./.env; set +a
  : "${RUNPOD_API_KEY:?}"; : "${B2_BUCKET:?}"
  : "${RCLONE_CONFIG_B2_TYPE:?}"; : "${RCLONE_CONFIG_B2_ACCOUNT:?}"; : "${RCLONE_CONFIG_B2_KEY:?}"
}

rp() {  # RunPod REST call: rp METHOD PATH [json-body]
  local method="$1" path="$2" body="${3:-}"
  curl -fsS -X "$method" "${RUNPOD_REST}${path}" \
    -H "Authorization: Bearer ${RUNPOD_API_KEY}" \
    -H "Content-Type: application/json" \
    ${body:+-d "$body"}
}

# True iff object <marker> exists in jobs/<job>/. Note: `rclone lsf b2:.../FILE`
# lists FILE as if it were a directory and exits 0 even when it is absent, so we
# must list the parent dir and match the exact name instead.
b2_has() { rclone lsf "b2:$B2_BUCKET/jobs/$1/" 2>/dev/null | grep -qx "$2"; }

cmd_build() {
  # --platform linux/amd64 forces an x86_64 build even on an Apple-Silicon host
  # (RunPod CPU pods are x86_64). On arm64 hosts this runs under emulation, so it
  # is slower than a native build but produces the right architecture.
  echo "Building x86_64 Linux binary (glibc, matches the $POD_IMAGE pod base) in Docker..."
  docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w \
    -e CARGO_TARGET_DIR=/w/target/amd64 rust:latest bash -c "cargo build --release"
  ls -la "$BINARY"; file "$BINARY"
}

cmd_push_data() {
  load_env
  echo "Uploading binary + parquet to b2:$B2_BUCKET ..."
  rclone copyto "$BINARY" "b2:$B2_BUCKET/bin/gr_raytracer" -P
  rclone copyto "$PARQUET" "b2:$B2_BUCKET/$PARQUET" -P
  echo "Done. (Re-run only when the binary or catalogue changes.)"
}

cmd_types() {  # list the CPU flavors the API accepts (from its own OpenAPI spec)
  load_env
  echo "Valid cpuFlavorIds (set POD_CPU_FLAVOR to one; c=compute, g=general, m=memory):"
  curl -fsS "$RUNPOD_REST/openapi.json" -H "Authorization: Bearer $RUNPOD_API_KEY" \
    | jq -r '.components.schemas.PodCreateInput.properties.cpuFlavorIds.items.enum[]' \
    | sed 's/^/  /'
}

# The command the pod runs at boot. It ALWAYS self-terminates (trap on EXIT),
# so a crash cannot leave a pod billing. Everything it needs is passed as env.
pod_command() {
  cat <<'POD'
set -uo pipefail
terminate() { curl -fsS -X DELETE "https://rest.runpod.io/v1/pods/${RUNPOD_POD_ID}" \
                -H "Authorization: Bearer ${RUNPOD_API_KEY}" >/dev/null 2>&1 || true; }
trap terminate EXIT
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq curl unzip ca-certificates >/dev/null 2>&1
curl -fsSL https://rclone.org/install.sh | bash >/dev/null 2>&1
mkdir -p /work/data && cd /work
rclone copyto "b2:${B2_BUCKET}/bin/gr_raytracer" gr_raytracer && chmod +x gr_raytracer
rclone copyto "b2:${B2_BUCKET}/${PARQUET}" "${PARQUET}"
rclone copyto "b2:${B2_BUCKET}/jobs/${JOB}/scene.toml" scene.toml
echo "started $(date -u +%FT%TZ), args=${RENDER_ARGS}, $(nproc) vCPU, pod ${RUNPOD_POD_ID}" \
  | rclone rcat "b2:${B2_BUCKET}/jobs/${JOB}/STARTED"
# Heartbeat: overwrite jobs/$JOB/progress with elapsed seconds every 30s while
# rendering. The renderer's indicatif progress bar only draws on a TTY, so on a
# headless pod this is how `status` sees the job is alive and how long it has run.
: > render.log
( t0=$(date +%s); while :; do sleep 30; \
    line="$(grep 'progress:' render.log | tail -1)"; \
    echo "elapsed $(( $(date +%s) - t0 ))s${line:+ | ${line}}" \
    | rclone rcat "b2:${B2_BUCKET}/jobs/${JOB}/progress"; done ) & HB=$!
set +e
RAYON_NUM_THREADS="$(nproc)" RUST_LOG=off ./gr_raytracer ${RENDER_ARGS} --config-file scene.toml render --filename "${OUT}" ${RENDER_SUBARGS} 2> render.log
rc=$?
set -e
kill "$HB" 2>/dev/null || true
if [ $rc -eq 0 ] && [ -f "${OUT}" ]; then
  rclone copyto "${OUT}" "b2:${B2_BUCKET}/jobs/${JOB}/${OUT}"
  echo "done $(date -u +%FT%TZ)" | rclone rcat "b2:${B2_BUCKET}/jobs/${JOB}/DONE"
else
  echo "render failed rc=$rc" | rclone rcat "b2:${B2_BUCKET}/jobs/${JOB}/FAILED"
fi
POD
}

cmd_render() {  # render <job> <scene.toml> <gr_raytracer args...>
  load_env
  local job="$1" scene="$2"; shift 2
  local render_args="$*"
  [ -n "$render_args" ] || { echo "need gr_raytracer args, e.g. --width=1280 --height=720 --exposure=2 --camera-position=-17,0,1.5 --theta=-3.14159 --psi=0 --phi=0"; exit 1; }
  echo "Uploading scene to b2:$B2_BUCKET/jobs/$job/ ..."
  rclone copyto "$scene" "b2:$B2_BUCKET/jobs/$job/scene.toml"
  for m in DONE FAILED STARTED progress; do
    rclone deletefile "b2:$B2_BUCKET/jobs/$job/$m" 2>/dev/null || true
  done

  # env handed to the pod (includes the B2 creds + api key so it can pull/push/self-terminate)
  local env_json
  env_json=$(jq -n \
    --arg JOB "$job" --arg RENDER_ARGS "$render_args" --arg B2_BUCKET "$B2_BUCKET" \
    --arg PARQUET "$PARQUET" --arg RUNPOD_API_KEY "$RUNPOD_API_KEY" \
    --arg OUT "${OUT:-out.hdr}" --arg RENDER_SUBARGS "${RENDER_SUBARGS:-}" \
    --arg T "$RCLONE_CONFIG_B2_TYPE" --arg A "$RCLONE_CONFIG_B2_ACCOUNT" --arg K "$RCLONE_CONFIG_B2_KEY" \
    '{JOB:$JOB, RENDER_ARGS:$RENDER_ARGS, B2_BUCKET:$B2_BUCKET, PARQUET:$PARQUET, OUT:$OUT, RENDER_SUBARGS:$RENDER_SUBARGS,
      RUNPOD_API_KEY:$RUNPOD_API_KEY,
      RCLONE_CONFIG_B2_TYPE:$T, RCLONE_CONFIG_B2_ACCOUNT:$A, RCLONE_CONFIG_B2_KEY:$K}')

  # CPU-pod create payload per the RunPod REST OpenAPI (rest.runpod.io/v1/openapi.json):
  # computeType=CPU, cpuFlavorIds (enum) + vcpuCount, cloudType SECURE/COMMUNITY.
  local body
  body=$(jq -n \
    --arg name "gr-render-$job" --arg image "$POD_IMAGE" --arg flavor "$POD_CPU_FLAVOR" \
    --arg cloud "$POD_CLOUD" --argjson vcpu "$POD_VCPU" \
    --argjson disk "$POD_DISK_GB" --argjson env "$env_json" \
    --arg cmd "$(pod_command)" \
    '{name:$name, imageName:$image, computeType:"CPU", cloudType:$cloud,
      cpuFlavorIds:[$flavor], cpuFlavorPriority:"availability", vcpuCount:$vcpu,
      containerDiskInGb:$disk, env:$env,
      dockerStartCmd:["bash","-lc",$cmd]}')

  echo "Creating CPU pod..."
  local resp; resp=$(rp POST "/pods" "$body")
  echo "$resp" | jq . 2>/dev/null || echo "$resp"
  echo "Pod launching. It will render, push to b2:$B2_BUCKET/jobs/$job/${OUT:-out.hdr}, and self-terminate."
  echo "Fetch when ready:  $0 fetch $job"
}

cmd_fetch() {  # fetch <job> [local-name] : downloads whatever out.* the pod wrote
  load_env
  local job="$1"
  if b2_has "$job" FAILED; then
    echo "job $job FAILED on the pod; check its logs. FAILED marker present."; exit 1
  fi
  if ! b2_has "$job" DONE; then
    echo "not done yet (no DONE marker). Re-run later."; exit 2
  fi
  local remote; remote=$(rclone lsf "b2:$B2_BUCKET/jobs/$job/" 2>/dev/null | grep -E '^out\.' | head -1)
  [ -n "$remote" ] || { echo "DONE but no out.* found in jobs/$job/"; exit 1; }
  local out="${2:-${job}_${remote}}"
  rclone copyto "b2:$B2_BUCKET/jobs/$job/$remote" "$out" -P
  echo "Fetched -> $out"
}

cmd_status() {  # status [job] : how many pods are running + per-job progress
  load_env
  local pods; pods=$(rp GET "/pods")
  local n; n=$(echo "$pods" | jq 'length')
  echo "RunPod pods live: $n"
  if [ "$n" -gt 0 ]; then
    # id, name, state, cpu flavor x vcpu, $/hr, when it entered its current state
    echo "$pods" | jq -r '.[] | "  \(.id)  \(.name)  \(.desiredStatus)  \(.cpuFlavorId // "?")x\(.vcpuCount // "?")vcpu  $\(.costPerHr // .adjustedCostPerHr // 0)/hr  since \(.lastStatusChange // "?")"'
  fi
  echo
  echo "Jobs in b2:$B2_BUCKET/jobs/ :"
  local jobs; jobs=$(rclone lsf "b2:$B2_BUCKET/jobs/" 2>/dev/null | sed 's#/$##')
  [ -n "$jobs" ] || { echo "  (none)"; return 0; }
  local j state prog
  while IFS= read -r j; do
    [ -n "$j" ] || continue
    if   b2_has "$j" DONE;    then state="DONE";
    elif b2_has "$j" FAILED;  then state="FAILED";
    elif b2_has "$j" STARTED; then state="RUNNING";
    else state="queued"; fi
    prog=""
    if [ "$state" = "RUNNING" ] && b2_has "$j" progress; then
      prog=" ($(rclone cat "b2:$B2_BUCKET/jobs/$j/progress" 2>/dev/null))"
    fi
    echo "  $j: $state$prog"
  done <<< "$jobs"
}

cmd_kill() {  # kill <podId> : delete a pod by hand (backstop if self-terminate didn't fire)
  load_env
  local id="${1:?usage: $0 kill <podId>  (get ids from: $0 status)}"
  rp DELETE "/pods/$id" >/dev/null && echo "deleted pod $id"
}

case "${1:-}" in
  build)      cmd_build ;;
  push-data)  cmd_push_data ;;
  types)      cmd_types ;;
  render)     shift; cmd_render "$@" ;;
  fetch)      shift; cmd_fetch "$@" ;;
  status|ps)  shift; cmd_status "$@" ;;
  kill)       shift; cmd_kill "$@" ;;
  *) cat <<EOF
usage: $0 <command>
  build                          build the x86_64 Linux binary in Docker
  push-data                      upload binary + parquet to B2 (re-run on change)
  types                          list valid RunPod cpuFlavorIds (set POD_CPU_FLAVOR)
  render <job> <scene.toml> <gr_raytracer args...>
                                 e.g. render hero scene.toml --width=1280 --height=720 \\
                                      --exposure=2 --camera-position=-17,0,1.5 \\
                                      --theta=-3.14159 --psi=0 --phi=0
  status | ps [job]              how many pods are live + per-job progress
  kill <podId>                   delete a pod by hand (ids from status)
  fetch <job> [out.hdr]          download the result once DONE
EOF
  ;;
esac
