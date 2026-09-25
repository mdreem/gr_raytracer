# Baked render image for scripts/render-remote.sh: rclone + the x86_64
# gr_raytracer binary + the Gaia parquet, so a pod boots in seconds with no B2
# download (which is the flaky window). Built via `render-remote.sh build-image`,
# which assembles a minimal build context (bin/, data/) around this file.
FROM debian:bookworm-slim
RUN apt-get update -qq \
 && apt-get install -y -qq curl unzip ca-certificates >/dev/null \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL https://rclone.org/install.sh | bash >/dev/null
WORKDIR /work
COPY bin/gr_raytracer /work/gr_raytracer
COPY data/ /work/data/
RUN chmod +x /work/gr_raytracer
