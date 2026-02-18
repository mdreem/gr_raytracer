#[path = "../src/configuration.rs"]
mod configuration;
#[path = "../src/rendering/color.rs"]
pub mod rendering_color;
mod rendering {
    pub use crate::rendering_color as color;
}

use configuration::RenderConfig;
use schemars::schema_for;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use std::error::Error;
use std::fs;
use std::path::Path;

fn main() -> Result<(), Box<dyn Error>> {
    let schema = schema_for!(RenderConfig);
    let schema_value = serde_json::to_value(&schema)?;
    let docs_dir = Path::new("docs");
    fs::create_dir_all(docs_dir)?;

    let schema_json_path = docs_dir.join("config.schema.json");
    let schema_markdown_path = docs_dir.join("config-schema.md");

    fs::write(
        &schema_json_path,
        format!("{}\n", serde_json::to_string_pretty(&schema_value)?),
    )?;
    fs::write(&schema_markdown_path, render_markdown(&schema_value))?;

    println!(
        "Wrote {} and {}",
        schema_json_path.display(),
        schema_markdown_path.display()
    );

    Ok(())
}

fn render_markdown(root: &Value) -> String {
    let mut out = String::new();
    out.push_str("# Render Config Schema\n\n");
    out.push_str("Generated from Rust types in `src/configuration.rs`.\n\n");

    let defs = root
        .get("$defs")
        .and_then(Value::as_object)
        .cloned()
        .unwrap_or_default();

    render_schema_section(
        "RenderConfig",
        root,
        root.get("required"),
        &defs,
        &mut out,
        true,
    );

    let mut rendered = BTreeSet::new();
    rendered.insert(String::from("RenderConfig"));
    for key in sorted_keys(&defs) {
        if rendered.contains(key.as_str()) {
            continue;
        }
        if let Some(schema) = defs.get(key.as_str()) {
            render_schema_section(
                key.as_str(),
                schema,
                schema.get("required"),
                &defs,
                &mut out,
                false,
            );
            rendered.insert(key);
        }
    }

    out
}

fn render_schema_section(
    title: &str,
    schema: &Value,
    required: Option<&Value>,
    defs: &Map<String, Value>,
    out: &mut String,
    include_title: bool,
) {
    if include_title {
        out.push_str("## RenderConfig\n\n");
    } else {
        out.push_str(&format!("## {}\n\n", title));
    }

    if let Some(desc) = schema.get("description").and_then(Value::as_str) {
        out.push_str(desc);
        out.push_str("\n\n");
    }

    if let Some(enum_values) = schema.get("enum").and_then(Value::as_array) {
        out.push_str("Allowed values: ");
        out.push_str(&format_enum_values(enum_values));
        out.push_str("\n\n");
    }

    if let Some(one_of) = schema.get("oneOf").and_then(Value::as_array) {
        out.push_str("Variants:\n\n");
        for variant in one_of {
            render_variant(variant, defs, out);
        }
        out.push('\n');
    }

    if let Some(properties) = schema.get("properties").and_then(Value::as_object) {
        let required_set = required_set(required);
        out.push_str("| Field | Type | Required | Description |\n");
        out.push_str("| --- | --- | --- | --- |\n");
        for name in sorted_keys(properties) {
            if let Some(property_schema) = properties.get(name.as_str()) {
                let ty = summarize_type(property_schema, defs)
                    .unwrap_or_else(|| String::from("unknown"));
                let req = if required_set.contains(name.as_str()) {
                    "yes"
                } else {
                    "no"
                };
                let desc = property_schema
                    .get("description")
                    .and_then(Value::as_str)
                    .unwrap_or("");
                out.push_str(&format!("| `{}` | `{}` | {} | {} |\n", name, ty, req, desc));
            }
        }
        out.push('\n');
    }
}

fn sorted_keys(map: &Map<String, Value>) -> Vec<String> {
    let mut keys: Vec<String> = map.keys().cloned().collect();
    keys.sort();
    keys
}

fn required_set(required: Option<&Value>) -> BTreeSet<String> {
    required
        .and_then(Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(Value::as_str)
                .map(String::from)
                .collect()
        })
        .unwrap_or_default()
}

fn summarize_type(schema: &Value, defs: &Map<String, Value>) -> Option<String> {
    if let Some(constant) = schema.get("const") {
        return Some(match constant {
            Value::String(s) => s.clone(),
            _ => constant.to_string(),
        });
    }

    if let Some(reference) = schema.get("$ref").and_then(Value::as_str) {
        return Some(reference_name(reference));
    }

    if let Some(type_name) = schema.get("type").and_then(Value::as_str) {
        if type_name == "array" {
            if let Some(items) = schema.get("items") {
                return Some(format!(
                    "array<{}>",
                    summarize_type(items, defs).unwrap_or_else(|| String::from("unknown"))
                ));
            }
        }
        return Some(type_name.to_string());
    }

    if let Some(type_names) = schema.get("type").and_then(Value::as_array) {
        let parts: Vec<String> = type_names
            .iter()
            .filter_map(Value::as_str)
            .map(String::from)
            .collect();
        if !parts.is_empty() {
            return Some(parts.join(" | "));
        }
    }

    if let Some(enum_values) = schema.get("enum").and_then(Value::as_array) {
        return Some(format!("enum {}", format_enum_values(enum_values)));
    }

    if let Some(any_of) = schema.get("anyOf").and_then(Value::as_array) {
        let parts: Vec<String> = any_of
            .iter()
            .filter_map(|schema| summarize_type(schema, defs))
            .collect();
        if !parts.is_empty() {
            return Some(parts.join(" | "));
        }
    }

    if let Some(one_of) = schema.get("oneOf").and_then(Value::as_array) {
        let parts: Vec<String> = one_of
            .iter()
            .filter_map(|schema| summarize_type(schema, defs))
            .collect();
        if !parts.is_empty() {
            return Some(parts.join(" | "));
        }
    }

    if let Some(properties) = schema.get("properties").and_then(Value::as_object) {
        if properties.len() == 1 {
            if let Some((name, value)) = properties.iter().next() {
                if let Some(reference) = value.get("$ref").and_then(Value::as_str) {
                    return Some(format!("{} ({})", name, reference_name(reference)));
                }
                return Some(name.to_string());
            }
        }
        return Some(String::from("object"));
    }

    if let Some(reference) = schema.get("$ref").and_then(Value::as_str) {
        let def_name = reference_name(reference);
        if let Some(def_schema) = defs.get(def_name.as_str()) {
            return summarize_type(def_schema, defs);
        }
        return Some(def_name);
    }

    None
}

fn render_variant(variant: &Value, defs: &Map<String, Value>, out: &mut String) {
    if let Some(constant) = variant.get("const").and_then(Value::as_str) {
        out.push_str(&format!("### `{}`\n\n", constant));
        if let Some(desc) = variant.get("description").and_then(Value::as_str) {
            out.push_str(desc);
            out.push('\n');
        }
        out.push('\n');
        return;
    }

    if let Some(properties) = variant.get("properties").and_then(Value::as_object) {
        if properties.len() == 1 {
            if let Some((variant_name, payload_schema)) = properties.iter().next() {
                out.push_str(&format!("### `{}`\n\n", variant_name));
                if let Some(desc) = variant.get("description").and_then(Value::as_str) {
                    out.push_str(desc);
                    out.push('\n');
                }
                out.push('\n');

                if let Some(payload_properties) =
                    payload_schema.get("properties").and_then(Value::as_object)
                {
                    let required = required_set(payload_schema.get("required"));
                    out.push_str("| Field | Type | Required | Description |\n");
                    out.push_str("| --- | --- | --- | --- |\n");
                    for name in sorted_keys(payload_properties) {
                        if let Some(property_schema) = payload_properties.get(name.as_str()) {
                            let ty = summarize_type(property_schema, defs)
                                .unwrap_or_else(|| String::from("unknown"));
                            let req = if required.contains(name.as_str()) {
                                "yes"
                            } else {
                                "no"
                            };
                            let desc = property_schema
                                .get("description")
                                .and_then(Value::as_str)
                                .unwrap_or("");
                            out.push_str(&format!(
                                "| `{}` | `{}` | {} | {} |\n",
                                name, ty, req, desc
                            ));
                        }
                    }
                    out.push('\n');
                }
                return;
            }
        }
    }

    out.push_str(&format!(
        "### `{}`\n\n",
        summarize_type(variant, defs).unwrap_or_else(|| String::from("variant"))
    ));
    if let Some(desc) = variant.get("description").and_then(Value::as_str) {
        out.push_str(desc);
        out.push('\n');
    }
    out.push('\n');
}

fn reference_name(reference: &str) -> String {
    reference
        .rsplit('/')
        .next()
        .unwrap_or(reference)
        .to_string()
}

fn format_enum_values(values: &[Value]) -> String {
    values
        .iter()
        .map(|value| match value {
            Value::String(s) => format!("`{}`", s),
            _ => format!("`{}`", value),
        })
        .collect::<Vec<_>>()
        .join(", ")
}
