//! gc-schema: versionierte JSON-Ausgaben aller Genome-Werkzeuge.
//!
//! Jede Ausgabe trägt `schema` (z. B. `genome.fastq-check.v1`), Werkzeug und Version,
//! damit Ergebnisse auch Jahre später eindeutig zuzuordnen sind.

use serde::Serialize;

/// Hülle um jedes Werkzeug-Ergebnis.
#[derive(Debug, Serialize)]
pub struct Envelope<T: Serialize> {
    pub schema: &'static str,
    pub tool: &'static str,
    pub version: &'static str,
    pub result: T,
}

impl<T: Serialize> Envelope<T> {
    pub fn new(schema: &'static str, result: T) -> Self {
        Self {
            schema,
            tool: "genome",
            version: env!("CARGO_PKG_VERSION"),
            result,
        }
    }

    /// Für Programme neben `genome` (z. B. `genome-validate`).
    pub fn for_tool(tool: &'static str, schema: &'static str, result: T) -> Self {
        Self {
            schema,
            tool,
            version: env!("CARGO_PKG_VERSION"),
            result,
        }
    }

    pub fn to_json(&self) -> serde_json::Result<String> {
        serde_json::to_string_pretty(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn envelope_carries_schema_and_version() {
        let j = Envelope::new("genome.test.v1", 42).to_json().unwrap();
        let v: serde_json::Value = serde_json::from_str(&j).unwrap();
        assert_eq!(v["schema"], "genome.test.v1");
        assert_eq!(v["result"], 42);
        assert_eq!(v["version"], env!("CARGO_PKG_VERSION"));
    }
}
