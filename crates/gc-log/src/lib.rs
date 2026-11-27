//! `gc-log`: ein Ereignisstrom je Lauf (Schema `genome.log.v1`).
//!
//! Regeln (Plan Abschnitt 15):
//! - Jedes Ereignis ist eine JSON-Zeile in `run.log.jsonl` und zugleich eine lesbare Zeile auf stderr.
//! - Jede Zeile wird sofort geschrieben, damit nach einem Abbruch alles bis zur letzten Zeile da ist.
//! - Fehler tragen einen festen Code, Ort, Ursache und Handlungsempfehlung.
//! - Protokolle enthalten keine Sequenzen, nur Zahlen, Pfade und Namen.

use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use serde_json::{Map, Value, json};

pub const SCHEMA: &str = "genome.log.v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl Level {
    pub fn as_str(self) -> &'static str {
        match self {
            Level::Error => "error",
            Level::Warn => "warn",
            Level::Info => "info",
            Level::Debug => "debug",
            Level::Trace => "trace",
        }
    }

    pub fn parse(s: &str) -> Option<Level> {
        Some(match s.to_ascii_lowercase().as_str() {
            "error" => Level::Error,
            "warn" => Level::Warn,
            "info" => Level::Info,
            "debug" => Level::Debug,
            "trace" => Level::Trace,
            _ => return None,
        })
    }
}

/// Einstellungen für einen Lauf.
pub struct Config {
    /// Name des Werkzeugs, z. B. `index.build`.
    pub tool: String,
    /// Programmversion.
    pub version: String,
    /// Ab dieser Stufe wird geschrieben.
    pub level: Level,
    /// Zieldatei für die JSON-Zeilen (wird angelegt bzw. überschrieben).
    pub jsonl: Option<PathBuf>,
    /// Lesbare Zeilen auf stderr.
    pub stderr: bool,
}

struct Sink {
    jsonl: Option<File>,
    stderr: bool,
    capture: Option<Vec<String>>,
}

struct Logger {
    tool: String,
    version: String,
    run_id: String,
    level: Level,
    started: Instant,
    sink: Mutex<Sink>,
}

static LOGGER: OnceLock<Logger> = OnceLock::new();

/// Zeit in UTC nach RFC 3339 mit Millisekunden, ohne fremde Crate.
fn rfc3339(ms: u128) -> String {
    let secs = (ms / 1000) as i64;
    let millis = (ms % 1000) as u32;
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    // Zivildatum aus Tagen seit 1970-01-01 (Verfahren nach Howard Hinnant).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}.{millis:03}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

fn now_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

/// Startet das Protokoll des Laufs. Ein zweiter Aufruf ändert nichts und liefert `false`.
pub fn init(config: Config) -> io::Result<bool> {
    let jsonl = match &config.jsonl {
        Some(p) => Some(File::create(p).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Protokolldatei {} nicht anlegbar: {e}", p.display()),
            )
        })?),
        None => None,
    };
    let run_id = format!("{:x}-{:x}", now_ms(), std::process::id());
    let logger = Logger {
        tool: config.tool,
        version: config.version,
        run_id,
        level: config.level,
        started: Instant::now(),
        sink: Mutex::new(Sink {
            jsonl,
            stderr: config.stderr,
            capture: None,
        }),
    };
    Ok(LOGGER.set(logger).is_ok())
}

/// Nur für Tests: Protokoll in den Speicher statt in eine Datei.
pub fn init_capture(tool: &str, level: Level) -> bool {
    let logger = Logger {
        tool: tool.to_string(),
        version: "test".to_string(),
        run_id: "test".to_string(),
        level,
        started: Instant::now(),
        sink: Mutex::new(Sink {
            jsonl: None,
            stderr: false,
            capture: Some(Vec::new()),
        }),
    };
    LOGGER.set(logger).is_ok()
}

/// Nur für Tests: bisher geschriebene JSON-Zeilen.
pub fn captured() -> Vec<String> {
    LOGGER
        .get()
        .and_then(|l| l.sink.lock().ok().and_then(|s| s.capture.clone()))
        .unwrap_or_default()
}

pub fn run_id() -> Option<&'static str> {
    LOGGER.get().map(|l| l.run_id.as_str())
}

pub fn enabled(level: Level) -> bool {
    LOGGER.get().is_some_and(|l| level <= l.level)
}

fn human_fields(fields: &Value) -> String {
    let mut out = String::new();
    if let Value::Object(map) = fields {
        for (k, v) in map {
            out.push(' ');
            out.push_str(k);
            out.push('=');
            match v {
                Value::String(s) => out.push_str(s),
                other => out.push_str(&other.to_string()),
            }
        }
    }
    out
}

/// Schreibt ein Ereignis. `code` ist ein fester, dokumentierter Name wie `INDEX.SA.FERTIG`.
/// `fields` ist ein JSON-Objekt mit Zahlen, Pfaden und Namen (keine Sequenzen).
pub fn event(level: Level, code: &str, message: &str, fields: Value) {
    let Some(l) = LOGGER.get() else { return };
    if level > l.level {
        return;
    }
    let t_ms = l.started.elapsed().as_millis() as u64;
    let mut obj = Map::new();
    obj.insert("schema".into(), json!(SCHEMA));
    obj.insert("time".into(), json!(rfc3339(now_ms())));
    obj.insert("t_ms".into(), json!(t_ms));
    obj.insert("run".into(), json!(l.run_id));
    obj.insert("tool".into(), json!(l.tool));
    obj.insert("version".into(), json!(l.version));
    obj.insert("level".into(), json!(level.as_str()));
    obj.insert("code".into(), json!(code));
    obj.insert("msg".into(), json!(message));
    if fields.as_object().is_some_and(|m| !m.is_empty()) {
        obj.insert("fields".into(), fields.clone());
    }
    let line = Value::Object(obj).to_string();
    let Ok(mut sink) = l.sink.lock() else { return };
    if let Some(f) = sink.jsonl.as_mut() {
        // Jede Zeile sofort schreiben: nach einem Abbruch ist alles bis hierher lesbar.
        let _ = f.write_all(line.as_bytes());
        let _ = f.write_all(b"\n");
        let _ = f.flush();
    }
    if let Some(c) = sink.capture.as_mut() {
        c.push(line);
    }
    if sink.stderr {
        let secs = t_ms as f64 / 1000.0;
        eprintln!(
            "[{secs:9.3}s] {:<5} {code}: {message}{}",
            level.as_str().to_ascii_uppercase(),
            human_fields(&fields)
        );
    }
}

pub fn info(code: &str, message: &str, fields: Value) {
    event(Level::Info, code, message, fields);
}

pub fn warn(code: &str, message: &str, fields: Value) {
    event(Level::Warn, code, message, fields);
}

pub fn debug(code: &str, message: &str, fields: Value) {
    event(Level::Debug, code, message, fields);
}

/// Eine automatische Wahl mit Grund festhalten ("kompakter Index, weil 16 GB RAM").
pub fn decision(code: &str, choice: &str, reason: &str, mut fields: Value) {
    if let Value::Object(m) = &mut fields {
        m.insert("wahl".into(), json!(choice));
        m.insert("grund".into(), json!(reason));
    } else {
        fields = json!({ "wahl": choice, "grund": reason });
    }
    event(
        Level::Info,
        code,
        &format!("{choice} (Grund: {reason})"),
        fields,
    );
}

/// Ein Fehler, der sich selbst erklärt: was, wo, wahrscheinliche Ursache, was zu tun ist.
pub struct Failure<'a> {
    pub code: &'a str,
    pub what: &'a str,
    pub location: Value,
    pub cause: &'a str,
    pub action: &'a str,
}

pub fn failure(f: Failure<'_>) {
    event(
        Level::Error,
        f.code,
        f.what,
        json!({ "wo": f.location, "ursache": f.cause, "was_tun": f.action }),
    );
}

/// Laufkopf: Befehlszeile, Parameter, System. Als erstes Ereignis jedes Laufs schreiben.
pub fn run_header(parameters: Value) {
    let args: Vec<String> = std::env::args().collect();
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    info(
        "LAUF.START",
        "Lauf beginnt",
        json!({
            "befehl": args,
            "parameter": parameters,
            "system": std::env::consts::OS,
            "architektur": std::env::consts::ARCH,
            "kerne": cores,
            "arbeitsordner": std::env::current_dir().ok().map(|p| p.display().to_string()),
        }),
    );
}

/// Eingabedatei mit Größe und Änderungszeit festhalten.
pub fn input_file(role: &str, path: &Path) {
    let meta = std::fs::metadata(path).ok();
    let mtime = meta
        .as_ref()
        .and_then(|m| m.modified().ok())
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| rfc3339(d.as_millis()));
    info(
        "LAUF.EINGABE",
        "Eingabedatei",
        json!({
            "rolle": role,
            "pfad": path.display().to_string(),
            "bytes": meta.as_ref().map(|m| m.len()),
            "geaendert": mtime,
        }),
    );
}

/// Abschluss des Laufs. `ok = false` nach einem Fehler.
pub fn run_end(ok: bool, summary: Value) {
    let (level, code, msg) = if ok {
        (Level::Info, "LAUF.FERTIG", "Lauf fertig")
    } else {
        (Level::Error, "LAUF.FEHLGESCHLAGEN", "Lauf fehlgeschlagen")
    };
    event(level, code, msg, summary);
}

/// Eine Stufe mit Anfang, Ende und Dauer. Ohne `finish` gilt die Stufe als abgebrochen
/// und wird beim Verlassen so protokolliert.
pub struct Stage {
    name: String,
    started: Instant,
    last_beat: Instant,
    done: bool,
}

pub fn stage(name: &str, fields: Value) -> Stage {
    event(Level::Info, "STUFE.START", name, fields);
    let now = Instant::now();
    Stage {
        name: name.to_string(),
        started: now,
        last_beat: now,
        done: false,
    }
}

impl Stage {
    /// Lebenszeichen, höchstens alle `every_secs` Sekunden. Billig genug für innere Schleifen,
    /// wenn es nur alle paar tausend Durchläufe gerufen wird.
    pub fn heartbeat(&mut self, every_secs: u64, fields: impl FnOnce() -> Value) {
        if self.last_beat.elapsed().as_secs() >= every_secs {
            self.last_beat = Instant::now();
            let mut f = fields();
            if let Value::Object(m) = &mut f {
                m.insert(
                    "sekunden".into(),
                    json!(self.started.elapsed().as_secs_f64()),
                );
                if let Some(b) = gc_sys::peak_memory_bytes() {
                    m.insert("speicher_spitze_bytes".into(), json!(b));
                }
            }
            event(Level::Info, "STUFE.LEBT", &self.name, f);
        }
    }

    pub fn finish(mut self, mut fields: Value) {
        self.done = true;
        let secs = self.started.elapsed().as_secs_f64();
        if let Value::Object(m) = &mut fields {
            m.insert("sekunden".into(), json!(secs));
        } else {
            fields = json!({ "sekunden": secs });
        }
        if let (Value::Object(m), Some(b)) = (&mut fields, gc_sys::peak_memory_bytes()) {
            m.insert("speicher_spitze_bytes".into(), json!(b));
        }
        event(Level::Info, "STUFE.ENDE", &self.name, fields);
    }
}

impl Drop for Stage {
    fn drop(&mut self) {
        if !self.done {
            event(
                Level::Error,
                "STUFE.ABGEBROCHEN",
                &self.name,
                json!({ "sekunden": self.started.elapsed().as_secs_f64() }),
            );
        }
    }
}

/// Programmabstürze (panic) mit Ort ins Protokoll schreiben, bevor das Programm endet.
pub fn install_panic_hook() {
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let location = info
            .location()
            .map(|l| format!("{}:{}", l.file(), l.line()));
        let text = info
            .payload()
            .downcast_ref::<&str>()
            .map(|s| s.to_string())
            .or_else(|| info.payload().downcast_ref::<String>().cloned())
            .unwrap_or_else(|| "unbekannt".to_string());
        failure(Failure {
            code: "PROGRAMM.ABSTURZ",
            what: &format!("Programmfehler: {text}"),
            location: json!({ "quelltext": location }),
            cause: "Fehler im Programm, nicht in den Eingabedaten",
            action: "Protokoll und Befehlszeile aufbewahren und als Fehler melden",
        });
        previous(info);
    }));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rfc3339_known_dates() {
        assert_eq!(rfc3339(0), "1970-01-01T00:00:00.000Z");
        assert_eq!(rfc3339(951_782_400_000), "2000-02-29T00:00:00.000Z");
        assert_eq!(rfc3339(1_795_010_400_123), "2026-11-18T14:00:00.123Z");
    }

    #[test]
    fn events_are_valid_schema_lines() {
        init_capture("test", Level::Info);
        run_header(json!({ "x": 1 }));
        decision("TEST.WAHL", "kompakt", "wenig RAM", json!({ "ram_gb": 16 }));
        let mut st = stage("zaehlen", json!({}));
        st.heartbeat(0, || json!({ "fortschritt": 0.5 }));
        st.finish(json!({ "datensaetze": 3 }));
        {
            let _abgebrochen = stage("bricht ab", json!({}));
        }
        debug(
            "TEST.DEBUG",
            "darf bei Stufe info nicht erscheinen",
            json!({}),
        );
        failure(Failure {
            code: "TEST.FEHLER",
            what: "Datei abgeschnitten",
            location: json!({ "datei": "a.fq", "zeile": 7 }),
            cause: "Kopie unvollständig",
            action: "neu kopieren",
        });
        run_end(false, json!({}));

        let lines = captured();
        let codes: Vec<String> = lines
            .iter()
            .map(|l| {
                let v: Value = serde_json::from_str(l).expect("jede Zeile ist gültiges JSON");
                assert_eq!(v["schema"], SCHEMA);
                for key in [
                    "time", "t_ms", "run", "tool", "version", "level", "code", "msg",
                ] {
                    assert!(!v[key].is_null(), "Pflichtfeld {key} fehlt in {l}");
                }
                v["code"].as_str().unwrap().to_string()
            })
            .collect();
        assert_eq!(
            codes,
            [
                "LAUF.START",
                "TEST.WAHL",
                "STUFE.START",
                "STUFE.LEBT",
                "STUFE.ENDE",
                "STUFE.START",
                "STUFE.ABGEBROCHEN",
                "TEST.FEHLER",
                "LAUF.FEHLGESCHLAGEN"
            ]
        );
        let err: Value = serde_json::from_str(&lines[7]).unwrap();
        assert_eq!(err["fields"]["wo"]["zeile"], 7);
        assert_eq!(err["fields"]["was_tun"], "neu kopieren");
    }
}
