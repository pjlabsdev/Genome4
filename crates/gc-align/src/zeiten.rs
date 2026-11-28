//! Zeitmessung je Paarart für die Forschung am Tempo (`GC_ZEITEN=1`). Ohne Schalter kostet sie
//! nichts und ändert nie ein Ergebnis. Arten: 0 = beide Reads über Schnell- oder Strobeweg,
//! 1 = einer über den vollen Weg, 2 = beide über den vollen Weg.

use std::sync::atomic::{AtomicU64, Ordering};

pub fn on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_ZEITEN").is_ok_and(|v| v == "1"))
}

#[allow(clippy::declare_interior_mutable_const)] // nur als Anfangswert der Felder unten
const Z: AtomicU64 = AtomicU64::new(0);
/// Nanosekunden und Anzahl: schneller Versuch je Read, voller Weg je Read.
pub static QUICK: [AtomicU64; 2] = [Z, Z];
pub static FULL: [AtomicU64; 2] = [Z, Z];
/// Je Paarart: Nanosekunden und Zahl der Fenster der Partnersuche.
pub static RESCUE: [[AtomicU64; 2]; 3] = [[Z, Z], [Z, Z], [Z, Z]];
/// Je Paarart: Nanosekunden und Zahl der Paare für Paarbewertung und Zeilen.
pub static LINES: [[AtomicU64; 2]; 3] = [[Z, Z], [Z, Z], [Z, Z]];
/// Paare je Art.
pub static PAIRS: [AtomicU64; 3] = [Z, Z, Z];

pub fn add(c: &[AtomicU64; 2], ns: u64, n: u64) {
    c[0].fetch_add(ns, Ordering::Relaxed);
    c[1].fetch_add(n, Ordering::Relaxed);
}

/// Für das Protokoll: Name, Wert.
pub fn counts() -> Vec<(String, u64)> {
    let g = |a: &AtomicU64| a.load(Ordering::Relaxed);
    let mut v = vec![
        ("schnell_ns".into(), g(&QUICK[0])),
        ("schnell_reads".into(), g(&QUICK[1])),
        ("voll_ns".into(), g(&FULL[0])),
        ("voll_reads".into(), g(&FULL[1])),
    ];
    for k in 0..3 {
        v.push((format!("paare_{k}"), g(&PAIRS[k])));
        v.push((format!("partnersuche_{k}_ns"), g(&RESCUE[k][0])));
        v.push((format!("partnersuche_{k}_fenster"), g(&RESCUE[k][1])));
        v.push((format!("zeilen_{k}_ns"), g(&LINES[k][0])));
    }
    v
}
