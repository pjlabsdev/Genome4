//! Das lokale Alignment muss in jedem Wert dasselbe liefern wie `ksw_align2` aus bwa.
//! Die Fälle samt Ergebnissen von bwa stehen in `daten/lokal_bwa.tsv`, erzeugt mit
//! `bench/m3/make_local_cases.py`.

use gc_align::{Scoring, local_align};

#[test]
fn alle_faelle_wie_bwa() {
    let text = include_str!("daten/lokal_bwa.tsv");
    let sc = Scoring::default();
    let mut wrong = Vec::new();
    let mut n = 0;
    for (nr, line) in text.lines().enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let codes = |s: &str| -> Vec<u8> { s.bytes().map(|b| b - b'0').collect() };
        let (q, t) = (codes(f[0]), codes(f[1]));
        let want: Vec<i32> = f[3..10].iter().map(|x| x.parse().unwrap()).collect();
        let r = local_align(&q, &t, &sc, f[2] == "1", 19, true);
        let got = vec![r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb];
        n += 1;
        if got != want {
            wrong.push(format!("Fall {}: bwa {want:?}, wir {got:?}", nr + 1));
        }
    }
    assert!(n >= 1000);
    assert!(
        wrong.is_empty(),
        "{} von {n} Fällen weichen ab:\n{}",
        wrong.len(),
        wrong[..wrong.len().min(15)].join("\n")
    );
}
