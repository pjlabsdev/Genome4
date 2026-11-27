//! Trimmen von Read-Paaren mit denselben Regeln wie fastp 1.4.0 in der Einstellung der alten
//! Kette (`--detect_adapter_for_pe`, sonst Vorgaben). Jede Regel, die Quelle in fastp und die
//! bewusst nachgebauten Eigenheiten stehen im Handbuch unter `verfahren/trimmen.md`.
//!
//! Ablauf: Vor dem Lauf entscheidet [`two_color`] über das Kürzen von PolyG und
//! [`detect_adapter`] sucht je Read-Seite die Adapterfolge im Anfang der Datei. Dann geht jedes
//! Paar durch [`Trimmer::trim_pair`], das die Reads kürzt und sagt, ob das Paar bleibt und warum.

use crate::fastq::Record;
use crate::known_adapters::KNOWN_ADAPTERS;

/// Einstellungen, Vorgaben wie fastp 1.4.0.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Settings {
    /// Kürzeste Überlappung der beiden Reads (`overlap_len_require`).
    pub overlap_require: usize,
    /// Höchstens so viele Unterschiede in der Überlappung (`overlap_diff_limit`).
    pub overlap_diff_limit: usize,
    /// Höchstens so viel Prozent Unterschiede in der Überlappung (`overlap_diff_percent_limit`).
    pub overlap_diff_percent: u32,
    /// Paare, die nach dem Entfernen von Adaptern höchstens so lang sind, gelten als Dimer.
    pub dimer_max_len: usize,
    /// Mindestlänge eines PolyG-Endes (`poly_g_min_len`).
    pub poly_g_min_len: usize,
    /// Base gilt ab dieser Phred-Qualität als gut (`qualified_quality_phred`).
    pub qualified_phred: u8,
    /// Höchstens so viel Prozent schlechte Basen je Read (`unqualified_percent_limit`).
    pub unqualified_percent: u32,
    /// Höchstens so viele N je Read (`n_base_limit`).
    pub n_base_limit: usize,
    /// Kürzeste erlaubte Read-Länge nach dem Kürzen (`length_required`).
    pub length_required: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            overlap_require: 30,
            overlap_diff_limit: 5,
            overlap_diff_percent: 20,
            dimer_max_len: 2,
            poly_g_min_len: 10,
            qualified_phred: 15,
            unqualified_percent: 40,
            n_base_limit: 5,
            length_required: 15,
        }
    }
}

/// Gerätekennungen mit zweifarbiger Chemie (`Evaluator::isTwoColorSystem`), ohne `@`.
const TWO_COLOR_PREFIXES: [&str; 11] = [
    "FS", "MN", "SH", "NS", "NB", "NDX", "VL", "VH", "A", "NA", "LH",
];

/// PolyG kürzen? fastp entscheidet allein am Namen des ersten Reads von R1 (ohne `@`).
pub fn two_color(first_name: &[u8]) -> Option<&'static str> {
    TWO_COLOR_PREFIXES
        .iter()
        .find(|p| first_name.starts_with(p.as_bytes()))
        .copied()
}

/// So viele Reads liest fastp höchstens für die Adaptersuche.
pub const DETECT_READS: usize = 256 * 1024;
/// So viele Basen liest fastp höchstens für die Adaptersuche.
pub const DETECT_BASES: usize = 151 * DETECT_READS;
/// Unter so vielen Reads sucht fastp keinen Adapter.
const DETECT_MIN_READS: usize = 10_000;

/// Wie eine Adapterfolge gefunden wurde.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Found {
    /// Ein bekannter Adapter traf genug Reads; Zahl der Treffer und der geprüften Reads.
    Known { hits: u64, checked: u64 },
    /// Aus der häufigsten Teilfolge verlängert; passt zu einem bekannten Adapter, wenn `known`.
    Seed { seed: String, known: bool },
}

/// Ergebnis der Adaptersuche für eine Read-Seite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    pub seq: Vec<u8>,
    /// Name aus der Liste bekannter Adapter, falls es einer ist.
    pub name: Option<&'static str>,
    pub found: Found,
}

/// Sucht die Adapterfolge in den ersten Reads einer Seite (`Evaluator::evalAdapterAndReadNum`).
/// `reads` sind die Sequenzen, so wie der Aufrufer sie bis [`DETECT_READS`] oder
/// [`DETECT_BASES`] gelesen hat (fastp hört auf, sobald eine der Grenzen erreicht ist).
/// `Err` mit dem Grund, wenn keiner gilt.
pub fn detect_adapter(reads: &[Vec<u8>]) -> Result<Adapter, String> {
    if reads.len() < DETECT_MIN_READS {
        return Err(format!(
            "nur {} Reads, fastp sucht erst ab {DETECT_MIN_READS}",
            reads.len()
        ));
    }
    let adapter = match check_known_adapters(reads) {
        Some((seq, hits, checked)) if seq.len() > 8 => Some(Adapter {
            name: known_name(seq),
            seq: seq.as_bytes().to_vec(),
            found: Found::Known { hits, checked },
        }),
        _ => search_by_seed(reads),
    };
    match adapter {
        // fastp kürzt längere Folgen mit `resize(0, 60)` auf die Länge 0 statt auf 60: Eine
        // erkannte Folge über 60 Basen gilt dort als "kein Adapter". Genau so nachgebaut.
        Some(a) if a.seq.len() > 60 => Err(format!(
            "erkannte Folge hat {} Basen; fastp verwirft Folgen über 60 Basen",
            a.seq.len()
        )),
        Some(a) => Ok(a),
        None => Err("kein Adapter erkannt".to_string()),
    }
}

fn known_name(seq: &str) -> Option<&'static str> {
    KNOWN_ADAPTERS
        .iter()
        .find(|(s, _)| *s == seq)
        .map(|(_, n)| n.trim_start_matches('>'))
}

/// Zahl der Bytes, in denen sich `head` und die 8 Bytes von `window` unterscheiden.
fn differing_bytes(head: u64, window: &[u8]) -> usize {
    let w = u64::from_le_bytes(window.try_into().expect("8 Bytes"));
    let x = head ^ w;
    // Je Byte das oberste Bit setzen, wenn das Byte nicht null ist.
    const LOW7: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    let nonzero = ((x & LOW7).wrapping_add(LOW7) | x) & !LOW7;
    nonzero.count_ones() as usize
}

/// `Evaluator::checkKnownAdapters`: Folge, Treffer, geprüfte Reads.
fn check_known_adapters(reads: &[Vec<u8>]) -> Option<(&'static str, u64, u64)> {
    const MAX_CHECK_READS: u64 = 100_000;
    const MAX_CHECK_BASES: u64 = MAX_CHECK_READS * 1000;
    const MAX_HIT: u64 = 1000;
    const MATCH_REQ: usize = 8;
    const ONE_MISMATCH_FOR_EACH: usize = 16;
    let heads: Vec<u64> = KNOWN_ADAPTERS
        .iter()
        .map(|(a, _)| {
            let b = a.as_bytes();
            let mut h = [0u8; 8];
            let n = b.len().min(8);
            h[..n].copy_from_slice(&b[..n]);
            u64::from_le_bytes(h)
        })
        .collect();
    // Ob ein Adapter einen Read trifft, hängt nur von beiden ab. Das wird blockweise auf allen
    // Kernen vorausberechnet; danach laufen die Regeln von fastp (Abbruch, Überspringen) der
    // Reihe nach über die fertigen Ergebnisse. Das Ergebnis ist dasselbe wie nacheinander.
    const NO_HIT: u8 = u8::MAX;
    let hit = |read: &[u8], a: usize| -> u8 {
        let adata = KNOWN_ADAPTERS[a].0.as_bytes();
        let (rlen, alen) = (read.len(), adata.len());
        if alen >= rlen {
            return NO_HIT;
        }
        for pos in 0..rlen.saturating_sub(MATCH_REQ) {
            let cmplen = (rlen - pos).min(alen);
            let allowed = cmplen / ONE_MISMATCH_FOR_EACH;
            // Vorprüfung, ändert kein Ergebnis: Die ersten 8 Basen liegen immer im Vergleich
            // (Adapter ab 19 Basen, pos + 8 < rlen). Haben sie schon mehr Unterschiede als
            // erlaubt, passt der Adapter hier nicht.
            if alen >= 8 && differing_bytes(heads[a], &read[pos..pos + 8]) > allowed {
                continue;
            }
            let mut mismatch = 0;
            let mut matched = true;
            for i in 0..cmplen {
                if adata[i] != read[i + pos] {
                    mismatch += 1;
                    if mismatch > allowed {
                        matched = false;
                        break;
                    }
                }
            }
            if matched {
                return mismatch.min(NO_HIT as usize - 1) as u8;
            }
        }
        NO_HIT
    };
    let n_adapters = KNOWN_ADAPTERS.len();
    let threads = std::thread::available_parallelism().map_or(1, usize::from);
    const CHUNK: usize = 8192;
    let mut counts = vec![0u64; n_adapters];
    let mut mismatches = vec![0u64; n_adapters];
    let (mut checked_reads, mut checked_bases, mut cur_max) = (0u64, 0u64, 0u64);
    let mut table = vec![NO_HIT; CHUNK * n_adapters];
    'chunks: for chunk in reads.chunks(CHUNK) {
        // Ein übersprungener Adapter bleibt übersprungen: Seine Zahl steht still, `cur_max`
        // wächst nur. Daher genügt es, die zu Beginn des Blocks noch aktiven zu berechnen.
        let active: Vec<usize> = (0..n_adapters)
            .filter(|&a| !(cur_max > 20 && counts[a] < cur_max / 10))
            .collect();
        let active = &active;
        let per = chunk.len().div_ceil(threads);
        std::thread::scope(|scope| {
            for (reads_part, out_part) in chunk.chunks(per).zip(table.chunks_mut(per * n_adapters))
            {
                scope.spawn(move || {
                    for (read, row) in reads_part.iter().zip(out_part.chunks_mut(n_adapters)) {
                        for &a in active {
                            row[a] = hit(read, a);
                        }
                    }
                });
            }
        });
        for (i, read) in chunk.iter().enumerate() {
            checked_reads += 1;
            checked_bases += read.len() as u64;
            if checked_reads > MAX_CHECK_READS || checked_bases > MAX_CHECK_BASES {
                break 'chunks;
            }
            if cur_max > MAX_HIT {
                break 'chunks;
            }
            let row = &table[i * n_adapters..(i + 1) * n_adapters];
            for a in 0..n_adapters {
                if cur_max > 20 && counts[a] < cur_max / 10 {
                    continue;
                }
                if row[a] != NO_HIT {
                    counts[a] += 1;
                    cur_max = cur_max.max(counts[a]);
                    mismatches[a] += u64::from(row[a]);
                }
            }
        }
    }
    // Größte Zahl, bei Gleichstand der erste in der Reihenfolge der Liste.
    let mut best: Option<usize> = None;
    for (a, &c) in counts.iter().enumerate() {
        if c > best.map_or(0, |b| counts[b]) {
            best = Some(a);
        }
    }
    let a = best?;
    let max = counts[a];
    if max > checked_reads / 50 || (max > checked_reads / 200 && mismatches[a] < checked_reads) {
        Some((KNOWN_ADAPTERS[a].0, max, checked_reads))
    } else {
        None
    }
}

const KEY_LEN: usize = 10;
/// fastp lässt die letzte Base jedes Reads bei der Suche weg (`shiftTail`, mindestens 1).
const SHIFT_TAIL: usize = 1;

/// Code einer Base wie in fastp: A 0, T 1, C 2, G 3, sonst keiner.
#[inline]
fn base_code(b: u8) -> Option<u32> {
    match b {
        b'A' => Some(0),
        b'T' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// `Evaluator::seq2int`: Schlüssel der 10 Basen ab `pos`, aus dem vorigen fortgeschrieben,
/// wenn es einen gibt.
#[inline]
fn seq2int(seq: &[u8], pos: usize, last: Option<u32>) -> Option<u32> {
    match last {
        Some(prev) => {
            let mask = (1u32 << (KEY_LEN * 2)) - 1;
            base_code(seq[pos + KEY_LEN - 1]).map(|c| ((prev << 2) & mask) + c)
        }
        None => seq[pos..pos + KEY_LEN]
            .iter()
            .try_fold(0u32, |key, &b| base_code(b).map(|c| (key << 2) + c)),
    }
}

fn int2seq(mut val: u32, len: usize) -> String {
    let mut out = vec![b'N'; len];
    for i in 0..len {
        out[len - i - 1] = b"ATCG"[(val & 3) as usize];
        val >>= 2;
    }
    String::from_utf8(out).expect("ASCII")
}

/// Ruft `f(pos, key)` für jede Stelle ab 20 bis zur letzten gewerteten Base, wie fastp.
fn for_each_key(read: &[u8], limit: usize, mut f: impl FnMut(usize, u32)) {
    let mut key = None;
    let last = read.len() as isize - KEY_LEN as isize - SHIFT_TAIL as isize;
    let mut pos = 20usize;
    while (pos as isize) <= last && pos < limit {
        key = seq2int(read, pos, key);
        if let Some(k) = key {
            f(pos, k);
        }
        pos += 1;
    }
}

/// Zweiter Weg in `evalAdapterAndReadNum`: häufigste Teilfolgen, dann `getAdapterWithSeed`.
fn search_by_seed(reads: &[Vec<u8>]) -> Option<Adapter> {
    let size = 1usize << (KEY_LEN * 2);
    let mut counts = vec![0u32; size];
    for read in reads {
        for_each_key(read, usize::MAX, |_, k| counts[k as usize] += 1);
    }
    counts[0] = 0;
    const TOP: usize = 10;
    let mut top = [0usize; TOP];
    let mut total: u64 = 0;
    for k in 0..size {
        let mut atcg = [0usize; 4];
        for i in 0..KEY_LEN {
            atcg[(k >> (i * 2)) & 3] += 1;
        }
        if atcg.iter().any(|&n| n >= KEY_LEN - 4) {
            continue;
        }
        if atcg[2] + atcg[3] >= KEY_LEN - 2 {
            continue;
        }
        if k >> 12 == 0xff {
            continue;
        }
        let val = counts[k];
        total += u64::from(val);
        // Einfügen in die Bestenliste genau wie fastp: absteigend, bei Gleichstand kommt der
        // neue Schlüssel vor die vorhandenen.
        let mut t = TOP - 1;
        loop {
            if val < counts[top[t]] {
                if t < TOP - 1 {
                    let mut m = TOP - 1;
                    while m > t + 1 {
                        top[m] = top[m - 1];
                        m -= 1;
                    }
                    top[t + 1] = k;
                }
                break;
            } else if t == 0 {
                let mut m = TOP - 1;
                while m > 0 {
                    top[m] = top[m - 1];
                    m -= 1;
                }
                top[0] = k;
                break;
            }
            t -= 1;
        }
    }
    const FOLD: u64 = 20;
    let mut candidates = Vec::new();
    for &key in &top {
        if key == 0 {
            continue;
        }
        let count = u64::from(counts[key]);
        if count < 10 || count * (size as u64) < total * FOLD {
            break;
        }
        let seq = int2seq(key as u32, KEY_LEN);
        let s = seq.as_bytes();
        let diff = s.windows(2).filter(|w| w[0] != w[1]).count();
        if diff < 3 {
            continue;
        }
        candidates.push(key as u32);
    }
    // fastp prüft die Keime der Reihe nach und nimmt den ersten, der eine Folge ergibt. Jede
    // Prüfung hängt nur vom Keim ab; sie laufen gleichzeitig, gewählt wird wie bei fastp.
    let found: Vec<Option<Adapter>> = std::thread::scope(|scope| {
        let handles: Vec<_> = candidates
            .iter()
            .map(|&key| scope.spawn(move || adapter_with_seed(key, reads)))
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("Keimsuche"))
            .collect()
    });
    found.into_iter().flatten().next()
}

/// `NucleotideTree`: Präfixbaum der Fortsetzungen mit Zählern.
#[derive(Default)]
struct Tree {
    nodes: Vec<([u32; 8], u32, u8)>, // Kinder (Index 0 = keins), Zähler, Base
}

impl Tree {
    fn new() -> Self {
        Self {
            nodes: vec![([0; 8], 0, b'N')],
        }
    }

    fn add(&mut self, seq: impl Iterator<Item = u8>) {
        let mut cur = 0usize;
        for b in seq {
            if b == b'N' {
                break;
            }
            let slot = (b & 7) as usize;
            let child = self.nodes[cur].0[slot];
            let child = if child == 0 {
                self.nodes.push(([0; 8], 0, b));
                let id = (self.nodes.len() - 1) as u32;
                self.nodes[cur].0[slot] = id;
                id
            } else {
                child
            };
            self.nodes[child as usize].1 += 1;
            cur = child as usize;
        }
    }

    /// `getDominantPath`: folgt dem Kind mit mindestens 95 % der Zähler, solange es mindestens
    /// 50 sind. `reached_leaf` wird falsch, wenn kein Kind überwiegt.
    fn dominant_path(&self, reached_leaf: &mut bool) -> Vec<u8> {
        let mut out = Vec::new();
        let mut cur = 0usize;
        loop {
            let kids = self.nodes[cur].0;
            let total: u32 = kids
                .iter()
                .filter(|&&c| c != 0)
                .map(|&c| self.nodes[c as usize].1)
                .sum();
            if total < 50 {
                break;
            }
            let next = kids
                .iter()
                .filter(|&&c| c != 0)
                .find(|&&c| f64::from(self.nodes[c as usize].1) / f64::from(total) >= 0.95);
            match next {
                Some(&c) => {
                    out.push(self.nodes[c as usize].2);
                    cur = c as usize;
                }
                None => {
                    *reached_leaf = false;
                    break;
                }
            }
        }
        out
    }
}

/// `Evaluator::getAdapterWithSeed`.
fn adapter_with_seed(seed: u32, reads: &[Vec<u8>]) -> Option<Adapter> {
    const MAX_SEARCH: usize = 500;
    let mut forward = Tree::new();
    let mut backward = Tree::new();
    for read in reads {
        let end = read.len().saturating_sub(SHIFT_TAIL);
        for_each_key(read, MAX_SEARCH, |pos, k| {
            if k == seed {
                forward.add(read[pos + KEY_LEN..end].iter().copied());
            }
        });
    }
    for read in reads {
        for_each_key(read, MAX_SEARCH, |pos, k| {
            if k == seed {
                backward.add(read[..pos].iter().rev().copied());
            }
        });
    }
    let mut reached_leaf = true;
    let fwd = forward.dominant_path(&mut reached_leaf);
    let bwd = backward.dominant_path(&mut reached_leaf);
    let seed_seq = int2seq(seed, KEY_LEN);
    let mut adapter: Vec<u8> = bwd.into_iter().rev().collect();
    adapter.extend_from_slice(seed_seq.as_bytes());
    adapter.extend_from_slice(&fwd);
    adapter.truncate(60);
    // `matchKnownAdapter`: ein bekannter Adapter, der Anfang der gefundenen Folge ist.
    if let Some((known, name)) = KNOWN_ADAPTERS
        .iter()
        .find(|(k, _)| adapter.len() >= k.len() && adapter.starts_with(k.as_bytes()))
    {
        return Some(Adapter {
            seq: known.as_bytes().to_vec(),
            name: Some(name.trim_start_matches('>')),
            found: Found::Seed {
                seed: seed_seq,
                known: true,
            },
        });
    }
    reached_leaf.then_some(Adapter {
        seq: adapter,
        name: None,
        found: Found::Seed {
            seed: seed_seq,
            known: false,
        },
    })
}

/// Warum ein Read den Filter nicht besteht; Zahlen wie in fastp (`common.h`), das Paar zählt
/// unter dem größeren Wert der beiden Reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Verdict {
    Pass = 0,
    NBase = 12,
    Length = 16,
    Quality = 20,
    Dimer = 28,
}

impl Verdict {
    pub fn name(self) -> &'static str {
        match self {
            Verdict::Pass => "behalten",
            Verdict::NBase => "zu_viele_n",
            Verdict::Length => "zu_kurz",
            Verdict::Quality => "schlechte_qualitaet",
            Verdict::Dimer => "adapter_dimer",
        }
    }
}

/// Was mit einem Paar geschah.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PairReport {
    /// Basen, die PolyG je Read entfernt hat.
    pub poly_g: [usize; 2],
    /// Über die Überlappung der beiden Reads gekürzt.
    pub by_overlap: bool,
    /// Über die Adapterfolge gekürzt, je Read.
    pub by_sequence: [bool; 2],
    /// Basen, die Adapter je Read entfernt haben.
    pub adapter_bases: [usize; 2],
    pub verdict: [Option<Verdict>; 2],
}

impl PairReport {
    /// Gilt für das Paar: der größere Grund der beiden Reads (wie fastp).
    pub fn pair_verdict(&self) -> Verdict {
        self.verdict[0]
            .unwrap_or(Verdict::Pass)
            .max(self.verdict[1].unwrap_or(Verdict::Pass))
    }

    pub fn kept(&self) -> bool {
        self.pair_verdict() == Verdict::Pass
    }
}

/// Trimmt Paare nach festen Entscheidungen, die vor dem Lauf getroffen wurden.
#[derive(Debug, Clone)]
pub struct Trimmer {
    pub settings: Settings,
    pub poly_g: bool,
    pub adapters: [Option<Vec<u8>>; 2],
}

impl Trimmer {
    /// Kürzt beide Reads an Ort und Stelle und sagt, ob das Paar bleibt.
    pub fn trim_pair(&self, r1: &mut Record, r2: &mut Record) -> PairReport {
        let s = &self.settings;
        let mut report = PairReport::default();
        if self.poly_g {
            report.poly_g[0] = trim_poly_g(r1, s.poly_g_min_len);
            report.poly_g[1] = trim_poly_g(r2, s.poly_g_min_len);
        }
        let before = [r1.seq.len(), r2.seq.len()];
        let ov = overlap(
            &r1.seq,
            &r2.seq,
            s.overlap_diff_limit,
            s.overlap_require,
            f64::from(s.overlap_diff_percent) / 100.0,
        );
        let mut trimmed = [false; 2];
        if let Some((offset, len)) = ov
            && offset < 0
        {
            let len1 = r1.seq.len().min(len);
            let len2 = r2.seq.len().min(len);
            resize(r1, len1);
            resize(r2, len2);
            report.by_overlap = true;
            trimmed = [true, true];
        } else {
            for (i, r) in [&mut *r1, &mut *r2].into_iter().enumerate() {
                if let Some(a) = &self.adapters[i]
                    && trim_by_sequence(r, a, 4)
                {
                    report.by_sequence[i] = true;
                    trimmed[i] = true;
                }
            }
        }
        report.adapter_bases = [before[0] - r1.seq.len(), before[1] - r2.seq.len()];
        let dimer = (trimmed[0] || trimmed[1])
            && r1.seq.len() <= s.dimer_max_len
            && r2.seq.len() <= s.dimer_max_len;
        report.verdict = if dimer {
            [Some(Verdict::Dimer), Some(Verdict::Dimer)]
        } else {
            [Some(self.filter(r1)), Some(self.filter(r2))]
        };
        report
    }

    /// `Filter::passFilter` mit den Vorgaben.
    fn filter(&self, r: &Record) -> Verdict {
        let s = &self.settings;
        let len = r.seq.len();
        if len == 0 {
            return Verdict::Length;
        }
        let threshold = s.qualified_phred + 33;
        let low = r.qual.iter().filter(|&&q| q < threshold).count();
        let n = r.seq.iter().filter(|&&b| b == b'N').count();
        if low as f64 > f64::from(s.unqualified_percent) * len as f64 / 100.0 {
            Verdict::Quality
        } else if n > s.n_base_limit {
            Verdict::NBase
        } else if len < s.length_required {
            Verdict::Length
        } else {
            Verdict::Pass
        }
    }
}

fn resize(r: &mut Record, len: usize) {
    if len <= r.seq.len() {
        r.seq.truncate(len);
        r.qual.truncate(len);
    }
}

/// `PolyX::trimPolyG`; gibt die Zahl der entfernten Basen zurück.
fn trim_poly_g(r: &mut Record, min_len: usize) -> usize {
    const ONE_MISMATCH_FOR_EACH: usize = 8;
    const MAX_MISMATCH: usize = 5;
    let rlen = r.seq.len();
    if rlen == 0 {
        return 0;
    }
    let mut mismatch = 0;
    let mut first_g = rlen - 1;
    let mut i = 0;
    while i < rlen {
        if r.seq[rlen - i - 1] != b'G' {
            mismatch += 1;
        } else {
            first_g = rlen - i - 1;
        }
        let allowed = (i + 1) / ONE_MISMATCH_FOR_EACH;
        if mismatch > MAX_MISMATCH || (mismatch > allowed && i + 1 >= min_len) {
            break;
        }
        i += 1;
    }
    if i >= min_len {
        resize(r, first_g);
        rlen - r.seq.len()
    } else {
        0
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

fn count_mismatches(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

/// Haben `a` und `b` (gleich lang) höchstens `limit` Unterschiede? Vergleicht 8 Bytes auf
/// einmal und hört auf, sobald die Grenze überschritten ist; das Ergebnis ist dasselbe wie
/// `count_mismatches(a, b) <= limit`.
fn within_mismatches(a: &[u8], b: &[u8], limit: usize) -> bool {
    debug_assert_eq!(a.len(), b.len());
    let mut diff = 0;
    let (ca, ra) = a.as_chunks::<8>();
    let (cb, rb) = b.as_chunks::<8>();
    for (x, y) in ca.iter().zip(cb) {
        diff += differing_bytes(u64::from_le_bytes(*x), y);
        if diff > limit {
            return false;
        }
    }
    diff + count_mismatches(ra, rb) <= limit
}

/// `OverlapAnalysis::analyze` ohne Lücke: `(offset, overlap_len)` der ersten passenden
/// Verschiebung.
fn overlap(
    r1: &[u8],
    r2: &[u8],
    diff_limit: usize,
    require: usize,
    percent: f64,
) -> Option<(isize, usize)> {
    const COMPLETE_COMPARE: usize = 50;
    let rc2: Vec<u8> = r2.iter().rev().map(|&b| complement(b)).collect();
    let (len1, len2) = (r1.len() as isize, rc2.len() as isize);
    let accept = |a: &[u8], b: &[u8], len: usize, limit: usize| {
        // Nur die ersten 50 Basen entscheiden; dahinter wird nur gezählt (fastp).
        let head = len.min(COMPLETE_COMPARE);
        within_mismatches(&a[..head], &b[..head], limit)
    };
    let limit_for = |len: usize| diff_limit.min((len as f64 * percent) as usize);
    let require = require as isize;
    let mut offset = 0isize;
    while offset < len1 - require {
        let len = (len1 - offset).min(len2) as usize;
        if accept(&r1[offset as usize..], &rc2, len, limit_for(len)) {
            return Some((offset, len));
        }
        offset += 1;
    }
    offset = 0;
    while offset > -(len2 - require) {
        let len = len1.min(len2 - offset.abs()) as usize;
        if accept(r1, &rc2[(-offset) as usize..], len, limit_for(len)) {
            return Some((offset, len));
        }
        offset -= 1;
    }
    None
}

/// `Matcher::matchWithOneInsertion`: passt `ins[..=cmplen]` mit einer eingefügten Base auf
/// `normal[..cmplen]`, mit höchstens `limit` Unterschieden?
fn match_with_one_insertion(ins: &[u8], normal: &[u8], cmplen: usize, limit: isize) -> bool {
    if cmplen == 0 {
        return false;
    }
    // Ohne Speicher anzulegen oder zu nullen (die Suche ruft das je Read Dutzende Male): Felder je
    // Faden wiederverwenden. Jeder gelesene Eintrag wird vorher gesetzt: `left` nur bis zum
    // Abbruch, dahinter liest die Schlussschleife nicht mehr (sie endet dort mit `false`).
    thread_local! {
        static BUF: std::cell::RefCell<(Vec<isize>, Vec<isize>)> =
            const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
    }
    BUF.with(|b| {
        let (left, right) = &mut *b.borrow_mut();
        if left.len() < cmplen {
            left.resize(cmplen, 0);
            right.resize(cmplen, 0);
        }
        one_insertion_with(
            ins,
            normal,
            cmplen,
            limit,
            &mut left[..cmplen],
            &mut right[..cmplen],
        )
    })
}

fn one_insertion_with(
    ins: &[u8],
    normal: &[u8],
    cmplen: usize,
    limit: isize,
    left: &mut [isize],
    right: &mut [isize],
) -> bool {
    left[0] = isize::from(ins[0] != normal[0]);
    right[cmplen - 1] = isize::from(ins[cmplen] != normal[cmplen - 1]);
    for i in 1..cmplen {
        left[i] = left[i - 1] + isize::from(ins[i] != normal[i]);
        if left[i] + right[cmplen - 1] > limit {
            break;
        }
    }
    let mut i = cmplen as isize - 2;
    while i >= 0 {
        let u = i as usize;
        right[u] = right[u + 1] + isize::from(ins[u + 1] != normal[u]);
        if right[u] + left[0] > limit {
            for p in right.iter_mut().take(u) {
                *p = limit + 1;
            }
            break;
        }
        i -= 1;
    }
    for i in 1..cmplen {
        if left[i - 1] + right[cmplen - 1] > limit {
            return false;
        }
        if left[i - 1] + right[i] <= limit {
            return true;
        }
    }
    false
}

/// Gleich lange Vergleiche der Suche mit Lücke nur einmal rechnen (Schalter zurück:
/// `GC_TRIM_GLEICHE_LAENGE=0`).
fn same_len_skip() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_TRIM_GLEICHE_LAENGE").map_or(true, |v| v != "0"))
}

/// `AdapterTrimmer::trimBySequence`; `true`, wenn gekürzt wurde.
fn trim_by_sequence(r: &mut Record, adapter: &[u8], match_req: usize) -> bool {
    const ONE_MISMATCH_FOR_EACH: isize = 8;
    let rlen = r.seq.len() as isize;
    let alen = adapter.len() as isize;
    if alen < match_req as isize {
        return false;
    }
    let req = match_req as isize;
    let start: isize = if alen >= 16 {
        -4
    } else if alen >= 12 {
        -3
    } else if alen >= 8 {
        -2
    } else {
        0
    };
    let seq = &r.seq;
    let mut found: Option<isize> = None;
    let mut pos = start;
    while pos < rlen - req {
        let cmplen = (rlen - pos).min(alen);
        let allowed = (cmplen / ONE_MISMATCH_FOR_EACH) as usize;
        let skip = (-pos).max(0);
        let n = (cmplen - skip) as usize;
        let a = &adapter[skip as usize..skip as usize + n];
        let b = &seq[(skip + pos) as usize..(skip + pos) as usize + n];
        if within_mismatches(a, b, allowed) {
            found = Some(pos);
            break;
        }
        pos += 1;
    }
    // Eigenheit von fastp 1.4.0: Beide Suchen mit Lücke vergleichen immer ab dem Anfang des
    // Reads, nicht ab `pos`; nur die Vergleichslänge hängt von `pos` ab.
    if found.is_none() {
        // Das Ergebnis hängt nur von `cmplen` ab: gleiche Länge wie eben heißt gleiches Nein.
        let mut pos = 0isize;
        let mut last = -1isize;
        while pos < rlen - req - 1 {
            let cmplen = (rlen - pos - 1).min(alen);
            if cmplen != last || !same_len_skip() {
                last = cmplen;
                let allowed = cmplen / ONE_MISMATCH_FOR_EACH - 1;
                if match_with_one_insertion(seq, adapter, cmplen as usize, allowed) {
                    found = Some(pos);
                    break;
                }
            }
            pos += 1;
        }
    }
    if found.is_none() {
        let mut pos = 0isize;
        let mut last = -1isize;
        while pos < rlen - req {
            let cmplen = (rlen - pos).min(alen - 1);
            if cmplen != last || !same_len_skip() {
                last = cmplen;
                let allowed = cmplen / ONE_MISMATCH_FOR_EACH - 1;
                if match_with_one_insertion(adapter, seq, cmplen as usize, allowed) {
                    found = Some(pos);
                    break;
                }
            }
            pos += 1;
        }
    }
    match found {
        Some(p) if p < 0 => {
            resize(r, 0);
            true
        }
        Some(p) => {
            resize(r, p as usize);
            true
        }
        None => false,
    }
}

/// Zähler über alle Paare eines Laufs; die Summen gehen auf (siehe [`Counts::balanced`]).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Counts {
    pub pairs: u64,
    pub kept: u64,
    pub dropped_quality: u64,
    pub dropped_length: u64,
    pub dropped_n: u64,
    pub dropped_dimer: u64,
    pub poly_g_reads: u64,
    pub poly_g_bases: u64,
    pub overlap_pairs: u64,
    pub sequence_reads: [u64; 2],
    pub adapter_bases: u64,
}

impl Counts {
    pub fn add(&mut self, r: &PairReport) {
        self.pairs += 1;
        match r.pair_verdict() {
            Verdict::Pass => self.kept += 1,
            Verdict::Quality => self.dropped_quality += 1,
            Verdict::Length => self.dropped_length += 1,
            Verdict::NBase => self.dropped_n += 1,
            Verdict::Dimer => self.dropped_dimer += 1,
        }
        for i in 0..2 {
            if r.poly_g[i] > 0 {
                self.poly_g_reads += 1;
                self.poly_g_bases += r.poly_g[i] as u64;
            }
            if r.by_sequence[i] {
                self.sequence_reads[i] += 1;
            }
            self.adapter_bases += r.adapter_bases[i] as u64;
        }
        if r.by_overlap {
            self.overlap_pairs += 1;
        }
    }

    pub fn merge(&mut self, o: &Counts) {
        self.pairs += o.pairs;
        self.kept += o.kept;
        self.dropped_quality += o.dropped_quality;
        self.dropped_length += o.dropped_length;
        self.dropped_n += o.dropped_n;
        self.dropped_dimer += o.dropped_dimer;
        self.poly_g_reads += o.poly_g_reads;
        self.poly_g_bases += o.poly_g_bases;
        self.overlap_pairs += o.overlap_pairs;
        self.sequence_reads[0] += o.sequence_reads[0];
        self.sequence_reads[1] += o.sequence_reads[1];
        self.adapter_bases += o.adapter_bases;
    }

    /// Jedes Paar ist genau einmal behalten oder aus einem Grund verworfen.
    pub fn balanced(&self) -> bool {
        self.pairs
            == self.kept
                + self.dropped_quality
                + self.dropped_length
                + self.dropped_n
                + self.dropped_dimer
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Der frühere Ablauf Read für Read, als Maßstab für die Fassung auf mehreren Kernen.
    fn check_known_adapters_sequential(reads: &[Vec<u8>]) -> Option<(&'static str, u64, u64)> {
        const MAX_CHECK_READS: u64 = 100_000;
        const MAX_CHECK_BASES: u64 = MAX_CHECK_READS * 1000;
        const MAX_HIT: u64 = 1000;
        const MATCH_REQ: usize = 8;
        const ONE_MISMATCH_FOR_EACH: usize = 16;
        let heads: Vec<u64> = KNOWN_ADAPTERS
            .iter()
            .map(|(a, _)| {
                let b = a.as_bytes();
                let mut h = [0u8; 8];
                let n = b.len().min(8);
                h[..n].copy_from_slice(&b[..n]);
                u64::from_le_bytes(h)
            })
            .collect();
        let mut counts = vec![0u64; KNOWN_ADAPTERS.len()];
        let mut mismatches = vec![0u64; KNOWN_ADAPTERS.len()];
        let (mut checked_reads, mut checked_bases, mut cur_max) = (0u64, 0u64, 0u64);
        for read in reads {
            let rlen = read.len();
            checked_reads += 1;
            checked_bases += rlen as u64;
            if checked_reads > MAX_CHECK_READS || checked_bases > MAX_CHECK_BASES {
                break;
            }
            if cur_max > MAX_HIT {
                break;
            }
            for (a, (adapter, _)) in KNOWN_ADAPTERS.iter().enumerate() {
                let adata = adapter.as_bytes();
                let alen = adata.len();
                if alen >= rlen {
                    continue;
                }
                if cur_max > 20 && counts[a] < cur_max / 10 {
                    continue;
                }
                for pos in 0..rlen.saturating_sub(MATCH_REQ) {
                    let cmplen = (rlen - pos).min(alen);
                    let allowed = cmplen / ONE_MISMATCH_FOR_EACH;
                    // Vorprüfung, ändert kein Ergebnis: Die ersten 8 Basen liegen immer im
                    // Vergleich (Adapter ab 19 Basen, pos + 8 < rlen). Haben sie schon mehr
                    // Unterschiede als erlaubt, passt der Adapter hier nicht.
                    if alen >= 8 && differing_bytes(heads[a], &read[pos..pos + 8]) > allowed {
                        continue;
                    }
                    let mut mismatch = 0;
                    let mut matched = true;
                    for i in 0..cmplen {
                        if adata[i] != read[i + pos] {
                            mismatch += 1;
                            if mismatch > allowed {
                                matched = false;
                                break;
                            }
                        }
                    }
                    if matched {
                        counts[a] += 1;
                        cur_max = cur_max.max(counts[a]);
                        mismatches[a] += mismatch as u64;
                        break;
                    }
                }
            }
        }
        // Größte Zahl, bei Gleichstand der erste in der Reihenfolge der Liste.
        let mut best: Option<usize> = None;
        for (a, &c) in counts.iter().enumerate() {
            if c > best.map_or(0, |b| counts[b]) {
                best = Some(a);
            }
        }
        let a = best?;
        let max = counts[a];
        if max > checked_reads / 50 || (max > checked_reads / 200 && mismatches[a] < checked_reads)
        {
            Some((KNOWN_ADAPTERS[a].0, max, checked_reads))
        } else {
            None
        }
    }

    #[test]
    fn known_adapters_parallel_equals_sequential() {
        let mut x: u64 = 0x1234_5678_9abc_def1;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        // Ohne Adapter, mit wenig Adapter, mit viel Adapter (früher Abbruch).
        for share in [0u64, 3, 40] {
            let adapter = KNOWN_ADAPTERS[17].0.as_bytes();
            let reads: Vec<Vec<u8>> = (0..30_000)
                .map(|_| {
                    let mut r: Vec<u8> = (0..150).map(|_| b"ACGT"[(next() & 3) as usize]).collect();
                    if next() % 100 < share {
                        let at = 40 + (next() % 100) as usize;
                        for (i, &b) in adapter.iter().enumerate() {
                            if at + i < r.len() {
                                r[at + i] = b;
                            }
                        }
                    }
                    r
                })
                .collect();
            assert_eq!(
                check_known_adapters(&reads),
                check_known_adapters_sequential(&reads)
            );
        }
    }

    #[test]
    fn within_mismatches_matches_count() {
        let mut x: u64 = 0x2545_f491_4f6c_dd1d;
        for _ in 0..5_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let len = (x % 70) as usize;
            let a: Vec<u8> = (0..len)
                .map(|i| b"ACGT"[((x >> (i % 60)) & 3) as usize])
                .collect();
            let mut b = a.clone();
            for (i, c) in b.iter_mut().enumerate() {
                if (x.rotate_left(i as u32 * 7) & 15) == 0 {
                    *c = b'N';
                }
            }
            let n = count_mismatches(&a, &b);
            for limit in 0..8 {
                assert_eq!(within_mismatches(&a, &b, limit), n <= limit);
            }
        }
    }

    #[test]
    fn differing_bytes_counts_like_plain_loop() {
        let mut x: u64 = 0x9e37_79b9_7f4a_7c15;
        for _ in 0..20_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let a = x.to_le_bytes();
            let mut b = a;
            for (i, byte) in b.iter_mut().enumerate() {
                if (x >> (i * 3)) & 3 == 0 {
                    *byte = byte.wrapping_add(((x >> 40) as u8) | 1);
                }
                if (x >> (i * 5 + 1)) & 7 == 0 {
                    *byte ^= 0x80;
                }
            }
            let plain = a.iter().zip(b.iter()).filter(|(p, q)| p != q).count();
            assert_eq!(differing_bytes(u64::from_le_bytes(a), &b), plain);
        }
    }

    fn rec(seq: &str, qual: u8) -> Record {
        Record {
            name: b"r".to_vec(),
            comment: Vec::new(),
            seq: seq.as_bytes().to_vec(),
            qual: vec![qual; seq.len()],
        }
    }

    fn rc(s: &str) -> String {
        s.bytes().rev().map(|b| complement(b) as char).collect()
    }

    const TRUSEQ1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    const TRUSEQ2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

    #[test]
    fn known_list_is_sorted_like_std_map() {
        assert!(KNOWN_ADAPTERS.windows(2).all(|w| w[0].0 < w[1].0));
        assert_eq!(KNOWN_ADAPTERS.len(), 234);
    }

    #[test]
    fn two_color_by_first_name() {
        assert_eq!(two_color(b"A00123:8:H7:1:1101:1000:1000"), Some("A"));
        assert_eq!(two_color(b"NB501234:1:X:1:1:1:1"), Some("NB"));
        assert_eq!(two_color(b"V350012345L1C001R0010000001"), None);
        assert_eq!(two_color(b"SRR123.1"), None);
    }

    #[test]
    fn poly_g_tail_is_cut() {
        let mut r = rec("ACGTACGTACGTACGTACGTGGGGGGGGGGGGGGG", b'I');
        // Wie fastp: Die Suche läuft über das T hinaus bis zum G von "ACGT" (zwei Fehler auf 18
        // Basen sind erlaubt) und schneidet dort, also zwei Basen mehr als das reine PolyG.
        assert_eq!(trim_poly_g(&mut r, 10), 17);
        assert_eq!(r.seq, b"ACGTACGTACGTACGTAC");
        // Zu kurz: bleibt.
        let mut r = rec("ACGTACGTACGTACGTACGTGGGGG", b'I');
        assert_eq!(trim_poly_g(&mut r, 10), 0);
        // Kein G am Ende, aber die Suche erreicht 10 Basen nicht: bleibt.
        let mut r = rec("ACGTACGTACGTACGTACGTACGT", b'I');
        assert_eq!(trim_poly_g(&mut r, 10), 0);
    }

    #[test]
    fn short_insert_is_cut_by_overlap() {
        // Insert 60 Basen, Reads 100 Basen: Rest ist Adapter.
        let insert = "TTGACCATGCAGTTCAGGATCCATGGACTTAGCCATGACTGCAATGGTCCAAGTTCCAGA";
        let s1 = format!("{insert}{TRUSEQ1}GGGGGGG");
        let s2 = format!("{}{TRUSEQ2}CCCCCCC", rc(insert));
        let (mut r1, mut r2) = (rec(&s1[..100], b'I'), rec(&s2[..100], b'I'));
        let t = Trimmer {
            settings: Settings::default(),
            poly_g: false,
            adapters: [None, None],
        };
        let rep = t.trim_pair(&mut r1, &mut r2);
        assert!(rep.by_overlap);
        assert_eq!(r1.seq, insert.as_bytes());
        assert_eq!(r2.seq, rc(insert).as_bytes());
        assert!(rep.kept());
    }

    #[test]
    fn adapter_by_sequence_and_filters() {
        let mut r = rec(&format!("ACGTTGCAAGCTTGCAACGTTGCATGCA{TRUSEQ1}"), b'I');
        assert!(trim_by_sequence(&mut r, TRUSEQ1.as_bytes(), 4));
        assert_eq!(r.seq, b"ACGTTGCAAGCTTGCAACGTTGCATGCA");
        let t = Trimmer {
            settings: Settings::default(),
            poly_g: false,
            adapters: [None, None],
        };
        assert_eq!(t.filter(&rec("ACGTACGTACGTAC", b'I')), Verdict::Length);
        assert_eq!(
            t.filter(&rec("ACGTACGTACGTACGTACGT", b'#')),
            Verdict::Quality
        );
        assert_eq!(
            t.filter(&rec("ACGTACNNNNNNGTACGTACGT", b'I')),
            Verdict::NBase
        );
        assert_eq!(t.filter(&rec("ACGTACGTACGTACGTACGT", b'I')), Verdict::Pass);
    }

    #[test]
    fn unknown_adapter_is_found_by_seed() {
        let mut x: u64 = 0x0bad_cafe_dead_beef;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        // Eine Folge, die in keiner Liste steht, als Adapter-Dimer am Anfang von 30 % der Reads,
        // danach immer dieselbe Fortsetzung. fastp verlangt, dass die Verlängerung in beide
        // Richtungen eindeutig bis zum Ende läuft; das ist bei Dimeren der Fall.
        let adapter = b"GTCAGTTCAAGCTGACATCGTAGCATCGAGTCTTGCAGACCTAGTCAGCA";
        assert!(
            KNOWN_ADAPTERS
                .iter()
                .all(|(k, _)| !k.as_bytes().starts_with(&adapter[..12]))
        );
        let reads: Vec<Vec<u8>> = (0..20_000)
            .map(|_| {
                let mut r: Vec<u8> = (0..150).map(|_| b"ACGT"[(next() & 3) as usize]).collect();
                if next() % 10 < 3 {
                    r[..adapter.len()].copy_from_slice(adapter);
                    for (i, b) in r[adapter.len()..].iter_mut().enumerate() {
                        *b = b"TGCA"[(i * 7 + i / 3) % 4];
                    }
                }
                r
            })
            .collect();
        let found = detect_adapter(&reads).expect("Folge erkannt");
        assert!(
            matches!(found.found, Found::Seed { .. }),
            "{:?}",
            found.found
        );
        let seq = String::from_utf8(found.seq).expect("ASCII");
        assert!(seq.as_bytes().starts_with(adapter), "gefunden {seq}");
    }

    #[test]
    fn known_adapter_is_detected() {
        // Zufällige Inserts, ein Teil kürzer als der Read mit TruSeq dahinter.
        let mut x = 7u64;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let reads: Vec<Vec<u8>> = (0..20_000)
            .map(|i| {
                let ins = if i % 10 == 0 { 60 } else { 150 };
                let mut r: Vec<u8> = (0..ins).map(|_| b"ACGT"[(next() & 3) as usize]).collect();
                r.extend_from_slice(TRUSEQ1.as_bytes());
                r.extend(std::iter::repeat_n(b'A', 100));
                r.truncate(150);
                r
            })
            .collect();
        let a = detect_adapter(&reads).expect("Adapter");
        assert_eq!(a.seq, TRUSEQ1.as_bytes());
        assert!(matches!(a.found, Found::Known { .. }));
        assert!(detect_adapter(&reads[..9_999]).is_err());
    }

    #[test]
    fn insertion_matcher_finds_one_extra_base() {
        let normal = b"ACGTACGTTGCAACGT";
        let ins = b"ACGTACGTXTGCAACGT";
        assert!(match_with_one_insertion(ins, normal, 16, 0));
        assert!(!match_with_one_insertion(
            b"TTTTTTTTTTTTTTTTT",
            normal,
            16,
            1
        ));
    }

    #[test]
    fn counts_balance() {
        let mut c = Counts::default();
        let mut r = PairReport {
            verdict: [Some(Verdict::Pass), Some(Verdict::Length)],
            ..Default::default()
        };
        c.add(&r);
        r.verdict = [Some(Verdict::Pass), Some(Verdict::Pass)];
        c.add(&r);
        assert!(c.balanced());
        assert_eq!((c.kept, c.dropped_length), (1, 1));
    }
}
