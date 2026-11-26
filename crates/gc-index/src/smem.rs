//! Seed-Suche: Suche in beide Richtungen auf dem FM-Index und längste exakte Treffer (SMEM).
//!
//! Der Indextext ist der Vorwärtsstrang gefolgt von seinem reversen Komplement. Er ist damit
//! sein eigenes reverses Komplement, und jedes Vorkommen eines Musters `P` an Textposition `i`
//! hat ein Vorkommen von `rk(P)` an Position `2L - i - |P|`. Deshalb reicht eine einzige BWT, um
//! ein Muster nach links und nach rechts zu verlängern (FMD-Index, Li 2012): Man führt die
//! Zeilenbereiche von `P` und von `rk(P)` gemeinsam.

use crate::fm::FmIndex;

/// Zeilenbereiche eines Musters `P` und seines reversen Komplements, beide gleich groß.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BiInterval {
    /// Erste Zeile der Suffixe, die mit `P` beginnen.
    pub k: u64,
    /// Erste Zeile der Suffixe, die mit `rk(P)` beginnen.
    pub kr: u64,
    /// Anzahl der Zeilen, also Vorkommen von `P` im Indextext (beide Stränge).
    pub size: u64,
}

/// Ein längster exakter Treffer eines Reads: `read[start..end]` kommt `interval.size`-mal vor
/// und lässt sich in keine Richtung verlängern, ohne unter die verlangte Häufigkeit zu fallen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Smem {
    pub start: u32,
    pub end: u32,
    pub interval: BiInterval,
}

impl Smem {
    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }
}

/// Längstes Muster, dessen Bereich in der Tabelle steht.
const SHORT_MAX: usize = 10;

/// Die Bereiche aller Muster bis zu einer festen Länge, fertig berechnet.
///
/// Die ersten Schritte jeder Suche laufen auf sehr breiten Bereichen: Die beiden Enden liegen in
/// verschiedenen Zählblöcken weit auseinander, jeder Schritt kostet zwei volle Zählungen und zwei
/// Zugriffe auf fernen Speicher. Für kurze Muster gibt es aber nur wenige Möglichkeiten (4^n bei
/// Länge n). Ihre Bereiche werden einmal beim Anlegen des Index mit genau denselben Schritten
/// berechnet, die sonst die Suche ginge, und danach nur noch nachgeschlagen. Das Ergebnis der
/// Suche ändert sich dadurch nicht.
///
/// Ablage: erst die 4 Muster der Länge 1, dann die 16 der Länge 2 und so weiter. Der Code eines
/// Musters ist seine Basenfolge als Zahl zur Basis 4, erste Base an höchster Stelle.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(crate) struct ShortTable {
    /// Längste enthaltene Musterlänge, 0 = keine Tabelle.
    longest: usize,
    entries: Vec<BiInterval>,
}

impl ShortTable {
    /// Längste enthaltene Musterlänge, 0 = keine Tabelle.
    #[inline]
    pub(crate) fn longest(&self) -> usize {
        self.longest
    }

    /// Stelle des ersten Musters der Länge `len` (ab 1) in `entries`: 4 + 16 + … davor.
    #[inline]
    fn offset(len: usize) -> usize {
        ((1usize << (2 * len)) - 4) / 3
    }

    /// Längste Musterlänge für einen Index mit `rows` Zeilen: höchstens `SHORT_MAX`, und nicht
    /// mehr Muster als Zeilen, damit die Tabelle bei kleinen Referenzen klein bleibt.
    fn longest_for(rows: u64) -> usize {
        (1..=SHORT_MAX)
            .take_while(|&n| 1u64 << (2 * n) <= rows)
            .last()
            .unwrap_or(0)
    }

    pub(crate) fn build(index: &FmIndex) -> Self {
        let longest = Self::longest_for(index.rows());
        if longest == 0 {
            return Self::default();
        }
        let mut entries = Vec::with_capacity(Self::offset(longest + 1));
        for c in 0..4u8 {
            entries.push(index.bi_start(c));
        }
        for len in 1..longest {
            for code in 0..1usize << (2 * len) {
                let iv = entries[Self::offset(len) + code];
                for c in 0..4u8 {
                    entries.push(index.extend_right(iv, c));
                }
            }
        }
        Self { longest, entries }
    }

    /// Bereich von `pattern` (Codes 0..3), wenn es kurz genug für die Tabelle ist.
    #[inline]
    pub(crate) fn get(&self, pattern: &[u8]) -> Option<BiInterval> {
        let len = pattern.len();
        if len == 0 || len > self.longest {
            return None;
        }
        let code = pattern
            .iter()
            .fold(0usize, |code, &c| code << 2 | usize::from(c));
        Some(self.entries[Self::offset(len) + code])
    }
}

impl FmIndex {
    /// `extend_right` für das Muster `read[start..end]`, dessen Bereich `iv` ist, um die Base
    /// `read[end]`: nachgeschlagen, wenn das neue Muster in der Tabelle steht.
    #[inline]
    fn grow_right(&self, read: &[u8], start: usize, end: usize, iv: BiInterval) -> BiInterval {
        match self.short.get(&read[start..=end]) {
            Some(found) => found,
            None => self.extend_right(iv, read[end]),
        }
    }

    /// `extend_left` für das Muster `read[start..end]`, dessen Bereich `iv` ist, um die Base
    /// `read[start - 1]`: nachgeschlagen, wenn das neue Muster in der Tabelle steht.
    #[inline]
    fn grow_left(&self, read: &[u8], start: usize, end: usize, iv: BiInterval) -> BiInterval {
        match self.short.get(&read[start - 1..end]) {
            Some(found) => found,
            None => self.extend_left(iv, read[start - 1]),
        }
    }

    /// Bereich eines einzelnen Zeichens (Code 0..3).
    pub fn bi_start(&self, c: u8) -> BiInterval {
        let end = |c: u8| {
            if c == 3 {
                self.rows()
            } else {
                self.first_row(c + 1)
            }
        };
        BiInterval {
            k: self.first_row(c),
            kr: self.first_row(3 - c),
            size: end(c) - self.first_row(c),
        }
    }

    /// Verlängert das Muster nach links um jede der vier Basen: Ergebnis `[cP für c in ACGT]`.
    pub fn extend_left_all(&self, iv: BiInterval) -> [BiInterval; 4] {
        let mut out = [BiInterval {
            k: 0,
            kr: 0,
            size: 0,
        }; 4];
        let mut total = 0;
        let (lo, hi) = self.occ4_pair(iv.k, iv.k + iv.size);
        for c in 0..4usize {
            out[c].k = self.first_row(c as u8) + lo[c];
            out[c].size = hi[c] - lo[c];
            total += hi[c] - lo[c];
        }
        // Im Bereich von rk(P) stehen die Fortsetzungen in der Reihenfolge Endzeichen, A, C, G, T.
        // rk(P)·x kommt so oft vor wie komp(x)·P. Das Endzeichen folgt höchstens einmal.
        let mut kr = iv.kr + (iv.size - total);
        for c in (0..4usize).rev() {
            out[c].kr = kr;
            kr += out[c].size;
        }
        out
    }

    /// Verlängert das Muster nach links um die Base `c`.
    pub fn extend_left(&self, iv: BiInterval, c: u8) -> BiInterval {
        self.extend_left_all(iv)[c as usize]
    }

    /// Verlängert das Muster nach rechts um die Base `c`: Linksverlängerung von `rk(P)` um
    /// das Komplement.
    pub fn extend_right(&self, iv: BiInterval, c: u8) -> BiInterval {
        let swapped = BiInterval {
            k: iv.kr,
            kr: iv.k,
            size: iv.size,
        };
        let r = self.extend_left(swapped, 3 - c);
        BiInterval {
            k: r.kr,
            kr: r.k,
            size: r.size,
        }
    }

    /// Alle längsten Treffer, die die Read-Position `x` überdecken und mindestens `min_size`-mal
    /// vorkommen. Hängt sie an `out` an, längster Anfang zuletzt, und gibt die Position zurück,
    /// an der der am weitesten nach rechts reichende Treffer endet. `read` enthält Codes 0..3,
    /// jeder andere Wert gilt als N und beendet jeden Treffer.
    ///
    /// Ablauf wie in BWA-MEM (Li 2013): erst von `x` nach rechts verlängern und jeden Bereich
    /// merken, bevor er kleiner wird; dann alle gemerkten Bereiche gemeinsam nach links
    /// verlängern, die weitest reichenden zuerst. Dieselben Schritte wie [`SmemSearch`], nur am
    /// Stück.
    pub fn smems_at(&self, read: &[u8], x: usize, min_size: u64, out: &mut Vec<Smem>) -> usize {
        let mut search = SmemSearch::new(self, read, x, min_size, out);
        loop {
            search.prefetch(self, read);
            if let Some(reach) = search.step(self, read, out) {
                return reach;
            }
        }
    }

    /// Dritte Seed-Runde von BWA-MEM (`bwt_seed_strategy1`): ab `x` nach rechts verlängern, bis
    /// das Muster mehr als `min_len` Basen hat und seltener als `max_size`-mal vorkommt. Gibt die
    /// Position zurück, an der die nächste Suche beginnt, und den Treffer, falls es einen gibt.
    /// Anders als ein SMEM ist dieser Treffer nicht nach links verlängert und kann in einem
    /// längeren Treffer liegen; er findet Orte, die ein längerer, anderswo passender Treffer
    /// verdeckt. Dieselben Schritte wie [`ForwardSeed`], nur am Stück.
    pub fn seed_forward(
        &self,
        read: &[u8],
        x: usize,
        min_len: u32,
        max_size: u64,
    ) -> (usize, Option<Smem>) {
        let mut search = ForwardSeed::new(self, read, x, min_len, max_size);
        loop {
            search.prefetch(self, read);
            if let Some(done) = search.step(self, read) {
                return done;
            }
        }
    }

    /// Alle längsten exakten Treffer (SMEM) eines Reads mit mindestens `min_len` Basen,
    /// nach Anfang geordnet.
    pub fn smems(&self, read: &[u8], min_len: u32) -> Vec<Smem> {
        let mut out = Vec::new();
        let mut x = 0;
        while x < read.len() {
            x = self.smems_at(read, x, 1, &mut out);
        }
        out.retain(|m| m.len() >= min_len);
        out
    }
}

/// Leerer Bereich.
const NONE: BiInterval = BiInterval {
    k: 0,
    kr: 0,
    size: 0,
};

/// Eine Suche wie [`FmIndex::smems_at`], in Schritte zerlegt. Jeder Schritt braucht höchstens
/// einen Zugriff auf fernen Speicher je Bereich; [`SmemSearch::prefetch`] fordert genau diese
/// Speicherzeilen vorher an. Wer viele Suchen abwechselnd Schritt für Schritt rechnet (erst alle
/// anfordern, dann alle rechnen), lässt den Prozessor auf viele Zugriffe zugleich warten statt
/// auf einen nach dem anderen. Das Ergebnis jeder Suche hängt davon nicht ab.
#[derive(Debug, Clone)]
pub struct SmemSearch {
    x: usize,
    min_size: u64,
    phase: Phase,
    /// Nach rechts: Ende und Bereich des Musters `read[x..end]`.
    end: usize,
    iv: BiInterval,
    /// Nach rechts gemerkte Bereiche, danach die noch lebenden beim Verlängern nach links.
    alive: Vec<(usize, BiInterval)>,
    next_alive: Vec<(usize, BiInterval)>,
    /// Nach links: nächste Read-Position.
    i: isize,
    reach: usize,
    /// Erster Eintrag dieser Suche in `out`.
    first_new: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    Right,
    Left,
    Done(usize),
}

impl SmemSearch {
    /// Beginnt die Suche an Read-Position `x`; `out` ist die Liste, an die sie später anhängt.
    pub fn new(index: &FmIndex, read: &[u8], x: usize, min_size: u64, out: &[Smem]) -> Self {
        let min_size = min_size.max(1);
        let startable = x < read.len() && read[x] <= 3;
        SmemSearch {
            x,
            min_size,
            phase: if startable {
                Phase::Right
            } else {
                Phase::Done(x + 1)
            },
            end: x + 1,
            iv: if startable {
                index.bi_start(read[x])
            } else {
                NONE
            },
            alive: Vec::new(),
            next_alive: Vec::new(),
            i: x as isize - 1,
            reach: x + 1,
            first_new: out.len(),
        }
    }

    /// Fertig? Dann die Position, an der der am weitesten reichende Treffer endet.
    pub fn finished(&self) -> Option<usize> {
        match self.phase {
            Phase::Done(reach) => Some(reach),
            _ => None,
        }
    }

    /// Fordert die Speicherzeilen an, die der nächste Schritt braucht. Ändert nichts.
    pub fn prefetch(&self, index: &FmIndex, read: &[u8]) {
        match self.phase {
            Phase::Right => {
                if self.end < read.len()
                    && read[self.end] <= 3
                    && index.short.get(&read[self.x..=self.end]).is_none()
                {
                    index.prefetch_row(self.iv.kr);
                    index.prefetch_row(self.iv.kr + self.iv.size);
                }
            }
            Phase::Left => {
                if self.i >= 0 && read[self.i as usize] <= 3 {
                    for &(e, iv) in &self.alive {
                        if index.short.get(&read[self.i as usize..e]).is_none() {
                            index.prefetch_row(iv.k);
                            index.prefetch_row(iv.k + iv.size);
                        }
                    }
                }
            }
            Phase::Done(_) => {}
        }
    }

    /// Rechnet bis einschließlich des nächsten Zugriffs auf fernen Speicher. Liefert am Ende die
    /// Position, an der der am weitesten reichende Treffer endet; die Treffer stehen dann in `out`.
    pub fn step(&mut self, index: &FmIndex, read: &[u8], out: &mut Vec<Smem>) -> Option<usize> {
        match self.phase {
            Phase::Right => {
                self.step_right(index, read);
                None
            }
            Phase::Left => {
                self.step_left(index, read, out);
                self.finished()
            }
            Phase::Done(reach) => Some(reach),
        }
    }

    fn step_right(&mut self, index: &FmIndex, read: &[u8]) {
        let x = self.x;
        loop {
            if !(self.end < read.len() && read[self.end] <= 3) {
                return self.finish_right();
            }
            let far = index.short.get(&read[x..=self.end]).is_none();
            let next = index.grow_right(read, x, self.end, self.iv);
            if next.size != self.iv.size {
                self.alive.push((self.end, self.iv));
                if next.size < self.min_size {
                    return self.finish_right();
                }
            }
            self.iv = next;
            self.end += 1;
            if far {
                return;
            }
        }
    }

    fn finish_right(&mut self) {
        let end = self.end;
        if self.alive.last().is_none_or(|&(e, _)| e != end) && self.iv.size >= self.min_size {
            self.alive.push((end, self.iv));
        }
        if self.alive.is_empty() {
            self.phase = Phase::Done(self.x + 1);
            return;
        }
        self.reach = self.alive.last().map_or(self.x + 1, |&(e, _)| e);
        // Nach links, alle gemeinsam; nach fallendem Ende geordnet.
        self.alive.reverse();
        self.phase = Phase::Left;
    }

    fn step_left(&mut self, index: &FmIndex, read: &[u8], out: &mut Vec<Smem>) {
        let i = self.i;
        let c = if i >= 0 { read[i as usize] } else { 4 };
        self.next_alive.clear();
        for &(e, iv) in &self.alive {
            let ext = if c <= 3 {
                index.grow_left(read, (i + 1) as usize, e, iv)
            } else {
                NONE
            };
            if ext.size < self.min_size {
                // Nicht weiter verlängerbar. Ein Treffer ist es nur, wenn kein weiter nach rechts
                // reichender Bereich diesen Schritt überlebt hat und nicht schon einer mit
                // demselben Anfang gemeldet wurde.
                let start = (i + 1) as usize;
                let contained = out[self.first_new..]
                    .last()
                    .is_some_and(|m| m.start as usize <= start);
                if self.next_alive.is_empty() && !contained {
                    out.push(Smem {
                        start: start as u32,
                        end: e as u32,
                        interval: iv,
                    });
                }
            } else if self
                .next_alive
                .last()
                .is_none_or(|&(_, l)| l.size != ext.size)
            {
                self.next_alive.push((e, ext));
            }
        }
        if self.next_alive.is_empty() {
            out[self.first_new..].reverse();
            self.phase = Phase::Done(self.reach);
            return;
        }
        std::mem::swap(&mut self.alive, &mut self.next_alive);
        self.i -= 1;
    }
}

/// [`FmIndex::seed_forward`] in Schritte zerlegt, wie [`SmemSearch`].
#[derive(Debug, Clone)]
pub struct ForwardSeed {
    x: usize,
    min_len: u32,
    max_size: u64,
    iv: BiInterval,
    /// Nächste Read-Position.
    i: usize,
    done: Option<(usize, Option<Smem>)>,
}

impl ForwardSeed {
    pub fn new(index: &FmIndex, read: &[u8], x: usize, min_len: u32, max_size: u64) -> Self {
        let startable = x < read.len() && read[x] <= 3;
        ForwardSeed {
            x,
            min_len,
            max_size,
            iv: if startable {
                index.bi_start(read[x])
            } else {
                NONE
            },
            i: x + 1,
            done: (!startable).then_some((x + 1, None)),
        }
    }

    /// Fordert die Speicherzeilen an, die der nächste Schritt braucht. Ändert nichts.
    pub fn prefetch(&self, index: &FmIndex, read: &[u8]) {
        if self.done.is_none()
            && self.i < read.len()
            && read[self.i] <= 3
            && index.short.get(&read[self.x..=self.i]).is_none()
        {
            index.prefetch_row(self.iv.kr);
            index.prefetch_row(self.iv.kr + self.iv.size);
        }
    }

    /// Rechnet bis einschließlich des nächsten Zugriffs auf fernen Speicher. Liefert am Ende die
    /// nächste Startposition und den Treffer, falls es einen gibt.
    pub fn step(&mut self, index: &FmIndex, read: &[u8]) -> Option<(usize, Option<Smem>)> {
        if let Some(done) = self.done {
            return Some(done);
        }
        let x = self.x;
        while self.i < read.len() {
            let i = self.i;
            let c = read[i];
            if c > 3 {
                self.done = Some((i + 1, None));
                return self.done;
            }
            let far = index.short.get(&read[x..=i]).is_none();
            self.iv = index.grow_right(read, x, i, self.iv);
            if self.iv.size < self.max_size && i - x >= self.min_len as usize {
                let found = (self.iv.size > 0).then_some(Smem {
                    start: x as u32,
                    end: (i + 1) as u32,
                    interval: self.iv,
                });
                self.done = Some((i + 1, found));
                return self.done;
            }
            self.i += 1;
            if far {
                return None;
            }
        }
        self.done = Some((read.len(), None));
        self.done
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fm::{Rng, Variant};
    use crate::reference::Reference;

    /// Indextext ohne Endzeichen als Codes 0..3.
    fn text_of(index: &FmIndex) -> Vec<u8> {
        let t = index.reference().index_text();
        t[..t.len() - 1].iter().map(|c| c - 1).collect()
    }

    fn occurrences(text: &[u8], p: &[u8]) -> u64 {
        if p.is_empty() || p.len() > text.len() {
            return 0;
        }
        text.windows(p.len()).filter(|w| *w == p).count() as u64
    }

    fn genome(seed: u64, len: usize) -> FmIndex {
        let mut rng = Rng(seed);
        let mut s: Vec<u8> = (0..len)
            .map(|_| b"ACGT"[(rng.next() & 3) as usize])
            .collect();
        // Wiederholungen in beiden Richtungen, damit Bereiche mit mehreren Zeilen vorkommen.
        let copy = s[10..90].to_vec();
        s[150..230].copy_from_slice(&copy);
        let rc: Vec<u8> = s[300..360]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        s[400..460].copy_from_slice(&rc);
        for b in &mut s[500..520] {
            *b = b'N';
        }
        let mut r = Reference::new();
        r.push_contig("a", &s[..len / 2]);
        r.push_contig("b", &s[len / 2..]);
        FmIndex::build(r, Variant::Kompakt, 0).unwrap().0
    }

    /// Ein Read: Stück des Textes mit Fehlern und gelegentlich einem N.
    fn read_from(text: &[u8], rng: &mut Rng) -> Vec<u8> {
        let len = 20 + (rng.next() % 80) as usize;
        let start = (rng.next() as usize) % (text.len() - len);
        let mut r = text[start..start + len].to_vec();
        for _ in 0..rng.next() % 5 {
            let i = (rng.next() as usize) % len;
            r[i] = (rng.next() & 3) as u8;
        }
        if rng.next().is_multiple_of(4) {
            let i = (rng.next() as usize) % len;
            r[i] = 4;
        }
        r
    }

    /// Jeder Eintrag der Tabelle gegen die einfache Suche, und die Suche mit Tabelle gegen
    /// die Suche ohne.
    #[test]
    fn short_table_changes_nothing() {
        let index = genome(5, 1200);
        assert_eq!(index.short.longest, 5);
        for len in 1..=index.short.longest {
            for code in 0..1usize << (2 * len) {
                let p: Vec<u8> = (0..len)
                    .map(|i| (code >> (2 * (len - 1 - i)) & 3) as u8)
                    .collect();
                let rc: Vec<u8> = p.iter().rev().map(|c| 3 - c).collect();
                let iv = index.short.get(&p).unwrap();
                let (lo, hi) = index.interval(&p);
                let (rlo, rhi) = index.interval(&rc);
                assert_eq!((iv.size, iv.size), (hi - lo, rhi - rlo), "Muster {p:?}");
                if iv.size > 0 {
                    assert_eq!((iv.k, iv.kr), (lo, rlo), "Muster {p:?}");
                }
            }
        }
        assert!(index.short.get(&[0; 6]).is_none());
        assert!(index.short.get(&[]).is_none());

        let mut plain = index.clone();
        plain.short = ShortTable::default();
        let text = text_of(&index);
        let mut rng = Rng(23);
        for _ in 0..400 {
            let mut read = read_from(&text, &mut rng);
            if rng.next().is_multiple_of(4) {
                let i = (rng.next() as usize) % read.len();
                read[i] = 4;
            }
            for min_size in [1, 2, 5] {
                for x in 0..read.len() {
                    let (mut with, mut without) = (Vec::new(), Vec::new());
                    let a = index.smems_at(&read, x, min_size, &mut with);
                    let b = plain.smems_at(&read, x, min_size, &mut without);
                    assert_eq!((a, &with), (b, &without));
                    for (min_len, max_size) in [(3, 4), (8, 2), (19, 20)] {
                        assert_eq!(
                            index.seed_forward(&read, x, min_len, max_size),
                            plain.seed_forward(&read, x, min_len, max_size)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn extensions_match_plain_search() {
        let index = genome(3, 1200);
        let text = text_of(&index);
        let mut rng = Rng(11);
        for _ in 0..300 {
            let len = 1 + (rng.next() % 30) as usize;
            let start = (rng.next() as usize) % (text.len() - len);
            let mut p = text[start..start + len].to_vec();
            if rng.next().is_multiple_of(3) {
                let i = (rng.next() as usize) % len;
                p[i] = (rng.next() & 3) as u8;
            }
            let rc: Vec<u8> = p.iter().rev().map(|c| 3 - c).collect();
            let check = |iv: BiInterval, upto: &[u8]| {
                let rcu: Vec<u8> = upto.iter().rev().map(|c| 3 - c).collect();
                let (lo, hi) = index.interval(upto);
                let (rlo, rhi) = index.interval(&rcu);
                assert_eq!(iv.size, hi - lo);
                assert_eq!(iv.size, rhi - rlo);
                if iv.size > 0 {
                    assert_eq!((iv.k, iv.kr), (lo, rlo), "Muster {upto:?}");
                }
            };
            // nach rechts aufbauen
            let mut iv = index.bi_start(p[0]);
            check(iv, &p[..1]);
            for i in 1..len {
                iv = index.extend_right(iv, p[i]);
                check(iv, &p[..=i]);
                if iv.size == 0 {
                    break;
                }
            }
            // nach links aufbauen
            let mut iv = index.bi_start(p[len - 1]);
            for i in (0..len - 1).rev() {
                iv = index.extend_left(iv, p[i]);
                check(iv, &p[i..]);
                if iv.size == 0 {
                    break;
                }
            }
            let _ = rc;
        }
    }

    /// Unabhängige Wahrheit: alle Abschnitte, die x überdecken, mindestens `min_size`-mal
    /// vorkommen und in keine Richtung verlängerbar sind.
    fn naive_at(text: &[u8], read: &[u8], x: usize, min_size: u64) -> Vec<(usize, usize, u64)> {
        let ok = |a: usize, b: usize| -> u64 {
            if read[a..b].iter().any(|&c| c > 3) {
                0
            } else {
                occurrences(text, &read[a..b])
            }
        };
        let mut out = Vec::new();
        for a in 0..=x {
            for b in x + 1..=read.len() {
                let n = ok(a, b);
                if n < min_size {
                    continue;
                }
                let left = a > 0 && ok(a - 1, b) >= min_size;
                let right = b < read.len() && ok(a, b + 1) >= min_size;
                if !left && !right {
                    out.push((a, b, n));
                }
            }
        }
        out
    }

    #[test]
    fn smems_at_match_naive_search() {
        let index = genome(5, 1200);
        let text = text_of(&index);
        let mut rng = Rng(21);
        let mut found = 0;
        for round in 0..120 {
            let read = read_from(&text, &mut rng);
            let x = (rng.next() as usize) % read.len();
            let min_size = 1 + round % 3;
            let mut got = Vec::new();
            let reach = index.smems_at(&read, x, min_size, &mut got);
            let got: Vec<(usize, usize, u64)> = got
                .iter()
                .map(|m| (m.start as usize, m.end as usize, m.interval.size))
                .collect();
            let want = naive_at(&text, &read, x, min_size);
            assert_eq!(
                got, want,
                "Runde {round} x={x} min={min_size} Read {read:?}"
            );
            let far = want.iter().map(|w| w.1).max().unwrap_or(x + 1);
            assert_eq!(reach, far);
            found += want.len();
        }
        assert!(found > 150);
    }

    #[test]
    fn forward_seeds_match_naive_search() {
        let index = genome(13, 1200);
        let text = text_of(&index);
        let mut rng = Rng(41);
        let mut found = 0;
        for round in 0..200 {
            let read = read_from(&text, &mut rng);
            let x = (rng.next() as usize) % read.len();
            let (min_len, max_size) = (5 + (round % 7) as u32, 1 + (round % 4) as u64);
            let got = index.seed_forward(&read, x, min_len, max_size);
            // Vergleich: erstes Ende, an dem das Muster lang genug und selten genug ist.
            let mut want = (read.len(), None);
            if read[x] > 3 {
                want = (x + 1, None);
            } else {
                for i in x + 1..read.len() {
                    if read[i] > 3 {
                        want = (i + 1, None);
                        break;
                    }
                    let n = occurrences(&text, &read[x..=i]);
                    if n < max_size && i - x >= min_len as usize {
                        want = (i + 1, (n > 0).then_some((x, i + 1, n)));
                        break;
                    }
                }
            }
            let got = (
                got.0,
                got.1
                    .map(|m| (m.start as usize, m.end as usize, m.interval.size)),
            );
            assert_eq!(got, want, "Runde {round} x={x} Read {read:?}");
            found += usize::from(want.1.is_some());
        }
        assert!(found > 40);
    }

    #[test]
    fn smems_of_read_match_naive_search() {
        let index = genome(8, 1200);
        let text = text_of(&index);
        let mut rng = Rng(31);
        for round in 0..150 {
            let read = read_from(&text, &mut rng);
            // Längste Übereinstimmung ab jeder Position.
            let longest: Vec<usize> = (0..read.len())
                .map(|i| {
                    let mut l = 0;
                    while i + l < read.len()
                        && read[i + l] <= 3
                        && occurrences(&text, &read[i..=i + l]) > 0
                    {
                        l += 1;
                    }
                    l
                })
                .collect();
            let mut want = Vec::new();
            let mut prev_end = 0;
            for (i, &l) in longest.iter().enumerate() {
                if l > 0 && i + l > prev_end {
                    want.push((i, i + l));
                }
                prev_end = prev_end.max(i + l);
            }
            let got = index.smems(&read, 1);
            let pairs: Vec<(usize, usize)> = got
                .iter()
                .map(|m| (m.start as usize, m.end as usize))
                .collect();
            assert_eq!(pairs, want, "Runde {round} Read {read:?}");
            for m in &got {
                let p = &read[m.start as usize..m.end as usize];
                assert_eq!(
                    index.interval(p),
                    (m.interval.k, m.interval.k + m.interval.size)
                );
            }
            // Mindestlänge filtert nur.
            let long: Vec<Smem> = got.iter().copied().filter(|m| m.len() >= 12).collect();
            assert_eq!(index.smems(&read, 12), long);
        }
    }
}
