//! Sortieren von BAM-Datensätzen nach Ort, mit festem Speicher.
//!
//! Die Datensätze kommen fertig kodiert an (jeder mit seinem Längenfeld vorn). Sie sammeln sich
//! im Speicher; ist das Budget voll, wird der Inhalt sortiert und als Teil in den Temp-Ordner
//! geschrieben (komprimiert, mit Zahl, Länge und Prüfsumme). Am Ende werden alle Teile und der
//! Rest im Speicher zusammengeführt.
//!
//! Reihenfolge wie `samtools sort`: Contig (ohne Contig zuletzt), Position, Vorwärtsstrang vor
//! Gegenstrang; gleiche Datensätze bleiben in der Reihenfolge der Eingabe.

use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::bufread::DeflateDecoder;
use flate2::write::DeflateEncoder;

const MAGIC: &[u8; 8] = b"GCSORT1\n";
/// Feste Felder eines BAM-Satzes hinter dem Längenfeld (SAM-Spezifikation 4.2).
const FIXED: usize = 32;

/// Sortierschlüssel eines kodierten Satzes (`record` ohne Längenfeld): Contig, Position, Strang.
/// Sätze ohne Contig bekommen den größten Schlüsselanfang und stehen damit am Ende.
pub fn sort_key(record: &[u8]) -> u64 {
    let ref_id = i32::from_le_bytes([record[0], record[1], record[2], record[3]]);
    let pos = i32::from_le_bytes([record[4], record[5], record[6], record[7]]);
    let flag = u16::from_le_bytes([record[14], record[15]]);
    // 31 Bit Contig (ohne Contig: der größte Wert), 32 Bit Position + 1 (kein Ort: 0), 1 Bit Strang.
    let contig = if ref_id < 0 {
        0x7FFF_FFFF
    } else {
        ref_id as u64
    };
    let pos = (i64::from(pos) + 1).max(0) as u64;
    contig << 33 | pos << 1 | u64::from(flag & 0x10 != 0)
}

#[derive(Debug, Clone)]
pub struct SortOptions {
    /// So viele Bytes kodierter Sätze liegen höchstens im Speicher. Dazu kommen 16 Bytes je Satz
    /// für die Sortierung.
    pub memory_bytes: usize,
    /// Ordner für die Teile. Darin entsteht ein eigener Unterordner, der am Ende gelöscht wird.
    pub temp_dir: PathBuf,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SortStats {
    pub records: u64,
    /// Bytes aller Sätze mit Längenfeld.
    pub bytes: u64,
    /// Teile, die auf die Platte geschrieben wurden.
    pub runs: u64,
    /// Größe dieser Teile auf der Platte.
    pub run_bytes: u64,
}

pub struct Sorter {
    opt: SortOptions,
    buf: Vec<u8>,
    /// (Schlüssel, Anfang des Satzes in `buf`)
    index: Vec<(u64, usize)>,
    dir: Option<PathBuf>,
    runs: Vec<PathBuf>,
    stats: SortStats,
    /// Teil, der gerade im Hintergrund geschrieben wird; liefert seine Größe auf der Platte.
    pending: Option<std::thread::JoinHandle<io::Result<u64>>>,
    /// Teile beim Zusammenführen in eigenen Threads lesen und entpacken ([`Sorter::set_read_ahead`]);
    /// `None`: nur bei höchstens [`AHEAD_MAX_RUNS`] Teilen.
    read_ahead: Option<bool>,
    /// Immer den Lese-Pool nehmen (`GC_SORT_AHEAD=2`, Tests).
    force_pool: bool,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

impl Sorter {
    pub fn new(opt: SortOptions) -> Self {
        Sorter {
            opt,
            buf: Vec::new(),
            index: Vec::new(),
            dir: None,
            runs: Vec::new(),
            stats: SortStats::default(),
            pending: None,
            read_ahead: match std::env::var_os("GC_SORT_AHEAD") {
                Some(v) if v == "0" => Some(false),
                Some(v) if v == "1" => Some(true),
                _ => None,
            },
            force_pool: std::env::var_os("GC_SORT_AHEAD").is_some_and(|v| v == "2"),
        }
    }

    /// Beim Zusammenführen liest und entpackt je Teil ein eigener Thread vorab, in Stücken von
    /// ganzen Sätzen. Vorgabe: nur bei höchstens [`AHEAD_MAX_RUNS`] Teilen, weil hunderte
    /// Threads das Mischen bremsen (Mac, 418 Teile: 24,6 s mit, 13,9 s ohne). `GC_SORT_AHEAD=0`
    /// schaltet es immer ab, `GC_SORT_AHEAD=1` immer an. Gleiche Sätze in gleicher Reihenfolge.
    pub fn set_read_ahead(&mut self, on: bool) {
        self.read_ahead = Some(on);
    }

    /// Teile immer über den Lese-Pool lesen (sonst erst ab [`AHEAD_MAX_RUNS`] Teilen).
    pub fn set_read_pool(&mut self) {
        self.force_pool = true;
    }

    /// Nimmt kodierte Sätze hintereinander, jeden mit seinem Längenfeld. `records` muss mit
    /// einem ganzen Satz enden.
    pub fn push(&mut self, records: &[u8]) -> io::Result<()> {
        let mut at = 0;
        while at < records.len() {
            if records.len() - at < 4 {
                return Err(invalid("abgeschnittenes Längenfeld"));
            }
            let len = u32::from_le_bytes([
                records[at],
                records[at + 1],
                records[at + 2],
                records[at + 3],
            ]) as usize;
            let end = at + 4 + len;
            if len < FIXED || end > records.len() {
                return Err(invalid("abgeschnittener oder zu kurzer BAM-Satz"));
            }
            if !self.buf.is_empty() && self.buf.len() + 4 + len > self.opt.memory_bytes {
                self.spill()?;
            }
            self.index
                .push((sort_key(&records[at + 4..end]), self.buf.len()));
            self.buf.extend_from_slice(&records[at..end]);
            self.stats.records += 1;
            self.stats.bytes += (4 + len) as u64;
            at = end;
        }
        Ok(())
    }

    /// Zwischenstand: bisher angenommene Sätze und geschriebene Teile.
    pub fn stats(&self) -> SortStats {
        self.stats
    }

    fn record_at(&self, offset: usize) -> &[u8] {
        let b = &self.buf[offset..];
        let len = u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize;
        &b[..4 + len]
    }

    /// Sortiert den Speicherinhalt. Stabil: Bei gleichem Schlüssel entscheidet die Stelle im
    /// Puffer, also die Reihenfolge der Eingabe.
    fn sort_memory(&mut self) {
        self.index.sort_unstable();
    }

    /// Schreibt den Speicherinhalt als sortierten Teil auf die Platte.
    fn spill(&mut self) -> io::Result<()> {
        self.sort_memory();
        let dir = match &self.dir {
            Some(d) => d.clone(),
            None => {
                let nonce = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_nanos());
                let d = self
                    .opt
                    .temp_dir
                    .join(format!("genome-sort-{}-{nonce}", std::process::id()));
                std::fs::create_dir_all(&d)?;
                self.dir = Some(d.clone());
                d
            }
        };
        let path = dir.join(format!("teil-{:06}.gcsort", self.runs.len()));
        // Der vorige Teil muss auf der Platte sein, bevor der nächste entsteht.
        self.wait_for_pending()?;
        if std::env::var_os("GC_SORT_PARALLEL").is_some_and(|v| v == "0") {
            let mut file = BufWriter::with_capacity(1 << 20, File::create(&path)?);
            file.write_all(MAGIC)?;
            let (crc, bytes) = self.write_run_serial(&mut file)?;
            let file = finish_run(file, self.index.len() as u64, bytes, &crc)?;
            self.stats.runs += 1;
            self.stats.run_bytes += file.metadata()?.len();
        } else {
            // Komprimieren auf mehreren Kernen in den Speicher, Schreiben im Hintergrund: Die
            // Platte kann langsam sein (eine USB-SSD schreibt teils nur 20 MB/s), das Annehmen
            // neuer Sätze wartet darauf nicht. Zusätzlicher Speicher: der komprimierte Teil.
            let mut data = Vec::with_capacity(self.buf.len() / 3 + 1024);
            data.extend_from_slice(MAGIC);
            let (crc, bytes) = self.write_run_parallel(&mut data)?;
            let records = self.index.len() as u64;
            let target = path.clone();
            self.pending = Some(std::thread::spawn(move || -> io::Result<u64> {
                let mut file = BufWriter::with_capacity(1 << 20, File::create(&target)?);
                file.write_all(&data)?;
                drop(data);
                let file = finish_run(file, records, bytes, &crc)?;
                Ok(file.metadata()?.len())
            }));
            self.stats.runs += 1;
        }
        self.runs.push(path);
        self.buf.clear();
        self.index.clear();
        Ok(())
    }

    /// Wartet, bis der Teil im Hintergrund geschrieben ist, und meldet seinen Fehler.
    fn wait_for_pending(&mut self) -> io::Result<()> {
        if let Some(h) = self.pending.take() {
            let len = match h.join() {
                Ok(r) => r?,
                Err(panic) => std::panic::resume_unwind(panic),
            };
            self.stats.run_bytes += len;
        }
        Ok(())
    }

    /// Schreibt die sortierten Sätze als einen Deflate-Strom (bisheriger Weg, ein Kern).
    fn write_run_serial(&self, file: &mut impl Write) -> io::Result<(flate2::Crc, u64)> {
        let mut crc = flate2::Crc::new();
        let mut bytes = 0u64;
        let mut enc = DeflateEncoder::new(file, Compression::fast());
        for &(_, offset) in &self.index {
            let rec = self.record_at(offset);
            crc.update(rec);
            bytes += rec.len() as u64;
            enc.write_all(rec)?;
        }
        enc.finish()?;
        Ok((crc, bytes))
    }

    /// Wie `write_run_serial`, aber in Stücken von `SPILL_CHUNK` Bytes, die mehrere Kerne
    /// zugleich komprimieren. Jedes Stück beginnt ohne Vorgeschichte und endet mit einem
    /// Sync-Flush, nur das letzte schließt den Strom ab; aneinandergehängt ergeben sie einen
    /// gültigen Deflate-Strom (wie bei pigz). Das Lesen bleibt unverändert. Das Ausrichten wartet
    /// so nicht mehr viele Sekunden auf jeden Teil. Im Speicher liegen zusätzlich höchstens
    /// `Kerne × SPILL_CHUNK` Bytes Rohdaten und ihre komprimierte Form.
    fn write_run_parallel(&self, file: &mut impl Write) -> io::Result<(flate2::Crc, u64)> {
        const SPILL_CHUNK: usize = 8 << 20;
        // Stückgrenzen in `index`: jedes Stück fasst ganze Sätze bis rund SPILL_CHUNK Bytes.
        let mut bounds = vec![0usize];
        let mut size = 0usize;
        for (i, &(_, offset)) in self.index.iter().enumerate() {
            size += self.record_at(offset).len();
            if size >= SPILL_CHUNK {
                bounds.push(i + 1);
                size = 0;
            }
        }
        if *bounds.last().unwrap_or(&0) != self.index.len() || bounds.len() == 1 {
            bounds.push(self.index.len());
        }
        let chunks: Vec<(usize, usize)> = bounds.windows(2).map(|w| (w[0], w[1])).collect();
        let workers = std::thread::available_parallelism()
            .map_or(1, usize::from)
            .clamp(1, 8);
        let mut crc = flate2::Crc::new();
        let mut bytes = 0u64;
        for (r, round) in chunks.chunks(workers).enumerate() {
            let first = r * workers;
            let done: Vec<io::Result<(Vec<u8>, flate2::Crc, u64)>> = std::thread::scope(|scope| {
                let handles: Vec<_> = round
                    .iter()
                    .enumerate()
                    .map(|(k, &(from, to))| {
                        let last = first + k + 1 == chunks.len();
                        scope.spawn(move || self.compress_chunk(from, to, last))
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|h| match h.join() {
                        Ok(r) => r,
                        Err(panic) => std::panic::resume_unwind(panic),
                    })
                    .collect()
            });
            for part in done {
                let (data, part_crc, part_bytes) = part?;
                file.write_all(&data)?;
                crc.combine(&part_crc);
                bytes += part_bytes;
            }
        }
        Ok((crc, bytes))
    }

    /// Komprimiert die Sätze `index[from..to]` roh (Deflate ohne Kopf); `last` schließt den
    /// Strom ab, sonst endet das Stück mit einem Sync-Flush auf einer Bytegrenze.
    fn compress_chunk(
        &self,
        from: usize,
        to: usize,
        last: bool,
    ) -> io::Result<(Vec<u8>, flate2::Crc, u64)> {
        let mut raw = Vec::new();
        for &(_, offset) in &self.index[from..to] {
            raw.extend_from_slice(self.record_at(offset));
        }
        let mut crc = flate2::Crc::new();
        crc.update(&raw);
        let mut c = flate2::Compress::new(Compression::fast(), false);
        let flush = if last {
            flate2::FlushCompress::Finish
        } else {
            flate2::FlushCompress::Sync
        };
        let mut out = Vec::with_capacity(raw.len() / 2 + 1024);
        loop {
            let consumed = c.total_in() as usize;
            let status = c
                .compress_vec(&raw[consumed..], &mut out, flush)
                .map_err(|e| io::Error::other(format!("Komprimieren: {e}")))?;
            let all_in = c.total_in() as usize == raw.len();
            let room_left = out.len() < out.capacity();
            match status {
                flate2::Status::StreamEnd => break,
                _ if !last && all_in && room_left => break,
                _ => out.reserve(out.capacity().max(1 << 16)),
            }
        }
        Ok((out, crc, raw.len() as u64))
    }

    /// Gibt alle Sätze sortiert an `sink`, jeden mit seinem Längenfeld, und räumt die Teile weg.
    pub fn finish(
        mut self,
        mut sink: impl FnMut(&[u8]) -> io::Result<()>,
    ) -> io::Result<SortStats> {
        self.wait_for_pending()?;
        self.sort_memory();
        if self.runs.is_empty() {
            for &(_, offset) in &self.index {
                sink(self.record_at(offset))?;
            }
            return Ok(self.stats);
        }
        // Quellen in der Reihenfolge der Eingabe: erst die Teile, zuletzt der Speicher. Bei
        // gleichem Schlüssel gewinnt die frühere Quelle, so bleibt die Sortierung stabil.
        let mode = match self.read_ahead {
            Some(false) => Ahead::Off,
            Some(true) => Ahead::PerRun,
            None if self.runs.len() <= AHEAD_MAX_RUNS => Ahead::PerRun,
            None => Ahead::Pool,
        };
        let mode = if self.force_pool { Ahead::Pool } else { mode };
        let read_buffer = read_buffer_bytes(self.runs.len());
        let mut readers = Vec::with_capacity(self.runs.len());
        match mode {
            Ahead::Off | Ahead::PerRun => {
                for path in &self.runs {
                    let r = RunReader::open(path, read_buffer)?;
                    readers.push(if mode == Ahead::PerRun {
                        Source::Ahead(AheadReader::start(r))
                    } else {
                        Source::Plain(Box::new(r))
                    });
                }
            }
            Ahead::Pool => {
                let workers = std::env::var("GC_SORT_POOL_THREADS")
                    .ok()
                    .and_then(|v| v.parse::<usize>().ok())
                    .filter(|&n| n > 0)
                    .unwrap_or_else(|| {
                        std::thread::available_parallelism()
                            .map_or(1, |n| n.get())
                            .clamp(1, POOL_MAX_THREADS)
                    });
                // Geöffnet wird im Pool: Das erste Füllen der großen Puffer (bei 418 Teilen
                // und 8 MB über 3 GB) soll nicht nacheinander im Hauptthread geschehen.
                let mut groups: Vec<Vec<(PathBuf, ChunkTx)>> =
                    (0..workers).map(|_| Vec::new()).collect();
                for (i, path) in self.runs.iter().enumerate() {
                    let (tx, rx) = std::sync::mpsc::sync_channel(1);
                    groups[i % workers].push((path.clone(), tx));
                    readers.push(Source::Ahead(AheadReader::from_rx(rx)));
                }
                for g in groups {
                    std::thread::spawn(move || pool_worker(g, read_buffer, POOL_CHUNK));
                }
            }
        }
        let memory = readers.len();
        let mut heap = BinaryHeap::new();
        let mut current: Vec<Vec<u8>> = vec![Vec::new(); memory];
        for (i, r) in readers.iter_mut().enumerate() {
            if r.next(&mut current[i])? {
                heap.push(std::cmp::Reverse((sort_key(&current[i][4..]), i)));
            }
        }
        let mut next_in_memory = 0;
        if let Some(&(key, _)) = self.index.first() {
            heap.push(std::cmp::Reverse((key, memory)));
        }
        while let Some(std::cmp::Reverse((_, source))) = heap.pop() {
            if source == memory {
                sink(self.record_at(self.index[next_in_memory].1))?;
                next_in_memory += 1;
                if let Some(&(key, _)) = self.index.get(next_in_memory) {
                    heap.push(std::cmp::Reverse((key, memory)));
                }
            } else {
                sink(&current[source])?;
                if readers[source].next(&mut current[source])? {
                    heap.push(std::cmp::Reverse((sort_key(&current[source][4..]), source)));
                }
            }
        }
        Ok(self.stats)
    }
}

impl Drop for Sorter {
    fn drop(&mut self) {
        if let Some(h) = self.pending.take() {
            let _ = h.join();
        }
        if let Some(dir) = &self.dir {
            let _ = std::fs::remove_dir_all(dir);
        }
    }
}

/// Schluss eines Teils: Zahl der Sätze, Bytes, Prüfsumme. Daran erkennt das Lesen einen
/// abgeschnittenen oder beschädigten Teil. Danach liegt der Teil sicher auf der Platte.
fn finish_run(
    mut file: BufWriter<File>,
    records: u64,
    bytes: u64,
    crc: &flate2::Crc,
) -> io::Result<File> {
    file.write_all(&records.to_le_bytes())?;
    file.write_all(&bytes.to_le_bytes())?;
    file.write_all(&crc.sum().to_le_bytes())?;
    file.flush()?;
    let file = file.into_inner().map_err(|e| e.into_error())?;
    file.sync_all()?;
    Ok(file)
}

/// Ein Teil beim Zusammenführen: direkt gelesen oder von einem Thread vorab.
enum Source {
    Plain(Box<RunReader>),
    Ahead(AheadReader),
}

impl Source {
    fn next(&mut self, record: &mut Vec<u8>) -> io::Result<bool> {
        match self {
            Source::Plain(r) => r.next(record),
            Source::Ahead(r) => r.next(record),
        }
    }
}

/// Höchstzahl Teile, bis zu der je Teil ein eigener Lese-Thread läuft. Bei mehr Teilen liest
/// ein kleiner Pool reihum ([`pool_worker`]).
pub const AHEAD_MAX_RUNS: usize = 32;

/// Lesepuffer aller Teile zusammen beim Zusammenführen.
const READ_BUFFER_TOTAL: usize = 2 << 30;

/// Lesepuffer je Teil: 2 GiB geteilt durch die Zahl der Teile, mindestens 256 KB, höchstens
/// 8 MB (`GC_SORT_LESEPUFFER_KB` setzt ihn fest). Mac M5, USB-SSD, 418 Teile mit 256 KB:
/// die Lese-Threads hingen zu 70–89 % im Lesen, ~45 MB/s statt ~390 MB/s.
fn read_buffer_bytes(runs: usize) -> usize {
    if let Some(kb) = std::env::var("GC_SORT_LESEPUFFER_KB")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&kb| kb > 0)
    {
        return kb << 10;
    }
    (READ_BUFFER_TOTAL / runs.max(1)).clamp(1 << 18, 8 << 20)
}

/// Höchstzahl Threads im Lese-Pool.
const POOL_MAX_THREADS: usize = 4;

/// Bytes je Stück im Lese-Pool (kleiner als [`AHEAD_CHUNK`], weil hunderte Teile warten).
const POOL_CHUNK: usize = 1 << 18;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Ahead {
    /// Jeder Teil im Hauptthread.
    Off,
    /// Je Teil ein Thread.
    PerRun,
    /// Wenige Threads lesen reihum alle Teile.
    Pool,
}

type ChunkTx = std::sync::mpsc::SyncSender<io::Result<Vec<u8>>>;

/// Liest bis `size` Bytes ganzer Sätze aus `r`. Leer heißt: Ende, Prüfung bestanden.
fn read_chunk(r: &mut RunReader, record: &mut Vec<u8>, size: usize) -> io::Result<Vec<u8>> {
    let mut chunk = Vec::with_capacity(size + 4096);
    while chunk.len() < size {
        if !r.next(record)? {
            break;
        }
        chunk.extend_from_slice(record);
    }
    Ok(chunk)
}

/// Ein Thread des Lese-Pools: liest und entpackt seine Teile reihum. Jeder Teil hat höchstens
/// ein wartendes Stück im Kanal und eines hier; ist ein Kanal voll, geht es mit dem nächsten
/// Teil weiter (nie blockierend, sonst könnte der Hauptthread auf einen Teil warten, dessen
/// Thread an einem anderen hängt). Stücke enden immer an Satzgrenzen; ein leeres Stück meldet
/// das geprüfte Ende, ein Fehler beendet den Teil.
fn pool_worker(parts: Vec<(PathBuf, ChunkTx)>, buffer: usize, size: usize) {
    use std::sync::mpsc::TrySendError;
    let n = parts.len();
    let (paths, senders): (Vec<PathBuf>, Vec<ChunkTx>) = parts.into_iter().unzip();
    let mut open: Vec<Option<RunReader>> = (0..n).map(|_| None).collect();
    let mut pending: Vec<Option<io::Result<Vec<u8>>>> = (0..n).map(|_| None).collect();
    let mut finished = vec![false; n];
    let mut ended = vec![false; n];
    let mut record = Vec::new();
    loop {
        let mut progress = false;
        let mut alive = false;
        for i in 0..n {
            if finished[i] {
                continue;
            }
            alive = true;
            let item = match pending[i].take() {
                Some(item) => item,
                None if ended[i] => Ok(Vec::new()),
                None => {
                    if open[i].is_none() {
                        match RunReader::open(&paths[i], buffer) {
                            Ok(r) => open[i] = Some(r),
                            Err(e) => {
                                pending[i] = Some(Err(e));
                            }
                        }
                    }
                    match (pending[i].take(), open[i].as_mut()) {
                        (Some(err), _) => err,
                        (None, Some(r)) => {
                            let c = read_chunk(r, &mut record, size);
                            if matches!(&c, Ok(v) if v.len() < size) {
                                ended[i] = true;
                            }
                            c
                        }
                        (None, None) => Err(invalid("Teil nicht geöffnet")),
                    }
                }
            };
            let last = match &item {
                Ok(v) => v.is_empty(),
                Err(_) => true,
            };
            match senders[i].try_send(item) {
                Ok(()) => {
                    progress = true;
                    finished[i] = last;
                }
                Err(TrySendError::Full(item)) => pending[i] = Some(item),
                Err(TrySendError::Disconnected(_)) => finished[i] = true,
            }
        }
        if !alive {
            return;
        }
        if !progress {
            std::thread::sleep(std::time::Duration::from_micros(200));
        }
    }
}

/// Bytes je Stück, das der Lese-Thread weitergibt (ganze Sätze).
const AHEAD_CHUNK: usize = 1 << 19;

/// Liest einen Teil in einem eigenen Thread mit [`RunReader`] (dieselben Prüfungen) und gibt
/// ganze Sätze in Stücken weiter. Höchstens zwei Stücke warten, damit der Speicher begrenzt ist.
struct AheadReader {
    rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    chunk: Vec<u8>,
    at: usize,
    done: bool,
}

impl AheadReader {
    fn start(mut reader: RunReader) -> Self {
        let (tx, rx) = std::sync::mpsc::sync_channel(2);
        std::thread::spawn(move || {
            let mut record = Vec::new();
            loop {
                let mut chunk = Vec::with_capacity(AHEAD_CHUNK + 4096);
                let mut end = false;
                while chunk.len() < AHEAD_CHUNK {
                    match reader.next(&mut record) {
                        Ok(true) => chunk.extend_from_slice(&record),
                        Ok(false) => {
                            end = true;
                            break;
                        }
                        Err(e) => {
                            let _ = tx.send(Err(e));
                            return;
                        }
                    }
                }
                if !chunk.is_empty() && tx.send(Ok(chunk)).is_err() {
                    return;
                }
                if end {
                    // Leeres Stück: Ende, Prüfung am Schluss des Teils bestanden.
                    let _ = tx.send(Ok(Vec::new()));
                    return;
                }
            }
        });
        AheadReader::from_rx(rx)
    }

    fn from_rx(rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>) -> Self {
        AheadReader {
            rx,
            chunk: Vec::new(),
            at: 0,
            done: false,
        }
    }

    fn next(&mut self, record: &mut Vec<u8>) -> io::Result<bool> {
        if self.at >= self.chunk.len() {
            if self.done {
                return Ok(false);
            }
            match self.rx.recv() {
                Ok(Ok(chunk)) if chunk.is_empty() => {
                    self.done = true;
                    return Ok(false);
                }
                Ok(Ok(chunk)) => {
                    self.chunk = chunk;
                    self.at = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Err(invalid("Lese-Thread eines Teils ohne Ende abgebrochen")),
            }
        }
        let b = &self.chunk[self.at..];
        let len = u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize;
        record.clear();
        record.extend_from_slice(&b[..4 + len]);
        self.at += 4 + len;
        Ok(true)
    }
}

/// Liest einen Teil Satz für Satz und prüft am Ende Zahl, Länge und Prüfsumme.
struct RunReader {
    path: PathBuf,
    inner: DeflateDecoder<BufReader<File>>,
    crc: flate2::Crc,
    records: u64,
    bytes: u64,
    done: bool,
}

impl RunReader {
    /// `buffer`: Bytes je Lesezugriff. Große Zugriffe halten eine USB-SSD bei hunderten
    /// abwechselnd gelesenen Teilen schnell.
    fn open(path: &Path, buffer: usize) -> io::Result<Self> {
        let mut file = BufReader::with_capacity(buffer, File::open(path)?);
        let mut magic = [0u8; 8];
        file.read_exact(&mut magic)
            .map_err(|_| invalid(format!("{}: Teil ist abgeschnitten", path.display())))?;
        if &magic != MAGIC {
            return Err(invalid(format!("{}: kein Sortier-Teil", path.display())));
        }
        Ok(RunReader {
            path: path.to_path_buf(),
            inner: DeflateDecoder::new(file),
            crc: flate2::Crc::new(),
            records: 0,
            bytes: 0,
            done: false,
        })
    }

    fn damaged(&self, what: &str) -> io::Error {
        invalid(format!(
            "{}: Teil ist beschädigt ({what})",
            self.path.display()
        ))
    }

    /// Nächster Satz nach `record` (mit Längenfeld). `false` am Ende, nach bestandener Prüfung.
    fn next(&mut self, record: &mut Vec<u8>) -> io::Result<bool> {
        if self.done {
            return Ok(false);
        }
        let mut len = [0u8; 4];
        let mut got = 0;
        while got < 4 {
            let n = self
                .inner
                .read(&mut len[got..])
                .map_err(|_| self.damaged("Daten nicht lesbar"))?;
            if n == 0 {
                break;
            }
            got += n;
        }
        if got == 0 {
            self.done = true;
            let mut tail = [0u8; 20];
            let file = self.inner.get_mut();
            file.read_exact(&mut tail)
                .map_err(|_| self.damaged("Schluss fehlt"))?;
            let records = u64::from_le_bytes(tail[0..8].try_into().unwrap_or_default());
            let bytes = u64::from_le_bytes(tail[8..16].try_into().unwrap_or_default());
            let crc = u32::from_le_bytes(tail[16..20].try_into().unwrap_or_default());
            if records != self.records || bytes != self.bytes {
                return Err(self.damaged("Zahl der Sätze oder Länge stimmt nicht"));
            }
            if crc != self.crc.sum() {
                return Err(self.damaged("Prüfsumme stimmt nicht"));
            }
            return Ok(false);
        }
        if got < 4 {
            return Err(self.damaged("Längenfeld abgeschnitten"));
        }
        let n = u32::from_le_bytes(len) as usize;
        if n < FIXED {
            return Err(self.damaged("Satz zu kurz"));
        }
        record.clear();
        record.extend_from_slice(&len);
        record.resize(4 + n, 0);
        self.inner
            .read_exact(&mut record[4..])
            .map_err(|_| self.damaged("Satz abgeschnitten"))?;
        self.crc.update(record);
        self.records += 1;
        self.bytes += record.len() as u64;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    /// Ein Satz mit Längenfeld: Contig, Position, Flag und eine laufende Nummer als Inhalt.
    fn record(ref_id: i32, pos: i32, flag: u16, serial: u32, extra: usize) -> Vec<u8> {
        let mut body = vec![0u8; FIXED + 4 + extra];
        body[0..4].copy_from_slice(&ref_id.to_le_bytes());
        body[4..8].copy_from_slice(&pos.to_le_bytes());
        body[14..16].copy_from_slice(&flag.to_le_bytes());
        body[FIXED..FIXED + 4].copy_from_slice(&serial.to_le_bytes());
        let mut out = (body.len() as u32).to_le_bytes().to_vec();
        out.extend(body);
        out
    }

    fn random_records(n: u32, seed: u64) -> Vec<Vec<u8>> {
        let mut state = seed;
        (0..n)
            .map(|serial| {
                let r = xorshift(&mut state);
                let (ref_id, pos) = if r.is_multiple_of(20) {
                    (-1, -1)
                } else {
                    ((r % 3) as i32, ((r >> 8) % 50) as i32 - 1)
                };
                let flag = if r >> 40 & 1 == 1 { 0x10 } else { 0 };
                record(ref_id, pos, flag, serial, (r >> 20) as usize % 40)
            })
            .collect()
    }

    /// Erwartete Reihenfolge: stabil nach (Contig ohne Vorzeichen, Position, Strang).
    fn expected(records: &[Vec<u8>]) -> Vec<Vec<u8>> {
        let mut sorted = records.to_vec();
        sorted.sort_by_key(|r| {
            let b = &r[4..];
            (
                u32::from_le_bytes([b[0], b[1], b[2], b[3]]),
                i32::from_le_bytes([b[4], b[5], b[6], b[7]]),
                b[14] & 0x10,
            )
        });
        sorted
    }

    fn run(records: &[Vec<u8>], memory_bytes: usize, tag: &str) -> (Vec<Vec<u8>>, SortStats) {
        let ahead = run_with(records, memory_bytes, tag, Some(true));
        let plain = run_with(records, memory_bytes, &format!("{tag}-ohne"), Some(false));
        let pool = run_with(records, memory_bytes, &format!("{tag}-pool"), None);
        // Vorauslesen (je Teil ein Thread oder im Pool) ändert nichts am Ergebnis
        assert_eq!(ahead, plain);
        assert_eq!(pool.0, plain.0);
        ahead
    }

    fn run_with(
        records: &[Vec<u8>],
        memory_bytes: usize,
        tag: &str,
        read_ahead: Option<bool>,
    ) -> (Vec<Vec<u8>>, SortStats) {
        let temp_dir = std::env::temp_dir().join(format!("gcbam-{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&temp_dir).unwrap();
        let mut sorter = Sorter::new(SortOptions {
            memory_bytes,
            temp_dir: temp_dir.clone(),
        });
        match read_ahead {
            Some(on) => sorter.set_read_ahead(on),
            None => sorter.set_read_pool(),
        }
        // in ungleichen Stücken hineingeben
        for chunk in records.chunks(7) {
            sorter.push(&chunk.concat()).unwrap();
        }
        let mut out = Vec::new();
        let stats = sorter
            .finish(|r| {
                out.push(r.to_vec());
                Ok(())
            })
            .unwrap();
        // der Unterordner mit den Teilen ist wieder weg
        assert_eq!(std::fs::read_dir(&temp_dir).unwrap().count(), 0);
        std::fs::remove_dir_all(temp_dir).ok();
        (out, stats)
    }

    #[test]
    fn order_is_contig_position_strand_and_stable() {
        let records = random_records(5000, 0x1111_2222_3333_4444);
        let want = expected(&records);
        // alles im Speicher
        let (got, stats) = run(&records, usize::MAX, "mem");
        assert_eq!(got, want);
        assert_eq!((stats.records, stats.runs), (5000, 0));
        // viele kleine Teile, ein Rest im Speicher
        let (got, stats) = run(&records, 20_000, "teile");
        assert_eq!(got, want);
        assert!(stats.runs > 10, "{stats:?}");
        // Budget kleiner als ein Satz: jeder Satz ein eigener Teil
        let few = random_records(50, 7);
        let (got, stats) = run(&few, 1, "einzeln");
        assert_eq!(got, expected(&few));
        assert_eq!(stats.runs, 49);
        // nichts hinein, nichts heraus
        let (got, stats) = run(&[], 100, "leer");
        assert!(got.is_empty());
        assert_eq!(stats, SortStats::default());
    }

    #[test]
    fn unmapped_come_last_and_reverse_after_forward() {
        let records = vec![
            record(-1, -1, 4, 0, 0),
            record(1, 10, 0x10, 1, 0),
            record(1, 10, 0, 2, 0),
            record(0, 99, 0, 3, 0),
            record(1, 10, 0, 4, 0),
            record(1, 9, 0x10, 5, 0),
        ];
        let (got, _) = run(&records, usize::MAX, "ordnung");
        let serials: Vec<u32> = got
            .iter()
            .map(|r| u32::from_le_bytes(r[4 + FIXED..4 + FIXED + 4].try_into().unwrap()))
            .collect();
        assert_eq!(serials, [3, 5, 2, 4, 1, 0]);
    }

    #[test]
    fn truncated_input_is_rejected() {
        let mut sorter = Sorter::new(SortOptions {
            memory_bytes: 1000,
            temp_dir: std::env::temp_dir(),
        });
        let rec = record(0, 1, 0, 0, 0);
        assert!(sorter.push(&rec[..rec.len() - 1]).is_err());
        assert!(sorter.push(&rec[..3]).is_err());
        assert!(sorter.push(&[8, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8]).is_err());
    }

    #[test]
    fn damaged_run_is_detected() {
        let temp_dir = std::env::temp_dir().join(format!("gcbam-schaden-{}", std::process::id()));
        for case in 0..9 {
            std::fs::create_dir_all(&temp_dir).unwrap();
            let mut sorter = Sorter::new(SortOptions {
                memory_bytes: 3000,
                temp_dir: temp_dir.clone(),
            });
            match case / 3 {
                0 => sorter.set_read_ahead(true),
                1 => sorter.set_read_ahead(false),
                _ => sorter.set_read_pool(),
            }
            for r in random_records(400, 99) {
                sorter.push(&r).unwrap();
            }
            let path = sorter.runs[1].clone();
            let mut data = std::fs::read(&path).unwrap();
            match case % 3 {
                0 => {
                    let mid = data.len() / 2;
                    data[mid] ^= 0x40; // ein Bit im Inhalt
                }
                1 => data.truncate(data.len() - 30), // abgeschnitten
                _ => {
                    let last = data.len() - 1;
                    data[last] ^= 1; // Prüfsumme im Schluss
                }
            }
            std::fs::write(&path, data).unwrap();
            let err = sorter.finish(|_| Ok(())).unwrap_err();
            assert!(err.to_string().contains("beschädigt"), "{case}: {err}");
            std::fs::remove_dir_all(&temp_dir).ok();
        }
    }
}
