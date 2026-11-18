//! `gc-align`: Reads ausrichten nach dem Verfahren von BWA-MEM (Meilenstein M3).
//!
//! Bisher gebaut: Seeds sammeln (SMEM mit Nachsuche), im Genom orten, zu Ketten verbinden,
//! filtern, zu Bereichen verlängern, Bereiche bereinigen und zusammenfügen, primär bestimmen,
//! MAPQ, das fertige Alignment mit CIGAR, NM und MD, die Ausgabe als SAM und Read-Paare
//! (Insertgröße, Partnersuche, Paarbewertung). Es fehlen BAM und mehrere Kerne.

pub mod align;
pub mod batch;
pub mod chain;
pub mod extend;
pub mod finalize;
pub mod global;
pub mod kbtree;
pub mod ksort;
pub mod local;
pub mod pair;
pub mod region;
pub mod sam;
pub mod seed;

pub use align::{Alignment, region_to_alignment};
pub use batch::{ReadInput, align_pairs, align_reads, parallel_map};
pub use chain::{Chain, ChainOptions, FilterOptions, chain_seeds, filter_chains};
pub use extend::{Extension, Scoring, extend};
pub use finalize::{FinalizeOptions, mapq, mark_primary, sort_dedup};
pub use global::{Cigar, GlobalAlignment, Op, cigar_string, gen_cigar, global};
pub use local::{LocalResult, local_align};
pub use pair::{
    BestPair, DIRECTIONS, InsertStats, PairOptions, PairedRead, best_pair, estimate_insert,
    mate_rescue, pair_lines, write_sam_pair,
};
pub use region::{Region, RegionOptions, chain_to_regions};
pub use sam::{
    AlignOptions, Line, SamOptions, TagValue, align_read, find_regions, read_lines, sam_header,
    write_sam,
};
pub use seed::{
    Seed, SeedOptions, collect_intervals, collect_seeds, locate_seeds, repeat_fraction,
};
