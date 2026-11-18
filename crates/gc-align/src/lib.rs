//! `gc-align`: Reads ausrichten nach dem Verfahren von BWA-MEM (Meilenstein M3).
//!
//! Bisher gebaut: Seeds sammeln (SMEM mit Nachsuche), im Genom orten, zu Ketten verbinden,
//! filtern, zu Bereichen verlängern, Bereiche bereinigen und zusammenfügen, primär bestimmen,
//! MAPQ, und das globale Alignment mit CIGAR. Es fehlen CIGAR je Read, Paarung und Ausgabe.

pub mod chain;
pub mod extend;
pub mod finalize;
pub mod global;
pub mod kbtree;
pub mod ksort;
pub mod region;
pub mod seed;

pub use chain::{Chain, ChainOptions, FilterOptions, chain_seeds, filter_chains};
pub use extend::{Extension, Scoring, extend};
pub use finalize::{FinalizeOptions, mapq, mark_primary, sort_dedup};
pub use global::{Cigar, GlobalAlignment, Op, cigar_string, gen_cigar, global};
pub use region::{Region, RegionOptions, chain_to_regions};
pub use seed::{
    Seed, SeedOptions, collect_intervals, collect_seeds, locate_seeds, repeat_fraction,
};
