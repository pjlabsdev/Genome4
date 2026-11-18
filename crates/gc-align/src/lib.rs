//! `gc-align`: Reads ausrichten nach dem Verfahren von BWA-MEM (Meilenstein M3).
//!
//! Bisher gebaut: Seeds sammeln (SMEM mit Nachsuche), im Genom orten, zu Ketten verbinden,
//! filtern, zu Bereichen verlängern, Bereiche bereinigen und zusammenfügen, primär bestimmen,
//! MAPQ, das fertige Alignment mit CIGAR, NM und MD und die Ausgabe einzelner Reads als SAM.
//! Es fehlen Paarung, BAM und mehrere Kerne.

pub mod align;
pub mod chain;
pub mod extend;
pub mod finalize;
pub mod global;
pub mod kbtree;
pub mod ksort;
pub mod region;
pub mod sam;
pub mod seed;

pub use align::{Alignment, region_to_alignment};
pub use chain::{Chain, ChainOptions, FilterOptions, chain_seeds, filter_chains};
pub use extend::{Extension, Scoring, extend};
pub use finalize::{FinalizeOptions, mapq, mark_primary, sort_dedup};
pub use global::{Cigar, GlobalAlignment, Op, cigar_string, gen_cigar, global};
pub use region::{Region, RegionOptions, chain_to_regions};
pub use sam::{AlignOptions, SamOptions, align_read, sam_header, write_sam};
pub use seed::{
    Seed, SeedOptions, collect_intervals, collect_seeds, locate_seeds, repeat_fraction,
};
