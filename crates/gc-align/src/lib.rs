//! `gc-align`: Reads ausrichten nach dem Verfahren von BWA-MEM (Meilenstein M3).
//!
//! Bisher gebaut: Seeds sammeln (SMEM mit Nachsuche), im Genom orten, zu Ketten verbinden und
//! einen Seed mit Smith-Waterman verlängern.

pub mod chain;
pub mod extend;
pub mod region;
pub mod seed;

pub use chain::{Chain, ChainOptions, FilterOptions, chain_seeds, filter_chains};
pub use extend::{Extension, Scoring, extend};
pub use region::{Region, RegionOptions, chain_to_regions};
pub use seed::{Seed, SeedOptions, collect_seeds};
