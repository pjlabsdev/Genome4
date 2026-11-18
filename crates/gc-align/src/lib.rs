//! `gc-align`: Reads ausrichten nach dem Verfahren von BWA-MEM (Meilenstein M3).
//!
//! Bisher gebaut: Seeds sammeln (SMEM mit Nachsuche), im Genom orten, zu Ketten verbinden und
//! einen Seed mit Smith-Waterman verlängern.

pub mod chain;
pub mod extend;
pub mod seed;

pub use chain::{Chain, ChainOptions, chain_seeds};
pub use extend::{Extension, Scoring, extend};
pub use seed::{Seed, SeedOptions, collect_seeds};
