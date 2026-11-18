//! `gc-align`: Reads ausrichten nach dem Verfahren von BWA-MEM (Meilenstein M3).
//!
//! Bisher gebaut: Seeds sammeln (SMEM mit Nachsuche), im Genom orten und zu Ketten verbinden.

pub mod chain;
pub mod seed;

pub use chain::{Chain, ChainOptions, chain_seeds};
pub use seed::{Seed, SeedOptions, collect_seeds};
