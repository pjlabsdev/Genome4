//! `gc-index`: Suffix-Array und FM-Index über das Referenzgenom (Meilenstein M2).
//!
//! Der Index deckt wie bei BWA beide Stränge ab: der Text ist der Vorwärtsstrang aller Contigs,
//! gefolgt von seinem reversen Komplement.

pub mod blockwise;
pub mod fm;
pub mod reference;
pub mod sais;
pub mod smem;

pub use fm::{BuildError, CheckReport, FmIndex, Hit, Variant};
pub use reference::{Contig, Hole, Reference};
pub use smem::{BiInterval, Smem};
