//! BAM-Datensätze verarbeiten, ohne sie in Text zu wandeln.

pub mod header;
pub mod sort;

pub use header::sorted_header;
pub use sort::{SortOptions, SortStats, Sorter, sort_key};
