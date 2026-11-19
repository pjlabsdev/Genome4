//! BAM-Datensätze verarbeiten: sortieren, Duplikate markieren (später Index).

pub mod header;
pub mod markdup;
pub mod sort;

pub use header::{sorted_header, with_program_line};
pub use markdup::{DupMarker, DupOptions, DupStats};
pub use sort::{SortOptions, SortStats, Sorter, sort_key};
