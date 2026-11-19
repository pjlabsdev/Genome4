//! BAM-Datensätze verarbeiten: sortieren, Duplikate markieren, Index bauen.

pub mod flagstat;
pub mod header;
pub mod index;
pub mod markdup;
pub mod sort;

pub use flagstat::{Counts, FlagStats};
pub use header::{sorted_header, with_coordinate_order, with_program_line};
pub use index::{Bai, BaiBuilder};
pub use markdup::{DupMarker, DupOptions, DupStats};
pub use sort::{SortOptions, SortStats, Sorter, sort_key};
