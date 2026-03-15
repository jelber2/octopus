// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;

pub trait ReadFilter: Send + Sync {
    fn passes(&self, read: &AlignedRead) -> bool;
    fn name(&self) -> &str;
}

pub struct CompositeReadFilter {
    filters: Vec<Box<dyn ReadFilter>>,
}

impl CompositeReadFilter {
    pub fn new() -> Self { CompositeReadFilter { filters: Vec::new() } }

    pub fn add(&mut self, filter: Box<dyn ReadFilter>) {
        self.filters.push(filter);
    }
}

impl ReadFilter for CompositeReadFilter {
    fn name(&self) -> &str { "Composite" }

    fn passes(&self, read: &AlignedRead) -> bool {
        self.filters.iter().all(|f| f.passes(read))
    }
}

impl Default for CompositeReadFilter {
    fn default() -> Self { Self::new() }
}
