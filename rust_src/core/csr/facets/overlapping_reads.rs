// Converted from C++ to Rust

use std::any::Any;
use std::collections::HashMap;
use super::facet::Facet;
use crate::basics::aligned_read::AlignedRead;

pub struct OverlappingReadsFacet {
    sample_reads: HashMap<String, Vec<AlignedRead>>,
}

impl OverlappingReadsFacet {
    pub fn new(sample_reads: HashMap<String, Vec<AlignedRead>>) -> Self {
        OverlappingReadsFacet { sample_reads }
    }
    pub fn get(&self, sample: &str) -> Option<&[AlignedRead]> {
        self.sample_reads.get(sample).map(|v| v.as_slice())
    }
}

impl Facet for OverlappingReadsFacet {
    fn name(&self) -> &str { "OverlappingReads" }
    fn as_any(&self) -> &dyn Any { self }
}
