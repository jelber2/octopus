// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;
use super::mappable_flat_set::HasRegion;

pub struct MappableBlock<T: HasRegion> {
    elements: Vec<T>,
    region: GenomicRegion,
}

impl<T: HasRegion + Clone> MappableBlock<T> {
    pub fn new(region: GenomicRegion) -> Self {
        MappableBlock { elements: Vec::new(), region }
    }

    pub fn push(&mut self, element: T) {
        self.elements.push(element);
    }

    pub fn region(&self) -> &GenomicRegion { &self.region }
    pub fn len(&self) -> usize { self.elements.len() }
    pub fn is_empty(&self) -> bool { self.elements.is_empty() }
    pub fn iter(&self) -> impl Iterator<Item = &T> { self.elements.iter() }
}
