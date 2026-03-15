// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;
use super::mappable_flat_set::HasRegion;

pub struct MappableFlatMultiSet<T: HasRegion + Ord> {
    elements: Vec<T>,
}

impl<T: HasRegion + Ord + Clone> MappableFlatMultiSet<T> {
    pub fn new() -> Self { MappableFlatMultiSet { elements: Vec::new() } }

    pub fn insert(&mut self, element: T) {
        let pos = self.elements.partition_point(|e| e < &element);
        self.elements.insert(pos, element);
    }

    pub fn len(&self) -> usize { self.elements.len() }
    pub fn is_empty(&self) -> bool { self.elements.is_empty() }
    pub fn iter(&self) -> impl Iterator<Item = &T> { self.elements.iter() }

    pub fn overlap_range<'a, 'b>(&'a self, region: &'b GenomicRegion) -> impl Iterator<Item = &'a T> + 'b
    where 'a: 'b
    {
        self.elements.iter().filter(move |e| e.region().overlaps(region))
    }
}

impl<T: HasRegion + Ord + Clone> Default for MappableFlatMultiSet<T> {
    fn default() -> Self { Self::new() }
}
