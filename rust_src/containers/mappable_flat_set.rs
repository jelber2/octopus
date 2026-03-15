// Converted from C++ to Rust
// Sorted flat set of Mappable types (sorted by GenomicRegion)

use crate::basics::genomic_region::GenomicRegion;

pub struct MappableFlatSet<T: HasRegion + Ord> {
    elements: Vec<T>,
}

pub trait HasRegion {
    fn region(&self) -> &GenomicRegion;
}

impl<T: HasRegion + Ord + Clone> MappableFlatSet<T> {
    pub fn new() -> Self { MappableFlatSet { elements: Vec::new() } }

    pub fn insert(&mut self, element: T) -> bool {
        let pos = self.elements.partition_point(|e| e < &element);
        if pos < self.elements.len() && self.elements[pos] == element {
            return false;
        }
        self.elements.insert(pos, element);
        true
    }

    pub fn remove(&mut self, element: &T) -> bool {
        if let Some(pos) = self.elements.iter().position(|e| e == element) {
            self.elements.remove(pos);
            true
        } else {
            false
        }
    }

    pub fn contains(&self, element: &T) -> bool {
        self.elements.binary_search(element).is_ok()
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

impl<T: HasRegion + Ord + Clone> Default for MappableFlatSet<T> {
    fn default() -> Self { Self::new() }
}
