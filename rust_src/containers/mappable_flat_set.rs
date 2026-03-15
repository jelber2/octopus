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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basics::contig_region::ContigRegion;

    /// Minimal test type: a ContigRegion tagged with a fixed contig name so
    /// we can implement HasRegion (which returns GenomicRegion).
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct R {
        region: GenomicRegion,
    }

    impl R {
        fn new(begin: u32, end: u32) -> Self {
            R {
                region: GenomicRegion::new("1", begin, end).unwrap(),
            }
        }
    }

    impl HasRegion for R {
        fn region(&self) -> &GenomicRegion { &self.region }
    }

    impl PartialOrd for R {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for R {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            let a = ContigRegion::new(self.region.begin(), self.region.end()).unwrap();
            let b = ContigRegion::new(other.region.begin(), other.region.end()).unwrap();
            a.cmp(&b)
        }
    }

    // ── insert / contains / len / is_empty ────────────────────────────────

    #[test]
    fn insert_works_and_rejects_duplicates() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();

        assert!(set.insert(R::new(0, 0)));
        assert_eq!(set.len(), 1);
        assert!(!set.insert(R::new(0, 0)));   // duplicate
        assert_eq!(set.len(), 1);

        assert!(set.insert(R::new(0, 1)));
        assert_eq!(set.len(), 2);
        assert!(!set.insert(R::new(0, 1)));   // duplicate
        assert_eq!(set.len(), 2);

        assert!(set.insert(R::new(0, 3)));
        assert!(set.insert(R::new(1, 1)));
        assert!(set.insert(R::new(2, 4)));
        assert!(set.insert(R::new(4, 5)));
        assert_eq!(set.len(), 6);
    }

    #[test]
    fn elements_are_kept_sorted_after_inserts() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        for (b, e) in [(3u32, 5u32), (0, 2), (1, 4), (0, 0)] {
            set.insert(R::new(b, e));
        }
        let items: Vec<_> = set.iter().collect();
        for w in items.windows(2) {
            assert!(w[0] <= w[1], "not sorted: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn is_empty_and_len() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        assert!(set.is_empty());
        set.insert(R::new(0, 1));
        assert!(!set.is_empty());
        assert_eq!(set.len(), 1);
    }

    // ── contains ─────────────────────────────────────────────────────────

    #[test]
    fn contains_finds_inserted_element() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        let r = R::new(5, 10);
        assert!(!set.contains(&r));
        set.insert(r.clone());
        assert!(set.contains(&r));
        assert!(!set.contains(&R::new(5, 11)));
        assert!(!set.contains(&R::new(4, 10)));
    }

    // ── remove ───────────────────────────────────────────────────────────

    #[test]
    fn remove_erases_existing_element() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        set.insert(R::new(0, 1));
        set.insert(R::new(1, 2));
        set.insert(R::new(2, 3));
        assert_eq!(set.len(), 3);

        assert!(set.remove(&R::new(1, 2)));
        assert_eq!(set.len(), 2);
        assert!(!set.contains(&R::new(1, 2)));

        assert!(!set.remove(&R::new(1, 2))); // already removed
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn remove_absent_element_returns_false() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        set.insert(R::new(0, 5));
        assert!(!set.remove(&R::new(1, 5)));
        assert_eq!(set.len(), 1);
    }

    // ── overlap_range ─────────────────────────────────────────────────────

    #[test]
    fn overlap_range_finds_overlapping_elements() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        set.insert(R::new(0, 10));
        set.insert(R::new(10, 20));
        set.insert(R::new(20, 30));
        set.insert(R::new(50, 60));

        let query = GenomicRegion::new("1", 8, 22).unwrap();
        let hits: Vec<_> = set.overlap_range(&query).collect();
        assert_eq!(hits.len(), 3); // [0,10), [10,20), [20,30) all touch [8,22)
    }

    #[test]
    fn overlap_range_empty_result_when_no_overlap() {
        let mut set: MappableFlatSet<R> = MappableFlatSet::new();
        set.insert(R::new(0, 10));
        set.insert(R::new(20, 30));

        let query = GenomicRegion::new("1", 11, 19).unwrap();
        let hits: Vec<_> = set.overlap_range(&query).collect();
        assert!(hits.is_empty());
    }

    // ── default ───────────────────────────────────────────────────────────

    #[test]
    fn default_creates_empty_set() {
        let set: MappableFlatSet<R> = MappableFlatSet::default();
        assert!(set.is_empty());
    }
}
