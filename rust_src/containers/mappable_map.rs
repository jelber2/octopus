// Converted from C++ to Rust

use std::collections::BTreeMap;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::contig_region::ContigRegion;

pub struct MappableMap<V> {
    map: BTreeMap<String, BTreeMap<ContigRegion, V>>,
}

impl<V: Clone> MappableMap<V> {
    pub fn new() -> Self {
        MappableMap { map: BTreeMap::new() }
    }

    pub fn insert(&mut self, region: GenomicRegion, value: V) {
        self.map.entry(region.contig_name().to_string())
            .or_default()
            .insert(*region.contig_region(), value);
    }

    pub fn get(&self, region: &GenomicRegion) -> Option<&V> {
        self.map.get(region.contig_name())?
            .get(region.contig_region())
    }

    pub fn contains(&self, region: &GenomicRegion) -> bool {
        self.map.get(region.contig_name())
            .map(|m| m.contains_key(region.contig_region()))
            .unwrap_or(false)
    }

    pub fn len(&self) -> usize {
        self.map.values().map(|m| m.len()).sum()
    }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn overlap_range(&self, region: &GenomicRegion) -> Vec<(&ContigRegion, &V)> {
        let mut result = Vec::new();
        if let Some(contig_map) = self.map.get(region.contig_name()) {
            for (cr, v) in contig_map {
                if cr.end() > region.begin() && cr.begin() < region.end() {
                    result.push((cr, v));
                }
            }
        }
        result
    }
}

impl<V: Clone> Default for MappableMap<V> {
    fn default() -> Self { Self::new() }
}
