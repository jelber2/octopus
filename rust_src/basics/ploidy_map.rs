// Converted from C++ to Rust

use std::collections::HashMap;
use super::genomic_region::GenomicRegion;

pub type Ploidy = u8;

#[derive(Debug, Clone, Default)]
pub struct PloidyMap {
    map: HashMap<String, Ploidy>,
    default_ploidy: Ploidy,
}

impl PloidyMap {
    pub fn new(default_ploidy: Ploidy) -> Self {
        PloidyMap { map: HashMap::new(), default_ploidy }
    }

    pub fn set(&mut self, contig: String, ploidy: Ploidy) {
        self.map.insert(contig, ploidy);
    }

    pub fn get(&self, region: &GenomicRegion) -> Ploidy {
        *self.map.get(region.contig_name()).unwrap_or(&self.default_ploidy)
    }

    pub fn get_by_contig(&self, contig: &str) -> Ploidy {
        *self.map.get(contig).unwrap_or(&self.default_ploidy)
    }
}
