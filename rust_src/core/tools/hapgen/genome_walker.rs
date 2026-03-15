// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;
use crate::io::reference::reference_genome::ReferenceGenome;
use crate::core::types::variant::Variant;

pub struct GenomeWalker {
    max_region_size: u32,
    max_variants: usize,
}

impl GenomeWalker {
    pub fn new(max_region_size: u32, max_variants: usize) -> Self {
        GenomeWalker { max_region_size, max_variants }
    }

    pub fn walk<'a>(&self, variants: &'a [Variant], contig: &str, contig_size: u32) -> Vec<(&'a [Variant], GenomicRegion)> {
        let mut result = Vec::new();
        let mut i = 0;

        while i < variants.len() {
            let start = variants[i].mapped_region().begin();
            let mut j = i;
            let mut end = start;

            while j < variants.len() && j - i < self.max_variants {
                let v_end = variants[j].mapped_region().end();
                if v_end - start > self.max_region_size { break; }
                end = v_end;
                j += 1;
            }

            if let Ok(region) = GenomicRegion::new(contig, start, end.max(start + 1)) {
                result.push((&variants[i..j], region));
            }
            i = j;
        }

        result
    }
}
