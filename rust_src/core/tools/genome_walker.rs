// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;
use crate::core::types::variant::Variant;
use crate::io::reference::reference_genome::ReferenceGenome;

pub struct GenomeWalker {
    max_alleles_per_region: usize,
    max_region_size: u32,
}

impl GenomeWalker {
    pub fn new(max_alleles_per_region: usize, max_region_size: u32) -> Self {
        GenomeWalker { max_alleles_per_region, max_region_size }
    }

    pub fn build_walk<'a>(
        &self,
        variants: &'a [Variant],
        reference: &ReferenceGenome,
    ) -> Vec<(GenomicRegion, &'a [Variant])> {
        if variants.is_empty() { return Vec::new(); }

        let mut result = Vec::new();
        let mut i = 0;

        while i < variants.len() {
            let start = variants[i].mapped_region().begin();
            let contig = variants[i].mapped_region().contig_name().to_string();
            let mut end = start;
            let mut j = i;

            while j < variants.len() && j - i < self.max_alleles_per_region {
                let v_end = variants[j].mapped_region().end();
                if v_end - start > self.max_region_size { break; }
                end = v_end;
                j += 1;
            }

            if let Ok(region) = GenomicRegion::new(&contig, start, end.max(start + 1)) {
                result.push((region, &variants[i..j]));
            }
            i = j;
        }

        result
    }
}
