// Converted from C++ to Rust

use std::collections::HashMap;
use crate::core::types::haplotype::Haplotype;
use crate::core::types::allele::{Allele, ContigAllele};
use crate::basics::genomic_region::GenomicRegion;
use crate::io::reference::reference_genome::ReferenceGenome;

pub struct HaplotypeTree {
    region: GenomicRegion,
    haplotypes: Vec<Haplotype>,
}

impl HaplotypeTree {
    pub fn new(region: GenomicRegion) -> Self {
        HaplotypeTree { region, haplotypes: Vec::new() }
    }

    pub fn extend(&mut self, alleles: &[Allele], reference: &ReferenceGenome) {
        for allele in alleles {
            if let Ok(seq) = reference.fetch_sequence(&allele.region) {
                let haplotype = Haplotype::new(self.region.clone(), seq);
                self.haplotypes.push(haplotype);
            }
        }
    }

    pub fn haplotypes(&self) -> &[Haplotype] { &self.haplotypes }
    pub fn num_haplotypes(&self) -> usize { self.haplotypes.len() }
    pub fn is_empty(&self) -> bool { self.haplotypes.is_empty() }
}
