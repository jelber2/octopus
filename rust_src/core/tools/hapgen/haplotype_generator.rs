// Converted from C++ to Rust

use crate::core::types::haplotype::Haplotype;
use crate::core::types::variant::Variant;
use crate::basics::genomic_region::GenomicRegion;
use crate::io::reference::reference_genome::ReferenceGenome;
use super::haplotype_tree::HaplotypeTree;

pub struct HaplotypeGenerator<'a> {
    reference: &'a ReferenceGenome,
    max_haplotypes: usize,
}

pub struct HaplotypePacket {
    pub haplotypes: Vec<Haplotype>,
    pub region: GenomicRegion,
}

impl<'a> HaplotypeGenerator<'a> {
    pub fn new(reference: &'a ReferenceGenome, max_haplotypes: usize) -> Self {
        HaplotypeGenerator { reference, max_haplotypes }
    }

    pub fn generate(&self, variants: &[Variant], region: &GenomicRegion) -> HaplotypePacket {
        let mut haplotypes = Vec::new();

        if let Ok(ref_seq) = self.reference.fetch_sequence(region) {
            let ref_haplotype = Haplotype::new(region.clone(), ref_seq);
            haplotypes.push(ref_haplotype);
        }

        for variant in variants {
            if haplotypes.len() >= self.max_haplotypes { break; }
            if let Ok(alt_seq) = self.build_haplotype_sequence(variant, region) {
                let haplotype = Haplotype::new(region.clone(), alt_seq);
                haplotypes.push(haplotype);
            }
        }

        HaplotypePacket { haplotypes, region: region.clone() }
    }

    fn build_haplotype_sequence(&self, variant: &Variant, region: &GenomicRegion) -> Result<Vec<u8>, String> {
        let mut ref_seq = self.reference.fetch_sequence(region)?;

        let var_begin = variant.mapped_region().begin();
        let var_end = variant.mapped_region().end();
        let region_begin = region.begin();

        if var_begin >= region_begin && var_end <= region.end() {
            let start = (var_begin - region_begin) as usize;
            let end = (var_end - region_begin) as usize;
            let alt = variant.alt_sequence();

            ref_seq.splice(start..end, alt.iter().cloned());
        }

        Ok(ref_seq)
    }
}
