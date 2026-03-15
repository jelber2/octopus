// Converted from C++ to Rust

use crate::core::types::haplotype::Haplotype;
use crate::basics::genomic_region::GenomicRegion;

pub fn filter_haplotypes(
    haplotypes: Vec<Haplotype>,
    max_haplotypes: usize,
    reference_haplotype: Option<&Haplotype>,
) -> Vec<Haplotype> {
    if haplotypes.len() <= max_haplotypes {
        return haplotypes;
    }

    let mut result: Vec<Haplotype>;

    if let Some(ref_hap) = reference_haplotype {
        result = vec![ref_hap.clone()];
        for h in haplotypes {
            if result.len() >= max_haplotypes { break; }
            if &h != ref_hap {
                result.push(h);
            }
        }
    } else {
        result = haplotypes;
        result.truncate(max_haplotypes);
    }

    result
}
