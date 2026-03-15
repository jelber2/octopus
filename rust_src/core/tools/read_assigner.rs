// Converted from C++ to Rust

use std::collections::HashMap;
use crate::core::types::haplotype::Haplotype;
use crate::basics::aligned_read::AlignedRead;
use crate::core::models::haplotype_likelihood::HaplotypeLikelihoodArray;

pub type HaplotypeReadAssignment = HashMap<usize, Vec<usize>>;

pub fn assign_reads_to_haplotypes(
    reads: &[AlignedRead],
    haplotypes: &[Haplotype],
    likelihoods: &HaplotypeLikelihoodArray,
) -> HaplotypeReadAssignment {
    let mut assignment: HaplotypeReadAssignment = HashMap::new();

    for (read_idx, _read) in reads.iter().enumerate() {
        let best_haplotype_idx = (0..haplotypes.len())
            .max_by(|&i, &j| {
                let li = likelihoods.log_likelihood(&haplotypes[i], read_idx).unwrap_or(f64::NEG_INFINITY);
                let lj = likelihoods.log_likelihood(&haplotypes[j], read_idx).unwrap_or(f64::NEG_INFINITY);
                li.partial_cmp(&lj).unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some(hap_idx) = best_haplotype_idx {
            assignment.entry(hap_idx).or_default().push(read_idx);
        }
    }

    assignment
}
