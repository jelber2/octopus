// Converted from C++ to Rust
// Haplotype likelihood array — stores log-likelihoods P(read | haplotype)

use std::collections::HashMap;
use crate::core::types::haplotype::Haplotype;
use crate::basics::aligned_read::AlignedRead;

pub type Probability = f64;
pub type LogProbability = f64;

#[derive(Debug, Clone)]
pub struct HaplotypeLikelihoodArray {
    haplotype_indices: HashMap<u64, usize>,
    likelihoods: Vec<Vec<LogProbability>>,
    num_reads: usize,
}

impl HaplotypeLikelihoodArray {
    pub fn new() -> Self {
        HaplotypeLikelihoodArray {
            haplotype_indices: HashMap::new(),
            likelihoods: Vec::new(),
            num_reads: 0,
        }
    }

    pub fn populate(
        &mut self,
        haplotypes: &[Haplotype],
        reads: &[AlignedRead],
        model: &dyn HaplotypeLikelihoodModel,
    ) {
        self.num_reads = reads.len();
        self.haplotype_indices.clear();
        self.likelihoods.clear();

        for (i, haplotype) in haplotypes.iter().enumerate() {
            self.haplotype_indices.insert(haplotype.get_hash(), i);
            let lls: Vec<LogProbability> = reads.iter()
                .map(|read| model.log_likelihood(read, haplotype))
                .collect();
            self.likelihoods.push(lls);
        }
    }

    pub fn num_haplotypes(&self) -> usize { self.likelihoods.len() }
    pub fn num_reads(&self) -> usize { self.num_reads }

    pub fn log_likelihood(&self, haplotype: &Haplotype, read_index: usize) -> Option<LogProbability> {
        let hap_idx = *self.haplotype_indices.get(&haplotype.get_hash())?;
        self.likelihoods.get(hap_idx)?.get(read_index).copied()
    }

    pub fn log_likelihoods(&self, haplotype: &Haplotype) -> Option<&[LogProbability]> {
        let hap_idx = *self.haplotype_indices.get(&haplotype.get_hash())?;
        self.likelihoods.get(hap_idx).map(|v| v.as_slice())
    }
}

impl Default for HaplotypeLikelihoodArray {
    fn default() -> Self { Self::new() }
}

pub trait HaplotypeLikelihoodModel: Send + Sync {
    fn log_likelihood(&self, read: &AlignedRead, haplotype: &Haplotype) -> LogProbability;
}

pub struct FlatHaplotypeLikelihoodModel {
    base_likelihood: LogProbability,
}

impl FlatHaplotypeLikelihoodModel {
    pub fn new(base_likelihood: LogProbability) -> Self {
        FlatHaplotypeLikelihoodModel { base_likelihood }
    }
}

impl HaplotypeLikelihoodModel for FlatHaplotypeLikelihoodModel {
    fn log_likelihood(&self, read: &AlignedRead, haplotype: &Haplotype) -> LogProbability {
        let mismatches = count_mismatches(read.sequence(), haplotype.sequence());
        let mismatch_penalty = -3.0 * mismatches as f64;
        self.base_likelihood + mismatch_penalty
    }
}

fn count_mismatches(seq1: &[u8], seq2: &[u8]) -> usize {
    seq1.iter().zip(seq2.iter())
        .filter(|(a, b)| a != b)
        .count()
}
