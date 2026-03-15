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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basics::genomic_region::GenomicRegion;
    use crate::basics::cigar_string::parse_cigar;

    fn make_haplotype(seq: &[u8]) -> Haplotype {
        let region = GenomicRegion::new("chr1", 0, seq.len() as u32).unwrap();
        Haplotype::new(region, seq.to_vec())
    }

    fn make_read(seq: &[u8]) -> AlignedRead {
        let region = GenomicRegion::new("chr1", 0, seq.len() as u32).unwrap();
        let cigar = parse_cigar(&format!("{}M", seq.len())).unwrap();
        let quals = vec![30u8; seq.len()];
        AlignedRead::new("read1".to_string(), region, seq.to_vec(), quals, cigar, 60, 0)
    }

    #[test]
    fn flat_model_perfect_match_gives_base_likelihood() {
        let haplotype = make_haplotype(b"ACGTACGT");
        let read = make_read(b"ACGTACGT");
        let model = FlatHaplotypeLikelihoodModel::new(-1.0);
        let ll = model.log_likelihood(&read, &haplotype);
        assert_eq!(ll, -1.0);
    }

    #[test]
    fn flat_model_one_mismatch_penalises() {
        let haplotype = make_haplotype(b"ACGTACGT");
        let read = make_read(b"ACGTACGA");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let ll = model.log_likelihood(&read, &haplotype);
        assert_eq!(ll, -3.0);
    }

    #[test]
    fn flat_model_two_mismatches() {
        let haplotype = make_haplotype(b"ACGTACGT");
        let read = make_read(b"ACGTTCGA");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let ll = model.log_likelihood(&read, &haplotype);
        assert_eq!(ll, -6.0);
    }

    #[test]
    fn flat_model_all_mismatches() {
        let haplotype = make_haplotype(b"AAAA");
        let read = make_read(b"TTTT");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let ll = model.log_likelihood(&read, &haplotype);
        assert_eq!(ll, -12.0);
    }

    #[test]
    fn array_new_is_empty() {
        let arr = HaplotypeLikelihoodArray::new();
        assert_eq!(arr.num_haplotypes(), 0);
        assert_eq!(arr.num_reads(), 0);
    }

    #[test]
    fn array_populate_fills_likelihoods() {
        let h1 = make_haplotype(b"ACGTACGT");
        let h2 = make_haplotype(b"ACGTACGA");
        let r1 = make_read(b"ACGTACGT");
        let r2 = make_read(b"ACGTACGA");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h1.clone(), h2.clone()], &[r1, r2], &model);

        assert_eq!(arr.num_haplotypes(), 2);
        assert_eq!(arr.num_reads(), 2);

        let h1_r1 = arr.log_likelihood(&h1, 0).unwrap();
        let h1_r2 = arr.log_likelihood(&h1, 1).unwrap();
        let h2_r1 = arr.log_likelihood(&h2, 0).unwrap();
        let h2_r2 = arr.log_likelihood(&h2, 1).unwrap();

        assert_eq!(h1_r1, 0.0);
        assert_eq!(h1_r2, -3.0);
        assert_eq!(h2_r1, -3.0);
        assert_eq!(h2_r2, 0.0);
    }

    #[test]
    fn array_log_likelihoods_slice() {
        let h = make_haplotype(b"ACGT");
        let r1 = make_read(b"ACGT");
        let r2 = make_read(b"ACGA");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h.clone()], &[r1, r2], &model);
        let lls = arr.log_likelihoods(&h).unwrap();
        assert_eq!(lls.len(), 2);
        assert_eq!(lls[0], 0.0);
        assert_eq!(lls[1], -3.0);
    }

    #[test]
    fn array_unknown_haplotype_returns_none() {
        let h = make_haplotype(b"ACGT");
        let other = make_haplotype(b"TTTT");
        let r = make_read(b"ACGT");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h], &[r], &model);
        assert!(arr.log_likelihood(&other, 0).is_none());
        assert!(arr.log_likelihoods(&other).is_none());
    }

    #[test]
    fn array_out_of_range_read_index_returns_none() {
        let h = make_haplotype(b"ACGT");
        let r = make_read(b"ACGT");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h.clone()], &[r], &model);
        assert!(arr.log_likelihood(&h, 99).is_none());
    }

    #[test]
    fn array_populate_clears_previous_data() {
        let h1 = make_haplotype(b"ACGT");
        let h2 = make_haplotype(b"TTTT");
        let r = make_read(b"ACGT");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h1.clone()], &[r.clone()], &model);
        assert_eq!(arr.num_haplotypes(), 1);
        arr.populate(&[h2.clone()], &[r.clone()], &model);
        assert_eq!(arr.num_haplotypes(), 1);
        assert!(arr.log_likelihood(&h1, 0).is_none());
        assert!(arr.log_likelihood(&h2, 0).is_some());
    }
}
