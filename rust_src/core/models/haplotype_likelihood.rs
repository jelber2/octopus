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

// ── FlatHaplotypeLikelihoodModel ──────────────────────────────────────────────

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

// ── BaseQualityLikelihoodModel ────────────────────────────────────────────────

/// Convert a Phred-scaled quality score to an error probability.
pub fn phred_to_prob(phred: u8) -> f64 {
    10_f64.powf(-(phred as f64) / 10.0)
}

/// Likelihood model that uses per-base quality scores.
///
/// For each position where the read overlaps the haplotype:
///   - match    → log(1 − ε)     where ε = 10^(−Q/10)
///   - mismatch → log(ε / 3)
///
/// Positions outside the overlap contribute 0.0 (neutral).
/// This correctly handles SNVs; for indels the comparison is positionally
/// offset past the indel, which reduces sensitivity but still distinguishes
/// ref vs alt haplotypes in most cases.
pub struct BaseQualityLikelihoodModel {
    min_quality: u8,
}

impl BaseQualityLikelihoodModel {
    pub fn new() -> Self {
        BaseQualityLikelihoodModel { min_quality: 1 }
    }
}

impl Default for BaseQualityLikelihoodModel {
    fn default() -> Self { Self::new() }
}

impl HaplotypeLikelihoodModel for BaseQualityLikelihoodModel {
    fn log_likelihood(&self, read: &AlignedRead, haplotype: &Haplotype) -> LogProbability {
        let hap_region  = haplotype.mapped_region();
        let read_region = read.mapped_region();

        if read_region.contig_name() != hap_region.contig_name() {
            return 0.0;
        }

        let hap_begin  = hap_region.begin();
        let hap_end    = hap_region.end();
        let read_begin = read_region.begin();
        let read_end   = read_region.end();

        let overlap_begin = hap_begin.max(read_begin);
        let overlap_end   = hap_end.min(read_end);
        if overlap_begin >= overlap_end {
            return 0.0;
        }

        let hap_seq  = haplotype.sequence();
        let read_seq = read.sequence();
        let quals    = read.qualities();

        let hap_off  = (overlap_begin - hap_begin)  as usize;
        let read_off = (overlap_begin - read_begin) as usize;
        let overlap_len = (overlap_end - overlap_begin) as usize;

        let usable = overlap_len
            .min(hap_seq.len().saturating_sub(hap_off))
            .min(read_seq.len().saturating_sub(read_off));

        let mut log_prob = 0.0_f64;
        for i in 0..usable {
            let hap_base  = hap_seq[hap_off + i];
            let read_base = read_seq[read_off + i];
            let raw_qual  = quals.get(read_off + i).copied().unwrap_or(20);
            let qual      = raw_qual.max(self.min_quality);

            let error_prob = phred_to_prob(qual).max(1e-30);

            if hap_base == read_base {
                log_prob += (1.0 - error_prob).ln();
            } else {
                log_prob += (error_prob / 3.0).ln();
            }
        }

        log_prob
    }
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
        let quals  = vec![30u8; seq.len()];
        let cigar  = parse_cigar(&format!("{}M", seq.len())).unwrap();
        AlignedRead::new(
            "r".to_string(), region, seq.to_vec(), quals, cigar, 60, 0x1,
        )
    }

    // ── FlatHaplotypeLikelihoodModel tests ────────────────────────────────────

    #[test]
    fn flat_perfect_match_returns_base() {
        let h = make_haplotype(b"ACGT");
        let r = make_read(b"ACGT");
        let m = FlatHaplotypeLikelihoodModel::new(0.0);
        assert_eq!(m.log_likelihood(&r, &h), 0.0);
    }

    #[test]
    fn flat_one_mismatch_penalty() {
        let h = make_haplotype(b"ACGT");
        let r = make_read(b"ACGA");
        let m = FlatHaplotypeLikelihoodModel::new(0.0);
        assert_eq!(m.log_likelihood(&r, &h), -3.0);
    }

    #[test]
    fn flat_base_likelihood_is_added() {
        let h = make_haplotype(b"ACGT");
        let r = make_read(b"ACGT");
        let m = FlatHaplotypeLikelihoodModel::new(-2.5);
        assert_eq!(m.log_likelihood(&r, &h), -2.5);
    }

    // ── BaseQualityLikelihoodModel tests ──────────────────────────────────────

    #[test]
    fn bq_perfect_match_positive_ll() {
        let h = make_haplotype(b"ACGT");
        let r = make_read(b"ACGT");
        let m = BaseQualityLikelihoodModel::new();
        let ll = m.log_likelihood(&r, &h);
        assert!(ll > -1.0 && ll <= 0.0, "expected small negative, got {}", ll);
    }

    #[test]
    fn bq_mismatch_lower_than_match() {
        let h  = make_haplotype(b"ACGT");
        let r_match    = make_read(b"ACGT");
        let r_mismatch = make_read(b"ACGA");
        let m = BaseQualityLikelihoodModel::new();
        assert!(m.log_likelihood(&r_match, &h) > m.log_likelihood(&r_mismatch, &h));
    }

    #[test]
    fn bq_different_contigs_returns_zero() {
        let hap_region  = GenomicRegion::new("chr1", 0, 4).unwrap();
        let read_region = GenomicRegion::new("chr2", 0, 4).unwrap();
        let h = Haplotype::new(hap_region, b"ACGT".to_vec());
        let cigar = parse_cigar("4M").unwrap();
        let r = AlignedRead::new("r".to_string(), read_region, b"ACGT".to_vec(), vec![30; 4], cigar, 60, 0x1);
        let m = BaseQualityLikelihoodModel::new();
        assert_eq!(m.log_likelihood(&r, &h), 0.0);
    }

    #[test]
    fn bq_no_overlap_returns_zero() {
        let hap_region  = GenomicRegion::new("chr1", 0, 4).unwrap();
        let read_region = GenomicRegion::new("chr1", 100, 104).unwrap();
        let h = Haplotype::new(hap_region, b"ACGT".to_vec());
        let cigar = parse_cigar("4M").unwrap();
        let r = AlignedRead::new("r".to_string(), read_region, b"ACGT".to_vec(), vec![30; 4], cigar, 60, 0x1);
        let m = BaseQualityLikelihoodModel::new();
        assert_eq!(m.log_likelihood(&r, &h), 0.0);
    }

    #[test]
    fn bq_partial_overlap_uses_only_overlap() {
        // Haplotype: chr1:0-8, read: chr1:4-8 (overlap = 4..8)
        let hap_region  = GenomicRegion::new("chr1", 0, 8).unwrap();
        let read_region = GenomicRegion::new("chr1", 4, 8).unwrap();
        let h = Haplotype::new(hap_region, b"ACGTTTTT".to_vec());
        let cigar = parse_cigar("4M").unwrap();
        let r = AlignedRead::new("r".to_string(), read_region, b"TTTT".to_vec(), vec![30; 4], cigar, 60, 0x1);
        let m = BaseQualityLikelihoodModel::new();
        let ll_match = m.log_likelihood(&r, &h);
        // All 4 bases in the overlap match (hap[4..8]=TTTT == read=TTTT)
        assert!(ll_match > -0.5, "all bases should match; got {}", ll_match);
    }

    // ── HaplotypeLikelihoodArray tests ────────────────────────────────────────

    #[test]
    fn array_log_likelihoods_slice() {
        let h  = make_haplotype(b"ACGT");
        let r1 = make_read(b"ACGT");
        let r2 = make_read(b"ACGA");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h.clone()], &[r1, r2], &model);
        let lls = arr.log_likelihoods(&h).unwrap();
        assert_eq!(lls.len(), 2);
        assert_eq!(lls[0],  0.0);
        assert_eq!(lls[1], -3.0);
    }

    #[test]
    fn array_unknown_haplotype_returns_none() {
        let h     = make_haplotype(b"ACGT");
        let other = make_haplotype(b"TTTT");
        let r     = make_read(b"ACGT");
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
        let r  = make_read(b"ACGT");
        let model = FlatHaplotypeLikelihoodModel::new(0.0);
        let mut arr = HaplotypeLikelihoodArray::new();
        arr.populate(&[h1.clone()], &[r.clone()], &model);
        assert_eq!(arr.num_haplotypes(), 1);
        arr.populate(&[h2.clone()], &[r.clone()], &model);
        assert_eq!(arr.num_haplotypes(), 1);
        assert!(arr.log_likelihood(&h1, 0).is_none());
        assert!(arr.log_likelihood(&h2, 0).is_some());
    }

    // ── phred_to_prob ─────────────────────────────────────────────────────────

    #[test]
    fn phred_to_prob_q30() {
        let p = phred_to_prob(30);
        assert!((p - 0.001).abs() < 1e-10, "Q30 should be 0.001, got {}", p);
    }

    #[test]
    fn phred_to_prob_q20() {
        let p = phred_to_prob(20);
        assert!((p - 0.01).abs() < 1e-10, "Q20 should be 0.01, got {}", p);
    }
}
