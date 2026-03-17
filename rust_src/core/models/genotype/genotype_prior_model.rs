// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;

pub type LogProbability = f64;

pub trait GenotypePriorModel: Send + Sync {
    fn log_prior(&self, genotype: &Genotype) -> LogProbability;
}

// ── UniformGenotypePriorModel ─────────────────────────────────────────────────

pub struct UniformGenotypePriorModel {
    num_genotypes: usize,
}

impl UniformGenotypePriorModel {
    pub fn new(num_genotypes: usize) -> Self {
        UniformGenotypePriorModel { num_genotypes }
    }
}

impl GenotypePriorModel for UniformGenotypePriorModel {
    fn log_prior(&self, _genotype: &Genotype) -> LogProbability {
        if self.num_genotypes > 0 {
            -(self.num_genotypes as f64).ln()
        } else {
            f64::NEG_INFINITY
        }
    }
}

// ── HardyWeinbergGenotypePriorModel ───────────────────────────────────────────

/// Hardy-Weinberg equilibrium genotype prior for a single-sample caller.
///
/// For a diploid genotype with `k` non-reference copies out of ploidy `n`:
///   log P(G) = log C(n,k) + k·log(θ) + (n−k)·log(1−θ)
///
/// `reference_hash` identifies the reference haplotype so that non-reference
/// allele counts can be determined by hash comparison alone.
pub struct HardyWeinbergGenotypePriorModel {
    log_snp_theta:      f64,
    log_one_minus_snp:  f64,
    reference_hash:     u64,
    /// Indel priors stored for future per-variant-type dispatch.
    #[allow(dead_code)]
    log_indel_theta:     f64,
    #[allow(dead_code)]
    log_one_minus_indel: f64,
}

impl HardyWeinbergGenotypePriorModel {
    pub fn new(snp_heterozygosity: f64, indel_heterozygosity: f64, reference_hash: u64) -> Self {
        let snp   = snp_heterozygosity.clamp(1e-30, 1.0 - 1e-30);
        let indel = indel_heterozygosity.clamp(1e-30, 1.0 - 1e-30);
        HardyWeinbergGenotypePriorModel {
            log_snp_theta:       snp.ln(),
            log_one_minus_snp:   (1.0 - snp).ln(),
            reference_hash,
            log_indel_theta:     indel.ln(),
            log_one_minus_indel: (1.0 - indel).ln(),
        }
    }
}

impl GenotypePriorModel for HardyWeinbergGenotypePriorModel {
    fn log_prior(&self, genotype: &Genotype) -> LogProbability {
        let ploidy = genotype.ploidy();
        if ploidy == 0 { return f64::NEG_INFINITY; }

        let num_alt = genotype.haplotypes().iter()
            .filter(|h| h.get_hash() != self.reference_hash)
            .count();

        let log_binom = log_binomial_coeff(ploidy, num_alt);
        let k   = num_alt as f64;
        let n_k = (ploidy - num_alt) as f64;

        log_binom
            + k   * self.log_snp_theta
            + n_k * self.log_one_minus_snp
    }
}

/// log C(n, k) computed without overflow.
fn log_binomial_coeff(n: usize, k: usize) -> f64 {
    if k > n { return f64::NEG_INFINITY; }
    if k == 0 || k == n { return 0.0; }
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

fn log_factorial(n: usize) -> f64 {
    (1..=n).map(|i| (i as f64).ln()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_prior_distributes_evenly() {
        use crate::core::types::haplotype::Haplotype;
        use crate::basics::genomic_region::GenomicRegion;
        use crate::core::types::genotype::Genotype;
        let region = GenomicRegion::new("chr1", 0, 10).unwrap();
        let h = Haplotype::new(region, b"ACGTACGTAC".to_vec());
        let g = Genotype::new(vec![h.clone(), h]);
        let model = UniformGenotypePriorModel::new(4);
        let lp = model.log_prior(&g);
        assert!((lp - (-(4_f64).ln())).abs() < 1e-10);
    }

    #[test]
    fn hw_hom_ref_dominates_at_low_het() {
        use crate::core::types::haplotype::Haplotype;
        use crate::basics::genomic_region::GenomicRegion;
        use crate::core::types::genotype::Genotype;
        let region = GenomicRegion::new("chr1", 0, 4).unwrap();
        let ref_h  = Haplotype::new(region.clone(), b"ACGT".to_vec());
        let alt_h  = Haplotype::new(region,          b"ACTT".to_vec());

        let model = HardyWeinbergGenotypePriorModel::new(0.001, 0.0001, ref_h.get_hash());

        let hom_ref = Genotype::new(vec![ref_h.clone(), ref_h.clone()]);
        let het     = Genotype::new(vec![ref_h.clone(), alt_h.clone()]);
        let hom_alt = Genotype::new(vec![alt_h.clone(), alt_h.clone()]);

        let p_hom_ref = model.log_prior(&hom_ref);
        let p_het     = model.log_prior(&het);
        let p_hom_alt = model.log_prior(&hom_alt);

        assert!(p_hom_ref > p_het,     "hom_ref ({}) > het ({}))",     p_hom_ref, p_het);
        assert!(p_het     > p_hom_alt, "het ({}) > hom_alt ({})",       p_het, p_hom_alt);
    }

    #[test]
    fn log_binomial_coeff_symmetry() {
        assert!((log_binomial_coeff(5, 2) - log_binomial_coeff(5, 3)).abs() < 1e-10);
    }

    #[test]
    fn log_binomial_coeff_known_values() {
        // C(4,2) = 6 → log(6)
        assert!((log_binomial_coeff(4, 2) - 6_f64.ln()).abs() < 1e-10);
        assert_eq!(log_binomial_coeff(5, 0), 0.0);
        assert_eq!(log_binomial_coeff(5, 5), 0.0);
    }
}
