// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;
use crate::core::types::haplotype::Haplotype;
use super::genotype_prior_model::GenotypePriorModel;

pub struct HardyWeinbergModel {
    allele_frequencies: Vec<f64>,
    ploidy: usize,
}

impl HardyWeinbergModel {
    pub fn new(allele_frequencies: Vec<f64>, ploidy: usize) -> Self {
        HardyWeinbergModel { allele_frequencies, ploidy }
    }
}

impl GenotypePriorModel for HardyWeinbergModel {
    fn log_prior(&self, genotype: &Genotype) -> f64 {
        if self.allele_frequencies.is_empty() {
            return -(genotype.ploidy() as f64).ln();
        }

        let n = self.allele_frequencies.len();
        let ploidy = genotype.ploidy();

        let freq_per_haplotype: f64 = genotype.haplotypes().iter()
            .enumerate()
            .map(|(i, _h)| {
                let idx = i % n;
                self.allele_frequencies[idx].max(1e-10).ln()
            })
            .sum();

        let multinomial_coef = log_multinomial_coefficient(genotype);
        freq_per_haplotype + multinomial_coef
    }
}

fn log_multinomial_coefficient(genotype: &Genotype) -> f64 {
    let ploidy = genotype.ploidy();
    if ploidy <= 1 { return 0.0; }

    let mut counts: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    for h in genotype.haplotypes() {
        *counts.entry(h.get_hash()).or_insert(0) += 1;
    }

    let log_ploidy_factorial = log_factorial(ploidy);
    let log_denom: f64 = counts.values().map(|&c| log_factorial(c)).sum();
    log_ploidy_factorial - log_denom
}

fn log_factorial(n: usize) -> f64 {
    (1..=n).map(|i| (i as f64).ln()).sum()
}
