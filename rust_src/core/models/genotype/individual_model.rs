// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;
use crate::core::types::haplotype::Haplotype;
use crate::core::models::haplotype_likelihood::{HaplotypeLikelihoodArray, LogProbability};
use super::genotype_prior_model::GenotypePriorModel;

pub struct IndividualModel<'a> {
    prior_model: &'a dyn GenotypePriorModel,
}

#[derive(Debug, Clone)]
pub struct GenotypeInferenceResult {
    pub posteriors: Vec<(Genotype, LogProbability)>,
    pub map_genotype_index: usize,
}

impl<'a> IndividualModel<'a> {
    pub fn new(prior_model: &'a dyn GenotypePriorModel) -> Self {
        IndividualModel { prior_model }
    }

    pub fn evaluate(
        &self,
        genotypes: &[Genotype],
        likelihoods: &HaplotypeLikelihoodArray,
    ) -> GenotypeInferenceResult {
        let mut log_posteriors: Vec<LogProbability> = genotypes.iter()
            .map(|g| self.log_joint(g, likelihoods))
            .collect();

        let log_norm = log_sum_exp(&log_posteriors);

        let posteriors: Vec<(Genotype, LogProbability)> = genotypes.iter().cloned()
            .zip(log_posteriors.iter().map(|&lp| (lp - log_norm).exp()))
            .collect();

        let map_genotype_index = log_posteriors.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);

        GenotypeInferenceResult { posteriors, map_genotype_index }
    }

    fn log_joint(&self, genotype: &Genotype, likelihoods: &HaplotypeLikelihoodArray) -> LogProbability {
        let log_prior = self.prior_model.log_prior(genotype);
        let log_likelihood = self.log_likelihood(genotype, likelihoods);
        log_prior + log_likelihood
    }

    fn log_likelihood(&self, genotype: &Genotype, likelihoods: &HaplotypeLikelihoodArray) -> LogProbability {
        if likelihoods.num_reads() == 0 { return 0.0; }

        let ploidy = genotype.ploidy();
        if ploidy == 0 { return f64::NEG_INFINITY; }
        let log_ploidy = (ploidy as f64).ln();

        let mut total = 0.0;
        for read_idx in 0..likelihoods.num_reads() {
            let mut read_lls: Vec<LogProbability> = genotype.haplotypes().iter()
                .filter_map(|h| likelihoods.log_likelihood(h, read_idx))
                .collect();

            if read_lls.is_empty() { continue; }

            let log_sum = log_sum_exp(&read_lls) - log_ploidy;
            total += log_sum;
        }
        total
    }
}

fn log_sum_exp(log_probs: &[f64]) -> f64 {
    if log_probs.is_empty() { return f64::NEG_INFINITY; }
    let max = log_probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max.is_infinite() { return max; }
    max + log_probs.iter().map(|&x| (x - max).exp()).sum::<f64>().ln()
}
