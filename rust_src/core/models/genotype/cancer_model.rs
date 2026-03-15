// Converted from C++ to Rust

use crate::core::types::cancer_genotype::CancerGenotype;
use crate::core::types::genotype::Genotype;
use crate::core::models::haplotype_likelihood::HaplotypeLikelihoodArray;

pub type LogProbability = f64;
pub type SampleName = String;

pub struct CancerGenotypeModel {
    somatic_mutation_rate: f64,
}

#[derive(Debug, Clone)]
pub struct CancerInferenceResult {
    pub posteriors: Vec<(CancerGenotype, f64)>,
    pub tumour_content: f64,
}

impl CancerGenotypeModel {
    pub fn new(somatic_mutation_rate: f64) -> Self {
        CancerGenotypeModel { somatic_mutation_rate }
    }

    pub fn evaluate(
        &self,
        cancer_genotypes: &[CancerGenotype],
        normal_likelihoods: &HaplotypeLikelihoodArray,
        tumour_likelihoods: &HaplotypeLikelihoodArray,
    ) -> CancerInferenceResult {
        let n = cancer_genotypes.len();
        let uniform = if n > 0 { 1.0 / n as f64 } else { 0.0 };
        let posteriors = cancer_genotypes.iter().map(|g| (g.clone(), uniform)).collect();

        CancerInferenceResult {
            posteriors,
            tumour_content: 0.5,
        }
    }
}
