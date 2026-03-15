// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;
use crate::core::models::haplotype_likelihood::HaplotypeLikelihoodArray;

pub type LogProbability = f64;

pub struct SubcloneModel {
    num_clones: usize,
    concentration: f64,
}

#[derive(Debug, Clone)]
pub struct SubcloneInferenceResult {
    pub genotype_posteriors: Vec<(Genotype, f64)>,
    pub clone_frequencies: Vec<f64>,
}

impl SubcloneModel {
    pub fn new(num_clones: usize, concentration: f64) -> Self {
        SubcloneModel { num_clones, concentration }
    }

    pub fn evaluate(
        &self,
        genotypes: &[Genotype],
        likelihoods: &HaplotypeLikelihoodArray,
    ) -> SubcloneInferenceResult {
        let n = genotypes.len();
        let uniform_genotype = if n > 0 { 1.0 / n as f64 } else { 0.0 };
        let uniform_clone = if self.num_clones > 0 { 1.0 / self.num_clones as f64 } else { 0.0 };

        SubcloneInferenceResult {
            genotype_posteriors: genotypes.iter().map(|g| (g.clone(), uniform_genotype)).collect(),
            clone_frequencies: vec![uniform_clone; self.num_clones],
        }
    }
}
