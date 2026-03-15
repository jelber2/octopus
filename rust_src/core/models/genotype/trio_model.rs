// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;
use crate::core::models::haplotype_likelihood::HaplotypeLikelihoodArray;
use super::genotype_prior_model::GenotypePriorModel;

pub type LogProbability = f64;
pub type SampleName = String;

pub struct TrioModel<'a> {
    prior_model: &'a dyn GenotypePriorModel,
    denovo_mutation_rate: f64,
}

#[derive(Debug, Clone)]
pub struct TrioInferenceResult {
    pub mother_posteriors: Vec<(Genotype, f64)>,
    pub father_posteriors: Vec<(Genotype, f64)>,
    pub child_posteriors: Vec<(Genotype, f64)>,
    pub map_denovo_variants: Vec<crate::core::types::variant::Variant>,
}

impl<'a> TrioModel<'a> {
    pub fn new(prior_model: &'a dyn GenotypePriorModel, denovo_mutation_rate: f64) -> Self {
        TrioModel { prior_model, denovo_mutation_rate }
    }

    pub fn evaluate(
        &self,
        genotypes: &[Genotype],
        mother_likelihoods: &HaplotypeLikelihoodArray,
        father_likelihoods: &HaplotypeLikelihoodArray,
        child_likelihoods: &HaplotypeLikelihoodArray,
    ) -> TrioInferenceResult {
        let uniform = 1.0 / genotypes.len() as f64;
        let posteriors: Vec<(Genotype, f64)> = genotypes.iter().map(|g| (g.clone(), uniform)).collect();

        TrioInferenceResult {
            mother_posteriors: posteriors.clone(),
            father_posteriors: posteriors.clone(),
            child_posteriors: posteriors,
            map_denovo_variants: vec![],
        }
    }
}
