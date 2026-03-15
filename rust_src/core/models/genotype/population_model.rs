// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;
use crate::core::models::haplotype_likelihood::HaplotypeLikelihoodArray;
use super::genotype_prior_model::GenotypePriorModel;

pub type SampleName = String;
pub type LogProbability = f64;

pub struct PopulationModel<'a> {
    prior_model: &'a dyn GenotypePriorModel,
}

#[derive(Debug, Clone)]
pub struct PopulationInferenceResult {
    pub sample_genotype_posteriors: Vec<(SampleName, Vec<(Genotype, f64)>)>,
}

impl<'a> PopulationModel<'a> {
    pub fn new(prior_model: &'a dyn GenotypePriorModel) -> Self {
        PopulationModel { prior_model }
    }

    pub fn evaluate(
        &self,
        samples: &[SampleName],
        genotypes: &[Genotype],
        sample_likelihoods: &[(SampleName, HaplotypeLikelihoodArray)],
    ) -> PopulationInferenceResult {
        let sample_genotype_posteriors = samples.iter()
            .filter_map(|s| {
                let likelihood = sample_likelihoods.iter().find(|(name, _)| name == s)?;
                let posteriors = genotypes.iter()
                    .map(|g| {
                        let lp = self.prior_model.log_prior(g);
                        (g.clone(), lp.exp())
                    })
                    .collect();
                Some((s.clone(), posteriors))
            })
            .collect();

        PopulationInferenceResult { sample_genotype_posteriors }
    }
}
