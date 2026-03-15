// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use crate::io::variant::vcf_record::VcfRecord;
use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};
use crate::core::types::genotype::make_genotypes;
use crate::core::types::haplotype::Haplotype;
use crate::core::models::haplotype_likelihood::{HaplotypeLikelihoodArray, FlatHaplotypeLikelihoodModel};
use crate::core::models::genotype::individual_model::IndividualModel;
use crate::core::models::genotype::genotype_prior_model::UniformGenotypePriorModel;

pub struct IndividualCaller {
    options: CallerOptions,
}

impl IndividualCaller {
    pub fn new(options: CallerOptions) -> Self {
        IndividualCaller { options }
    }
}

impl Caller for IndividualCaller {
    fn name(&self) -> &str { "individual" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        let mut calls = CallBuffer::new();

        for (sample, reads) in &env.reads {
            if reads.len() < self.options.min_read_depth {
                continue;
            }

            let ref_sequence = env.reference.fetch_sequence(&env.region)
                .map_err(|e| e.to_string())?;

            let ref_haplotype = Haplotype::new(env.region.clone(), ref_sequence.clone());
            let haplotypes = vec![ref_haplotype];

            let mut likelihood_array = HaplotypeLikelihoodArray::new();
            let model = FlatHaplotypeLikelihoodModel::new(-1.0);
            likelihood_array.populate(&haplotypes, reads, &model);

            let genotypes = make_genotypes(&haplotypes, self.options.ploidy);
            let prior_model = UniformGenotypePriorModel::new(genotypes.len());
            let inference = IndividualModel::new(&prior_model);
            let result = inference.evaluate(&genotypes, &likelihood_array);
        }

        Ok(calls)
    }
}
