// Converted from C++ to Rust

use crate::core::types::genotype::Genotype;

pub type LogProbability = f64;

pub trait GenotypePriorModel: Send + Sync {
    fn log_prior(&self, genotype: &Genotype) -> LogProbability;
}

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
