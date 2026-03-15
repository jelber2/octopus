// Converted from C++ to Rust

use std::any::Any;
use std::collections::HashMap;
use super::facet::Facet;
use crate::core::types::genotype::Genotype;

pub struct GenotypesFacet {
    sample_genotypes: HashMap<String, Genotype>,
}

impl GenotypesFacet {
    pub fn new(sample_genotypes: HashMap<String, Genotype>) -> Self {
        GenotypesFacet { sample_genotypes }
    }
    pub fn get(&self, sample: &str) -> Option<&Genotype> {
        self.sample_genotypes.get(sample)
    }
}

impl Facet for GenotypesFacet {
    fn name(&self) -> &str { "Genotypes" }
    fn as_any(&self) -> &dyn Any { self }
}
