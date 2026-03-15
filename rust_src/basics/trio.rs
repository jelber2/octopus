// Converted from C++ to Rust

use super::pedigree::Sample;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trio {
    mother: Sample,
    father: Sample,
    child: Sample,
}

impl Trio {
    pub fn new(mother: Sample, father: Sample, child: Sample) -> Self {
        Trio { mother, father, child }
    }

    pub fn mother(&self) -> &Sample { &self.mother }
    pub fn father(&self) -> &Sample { &self.father }
    pub fn child(&self) -> &Sample { &self.child }
}
