// Converted from C++ to Rust

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Sample(pub String);

impl Sample {
    pub fn new(name: impl Into<String>) -> Self {
        Sample(name.into())
    }
    pub fn name(&self) -> &str {
        &self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pedigree {
    mother: Option<Sample>,
    father: Option<Sample>,
    child: Sample,
}

impl Pedigree {
    pub fn new(child: Sample) -> Self {
        Pedigree { mother: None, father: None, child }
    }

    pub fn with_parents(child: Sample, mother: Sample, father: Sample) -> Self {
        Pedigree { mother: Some(mother), father: Some(father), child }
    }

    pub fn child(&self) -> &Sample { &self.child }
    pub fn mother(&self) -> Option<&Sample> { self.mother.as_ref() }
    pub fn father(&self) -> Option<&Sample> { self.father.as_ref() }
}
