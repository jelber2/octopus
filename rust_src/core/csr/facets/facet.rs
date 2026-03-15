// Converted from C++ to Rust

use std::any::Any;
use crate::io::variant::vcf_record::VcfRecord;

pub trait Facet: Any + Send + Sync {
    fn name(&self) -> &str;
    fn as_any(&self) -> &dyn Any;
}

pub struct FacetWrapper(pub Box<dyn Facet>);
