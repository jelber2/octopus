// Converted from C++ to Rust

use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;
use std::fmt;

#[derive(Debug, Clone)]
pub enum MeasureValue {
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(String),
    None,
    Vec(Vec<MeasureValue>),
}

impl fmt::Display for MeasureValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MeasureValue::Int(v) => write!(f, "{}", v),
            MeasureValue::Float(v) => write!(f, "{:.4}", v),
            MeasureValue::Bool(v) => write!(f, "{}", v),
            MeasureValue::Str(v) => write!(f, "{}", v),
            MeasureValue::None => write!(f, "."),
            MeasureValue::Vec(v) => write!(f, "[{}]", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
        }
    }
}

pub trait Measure: Send + Sync {
    fn name(&self) -> &str;
    fn evaluate(&self, record: &VcfRecord, facets: &[FacetWrapper]) -> MeasureValue;
}
