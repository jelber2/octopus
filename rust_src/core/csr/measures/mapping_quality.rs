// Converted from C++ to Rust

use super::measure::{Measure, MeasureValue};
use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;

pub struct MappingQualityMeasure;

impl Measure for MappingQualityMeasure {
    fn name(&self) -> &str { "MQ" }

    fn evaluate(&self, record: &VcfRecord, _facets: &[FacetWrapper]) -> MeasureValue {
        if let Some(mq_vals) = record.info_value("MQ") {
            if let Some(v) = mq_vals.first().and_then(|s| s.parse::<f64>().ok()) {
                return MeasureValue::Float(v);
            }
        }
        MeasureValue::None
    }
}
