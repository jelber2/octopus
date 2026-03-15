// Converted from C++ to Rust

use super::measure::{Measure, MeasureValue};
use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;

pub struct DepthMeasure;

impl Measure for DepthMeasure {
    fn name(&self) -> &str { "DP" }

    fn evaluate(&self, record: &VcfRecord, _facets: &[FacetWrapper]) -> MeasureValue {
        if let Some(dp_vals) = record.info_value("DP") {
            if let Some(v) = dp_vals.first().and_then(|s| s.parse::<i64>().ok()) {
                return MeasureValue::Int(v);
            }
        }
        MeasureValue::None
    }
}
