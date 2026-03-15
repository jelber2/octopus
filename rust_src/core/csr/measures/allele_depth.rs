// Converted from C++ to Rust

use super::measure::{Measure, MeasureValue};
use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;

pub struct AlleleDepth;

impl Measure for AlleleDepth {
    fn name(&self) -> &str { "AD" }

    fn evaluate(&self, record: &VcfRecord, _facets: &[FacetWrapper]) -> MeasureValue {
        for sample in record.samples() {
            if let Some(ad_vals) = record.get_sample_value(sample, "AD") {
                let values: Vec<MeasureValue> = ad_vals.iter()
                    .filter_map(|v| v.parse::<i64>().ok().map(MeasureValue::Int))
                    .collect();
                if !values.is_empty() {
                    return MeasureValue::Vec(values);
                }
            }
        }
        MeasureValue::None
    }
}
