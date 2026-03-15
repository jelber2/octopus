// Converted from C++ to Rust

use super::measure::{Measure, MeasureValue};
use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;

pub struct StrandBias;

impl Measure for StrandBias {
    fn name(&self) -> &str { "SB" }

    fn evaluate(&self, record: &VcfRecord, _facets: &[FacetWrapper]) -> MeasureValue {
        if let Some(sb_vals) = record.info_value("SB") {
            let values: Vec<MeasureValue> = sb_vals.iter()
                .filter_map(|v| v.parse::<i64>().ok().map(MeasureValue::Int))
                .collect();
            if !values.is_empty() {
                return MeasureValue::Vec(values);
            }
        }
        MeasureValue::None
    }
}
