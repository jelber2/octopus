// Converted from C++ to Rust

pub const VCF_VERSION: &str = "VCFv4.3";
pub const MISSING_VALUE: &str = ".";
pub const PASS_FILTER: &str = "PASS";
pub const INFO_SEP: char = ';';
pub const FORMAT_SEP: char = ':';
pub const ALLELE_SEP: char = ',';
pub const PHASED_ALLELE_SEP: char = '|';
pub const UNPHASED_ALLELE_SEP: char = '/';
pub const MISSING_ALLELE: &str = ".";

pub const GT_KEY: &str = "GT";
pub const DP_KEY: &str = "DP";
pub const GQ_KEY: &str = "GQ";
pub const AD_KEY: &str = "AD";
pub const PL_KEY: &str = "PL";
pub const PS_KEY: &str = "PS";
pub const AF_KEY: &str = "AF";

pub fn is_info_flag(number: &str) -> bool {
    number == "0"
}

pub fn is_missing(value: &str) -> bool {
    value == MISSING_VALUE
}
