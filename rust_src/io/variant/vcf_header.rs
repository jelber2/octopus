// Converted from C++ to Rust

use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct VcfHeader {
    file_format: String,
    contigs: Vec<ContigInfo>,
    filters: Vec<FilterInfo>,
    info_fields: Vec<InfoField>,
    format_fields: Vec<FormatField>,
    samples: Vec<String>,
    extra: HashMap<String, Vec<String>>,
}

#[derive(Debug, Clone)]
pub struct ContigInfo {
    pub id: String,
    pub length: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct FilterInfo {
    pub id: String,
    pub description: String,
}

#[derive(Debug, Clone)]
pub struct InfoField {
    pub id: String,
    pub number: String,
    pub type_: String,
    pub description: String,
}

#[derive(Debug, Clone)]
pub struct FormatField {
    pub id: String,
    pub number: String,
    pub type_: String,
    pub description: String,
}

impl VcfHeader {
    pub fn new(file_format: impl Into<String>) -> Self {
        VcfHeader {
            file_format: file_format.into(),
            ..Default::default()
        }
    }

    pub fn file_format(&self) -> &str { &self.file_format }
    pub fn samples(&self) -> &[String] { &self.samples }
    pub fn num_samples(&self) -> usize { self.samples.len() }

    pub fn add_sample(&mut self, sample: impl Into<String>) {
        self.samples.push(sample.into());
    }

    pub fn add_contig(&mut self, id: impl Into<String>, length: Option<u64>) {
        self.contigs.push(ContigInfo { id: id.into(), length });
    }

    pub fn add_filter(&mut self, id: impl Into<String>, description: impl Into<String>) {
        self.filters.push(FilterInfo { id: id.into(), description: description.into() });
    }

    pub fn add_info(&mut self, id: impl Into<String>, number: impl Into<String>, type_: impl Into<String>, description: impl Into<String>) {
        self.info_fields.push(InfoField { id: id.into(), number: number.into(), type_: type_.into(), description: description.into() });
    }

    pub fn add_format(&mut self, id: impl Into<String>, number: impl Into<String>, type_: impl Into<String>, description: impl Into<String>) {
        self.format_fields.push(FormatField { id: id.into(), number: number.into(), type_: type_.into(), description: description.into() });
    }

    pub fn has_contig(&self, id: &str) -> bool {
        self.contigs.iter().any(|c| c.id == id)
    }
}

impl std::fmt::Display for VcfHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "##fileformat={}", self.file_format)?;
        for c in &self.contigs {
            if let Some(len) = c.length {
                writeln!(f, "##contig=<ID={},length={}>", c.id, len)?;
            } else {
                writeln!(f, "##contig=<ID={}>", c.id)?;
            }
        }
        for fi in &self.filters {
            writeln!(f, "##FILTER=<ID={},Description=\"{}\">", fi.id, fi.description)?;
        }
        for info in &self.info_fields {
            writeln!(f, "##INFO=<ID={},Number={},Type={},Description=\"{}\">", info.id, info.number, info.type_, info.description)?;
        }
        for fmt in &self.format_fields {
            writeln!(f, "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">", fmt.id, fmt.number, fmt.type_, fmt.description)?;
        }
        write!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        if !self.samples.is_empty() {
            write!(f, "\tFORMAT")?;
            for s in &self.samples {
                write!(f, "\t{}", s)?;
            }
        }
        Ok(())
    }
}
