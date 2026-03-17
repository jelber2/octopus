// Converted from C++ to Rust
// Writes VCF records to a file or stdout, including FORMAT/sample columns.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use super::vcf_record::VcfRecord;
use super::vcf_header::VcfHeader;

pub struct VcfWriter {
    writer: Box<dyn Write>,
    header_written: bool,
    sample_names: Vec<String>,
}

impl VcfWriter {
    /// Write to a file.
    pub fn to_file(path: impl AsRef<Path>) -> Result<Self, String> {
        let file = File::create(path.as_ref()).map_err(|e| e.to_string())?;
        Ok(VcfWriter {
            writer: Box::new(BufWriter::new(file)),
            header_written: false,
            sample_names: Vec::new(),
        })
    }

    /// Write to stdout.
    pub fn to_stdout() -> Self {
        VcfWriter {
            writer: Box::new(BufWriter::new(io::stdout())),
            header_written: false,
            sample_names: Vec::new(),
        }
    }

    pub fn write_header(&mut self, header: &VcfHeader) -> Result<(), String> {
        self.sample_names = header.samples().to_vec();
        write!(self.writer, "{}", header).map_err(|e| e.to_string())?;
        writeln!(self.writer).map_err(|e| e.to_string())?;
        self.header_written = true;
        Ok(())
    }

    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), String> {
        if !self.header_written {
            return Err("Header not written before record".to_string());
        }

        let qual_str = record
            .qual()
            .map(|q| format!("{:.2}", q))
            .unwrap_or_else(|| ".".to_string());

        let filter_str = if record.filters().is_empty() {
            ".".to_string()
        } else {
            record.filters().join(";")
        };

        let info_str = if record.info_keys().is_empty() {
            ".".to_string()
        } else {
            record
                .info_keys()
                .iter()
                .map(|k| {
                    if let Some(vals) = record.info_value(k) {
                        if vals.len() == 1 && vals[0] == "true" {
                            k.to_string()
                        } else {
                            format!("{}={}", k, vals.join(","))
                        }
                    } else {
                        k.to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(";")
        };

        write!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.chrom(),
            record.pos(),
            record.id(),
            record.ref_allele(),
            record.alt().join(","),
            qual_str,
            filter_str,
            info_str,
        )
        .map_err(|e| e.to_string())?;

        // FORMAT + one column per sample (in header order)
        if !record.format().is_empty() && !self.sample_names.is_empty() {
            write!(self.writer, "\t{}", record.format().join(":"))
                .map_err(|e| e.to_string())?;

            for sample_name in &self.sample_names {
                write!(self.writer, "\t").map_err(|e| e.to_string())?;
                let has_data = record.has_sample(sample_name);
                let field_values: Vec<String> = record
                    .format()
                    .iter()
                    .map(|key| {
                        if has_data {
                            record
                                .get_sample_value(sample_name, key)
                                .map(|vals| vals.join(","))
                                .unwrap_or_else(|| ".".to_string())
                        } else {
                            ".".to_string()
                        }
                    })
                    .collect();
                write!(self.writer, "{}", field_values.join(":"))
                    .map_err(|e| e.to_string())?;
            }
        }

        writeln!(self.writer).map_err(|e| e.to_string())
    }

    pub fn flush(&mut self) -> Result<(), String> {
        self.writer.flush().map_err(|e| e.to_string())
    }
}
