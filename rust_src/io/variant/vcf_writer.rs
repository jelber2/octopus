// Converted from C++ to Rust

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use super::vcf_record::VcfRecord;
use super::vcf_header::VcfHeader;

pub struct VcfWriter {
    writer: BufWriter<File>,
    header_written: bool,
}

impl VcfWriter {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, String> {
        let file = File::create(path.as_ref()).map_err(|e| e.to_string())?;
        Ok(VcfWriter { writer: BufWriter::new(file), header_written: false })
    }

    pub fn write_header(&mut self, header: &VcfHeader) -> Result<(), String> {
        write!(self.writer, "{}", header).map_err(|e| e.to_string())?;
        writeln!(self.writer).map_err(|e| e.to_string())?;
        self.header_written = true;
        Ok(())
    }

    pub fn write_record(&mut self, record: &VcfRecord) -> Result<(), String> {
        if !self.header_written {
            return Err("Header not written before record".to_string());
        }

        let qual_str = record.qual()
            .map(|q| format!("{:.2}", q))
            .unwrap_or(".".to_string());

        let filter_str = if record.filters().is_empty() {
            ".".to_string()
        } else {
            record.filters().join(";")
        };

        let info_str = if record.info_keys().is_empty() {
            ".".to_string()
        } else {
            record.info_keys().iter()
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

        let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.chrom(), record.pos(), record.id(),
            record.ref_allele(), record.alt().join(","),
            qual_str, filter_str, info_str);

        writeln!(self.writer, "{}", line).map_err(|e| e.to_string())
    }

    pub fn flush(&mut self) -> Result<(), String> {
        self.writer.flush().map_err(|e| e.to_string())
    }
}
