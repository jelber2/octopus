// Converted from C++ to Rust
// FASTA file reader implementation

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use super::reference_genome::{ContigName, GeneticSequence, ReferenceReader};
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::contig_region::Size;

pub struct FastaReader {
    path: PathBuf,
    name: String,
    contig_offsets: HashMap<String, u64>,
    contig_sizes: HashMap<String, Size>,
    contig_names: Vec<String>,
    line_length: HashMap<String, usize>,
}

impl FastaReader {
    pub fn new(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref().to_path_buf();
        let name = path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        let mut reader = FastaReader {
            path: path.clone(),
            name,
            contig_offsets: HashMap::new(),
            contig_sizes: HashMap::new(),
            contig_names: Vec::new(),
            line_length: HashMap::new(),
        };

        reader.index()?;
        Ok(reader)
    }

    fn index(&mut self) -> Result<(), String> {
        let file = File::open(&self.path).map_err(|e| e.to_string())?;
        let mut reader = BufReader::new(file);
        let mut current_name = String::new();
        let mut current_offset = 0u64;
        let mut current_size = 0u32;
        let mut line = String::new();
        let mut offset = 0u64;

        loop {
            line.clear();
            let bytes = reader.read_line(&mut line).map_err(|e| e.to_string())?;
            if bytes == 0 { break; }

            if line.starts_with('>') {
                if !current_name.is_empty() {
                    self.contig_sizes.insert(current_name.clone(), current_size);
                    current_size = 0;
                }
                current_name = line[1..].trim().split_whitespace().next().unwrap_or("").to_string();
                self.contig_names.push(current_name.clone());
                self.contig_offsets.insert(current_name.clone(), offset + bytes as u64);
                current_offset = offset + bytes as u64;
            } else {
                let seq_len = line.trim().len() as u32;
                current_size += seq_len;
                if !current_name.is_empty() && !self.line_length.contains_key(&current_name) && seq_len > 0 {
                    self.line_length.insert(current_name.clone(), line.trim().len());
                }
            }
            offset += bytes as u64;
        }

        if !current_name.is_empty() {
            self.contig_sizes.insert(current_name, current_size);
        }

        Ok(())
    }
}

impl ReferenceReader for FastaReader {
    fn name(&self) -> &str {
        &self.name
    }

    fn contig_names(&self) -> Vec<ContigName> {
        self.contig_names.clone()
    }

    fn contig_size(&self, contig: &str) -> Option<Size> {
        self.contig_sizes.get(contig).copied()
    }

    fn fetch_sequence(&self, region: &GenomicRegion) -> Result<GeneticSequence, String> {
        let contig = region.contig_name();
        let begin = region.begin() as usize;
        let end = region.end() as usize;
        let length = end - begin;

        let offset = *self.contig_offsets.get(contig)
            .ok_or_else(|| format!("Contig {} not found", contig))?;
        let line_len = *self.line_length.get(contig).unwrap_or(&60);

        let mut file = File::open(&self.path).map_err(|e| e.to_string())?;

        let newlines_before = begin / line_len;
        let seek_offset = offset + begin as u64 + newlines_before as u64;
        file.seek(SeekFrom::Start(seek_offset)).map_err(|e| e.to_string())?;

        let mut result = Vec::with_capacity(length);
        let mut reader = BufReader::new(file);

        while result.len() < length {
            let mut line = String::new();
            let bytes = reader.read_line(&mut line).map_err(|e| e.to_string())?;
            if bytes == 0 { break; }
            if line.starts_with('>') { break; }
            let trimmed = line.trim().as_bytes();
            let remaining = length - result.len();
            let take = std::cmp::min(remaining, trimmed.len());
            result.extend_from_slice(&trimmed[..take]);
        }

        Ok(result.iter().map(|&b| b.to_ascii_uppercase()).collect())
    }
}
