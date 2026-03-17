// Converted from C++ to Rust

use std::path::PathBuf;
use std::collections::HashMap;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;
use super::bam_reader::BamReader;

pub type SampleName = String;
pub type ReadContainer = Vec<AlignedRead>;
pub type SampleReadMap = HashMap<SampleName, ReadContainer>;

pub trait ReadReaderImpl: Send + Sync {
    fn samples(&self) -> Vec<SampleName>;
    fn fetch(&self, sample: &str, region: &GenomicRegion) -> ReadContainer;
}

pub struct ReadManager {
    readers: Vec<Box<dyn ReadReaderImpl>>,
    samples: Vec<SampleName>,
    max_open_files: usize,
}

impl ReadManager {
    /// Open all BAM/CRAM files in `paths`.  Returns an error if any file
    /// cannot be opened or its header cannot be read.
    pub fn new(paths: Vec<PathBuf>, max_open_files: usize) -> Result<Self, String> {
        let mut readers: Vec<Box<dyn ReadReaderImpl>> = Vec::new();
        let mut samples: Vec<SampleName> = Vec::new();

        for path in &paths {
            let reader = BamReader::open(path)?;
            for sample in reader.samples() {
                if !samples.contains(&sample) {
                    samples.push(sample);
                }
            }
            readers.push(Box::new(reader));
        }

        Ok(ReadManager { readers, samples, max_open_files })
    }

    pub fn num_files(&self) -> usize { self.readers.len() }
    pub fn num_samples(&self) -> usize { self.samples.len() }
    pub fn samples(&self) -> &[SampleName] { &self.samples }

    pub fn fetch(&self, region: &GenomicRegion) -> SampleReadMap {
        let mut result = SampleReadMap::new();
        for sample in &self.samples {
            for reader in &self.readers {
                if reader.samples().iter().any(|s| s == sample) {
                    let reads = reader.fetch(sample, region);
                    result.entry(sample.clone()).or_default().extend(reads);
                }
            }
        }
        result
    }

    pub fn fetch_for_sample(&self, sample: &str, region: &GenomicRegion) -> ReadContainer {
        let mut result = ReadContainer::new();
        for reader in &self.readers {
            if reader.samples().iter().any(|s| s == sample) {
                result.extend(reader.fetch(sample, region));
            }
        }
        result
    }
}
