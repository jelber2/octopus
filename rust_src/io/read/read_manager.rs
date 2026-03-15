// Converted from C++ to Rust

use std::path::PathBuf;
use std::collections::HashMap;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;

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
    pub fn new(paths: Vec<PathBuf>, max_open_files: usize) -> Self {
        ReadManager {
            readers: Vec::new(),
            samples: Vec::new(),
            max_open_files,
        }
    }

    pub fn num_files(&self) -> usize { self.readers.len() }
    pub fn num_samples(&self) -> usize { self.samples.len() }
    pub fn samples(&self) -> &[SampleName] { &self.samples }

    pub fn fetch(&self, region: &GenomicRegion) -> SampleReadMap {
        let mut result = SampleReadMap::new();
        for reader in &self.readers {
            for sample in reader.samples() {
                let reads = reader.fetch(&sample, region);
                result.entry(sample).or_default().extend(reads);
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
