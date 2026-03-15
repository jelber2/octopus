// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::collections::HashMap;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;
use crate::io::read::read_manager::ReadManager;
use super::filtering::read_filter::ReadFilter;
use super::transformers::read_transformer::ReadTransformer;
use super::downsampling::downsampler::Downsampler;

pub type SampleName = String;
pub type SampleReadMap = HashMap<SampleName, Vec<AlignedRead>>;

pub struct ReadPipe {
    manager: ReadManager,
    filters: Vec<Box<dyn ReadFilter>>,
    transformers: Vec<Box<dyn ReadTransformer>>,
    downsampler: Option<Box<dyn Downsampler>>,
}

impl ReadPipe {
    pub fn new(manager: ReadManager) -> Self {
        ReadPipe {
            manager,
            filters: Vec::new(),
            transformers: Vec::new(),
            downsampler: None,
        }
    }

    pub fn add_filter(&mut self, filter: Box<dyn ReadFilter>) {
        self.filters.push(filter);
    }

    pub fn add_transformer(&mut self, transformer: Box<dyn ReadTransformer>) {
        self.transformers.push(transformer);
    }

    pub fn set_downsampler(&mut self, downsampler: Box<dyn Downsampler>) {
        self.downsampler = Some(downsampler);
    }

    pub fn fetch(&self, region: &GenomicRegion) -> SampleReadMap {
        let mut reads = self.manager.fetch(region);
        self.process(&mut reads, region);
        reads
    }

    pub fn fetch_for_sample(&self, sample: &str, region: &GenomicRegion) -> Vec<AlignedRead> {
        let raw = self.manager.fetch_for_sample(sample, region);
        let mut reads = HashMap::from([(sample.to_string(), raw)]);
        self.process(&mut reads, region);
        reads.remove(sample).unwrap_or_default()
    }

    fn process(&self, reads: &mut SampleReadMap, region: &GenomicRegion) {
        for (_, sample_reads) in reads.iter_mut() {
            for filter in &self.filters {
                sample_reads.retain(|r| filter.passes(r));
            }
            for transformer in &self.transformers {
                for read in sample_reads.iter_mut() {
                    transformer.transform(read);
                }
            }
        }

        if let Some(downsampler) = &self.downsampler {
            for (_, sample_reads) in reads.iter_mut() {
                downsampler.downsample(sample_reads, region);
            }
        }
    }
}
