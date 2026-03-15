// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;
use crate::basics::tandem_repeat::TandemRepeat;

pub fn find_repeats(region: &GenomicRegion, sequence: &[u8], min_period: usize, max_period: usize) -> Vec<TandemRepeat> {
    let mut repeats = Vec::new();
    let len = sequence.len();

    for period in min_period..=std::cmp::min(max_period, len / 2) {
        let mut pos = 0;
        while pos + period < len {
            let mut repeat_len = period;
            while pos + repeat_len + period <= len
                && &sequence[pos..pos + period] == &sequence[pos + repeat_len..pos + repeat_len + period]
            {
                repeat_len += period;
            }
            if repeat_len > period {
                let begin = region.begin() + pos as u32;
                let end = begin + repeat_len as u32;
                if let Ok(r) = GenomicRegion::new(region.contig_name(), begin, end) {
                    let seq_str = String::from_utf8_lossy(&sequence[pos..pos + period]).to_string();
                    repeats.push(TandemRepeat::new(r, period, seq_str));
                }
                pos += repeat_len;
            } else {
                pos += 1;
            }
        }
    }

    repeats
}
