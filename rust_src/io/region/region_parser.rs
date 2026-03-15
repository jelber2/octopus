// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;

pub fn parse_region(region_str: &str) -> Result<GenomicRegion, String> {
    if let Some(colon_pos) = region_str.find(':') {
        let contig = &region_str[..colon_pos];
        let range = &region_str[colon_pos + 1..];

        if let Some(dash_pos) = range.find('-') {
            let begin: u32 = range[..dash_pos].replace(',', "").parse()
                .map_err(|_| format!("Invalid begin position in region: {}", region_str))?;
            let end: u32 = range[dash_pos + 1..].replace(',', "").parse()
                .map_err(|_| format!("Invalid end position in region: {}", region_str))?;
            GenomicRegion::new(contig, begin.saturating_sub(1), end)
                .map_err(|e| e.to_string())
        } else {
            let pos: u32 = range.replace(',', "").parse()
                .map_err(|_| format!("Invalid position in region: {}", region_str))?;
            let begin = pos.saturating_sub(1);
            GenomicRegion::new(contig, begin, begin + 1)
                .map_err(|e| e.to_string())
        }
    } else {
        GenomicRegion::new(region_str, 0, u32::MAX)
            .map_err(|e| e.to_string())
    }
}

pub fn parse_regions(region_strs: &[String]) -> Result<Vec<GenomicRegion>, String> {
    region_strs.iter().map(|s| parse_region(s)).collect()
}
