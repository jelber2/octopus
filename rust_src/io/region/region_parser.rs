// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;

/// Parse a region string in any of these formats (all 0-based, half-open):
///   "contig"           → [0, u32::MAX)  (whole contig)
///   "contig:pos"       → [pos, pos+1)   (single position)
///   "contig:begin-end" → [begin, end)
///
/// Commas in numeric fields are stripped (thousands-separator convention).
pub fn parse_region(region_str: &str) -> Result<GenomicRegion, String> {
    if region_str.is_empty() {
        return Err("Empty region string".to_string());
    }

    if let Some(colon_pos) = region_str.find(':') {
        let contig = &region_str[..colon_pos];
        let range  = &region_str[colon_pos + 1..];

        if contig.is_empty() {
            return Err(format!("Empty contig name in region: '{}'", region_str));
        }
        if range.is_empty() {
            return Err(format!("Missing position after ':' in region: '{}'", region_str));
        }

        // Reject if there is a second colon inside the range part
        if range.contains(':') {
            return Err(format!("Invalid region format (extra ':'): '{}'", region_str));
        }

        if let Some(dash_pos) = range.find('-') {
            let begin_str = range[..dash_pos].replace(',', "");
            let end_str   = range[dash_pos + 1..].replace(',', "");

            if begin_str.is_empty() || end_str.is_empty() {
                return Err(format!("Invalid begin or end in region: '{}'", region_str));
            }

            let begin: u32 = begin_str.parse()
                .map_err(|_| format!("Invalid begin position in region: '{}'", region_str))?;
            let end: u32 = end_str.parse()
                .map_err(|_| format!("Invalid end position in region: '{}'", region_str))?;

            if begin > end {
                return Err(format!("Begin ({}) > end ({}) in region: '{}'", begin, end, region_str));
            }

            GenomicRegion::new(contig, begin, end)
                .map_err(|e| e.to_string())
        } else {
            // Single position: 0-based, expand to [pos, pos+1)
            let pos_str = range.replace(',', "");
            let pos: u32 = pos_str.parse()
                .map_err(|_| format!("Invalid position in region: '{}'", region_str))?;
            GenomicRegion::new(contig, pos, pos + 1)
                .map_err(|e| e.to_string())
        }
    } else {
        // No colon → treat the whole string as a contig name.
        // Reject strings that look like bare numeric ranges or lone dashes.
        if region_str.starts_with('-') || region_str.starts_with(':') {
            return Err(format!("Invalid region: '{}'", region_str));
        }
        GenomicRegion::new(region_str, 0, u32::MAX)
            .map_err(|e| e.to_string())
    }
}

pub fn parse_regions(region_strs: &[String]) -> Result<Vec<GenomicRegion>, String> {
    region_strs.iter().map(|s| parse_region(s)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Invalid inputs — ported from region_parser_tests.cpp ─────────────

    #[test]
    fn parse_region_rejects_empty_string() {
        assert!(parse_region("").is_err());
    }

    #[test]
    fn parse_region_rejects_dash_only() {
        assert!(parse_region("-").is_err());
    }

    #[test]
    fn parse_region_rejects_inverted_range() {
        assert!(parse_region("chr1:100-99").is_err());
        assert!(parse_region("5:200-100").is_err());
    }

    #[test]
    fn parse_region_rejects_non_numeric_positions() {
        assert!(parse_region("chr1:o-1").is_err());
        assert!(parse_region("chr1:0-1o").is_err());
        assert!(parse_region("chr1:0-1o0").is_err());
    }

    #[test]
    fn parse_region_rejects_double_colon() {
        assert!(parse_region("chr2::0-323").is_err());
    }

    #[test]
    fn parse_region_rejects_empty_contig() {
        assert!(parse_region(":0-10").is_err());
    }

    #[test]
    fn parse_region_rejects_missing_end_position() {
        assert!(parse_region("chr1:100-").is_err());
        assert!(parse_region("chr1:-200").is_err());
    }

    #[test]
    fn parse_region_rejects_colon_with_no_position() {
        assert!(parse_region("chr1:").is_err());
    }

    // ── Valid inputs — ported from region_parser_tests.cpp ───────────────

    #[test]
    fn parse_region_contig_only_spans_full_extent() {
        let r = parse_region("chr1").unwrap();
        assert_eq!(r.contig_name(), "chr1");
        assert_eq!(r.begin(), 0);
        assert_eq!(r.end(), u32::MAX);
    }

    #[test]
    fn parse_region_with_explicit_range() {
        let r = parse_region("chr1:100-200").unwrap();
        assert_eq!(r.contig_name(), "chr1");
        assert_eq!(r.begin(), 100);
        assert_eq!(r.end(), 200);
    }

    #[test]
    fn parse_region_single_position_becomes_unit_interval() {
        // "contig:pos" → 0-based position, expanded to [pos, pos+1)
        let r = parse_region("chr1:10").unwrap();
        assert_eq!(r.contig_name(), "chr1");
        assert_eq!(r.begin(), 10);
        assert_eq!(r.end(), 11);
    }

    #[test]
    fn parse_region_equal_begin_end_is_empty_region() {
        let r = parse_region("chr1:10-10").unwrap();
        assert_eq!(r.contig_name(), "chr1");
        assert_eq!(r.begin(), 10);
        assert_eq!(r.end(), 10);
    }

    #[test]
    fn parse_region_zero_range_is_valid() {
        let r = parse_region("chr1:0-0").unwrap();
        assert_eq!(r.begin(), 0);
        assert_eq!(r.end(), 0);
    }

    #[test]
    fn parse_region_comma_separated_numbers() {
        // Commas are thousands separators: "1,21" == 121, "2,91" == 291
        let r = parse_region("chr4:1,21-2,91").unwrap();
        assert_eq!(r.contig_name(), "chr4");
        assert_eq!(r.begin(), 121);
        assert_eq!(r.end(), 291);
    }

    #[test]
    fn parse_region_with_leading_zeros() {
        let r = parse_region("chr6:00-0100").unwrap();
        assert_eq!(r.contig_name(), "chr6");
        assert_eq!(r.begin(), 0);
        assert_eq!(r.end(), 100);
    }

    #[test]
    fn parse_region_numeric_contig_name() {
        let r = parse_region("5:3").unwrap();
        assert_eq!(r.contig_name(), "5");
        assert_eq!(r.begin(), 3);
        assert_eq!(r.end(), 4);
    }

    #[test]
    fn parse_regions_batch() {
        let inputs: Vec<String> = vec![
            "chr1:0-100".to_string(),
            "chr2:200-300".to_string(),
        ];
        let regions = parse_regions(&inputs).unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].contig_name(), "chr1");
        assert_eq!(regions[1].contig_name(), "chr2");
    }

    #[test]
    fn parse_regions_batch_fails_on_bad_entry() {
        let inputs: Vec<String> = vec![
            "chr1:0-100".to_string(),
            "chr1:50-10".to_string(), // inverted
        ];
        assert!(parse_regions(&inputs).is_err());
    }
}
