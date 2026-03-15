// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};
use super::contig_region::{ContigRegion, Position, Size, Distance};
use super::contig_region;

#[derive(Debug, Clone, Default, Eq, PartialEq)]
pub struct GenomicRegion {
    contig_name: String,
    contig_region: ContigRegion,
}

#[derive(Debug, Clone)]
pub struct BadRegionCompare {
    first: String,
    second: String,
}

impl fmt::Display for BadRegionCompare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "BadRegionCompare: comparing regions on different contigs: {} & {}", self.first, self.second)
    }
}

impl std::error::Error for BadRegionCompare {}

impl GenomicRegion {
    pub fn new(contig_name: impl Into<String>, begin: Position, end: Position) -> Result<Self, super::contig_region::BadRegion> {
        let cr = ContigRegion::new(begin, end)?;
        Ok(GenomicRegion { contig_name: contig_name.into(), contig_region: cr })
    }

    pub fn from_region(contig_name: impl Into<String>, region: ContigRegion) -> Self {
        GenomicRegion { contig_name: contig_name.into(), contig_region: region }
    }

    pub fn contig_name(&self) -> &str {
        &self.contig_name
    }

    pub fn contig_region(&self) -> &ContigRegion {
        &self.contig_region
    }

    pub fn begin(&self) -> Position {
        self.contig_region.begin()
    }

    pub fn end(&self) -> Position {
        self.contig_region.end()
    }

    pub fn overlaps(&self, other: &GenomicRegion) -> bool {
        overlaps(self, other)
    }

    pub fn contains_region(&self, other: &GenomicRegion) -> bool {
        contains(self, other)
    }
}

impl PartialOrd for GenomicRegion {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.contig_name != other.contig_name {
            return None;
        }
        self.contig_region.partial_cmp(&other.contig_region)
    }
}

impl Hash for GenomicRegion {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.contig_name.hash(state);
        self.contig_region.hash(state);
    }
}

impl fmt::Display for GenomicRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig_name, self.begin(), self.end())
    }
}

pub fn is_same_contig(lhs: &GenomicRegion, rhs: &GenomicRegion) -> bool {
    lhs.contig_name() == rhs.contig_name()
}

pub fn mapped_begin(region: &GenomicRegion) -> Position {
    region.begin()
}

pub fn mapped_end(region: &GenomicRegion) -> Position {
    region.end()
}

pub fn is_empty(region: &GenomicRegion) -> bool {
    contig_region::is_empty(region.contig_region())
}

pub fn size(region: &GenomicRegion) -> Size {
    contig_region::size(region.contig_region())
}

pub fn is_position(region: &GenomicRegion) -> bool {
    contig_region::is_position(region.contig_region())
}

fn check_same_contig(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<(), BadRegionCompare> {
    if !is_same_contig(lhs, rhs) {
        Err(BadRegionCompare { first: lhs.to_string(), second: rhs.to_string() })
    } else {
        Ok(())
    }
}

pub fn begins_equal(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<bool, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::begins_equal(lhs.contig_region(), rhs.contig_region()))
}

pub fn ends_equal(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<bool, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::ends_equal(lhs.contig_region(), rhs.contig_region()))
}

pub fn begins_before(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<bool, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::begins_before(lhs.contig_region(), rhs.contig_region()))
}

pub fn ends_before(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<bool, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::ends_before(lhs.contig_region(), rhs.contig_region()))
}

pub fn is_before(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<bool, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::is_before(lhs.contig_region(), rhs.contig_region()))
}

pub fn is_after(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<bool, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::is_after(lhs.contig_region(), rhs.contig_region()))
}

pub fn are_adjacent(lhs: &GenomicRegion, rhs: &GenomicRegion) -> bool {
    is_same_contig(lhs, rhs) && contig_region::are_adjacent(lhs.contig_region(), rhs.contig_region())
}

pub fn overlap_size(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Distance {
    if is_same_contig(lhs, rhs) {
        contig_region::overlap_size(lhs.contig_region(), rhs.contig_region())
    } else {
        0
    }
}

pub fn overlaps(lhs: &GenomicRegion, rhs: &GenomicRegion) -> bool {
    is_same_contig(lhs, rhs) && contig_region::overlaps(lhs.contig_region(), rhs.contig_region())
}

pub fn contains(lhs: &GenomicRegion, rhs: &GenomicRegion) -> bool {
    is_same_contig(lhs, rhs) && contig_region::contains(lhs.contig_region(), rhs.contig_region())
}

pub fn inner_distance(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<Distance, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::inner_distance(lhs.contig_region(), rhs.contig_region()))
}

pub fn outer_distance(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<Distance, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(contig_region::outer_distance(lhs.contig_region(), rhs.contig_region()))
}

pub fn shift(region: &GenomicRegion, n: Distance) -> Result<GenomicRegion, String> {
    let cr = contig_region::shift(region.contig_region(), n)?;
    Ok(GenomicRegion::from_region(region.contig_name(), cr))
}

pub fn next_position(region: &GenomicRegion) -> Result<GenomicRegion, super::contig_region::BadRegion> {
    let cr = contig_region::next_position(region.contig_region())?;
    Ok(GenomicRegion::from_region(region.contig_name(), cr))
}

pub fn expand_lhs(region: &GenomicRegion, n: Distance) -> Result<GenomicRegion, String> {
    let cr = contig_region::expand_lhs(region.contig_region(), n)?;
    Ok(GenomicRegion::from_region(region.contig_name(), cr))
}

pub fn expand_rhs(region: &GenomicRegion, n: Distance) -> Result<GenomicRegion, String> {
    let cr = contig_region::expand_rhs(region.contig_region(), n)?;
    Ok(GenomicRegion::from_region(region.contig_name(), cr))
}

pub fn expand(region: &GenomicRegion, n: Distance) -> Result<GenomicRegion, super::contig_region::BadRegion> {
    let cr = contig_region::expand(region.contig_region(), n)?;
    Ok(GenomicRegion::from_region(region.contig_name(), cr))
}

pub fn encompassing_region(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<GenomicRegion, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    let cr = contig_region::encompassing_region(lhs.contig_region(), rhs.contig_region())
        .map_err(|_| BadRegionCompare { first: lhs.to_string(), second: rhs.to_string() })?;
    Ok(GenomicRegion::from_region(lhs.contig_name(), cr))
}

pub fn intervening_region(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Option<GenomicRegion> {
    if !is_same_contig(lhs, rhs) { return None; }
    let cr = contig_region::intervening_region(lhs.contig_region(), rhs.contig_region())?;
    Some(GenomicRegion::from_region(lhs.contig_name(), cr))
}

pub fn intervening_region_size(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Size {
    if !is_same_contig(lhs, rhs) { return 0; }
    contig_region::intervening_region_size(lhs.contig_region(), rhs.contig_region())
}

pub fn overlapped_region(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Option<GenomicRegion> {
    if !overlaps(lhs, rhs) { return None; }
    let cr = contig_region::overlapped_region(lhs.contig_region(), rhs.contig_region())?;
    Some(GenomicRegion::from_region(lhs.contig_name(), cr))
}

pub fn left_overhang_size(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Size {
    if is_same_contig(lhs, rhs) {
        contig_region::left_overhang_size(lhs.contig_region(), rhs.contig_region())
    } else {
        0
    }
}

pub fn right_overhang_size(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Size {
    if is_same_contig(lhs, rhs) {
        contig_region::right_overhang_size(lhs.contig_region(), rhs.contig_region())
    } else {
        0
    }
}

pub fn left_overhang_region(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<GenomicRegion, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(GenomicRegion::from_region(lhs.contig_name(), contig_region::left_overhang_region(lhs.contig_region(), rhs.contig_region())))
}

pub fn right_overhang_region(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<GenomicRegion, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    Ok(GenomicRegion::from_region(lhs.contig_name(), contig_region::right_overhang_region(lhs.contig_region(), rhs.contig_region())))
}

pub fn closed_region(lhs: &GenomicRegion, rhs: &GenomicRegion) -> Result<GenomicRegion, BadRegionCompare> {
    check_same_contig(lhs, rhs)?;
    let cr = contig_region::closed_region(lhs.contig_region(), rhs.contig_region())
        .map_err(|_| BadRegionCompare { first: lhs.to_string(), second: rhs.to_string() })?;
    Ok(GenomicRegion::from_region(lhs.contig_name(), cr))
}

pub fn head_region(region: &GenomicRegion, n: Size) -> GenomicRegion {
    GenomicRegion::from_region(region.contig_name(), contig_region::head_region(region.contig_region(), n))
}

pub fn head_position(region: &GenomicRegion) -> GenomicRegion {
    GenomicRegion::from_region(region.contig_name(), contig_region::head_position(region.contig_region()))
}

pub fn tail_region(region: &GenomicRegion, n: Size) -> GenomicRegion {
    GenomicRegion::from_region(region.contig_name(), contig_region::tail_region(region.contig_region(), n))
}

pub fn tail_position(region: &GenomicRegion) -> GenomicRegion {
    GenomicRegion::from_region(region.contig_name(), contig_region::tail_position(region.contig_region()))
}

pub fn begin_distance(first: &GenomicRegion, second: &GenomicRegion) -> Result<Distance, BadRegionCompare> {
    check_same_contig(first, second)?;
    Ok(contig_region::begin_distance(first.contig_region(), second.contig_region()))
}

pub fn end_distance(first: &GenomicRegion, second: &GenomicRegion) -> Result<Distance, BadRegionCompare> {
    check_same_contig(first, second)?;
    Ok(contig_region::end_distance(first.contig_region(), second.contig_region()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(contig: &str, begin: u32, end: u32) -> GenomicRegion {
        GenomicRegion::new(contig, begin, end).unwrap()
    }

    #[test]
    fn new_valid_region() {
        let region = r("chr1", 100, 200);
        assert_eq!(region.contig_name(), "chr1");
        assert_eq!(region.begin(), 100);
        assert_eq!(region.end(), 200);
    }

    #[test]
    fn new_zero_length_region_is_ok() {
        let region = r("chr1", 50, 50);
        assert_eq!(region.begin(), 50);
        assert_eq!(region.end(), 50);
        assert!(is_empty(&region));
    }

    #[test]
    fn new_invalid_region_returns_err() {
        assert!(GenomicRegion::new("chr1", 200, 100).is_err());
    }

    #[test]
    fn display_format() {
        let region = r("chr1", 100, 200);
        assert_eq!(region.to_string(), "chr1:100-200");
    }

    #[test]
    fn size_is_correct() {
        assert_eq!(size(&r("chr1", 100, 200)), 100);
        assert_eq!(size(&r("chrX", 0, 0)), 0);
        assert_eq!(size(&r("chr2", 500, 501)), 1);
    }

    #[test]
    fn is_position_single_base() {
        assert!(is_position(&r("chr1", 10, 11)));
        assert!(!is_position(&r("chr1", 10, 12)));
        assert!(!is_position(&r("chr1", 10, 10)));
    }

    #[test]
    fn overlaps_same_contig() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 150, 250);
        let c = r("chr1", 200, 300);
        let d = r("chr1", 300, 400);
        assert!(overlaps(&a, &b));
        assert!(!overlaps(&a, &c));
        assert!(!overlaps(&a, &d));
    }

    #[test]
    fn overlaps_empty_regions_at_boundary() {
        let point = r("chr1", 100, 100);
        let spanning = r("chr1", 50, 150);
        assert!(overlaps(&point, &spanning));
    }

    #[test]
    fn overlaps_different_contig_returns_false() {
        let a = r("chr1", 100, 200);
        let b = r("chr2", 100, 200);
        assert!(!overlaps(&a, &b));
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn contains_correct() {
        let outer = r("chr1", 100, 300);
        let inner = r("chr1", 150, 250);
        let partial = r("chr1", 50, 150);
        assert!(contains(&outer, &inner));
        assert!(outer.contains_region(&inner));
        assert!(!contains(&outer, &partial));
        assert!(!contains(&inner, &outer));
    }

    #[test]
    fn is_same_contig_check() {
        assert!(is_same_contig(&r("chr1", 0, 10), &r("chr1", 5, 15)));
        assert!(!is_same_contig(&r("chr1", 0, 10), &r("chr2", 0, 10)));
    }

    #[test]
    fn is_before_and_after() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 200, 300);
        let c = r("chr1", 50, 99);
        assert!(is_before(&a, &b).unwrap());
        assert!(is_after(&b, &a).unwrap());
        assert!(is_before(&c, &a).unwrap());
    }

    #[test]
    fn before_after_different_contig_errors() {
        let a = r("chr1", 100, 200);
        let b = r("chr2", 300, 400);
        assert!(is_before(&a, &b).is_err());
        assert!(is_after(&a, &b).is_err());
    }

    #[test]
    fn are_adjacent_regions() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 200, 300);
        let c = r("chr1", 201, 300);
        assert!(are_adjacent(&a, &b));
        assert!(!are_adjacent(&a, &c));
        assert!(!are_adjacent(&r("chr1", 0, 10), &r("chr2", 10, 20)));
    }

    #[test]
    fn overlap_size_computed() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 150, 250);
        assert_eq!(overlap_size(&a, &b), 50);
        let non_overlapping = r("chr1", 300, 400);
        assert!(overlap_size(&a, &non_overlapping) <= 0);
    }

    #[test]
    fn encompassing_region_spans_both() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 150, 300);
        let enc = encompassing_region(&a, &b).unwrap();
        assert_eq!(enc.begin(), 100);
        assert_eq!(enc.end(), 300);
    }

    #[test]
    fn encompassing_region_different_contigs_errors() {
        assert!(encompassing_region(&r("chr1", 0, 10), &r("chr2", 0, 10)).is_err());
    }

    #[test]
    fn overlapped_region_intersection() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 150, 250);
        let ov = overlapped_region(&a, &b).unwrap();
        assert_eq!(ov.begin(), 150);
        assert_eq!(ov.end(), 200);
    }

    #[test]
    fn overlapped_region_no_overlap_returns_none() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 300, 400);
        assert!(overlapped_region(&a, &b).is_none());
    }

    #[test]
    fn intervening_region_between_non_adjacent() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 250, 350);
        let iv = intervening_region(&a, &b).unwrap();
        assert_eq!(iv.begin(), 200);
        assert_eq!(iv.end(), 250);
    }

    #[test]
    fn intervening_region_adjacent_is_empty_region() {
        // Adjacent regions [100,200) and [200,300) yield an empty intervening region at
        // position 200, not None — the gap exists but has zero size.
        let a = r("chr1", 100, 200);
        let b = r("chr1", 200, 300);
        let iv = intervening_region(&a, &b);
        assert!(iv.is_some());
        let iv = iv.unwrap();
        assert_eq!(iv.begin(), 200);
        assert_eq!(iv.end(), 200);
        assert!(is_empty(&iv));
    }

    #[test]
    fn shift_positive_and_negative() {
        let region = r("chr1", 100, 200);
        let shifted = shift(&region, 50).unwrap();
        assert_eq!(shifted.begin(), 150);
        assert_eq!(shifted.end(), 250);
        let shifted_back = shift(&region, -50).unwrap();
        assert_eq!(shifted_back.begin(), 50);
        assert_eq!(shifted_back.end(), 150);
    }

    #[test]
    fn shift_past_start_errors() {
        let region = r("chr1", 10, 20);
        assert!(shift(&region, -50).is_err());
    }

    #[test]
    fn expand_grows_both_sides() {
        let region = r("chr1", 100, 200);
        let expanded = expand(&region, 20).unwrap();
        assert_eq!(expanded.begin(), 80);
        assert_eq!(expanded.end(), 220);
    }

    #[test]
    fn head_and_tail_regions() {
        let region = r("chr1", 100, 200);
        let h = head_region(&region, 10);
        assert_eq!(h.begin(), 100);
        assert_eq!(h.end(), 110);
        let t = tail_region(&region, 10);
        assert_eq!(t.begin(), 190);
        assert_eq!(t.end(), 200);
    }

    #[test]
    fn inner_and_outer_distance() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 250, 350);
        assert_eq!(inner_distance(&a, &b).unwrap(), 50);
        let overlapping = r("chr1", 150, 250);
        assert_eq!(inner_distance(&a, &overlapping).unwrap(), 0);
    }

    #[test]
    fn partial_order_same_contig() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 200, 300);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn partial_order_different_contig_is_none() {
        let a = r("chr1", 100, 200);
        let b = r("chr2", 100, 200);
        assert_eq!(a.partial_cmp(&b), None);
    }

    #[test]
    fn equality() {
        let a = r("chr1", 100, 200);
        let b = r("chr1", 100, 200);
        let c = r("chr1", 100, 201);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn left_and_right_overhang() {
        let a = r("chr1", 100, 300);
        let b = r("chr1", 150, 350);
        assert_eq!(left_overhang_size(&a, &b), 50);
        assert_eq!(right_overhang_size(&a, &b), 0);
        assert_eq!(right_overhang_size(&b, &a), 50);
    }
}
