// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};

pub type Position = u32;
pub type Size = u32;
pub type Distance = i64;

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub struct ContigRegion {
    begin: Position,
    end: Position,
}

#[derive(Debug, Clone)]
pub struct BadRegion {
    begin: Position,
    end: Position,
}

impl fmt::Display for BadRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "BadRegion: begin={} end={}", self.begin, self.end)
    }
}

impl std::error::Error for BadRegion {}

impl ContigRegion {
    pub fn new(begin: Position, end: Position) -> Result<Self, BadRegion> {
        if end < begin {
            Err(BadRegion { begin, end })
        } else {
            Ok(ContigRegion { begin, end })
        }
    }

    pub fn begin(&self) -> Position {
        self.begin
    }

    pub fn end(&self) -> Position {
        self.end
    }
}

impl PartialOrd for ContigRegion {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ContigRegion {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.begin.cmp(&other.begin)
            .then(self.end.cmp(&other.end))
    }
}

impl Hash for ContigRegion {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.begin.hash(state);
        self.end.hash(state);
    }
}

impl fmt::Display for ContigRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}", self.begin, self.end)
    }
}

pub fn mapped_begin(region: &ContigRegion) -> Position {
    region.begin()
}

pub fn mapped_end(region: &ContigRegion) -> Position {
    region.end()
}

pub fn is_empty(region: &ContigRegion) -> bool {
    region.begin() == region.end()
}

pub fn size(region: &ContigRegion) -> Size {
    region.end() - region.begin()
}

pub fn is_position(region: &ContigRegion) -> bool {
    size(region) == 1
}

pub fn begins_equal(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.begin() == rhs.begin()
}

pub fn ends_equal(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.end() == rhs.end()
}

pub fn begins_before(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.begin() < rhs.begin()
}

pub fn ends_before(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.end() < rhs.end()
}

pub fn is_before(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.end() <= rhs.begin() && lhs != rhs
}

pub fn is_after(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    rhs.end() <= lhs.begin() && lhs != rhs
}

pub fn are_adjacent(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.begin() == rhs.end() || lhs.end() == rhs.begin()
}

pub fn overlap_size(lhs: &ContigRegion, rhs: &ContigRegion) -> Distance {
    (std::cmp::min(lhs.end(), rhs.end()) as Distance)
        - (std::cmp::max(lhs.begin(), rhs.begin()) as Distance)
}

pub fn overlaps(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    let overlapped = overlap_size(lhs, rhs);
    overlapped > 0 || (overlapped == 0 && (is_empty(lhs) || is_empty(rhs)))
}

pub fn contains(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    lhs.begin() <= rhs.begin() && rhs.end() <= lhs.end()
}

pub fn inner_distance(lhs: &ContigRegion, rhs: &ContigRegion) -> Distance {
    if overlaps(lhs, rhs) {
        0
    } else if begins_before(lhs, rhs) {
        (rhs.begin() - lhs.end()) as Distance
    } else {
        -((lhs.begin() - rhs.end()) as Distance)
    }
}

pub fn outer_distance(lhs: &ContigRegion, rhs: &ContigRegion) -> Distance {
    if contains(lhs, rhs) || contains(rhs, lhs) {
        return 0;
    }
    rhs.end() as Distance - lhs.begin() as Distance
}

pub fn shift(region: &ContigRegion, n: Distance) -> Result<ContigRegion, String> {
    if n < 0 && ((-n) as u32) > region.begin() {
        return Err("ContigRegion: shifted past contig start".to_string());
    }
    ContigRegion::new(
        (region.begin() as Distance + n) as Position,
        (region.end() as Distance + n) as Position,
    ).map_err(|e| e.to_string())
}

pub fn next_position(region: &ContigRegion) -> Result<ContigRegion, BadRegion> {
    ContigRegion::new(region.end(), region.end() + 1)
}

pub fn expand_lhs(region: &ContigRegion, n: Distance) -> Result<ContigRegion, String> {
    if n > 0 && (n as u32) > region.begin() {
        return Err("ContigRegion: expanded past contig start".to_string());
    }
    ContigRegion::new(
        (region.begin() as Distance - n) as Position,
        region.end(),
    ).map_err(|e| e.to_string())
}

pub fn expand_rhs(region: &ContigRegion, n: Distance) -> Result<ContigRegion, String> {
    if n < 0 && ((-n) as u32) > region.end() {
        return Err("ContigRegion: compressed past contig start".to_string());
    }
    ContigRegion::new(
        region.begin(),
        (region.end() as Distance + n) as Position,
    ).map_err(|e| e.to_string())
}

pub fn expand(region: &ContigRegion, n: Distance) -> Result<ContigRegion, BadRegion> {
    ContigRegion::new(
        std::cmp::max(0i64, region.begin() as Distance - n) as Position,
        (region.end() as Distance + n) as Position,
    )
}

pub fn expand2(region: &ContigRegion, lhs: Distance, rhs: Distance) -> Result<ContigRegion, String> {
    let r = expand_rhs(region, rhs)?;
    expand_lhs(&r, lhs)
}

pub fn overlapped_region(lhs: &ContigRegion, rhs: &ContigRegion) -> Option<ContigRegion> {
    if !overlaps(lhs, rhs) {
        return None;
    }
    ContigRegion::new(
        std::cmp::max(lhs.begin(), rhs.begin()),
        std::cmp::min(lhs.end(), rhs.end()),
    ).ok()
}

pub fn encompassing_region(lhs: &ContigRegion, rhs: &ContigRegion) -> Result<ContigRegion, BadRegion> {
    ContigRegion::new(
        std::cmp::min(lhs.begin(), rhs.begin()),
        std::cmp::max(lhs.end(), rhs.end()),
    )
}

pub fn intervening_region(lhs: &ContigRegion, rhs: &ContigRegion) -> Option<ContigRegion> {
    if begins_before(rhs, lhs) || overlaps(lhs, rhs) {
        return None;
    }
    ContigRegion::new(lhs.end(), rhs.begin()).ok()
}

pub fn intervening_region_size(lhs: &ContigRegion, rhs: &ContigRegion) -> Size {
    if begins_before(rhs, lhs) || overlaps(lhs, rhs) {
        0
    } else {
        rhs.begin() - lhs.end()
    }
}

pub fn left_overhang_size(lhs: &ContigRegion, rhs: &ContigRegion) -> Size {
    if begins_before(lhs, rhs) { rhs.begin() - lhs.begin() } else { 0 }
}

pub fn right_overhang_size(lhs: &ContigRegion, rhs: &ContigRegion) -> Size {
    if ends_before(lhs, rhs) { 0 } else { lhs.end() - rhs.end() }
}

pub fn left_overhangs(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    left_overhang_size(lhs, rhs) > 0
}

pub fn right_overhangs(lhs: &ContigRegion, rhs: &ContigRegion) -> bool {
    right_overhang_size(lhs, rhs) > 0
}

pub fn left_overhang_region(lhs: &ContigRegion, rhs: &ContigRegion) -> ContigRegion {
    if begins_before(rhs, lhs) {
        ContigRegion::new(lhs.begin(), lhs.begin()).unwrap()
    } else {
        ContigRegion::new(lhs.begin(), rhs.begin()).unwrap()
    }
}

pub fn right_overhang_region(lhs: &ContigRegion, rhs: &ContigRegion) -> ContigRegion {
    if ends_before(lhs, rhs) {
        ContigRegion::new(lhs.end(), lhs.end()).unwrap()
    } else {
        ContigRegion::new(rhs.end(), lhs.end()).unwrap()
    }
}

pub fn closed_region(lhs: &ContigRegion, rhs: &ContigRegion) -> Result<ContigRegion, BadRegion> {
    ContigRegion::new(lhs.begin(), rhs.end())
}

pub fn head_region(region: &ContigRegion, n: Size) -> ContigRegion {
    let begin = region.begin();
    ContigRegion::new(begin, std::cmp::min(begin + n, region.end())).unwrap()
}

pub fn head_position(region: &ContigRegion) -> ContigRegion {
    head_region(region, 1)
}

pub fn tail_region(region: &ContigRegion, n: Size) -> ContigRegion {
    let end = region.end();
    ContigRegion::new(if end >= n { end - n } else { 0 }, end).unwrap()
}

pub fn tail_position(region: &ContigRegion) -> ContigRegion {
    tail_region(region, 1)
}

pub fn begin_distance(first: &ContigRegion, second: &ContigRegion) -> Distance {
    second.begin() as Distance - first.begin() as Distance
}

pub fn end_distance(first: &ContigRegion, second: &ContigRegion) -> Distance {
    second.end() as Distance - first.end() as Distance
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(begin: u32, end: u32) -> ContigRegion {
        ContigRegion::new(begin, end).unwrap()
    }

    // ── Construction ─────────────────────────────────────────────────────

    #[test]
    fn constructing_negative_region_is_error() {
        assert!(ContigRegion::new(0, 0).is_ok());
        assert!(ContigRegion::new(0, 1).is_ok());
        assert!(ContigRegion::new(1, 0).is_err());
        assert!(ContigRegion::new(5, 3).is_err());
    }

    // ── Ordering ─────────────────────────────────────────────────────────

    #[test]
    fn ordering_is_by_begin_then_end() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(1, 1);
        let r4 = r(0, 2);

        assert_ne!(r1, r2);
        assert!(r1 < r2);

        assert_ne!(r2, r3);
        assert!(r2 < r3);

        assert_ne!(r1, r4);
        assert!(r1 < r4);

        assert_ne!(r2, r4);
        assert!(r2 < r4);

        assert_ne!(r3, r4);
        assert!(r4 < r3);
    }

    // ── is_before ────────────────────────────────────────────────────────

    #[test]
    fn is_before_is_consistent() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(1, 1);
        let r4 = r(0, 2);
        let r5 = r(2, 2);

        assert!(!is_before(&r1, &r1));
        assert!(!is_before(&r2, &r2));

        assert!(is_before(&r1, &r2));
        assert!(!is_before(&r2, &r1));

        assert!(is_before(&r1, &r3));
        assert!(!is_before(&r3, &r1));

        assert!(is_before(&r1, &r4));
        assert!(!is_before(&r4, &r1));

        assert!(is_before(&r4, &r5));
        assert!(!is_before(&r5, &r4));

        assert!(!is_before(&r3, &r4));
        assert!(!is_before(&r4, &r3));
    }

    // ── is_after ─────────────────────────────────────────────────────────

    #[test]
    fn is_after_is_consistent() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(1, 1);
        let r4 = r(0, 2);
        let r5 = r(2, 2);

        assert!(!is_after(&r1, &r1));
        assert!(!is_after(&r2, &r2));

        assert!(is_after(&r2, &r1));
        assert!(!is_after(&r1, &r2));

        assert!(is_after(&r3, &r1));
        assert!(!is_after(&r1, &r3));

        assert!(is_after(&r4, &r1));
        assert!(!is_after(&r1, &r4));

        assert!(is_after(&r5, &r2));
        assert!(!is_after(&r2, &r5));

        assert!(is_after(&r5, &r3));
        assert!(!is_after(&r3, &r5));

        assert!(!is_after(&r3, &r4));
        assert!(!is_after(&r4, &r3));
    }

    // ── overlap_size ─────────────────────────────────────────────────────

    #[test]
    fn overlap_size_returns_number_of_overlapped_positions() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(0, 2);
        let r4 = r(0, 4);

        assert_eq!(overlap_size(&r1, &r1), 0);
        assert_eq!(overlap_size(&r1, &r2), 0);
        assert_eq!(overlap_size(&r1, &r3), 0);
        assert_eq!(overlap_size(&r1, &r4), 0);
        assert_eq!(overlap_size(&r2, &r1), 0);
        assert_eq!(overlap_size(&r3, &r1), 0);
        assert_eq!(overlap_size(&r4, &r1), 0);

        assert_eq!(overlap_size(&r2, &r3), 1);
        assert_eq!(overlap_size(&r3, &r2), 1);

        assert_eq!(overlap_size(&r2, &r4), 1);
        assert_eq!(overlap_size(&r4, &r2), 1);

        assert_eq!(overlap_size(&r3, &r4), 2);
        assert_eq!(overlap_size(&r4, &r3), 2);
    }

    // ── overlaps ─────────────────────────────────────────────────────────

    #[test]
    fn overlaps_is_consistent() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(1, 1);
        let r4 = r(0, 2);
        let r5 = r(2, 2);

        // every region overlaps itself (including empty ones)
        assert!(overlaps(&r1, &r1));
        assert!(overlaps(&r2, &r2));
        assert!(overlaps(&r3, &r3));
        assert!(overlaps(&r4, &r4));
        assert!(overlaps(&r5, &r5));

        assert!(overlaps(&r1, &r2));
        assert!(overlaps(&r2, &r1));
        assert!(!overlaps(&r1, &r3));
        assert!(!overlaps(&r3, &r1));
        assert!(overlaps(&r2, &r3));
        assert!(overlaps(&r3, &r2));

        assert!(overlaps(&r1, &r4));
        assert!(overlaps(&r2, &r4));
        assert!(overlaps(&r3, &r4));
        assert!(overlaps(&r4, &r1));
        assert!(overlaps(&r4, &r2));
        assert!(overlaps(&r4, &r3));

        assert!(!overlaps(&r1, &r5));
        assert!(!overlaps(&r2, &r5));
        assert!(!overlaps(&r3, &r5));
        assert!(!overlaps(&r5, &r1));
        assert!(!overlaps(&r5, &r2));
        assert!(!overlaps(&r5, &r3));
    }

    // ── contains ─────────────────────────────────────────────────────────

    #[test]
    fn contains_is_consistent() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(1, 1);
        let r4 = r(0, 2);
        let r5 = r(2, 2);

        assert!(contains(&r1, &r1));
        assert!(contains(&r2, &r2));
        assert!(contains(&r3, &r3));
        assert!(contains(&r4, &r4));
        assert!(contains(&r5, &r5));

        assert!(contains(&r2, &r1));
        assert!(!contains(&r1, &r2));

        assert!(contains(&r2, &r3));
        assert!(!contains(&r3, &r2));

        assert!(contains(&r4, &r1));
        assert!(contains(&r4, &r2));
        assert!(contains(&r4, &r3));
        assert!(contains(&r4, &r5));
        assert!(!contains(&r1, &r4));
        assert!(!contains(&r2, &r4));
        assert!(!contains(&r3, &r4));
        assert!(!contains(&r5, &r4));
    }

    // ── are_adjacent ─────────────────────────────────────────────────────

    #[test]
    fn overlapping_empty_regions_are_adjacent() {
        let r1 = r(0, 0);
        let r2 = r(0, 1);
        let r3 = r(1, 1);
        let r4 = r(0, 2);
        let r5 = r(2, 2);

        // empty regions that overlap each other are adjacent
        assert!(are_adjacent(&r1, &r1));
        assert!(are_adjacent(&r3, &r3));
        assert!(are_adjacent(&r5, &r5));

        // non-empty regions are not self-adjacent
        assert!(!are_adjacent(&r2, &r2));
        assert!(!are_adjacent(&r4, &r4));
    }

    #[test]
    fn are_adjacent_between_distinct_regions() {
        let r1 = r(0, 5);
        let r2 = r(5, 10);
        let r3 = r(6, 10);
        let r4 = r(0, 6);

        assert!(are_adjacent(&r1, &r2));
        assert!(are_adjacent(&r2, &r1));
        assert!(!are_adjacent(&r1, &r3));
        assert!(!are_adjacent(&r1, &r4));
    }

    // ── size and is_empty ─────────────────────────────────────────────────

    #[test]
    fn size_and_is_empty() {
        assert_eq!(size(&r(0, 0)), 0);
        assert!(is_empty(&r(0, 0)));
        assert_eq!(size(&r(0, 1)), 1);
        assert!(!is_empty(&r(0, 1)));
        assert_eq!(size(&r(3, 10)), 7);
    }

    // ── inner_distance ────────────────────────────────────────────────────

    #[test]
    fn inner_distance_for_non_overlapping() {
        let r1 = r(0, 5);
        let r2 = r(8, 12);
        assert_eq!(inner_distance(&r1, &r2), 3);
        assert_eq!(inner_distance(&r2, &r1), -3);
    }

    #[test]
    fn inner_distance_for_overlapping() {
        let r1 = r(0, 10);
        let r2 = r(5, 15);
        assert_eq!(inner_distance(&r1, &r2), 0);
    }

    // ── head / tail regions ───────────────────────────────────────────────

    #[test]
    fn head_region_clips_to_region_end() {
        let region = r(10, 20);
        assert_eq!(head_region(&region, 5), r(10, 15));
        assert_eq!(head_region(&region, 15), r(10, 20)); // clamped
        assert_eq!(head_position(&region), r(10, 11));
    }

    #[test]
    fn tail_region_clips_to_region_begin() {
        let region = r(10, 20);
        assert_eq!(tail_region(&region, 5),  r(15, 20));
        // n=15 < end=20, so begin = end-n = 5 (not clamped to region.begin())
        assert_eq!(tail_region(&region, 15), r(5, 20));
        // n > end → clamp to 0
        assert_eq!(tail_region(&region, 25), r(0, 20));
        assert_eq!(tail_position(&region), r(19, 20));
    }
}
