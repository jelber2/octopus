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
