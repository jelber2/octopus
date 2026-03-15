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
