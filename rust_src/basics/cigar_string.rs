// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};

pub type CigarSize = u32;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CigarFlag {
    AlignmentMatch = b'M',
    SequenceMatch  = b'=',
    Substitution   = b'X',
    Insertion      = b'I',
    Deletion       = b'D',
    SoftClipped    = b'S',
    HardClipped    = b'H',
    Padding        = b'P',
    Skipped        = b'N',
}

impl CigarFlag {
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            'M' => Some(CigarFlag::AlignmentMatch),
            '=' => Some(CigarFlag::SequenceMatch),
            'X' => Some(CigarFlag::Substitution),
            'I' => Some(CigarFlag::Insertion),
            'D' => Some(CigarFlag::Deletion),
            'S' => Some(CigarFlag::SoftClipped),
            'H' => Some(CigarFlag::HardClipped),
            'P' => Some(CigarFlag::Padding),
            'N' => Some(CigarFlag::Skipped),
            _   => None,
        }
    }

    pub fn as_char(self) -> char {
        self as u8 as char
    }
}

impl fmt::Display for CigarFlag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

impl PartialOrd for CigarFlag {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        (*self as u8).partial_cmp(&(*other as u8))
    }
}

impl Ord for CigarFlag {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (*self as u8).cmp(&(*other as u8))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOperation {
    size: CigarSize,
    flag: CigarFlag,
}

impl CigarOperation {
    pub fn new(size: CigarSize, flag: CigarFlag) -> Self {
        CigarOperation { size, flag }
    }

    pub fn flag(&self) -> CigarFlag {
        self.flag
    }

    pub fn size(&self) -> CigarSize {
        self.size
    }

    pub fn set_flag(&mut self, flag: CigarFlag) {
        self.flag = flag;
    }

    pub fn set_size(&mut self, size: CigarSize) {
        self.size = size;
    }
}

impl PartialOrd for CigarOperation {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CigarOperation {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.flag == other.flag {
            self.size.cmp(&other.size)
        } else {
            self.flag.cmp(&other.flag)
        }
    }
}

impl Hash for CigarOperation {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (self.flag as u8).hash(state);
        self.size.hash(state);
    }
}

impl fmt::Display for CigarOperation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.size, self.flag)
    }
}

pub fn increment_size(op: &mut CigarOperation, n: CigarSize) {
    op.set_size(op.size() + n);
}

pub fn decrement_size(op: &mut CigarOperation, n: CigarSize) {
    op.set_size(op.size() - n);
}

pub fn advances_reference_flag(flag: CigarFlag) -> bool {
    !matches!(flag, CigarFlag::Insertion | CigarFlag::HardClipped | CigarFlag::Padding)
}

pub fn advances_reference(op: &CigarOperation) -> bool {
    advances_reference_flag(op.flag())
}

pub fn advances_sequence_flag(flag: CigarFlag) -> bool {
    !matches!(flag, CigarFlag::Deletion | CigarFlag::HardClipped)
}

pub fn advances_sequence(op: &CigarOperation) -> bool {
    advances_sequence_flag(op.flag())
}

pub fn is_match_flag(flag: CigarFlag) -> bool {
    matches!(flag, CigarFlag::AlignmentMatch | CigarFlag::SequenceMatch)
}

pub fn is_match(op: &CigarOperation) -> bool {
    is_match_flag(op.flag())
}

pub fn is_substitution_flag(flag: CigarFlag) -> bool {
    flag == CigarFlag::Substitution
}

pub fn is_substitution(op: &CigarOperation) -> bool {
    is_substitution_flag(op.flag())
}

pub fn is_match_or_substitution_flag(flag: CigarFlag) -> bool {
    is_match_flag(flag) || is_substitution_flag(flag)
}

pub fn is_match_or_substitution(op: &CigarOperation) -> bool {
    is_match_or_substitution_flag(op.flag())
}

pub fn is_insertion_flag(flag: CigarFlag) -> bool {
    flag == CigarFlag::Insertion
}

pub fn is_insertion(op: &CigarOperation) -> bool {
    is_insertion_flag(op.flag())
}

pub fn is_deletion_flag(flag: CigarFlag) -> bool {
    flag == CigarFlag::Deletion
}

pub fn is_deletion(op: &CigarOperation) -> bool {
    is_deletion_flag(op.flag())
}

pub fn is_indel_flag(flag: CigarFlag) -> bool {
    is_insertion_flag(flag) || is_deletion_flag(flag)
}

pub fn is_indel(op: &CigarOperation) -> bool {
    is_indel_flag(op.flag())
}

pub fn is_clipping_flag(flag: CigarFlag) -> bool {
    matches!(flag, CigarFlag::SoftClipped | CigarFlag::HardClipped)
}

pub fn is_clipping(op: &CigarOperation) -> bool {
    is_clipping_flag(op.flag())
}

pub type CigarString = Vec<CigarOperation>;

pub fn parse_cigar(cigar: &str) -> Result<CigarString, String> {
    let mut result = Vec::with_capacity(cigar.len() / 2);
    let mut digits = String::with_capacity(3);

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            let size: CigarSize = digits.parse()
                .map_err(|_| format!("parse_cigar: invalid size in {}", cigar))?;
            let flag = CigarFlag::from_char(c)
                .ok_or_else(|| format!("parse_cigar: unknown operation '{}' in {}", c, cigar))?;
            result.push(CigarOperation::new(size, flag));
            digits.clear();
        }
    }

    if !digits.is_empty() {
        return Err(format!("parse_cigar: unparsed characters in {}", cigar));
    }

    Ok(result)
}

pub fn is_valid_cigar(cigar: &CigarString) -> bool {
    !cigar.is_empty() && cigar.iter().all(|op| op.size() > 0)
}

pub fn is_minimal(cigar: &CigarString) -> bool {
    cigar.windows(2).all(|w| w[0].flag() != w[1].flag())
}

pub fn is_front_soft_clipped(cigar: &CigarString) -> bool {
    cigar.first().map_or(false, |op| op.flag() == CigarFlag::SoftClipped)
}

pub fn is_back_soft_clipped(cigar: &CigarString) -> bool {
    cigar.last().map_or(false, |op| op.flag() == CigarFlag::SoftClipped)
}

pub fn is_soft_clipped(cigar: &CigarString) -> bool {
    is_front_soft_clipped(cigar) || is_back_soft_clipped(cigar)
}

pub fn sum_matches(cigar: &CigarString) -> i32 {
    cigar.iter()
        .map(|op| if is_match(op) { op.size() as i32 } else { 0 })
        .sum()
}

pub fn sum_non_matches(cigar: &CigarString) -> i32 {
    cigar.iter()
        .map(|op| if is_match(op) { 0 } else { op.size() as i32 })
        .sum()
}

pub fn has_indel(cigar: &CigarString) -> bool {
    cigar.iter().any(|op| is_indel(op))
}

fn indel_size(op: &CigarOperation) -> i32 {
    if is_indel(op) {
        if is_insertion(op) { op.size() as i32 } else { -(op.size() as i32) }
    } else {
        0
    }
}

pub fn sum_indel_sizes(cigar: &CigarString) -> i32 {
    cigar.iter().map(|op| indel_size(op)).sum()
}

pub fn max_indel_size(cigar: &CigarString) -> i32 {
    cigar.iter()
        .map(|op| indel_size(op))
        .max_by_key(|&x| x.abs())
        .unwrap_or(0)
}

pub fn get_soft_clipped_sizes(cigar: &CigarString) -> (CigarSize, CigarSize) {
    let front = if is_front_soft_clipped(cigar) { cigar[0].size() } else { 0 };
    let back = if cigar.len() > 1 && is_back_soft_clipped(cigar) { cigar.last().unwrap().size() } else { 0 };
    (front, back)
}

pub fn sum_operation_sizes(cigar: &CigarString) -> CigarSize {
    cigar.iter().map(|op| op.size()).sum()
}

pub fn reference_size(cigar: &CigarString) -> CigarSize {
    cigar.iter()
        .map(|op| if advances_reference(op) { op.size() } else { 0 })
        .sum()
}

pub fn sequence_size(cigar: &CigarString) -> CigarSize {
    cigar.iter()
        .map(|op| if advances_sequence(op) { op.size() } else { 0 })
        .sum()
}

pub fn clipped_begin(cigar: &CigarString, unclipped_begin: CigarSize) -> CigarSize {
    if is_front_soft_clipped(cigar) {
        unclipped_begin.saturating_sub(cigar[0].size())
    } else {
        unclipped_begin
    }
}

pub fn clipped_end(cigar: &CigarString, unclipped_end: CigarSize) -> CigarSize {
    if is_back_soft_clipped(cigar) {
        unclipped_end + cigar.last().unwrap().size()
    } else {
        unclipped_end
    }
}

pub fn decompose(cigar: &CigarString) -> Vec<CigarFlag> {
    let mut result = Vec::with_capacity(sum_operation_sizes(cigar) as usize);
    for op in cigar {
        for _ in 0..op.size() {
            result.push(op.flag());
        }
    }
    result
}

pub fn collapse_matches(cigar: &CigarString) -> CigarString {
    let mut result = Vec::with_capacity(cigar.len());
    let mut i = 0;
    while i < cigar.len() {
        if is_match_or_substitution(&cigar[i]) {
            let start = i;
            let mut total_size = 0u32;
            while i < cigar.len() && is_match_or_substitution(&cigar[i]) {
                total_size += cigar[i].size();
                i += 1;
            }
            result.push(CigarOperation::new(total_size, CigarFlag::AlignmentMatch));
            let _ = start;
        } else {
            result.push(cigar[i]);
            i += 1;
        }
    }
    result
}

pub fn cigar_to_string(cigar: &CigarString) -> String {
    cigar.iter().map(|op| op.to_string()).collect()
}
