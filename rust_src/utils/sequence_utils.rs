// Converted from C++ to Rust

pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        b'N' | b'n' => b'N',
        _ => b'N',
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

pub fn is_valid_dna_base(base: u8) -> bool {
    matches!(base, b'A' | b'T' | b'G' | b'C' | b'N' | b'a' | b't' | b'g' | b'c' | b'n')
}

pub fn is_ambiguous(base: u8) -> bool {
    matches!(base, b'N' | b'n')
}

pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() { return 0.0; }
    let gc = seq.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count();
    gc as f64 / seq.len() as f64
}

pub fn count_mismatches(seq1: &[u8], seq2: &[u8]) -> usize {
    seq1.iter().zip(seq2.iter()).filter(|(&a, &b)| a != b).count()
}
