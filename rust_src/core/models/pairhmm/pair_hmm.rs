// Converted from C++ to Rust
// Pair HMM for computing read-haplotype alignment likelihoods

/// Computes the log-likelihood of a read given a haplotype using a pair HMM.
/// This is the core alignment model used by Octopus.
pub fn compute_log_likelihood(
    read_sequence: &[u8],
    read_qualities: &[u8],
    haplotype_sequence: &[u8],
    gap_open: f64,
    gap_extend: f64,
) -> f64 {
    let n = read_sequence.len();
    let m = haplotype_sequence.len();

    if n == 0 { return 0.0; }
    if m == 0 { return f64::NEG_INFINITY; }

    let go = gap_open;
    let ge = gap_extend;

    const MATCH: usize = 0;
    const INSERT: usize = 1;
    const DELETE: usize = 2;

    let mut dp = vec![vec![vec![f64::NEG_INFINITY; 3]; m + 1]; n + 1];
    dp[0][0][MATCH] = 0.0;

    for j in 1..=m {
        dp[0][j][DELETE] = dp[0][j-1][MATCH] + go + ge;
    }

    for i in 1..=n {
        let base_i = read_sequence[i - 1];
        let qual_i = read_qualities.get(i - 1).copied().unwrap_or(30);
        let mismatch_prob = phred_to_prob(qual_i);

        for j in 1..=m {
            let base_j = haplotype_sequence[j - 1];
            let match_emit = if base_i == base_j { 1.0 - mismatch_prob } else { mismatch_prob / 3.0 };
            let log_match_emit = match_emit.max(1e-300).ln();

            let from_match = dp[i-1][j-1][MATCH];
            let from_ins = dp[i-1][j-1][INSERT];
            let from_del = dp[i-1][j-1][DELETE];
            let best_from = log_sum_exp3(from_match, from_ins, from_del);
            dp[i][j][MATCH] = best_from + log_match_emit;

            let ins_from_match = dp[i-1][j][MATCH] + go;
            let ins_from_ins = dp[i-1][j][INSERT] + ge;
            dp[i][j][INSERT] = ins_from_match.max(ins_from_ins);

            let del_from_match = dp[i][j-1][MATCH] + go;
            let del_from_del = dp[i][j-1][DELETE] + ge;
            dp[i][j][DELETE] = del_from_match.max(del_from_del);
        }
    }

    let mut result = f64::NEG_INFINITY;
    for j in 0..=m {
        let score = log_sum_exp3(dp[n][j][MATCH], dp[n][j][INSERT], dp[n][j][DELETE]);
        if score > result { result = score; }
    }
    result
}

fn phred_to_prob(phred: u8) -> f64 {
    10.0_f64.powf(-(phred as f64) / 10.0)
}

fn log_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    let max = a.max(b).max(c);
    if max.is_infinite() { return max; }
    max + ((a - max).exp() + (b - max).exp() + (c - max).exp()).ln()
}
