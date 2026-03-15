// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Phred {
    score: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Probability(pub f64);

impl Phred {
    pub fn from_score(score: f64) -> Result<Self, String> {
        let score = if score == 0.0 { 0.0 } else { score };
        if score < 0.0 {
            Err(format!("Phred: negative score {}", score))
        } else {
            Ok(Phred { score })
        }
    }

    pub fn from_probability(error: Probability) -> Result<Self, String> {
        if error.0 < 0.0 {
            return Err(format!("Phred: negative error probability {}", error.0));
        }
        let clamped = error.0.max(f64::MIN_POSITIVE).min(1.0);
        let score = (-10.0 * clamped.log10()).abs();
        Ok(Phred { score })
    }

    pub fn score(&self) -> f64 {
        self.score
    }

    pub fn probability_true(&self) -> Probability {
        Probability(1.0 - self.probability_false().0)
    }

    pub fn probability_false(&self) -> Probability {
        Probability(10_f64.powf(-self.score / 10.0))
    }
}

impl Default for Phred {
    fn default() -> Self {
        Phred { score: 0.0 }
    }
}

impl fmt::Display for Phred {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.score)
    }
}

pub fn to_phred(probability: f64) -> f64 {
    (-10.0 * probability.max(f64::MIN_POSITIVE).log10()).abs()
}

pub fn to_probability(phred: f64) -> f64 {
    10_f64.powf(-phred / 10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < TOL
    }

    // ── Construction ─────────────────────────────────────────────────────

    #[test]
    fn phreds_must_be_non_negative() {
        assert!(Phred::from_score(-1.0).is_err());
        assert!(Phred::from_score(-0.001).is_err());
        assert!(Phred::from_score(0.0).is_ok());
        assert!(Phred::from_score(30.0).is_ok());
    }

    #[test]
    fn zero_score_is_valid() {
        let p = Phred::from_score(0.0).unwrap();
        assert_eq!(p.score(), 0.0);
    }

    // ── score round-trip ─────────────────────────────────────────────────

    #[test]
    fn phreds_can_be_converted_to_scores() {
        for x in 0..100u32 {
            let score = x as f64;
            let p = Phred::from_score(score).unwrap();
            assert!(approx_eq(p.score(), score),
                "score mismatch at {}: got {}", score, p.score());
        }
    }

    // ── probability conversion ────────────────────────────────────────────

    #[test]
    fn phred_10_is_error_probability_0_1() {
        let p = Phred::from_score(10.0).unwrap();
        assert!(approx_eq(p.probability_false().0, 0.1),
            "expected 0.1, got {}", p.probability_false().0);
    }

    #[test]
    fn phred_20_is_error_probability_0_01() {
        let p = Phred::from_score(20.0).unwrap();
        assert!(approx_eq(p.probability_false().0, 0.01),
            "expected 0.01, got {}", p.probability_false().0);
    }

    #[test]
    fn phred_30_is_error_probability_0_001() {
        let p = Phred::from_score(30.0).unwrap();
        assert!(approx_eq(p.probability_false().0, 0.001),
            "expected 0.001, got {}", p.probability_false().0);
    }

    #[test]
    fn probability_true_plus_false_equals_one() {
        for score in [0.0_f64, 10.0, 20.0, 30.0, 40.0] {
            let p = Phred::from_score(score).unwrap();
            let sum = p.probability_true().0 + p.probability_false().0;
            assert!(approx_eq(sum, 1.0),
                "sum was {} for score {}", sum, score);
        }
    }

    // ── from_probability round-trip ───────────────────────────────────────

    #[test]
    fn phred_can_be_constructed_from_probability() {
        // p=0.1 → phred=10, p=0.01 → phred=20, p=0.001 → phred=30
        let cases = [(0.1_f64, 10.0_f64), (0.01, 20.0), (0.001, 30.0)];
        for (prob, expected_score) in cases {
            let p = Phred::from_probability(Probability(prob)).unwrap();
            assert!(approx_eq(p.score(), expected_score),
                "from p={}: expected score {}, got {}", prob, expected_score, p.score());
        }
    }

    #[test]
    fn negative_probability_is_rejected() {
        assert!(Phred::from_probability(Probability(-0.1)).is_err());
    }

    // ── free functions ────────────────────────────────────────────────────

    #[test]
    fn to_phred_and_to_probability_are_inverses() {
        for score in [0.0_f64, 10.0, 20.0, 30.0, 40.0] {
            let prob = to_probability(score);
            let recovered = to_phred(prob);
            assert!(approx_eq(recovered, score),
                "round-trip failed for score {}: got {}", score, recovered);
        }
    }

    // ── ordering ─────────────────────────────────────────────────────────

    #[test]
    fn higher_phred_is_greater() {
        let low  = Phred::from_score(10.0).unwrap();
        let high = Phred::from_score(30.0).unwrap();
        assert!(high > low);
        assert!(low < high);
    }

    #[test]
    fn default_phred_is_zero() {
        let p = Phred::default();
        assert_eq!(p.score(), 0.0);
    }
}
