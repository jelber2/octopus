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
