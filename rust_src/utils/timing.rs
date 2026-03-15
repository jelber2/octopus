// Converted from C++ to Rust

use std::time::{Duration, Instant};
use std::fmt;

pub struct TimeInterval {
    duration: Duration,
}

impl TimeInterval {
    pub fn new(start: Instant, end: Instant) -> Self {
        TimeInterval { duration: end.duration_since(start) }
    }

    pub fn duration(&self) -> Duration {
        self.duration
    }
}

impl fmt::Display for TimeInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let secs = self.duration.as_secs();
        let millis = self.duration.subsec_millis();
        if secs >= 3600 {
            write!(f, "{}h {}m {}s", secs / 3600, (secs % 3600) / 60, secs % 60)
        } else if secs >= 60 {
            write!(f, "{}m {}s", secs / 60, secs % 60)
        } else if secs > 0 {
            write!(f, "{}.{}s", secs, millis / 100)
        } else {
            write!(f, "{}ms", self.duration.as_millis())
        }
    }
}
