// Converted from C++ to Rust

pub fn split(s: &str, delim: char) -> Vec<String> {
    s.split(delim).map(|p| p.to_string()).collect()
}

pub fn split_multi(s: &str, delims: &str) -> Vec<String> {
    s.split(|c| delims.contains(c)).map(|p| p.to_string()).collect()
}

pub fn join(strings: &[String], delim: &str) -> String {
    strings.join(delim)
}

pub fn join_char(strings: &[String], delim: char) -> String {
    strings.join(&delim.to_string())
}

pub fn is_prefix(prefix: &str, text: &str) -> bool {
    text.starts_with(prefix)
}

pub fn is_suffix(suffix: &str, text: &str) -> bool {
    text.ends_with(suffix)
}

pub fn capitalise(s: &str) -> String {
    s.to_uppercase()
}

pub fn capitalise_front(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}

pub fn to_lower(s: &str) -> String {
    s.to_lowercase()
}

pub fn strip_leading_zeroes(s: &str) -> String {
    let stripped = s.trim_start_matches('0');
    if stripped.is_empty() { "0".to_string() } else { stripped.to_string() }
}

pub fn is_vowel(c: char) -> bool {
    matches!(c, 'a' | 'e' | 'i' | 'o' | 'u' | 'A' | 'E' | 'I' | 'O' | 'U')
}

pub fn begins_with_vowel(s: &str) -> bool {
    s.chars().next().map_or(false, is_vowel)
}

pub fn format_with_commas(value: i64) -> String {
    let s = value.abs().to_string();
    let chars: Vec<char> = s.chars().rev().collect();
    let with_commas: String = chars.chunks(3)
        .map(|c| c.iter().collect::<String>())
        .collect::<Vec<_>>()
        .join(",");
    let result: String = with_commas.chars().rev().collect();
    if value < 0 { format!("-{}", result) } else { result }
}
