// Converted from C++ to Rust

use std::collections::HashMap;
use std::marker::PhantomData;

pub struct ProbabilityMatrix<R, C> {
    data: Vec<Vec<f64>>,
    row_index: HashMap<String, usize>,
    col_index: HashMap<String, usize>,
    row_keys: Vec<String>,
    col_keys: Vec<String>,
    _phantom: PhantomData<(R, C)>,
}

impl<R: std::fmt::Display + std::hash::Hash + Eq + Clone,
     C: std::fmt::Display + std::hash::Hash + Eq + Clone>
ProbabilityMatrix<R, C>
{
    pub fn new() -> Self {
        ProbabilityMatrix {
            data: Vec::new(),
            row_index: HashMap::new(),
            col_index: HashMap::new(),
            row_keys: Vec::new(),
            col_keys: Vec::new(),
            _phantom: PhantomData,
        }
    }

    pub fn set_rows(&mut self, rows: &[R]) {
        self.row_index.clear();
        self.row_keys.clear();
        for (i, r) in rows.iter().enumerate() {
            let key = r.to_string();
            self.row_index.insert(key.clone(), i);
            self.row_keys.push(key);
        }
        self.resize();
    }

    pub fn set_cols(&mut self, cols: &[C]) {
        self.col_index.clear();
        self.col_keys.clear();
        for (i, c) in cols.iter().enumerate() {
            let key = c.to_string();
            self.col_index.insert(key.clone(), i);
            self.col_keys.push(key);
        }
        self.resize();
    }

    fn resize(&mut self) {
        let nrows = self.row_keys.len();
        let ncols = self.col_keys.len();
        self.data = vec![vec![0.0; ncols]; nrows];
    }

    pub fn set(&mut self, row: &R, col: &C, value: f64) {
        if let (Some(&ri), Some(&ci)) = (self.row_index.get(&row.to_string()), self.col_index.get(&col.to_string())) {
            self.data[ri][ci] = value;
        }
    }

    pub fn get(&self, row: &R, col: &C) -> Option<f64> {
        let ri = *self.row_index.get(&row.to_string())?;
        let ci = *self.col_index.get(&col.to_string())?;
        Some(self.data[ri][ci])
    }

    pub fn num_rows(&self) -> usize { self.row_keys.len() }
    pub fn num_cols(&self) -> usize { self.col_keys.len() }
}

impl<R: std::fmt::Display + std::hash::Hash + Eq + Clone,
     C: std::fmt::Display + std::hash::Hash + Eq + Clone>
Default for ProbabilityMatrix<R, C>
{
    fn default() -> Self { Self::new() }
}
