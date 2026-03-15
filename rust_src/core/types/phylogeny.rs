// Converted from C++ to Rust

use std::collections::HashMap;

pub type CloneId = usize;

#[derive(Debug, Clone)]
pub struct Phylogeny {
    parent: HashMap<CloneId, CloneId>,
    children: HashMap<CloneId, Vec<CloneId>>,
    roots: Vec<CloneId>,
    next_id: CloneId,
}

impl Phylogeny {
    pub fn new() -> Self {
        Phylogeny {
            parent: HashMap::new(),
            children: HashMap::new(),
            roots: Vec::new(),
            next_id: 0,
        }
    }

    pub fn add_node(&mut self) -> CloneId {
        let id = self.next_id;
        self.next_id += 1;
        self.children.insert(id, Vec::new());
        self.roots.push(id);
        id
    }

    pub fn add_edge(&mut self, parent: CloneId, child: CloneId) {
        self.parent.insert(child, parent);
        self.children.entry(parent).or_default().push(child);
        self.roots.retain(|&r| r != child);
    }

    pub fn parent(&self, node: CloneId) -> Option<CloneId> {
        self.parent.get(&node).copied()
    }

    pub fn children(&self, node: CloneId) -> &[CloneId] {
        self.children.get(&node).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn roots(&self) -> &[CloneId] { &self.roots }
    pub fn num_nodes(&self) -> usize { self.next_id }
}

impl Default for Phylogeny {
    fn default() -> Self { Self::new() }
}
