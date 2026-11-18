//! B-Baum mit genau dem Verhalten von klib `kbtree.h` (MIT-Lizenz), soweit bwa ihn beim
//! Verketten benutzt: einfügen, "größter Schlüssel kleiner oder gleich" suchen, der Reihe nach
//! durchlaufen.
//!
//! Der Nachbau ist nötig, weil beim Verketten viele Ketten denselben Schlüssel haben (mehrere
//! Seeds verschiedener Länge an derselben Stelle). Welche davon die Suche liefert, hängt bei
//! klib vom Aufbau des Baums ab: Ein gleicher Schlüssel in einem inneren Knoten wird vor denen
//! in den Blättern gefunden. Eine gewöhnliche sortierte Liste liefert in diesen Fällen eine
//! andere Kette und am Ende andere Ketten als bwa.

/// Mindestgrad. bwa legt Knoten von 512 Byte an; mit 40 Byte je Kette (64-Bit-Systeme) ergibt
/// das `((512 - 4 - 8) / (8 + 40) + 1) >> 1 = 5`, also höchstens 9 Schlüssel je Knoten.
const T: usize = 5;
const MAX_KEYS: usize = 2 * T - 1;

struct Node<V> {
    keys: Vec<(u64, V)>,
    /// Leer bei Blättern, sonst ein Kind mehr als Schlüssel.
    children: Vec<usize>,
}

pub struct KbTree<V> {
    nodes: Vec<Node<V>>,
    root: usize,
    len: usize,
}

impl<V> Default for KbTree<V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<V> KbTree<V> {
    pub fn new() -> Self {
        Self {
            nodes: vec![Node {
                keys: Vec::new(),
                children: Vec::new(),
            }],
            root: 0,
            len: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// `__kb_getp_aux`: Stelle im Knoten und Vergleich. Liefert (Index, Vergleich von `k` mit
    /// dem Schlüssel dort); Index -1, wenn alle Schlüssel größer sind.
    fn search(node: &Node<V>, k: u64) -> (isize, i32) {
        let n = node.keys.len();
        if n == 0 {
            return (-1, 0);
        }
        let (mut begin, mut end) = (0usize, n);
        while begin < end {
            let mid = (begin + end) >> 1;
            if node.keys[mid].0 < k {
                begin = mid + 1;
            } else {
                end = mid;
            }
        }
        if begin == n {
            return (n as isize - 1, 1);
        }
        let key = node.keys[begin].0;
        let r = i32::from(key < k) - i32::from(k < key);
        if r < 0 {
            (begin as isize - 1, r)
        } else {
            (begin as isize, r)
        }
    }

    /// Der Eintrag mit dem größten Schlüssel kleiner oder gleich `k`, so wie ihn
    /// `kb_intervalp` als `lower` liefert.
    pub fn lower_mut(&mut self, k: u64) -> Option<&mut V> {
        let mut found: Option<(usize, usize)> = None;
        let mut x = self.root;
        loop {
            let node = &self.nodes[x];
            let (i, r) = Self::search(node, k);
            if i >= 0 {
                found = Some((x, i as usize));
                if r == 0 {
                    break;
                }
            }
            if node.children.is_empty() {
                break;
            }
            x = node.children[(i + 1) as usize];
        }
        found.map(|(x, i)| &mut self.nodes[x].keys[i].1)
    }

    /// `__kb_split`: teilt das volle Kind `i` von `x`.
    fn split(&mut self, x: usize, i: usize) {
        let y = self.nodes[x].children[i];
        let z_keys = self.nodes[y].keys.split_off(T);
        let z_children = if self.nodes[y].children.is_empty() {
            Vec::new()
        } else {
            self.nodes[y].children.split_off(T)
        };
        let median = self.nodes[y].keys.pop().expect("volles Kind");
        let z = self.nodes.len();
        self.nodes.push(Node {
            keys: z_keys,
            children: z_children,
        });
        self.nodes[x].children.insert(i + 1, z);
        self.nodes[x].keys.insert(i, median);
    }

    /// `kb_putp`: fügt ein, auch wenn der Schlüssel schon vorkommt.
    pub fn insert(&mut self, k: u64, value: V) {
        self.len += 1;
        if self.nodes[self.root].keys.len() == MAX_KEYS {
            let s = self.nodes.len();
            self.nodes.push(Node {
                keys: Vec::new(),
                children: vec![self.root],
            });
            self.root = s;
            self.split(s, 0);
        }
        let mut x = self.root;
        loop {
            let (i, _) = Self::search(&self.nodes[x], k);
            if self.nodes[x].children.is_empty() {
                self.nodes[x].keys.insert((i + 1) as usize, (k, value));
                return;
            }
            let mut c = (i + 1) as usize;
            if self.nodes[self.nodes[x].children[c]].keys.len() == MAX_KEYS {
                self.split(x, c);
                if k > self.nodes[x].keys[c].0 {
                    c += 1;
                }
            }
            x = self.nodes[x].children[c];
        }
    }

    /// Alle Einträge der Reihe nach (`__kb_traverse`).
    pub fn into_values(mut self) -> Vec<V> {
        fn walk<V>(nodes: &mut [Node<V>], x: usize, out: &mut Vec<V>) {
            let keys = std::mem::take(&mut nodes[x].keys);
            let children = std::mem::take(&mut nodes[x].children);
            for (i, (_, v)) in keys.into_iter().enumerate() {
                if let Some(&c) = children.get(i) {
                    walk(nodes, c, out);
                }
                out.push(v);
            }
            if let Some(&c) = children.last() {
                walk(nodes, c, out);
            }
        }
        let mut out = Vec::with_capacity(self.len);
        let root = self.root;
        walk(&mut self.nodes, root, &mut out);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_all_entries_in_key_order_and_finds_the_largest_key_not_above() {
        let mut state = 0x1234_5678_9ABC_DEF1u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut tree = KbTree::new();
        let mut keys = Vec::new();
        for n in 0..3000u32 {
            let k = next() % 500; // viele gleiche Schlüssel
            if n % 3 == 0 {
                let want = keys.iter().copied().filter(|&x| x <= k).max();
                let got = tree.lower_mut(k).map(|v: &mut (u64, u32)| v.0);
                assert_eq!(got, want);
            }
            tree.insert(k, (k, n));
            keys.push(k);
        }
        assert_eq!(tree.len(), 3000);
        let values = tree.into_values();
        assert_eq!(values.len(), 3000);
        assert!(values.windows(2).all(|w| w[0].0 <= w[1].0));
        let mut ids: Vec<u32> = values.iter().map(|v| v.1).collect();
        ids.sort_unstable();
        assert_eq!(ids, (0..3000).collect::<Vec<_>>());
    }

    /// Bei gleichen Schlüsseln in einem Blatt liefert die Suche den ersten, und ein weiterer
    /// gleicher wird direkt hinter dem ersten eingefügt, nicht am Ende.
    #[test]
    fn equal_keys_follow_the_order_of_klib() {
        let mut tree = KbTree::new();
        for name in ['a', 'b', 'c'] {
            tree.insert(7, name);
        }
        assert_eq!(tree.lower_mut(7).copied(), Some('a'));
        // größerer Suchschlüssel: der letzte im Blatt, und das ist b (Reihenfolge a, c, b)
        assert_eq!(tree.lower_mut(9).copied(), Some('b'));
        assert_eq!(tree.into_values(), vec!['a', 'c', 'b']);
    }
}
