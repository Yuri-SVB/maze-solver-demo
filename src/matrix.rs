//! Boolean reachability matrices with "or/and" matrix arithmetic.
//!
//! Over the boolean semiring, matrix `+` is element-wise OR and matrix `*` is
//! the usual product with `+` replaced by OR and `*` by AND:
//!
//! ```text
//! (A * B)[i][j] = OR over k of (A[i][k] AND B[k][j])
//! ```
//!
//! Every vertex can always "stay put", so the diagonal is always true. With
//! that, `A^d[i][j]` is true exactly when `j` is reachable from `i` in *at
//! most* `d` steps, for every `d`: `A^0` is the identity and `A^1` is `A`.

use std::fmt;
use std::ops::{Add, Mul};

/// A square boolean matrix: `get(i, j)` is true when `j` is reachable from
/// `i` in at most one step. Invariant: the diagonal is always true.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReachabilityMatrix {
    size: usize,
    rows: Vec<Vec<bool>>,
}

impl ReachabilityMatrix {
    /// A `size x size` matrix with no edges: each vertex reaches only
    /// itself. This is the identity.
    pub fn new(size: usize) -> Self {
        let mut rows = vec![vec![false; size]; size];
        for (i, row) in rows.iter_mut().enumerate() {
            row[i] = true;
        }
        ReachabilityMatrix { size, rows }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        self.rows[i][j]
    }

    /// Add or remove the edge `i -> j`. The diagonal can't be cleared.
    pub fn set(&mut self, i: usize, j: usize, value: bool) {
        assert!(i != j || value, "the diagonal of a reachability matrix is always true");
        self.rows[i][j] = value;
    }

    /// Add an undirected edge: `i -> j` and `j -> i`.
    pub fn connect(&mut self, i: usize, j: usize) {
        self.rows[i][j] = true;
        self.rows[j][i] = true;
    }

    /// `self^d`, by repeated squaring. `self^0` is the identity.
    pub fn pow(&self, mut d: u32) -> Self {
        let mut result = Self::new(self.size);
        let mut base = self.clone();
        while d > 0 {
            if d & 1 == 1 {
                result = &result * &base;
            }
            base = &base * &base;
            d >>= 1;
        }
        result
    }

    /// Shortest path from `a` to `b` as a list of vertices (both ends
    /// included), or `None` if `b` is unreachable.
    ///
    /// 1. Compute `A^0, A^1, A^2, ...` until the first power with `[a][b]`
    ///    set; that `d` is the shortest distance. Powers only ever gain
    ///    entries. If one equals the previous power (it's idempotent) with
    ///    `[a][b]` still false, nothing new can ever become reachable, so
    ///    there's no path. This happens within `size` steps.
    /// 2. Walk backwards. If the current vertex is at distance `s + 1`, its
    ///    predecessor is any `v` with `A[v][current]` and `A^s[a][v]` set:
    ///    `v` is at most `s` steps from `a`, so exactly `s`. The current
    ///    vertex itself never qualifies, because `A^s[a][current]` is false.
    pub fn shortest_path(&self, a: usize, b: usize) -> Option<Vec<usize>> {
        assert!(a < self.size && b < self.size, "vertex out of range");
        let mut powers = vec![Self::new(self.size)];
        loop {
            let last = powers.last().expect("powers starts non-empty");
            if last.get(a, b) {
                break;
            }
            // self is sparse, so multiply it on the left: the product skips
            // the left operand's false entries.
            let next = self * last;
            if &next == last {
                return None;
            }
            powers.push(next);
        }

        let distance = powers.len() - 1;
        let mut path = vec![b];
        let mut current = b;
        for s in (0..distance).rev() {
            current = (0..self.size)
                .find(|&v| powers[s].get(a, v) && self.get(v, current))
                .expect("a vertex at distance s + 1 has a neighbour at distance s");
            path.push(current);
        }
        path.reverse();
        Some(path)
    }
}

/// Boolean matrix product: OR over k of (self[i][k] AND rhs[k][j]). The
/// product of two matrices with a true diagonal has a true diagonal.
impl Mul for &ReachabilityMatrix {
    type Output = ReachabilityMatrix;

    fn mul(self, rhs: &ReachabilityMatrix) -> ReachabilityMatrix {
        assert_eq!(self.size, rhs.size, "matrix sizes differ");
        let mut out = ReachabilityMatrix::new(self.size);
        for (out_row, self_row) in out.rows.iter_mut().zip(&self.rows) {
            // Row i of the product is the OR of rhs's rows k for every
            // self[i][k] that is true.
            for (k, _) in self_row.iter().enumerate().filter(|(_, &edge)| edge) {
                for (o, &r) in out_row.iter_mut().zip(&rhs.rows[k]) {
                    *o |= r;
                }
            }
        }
        out
    }
}

impl Mul for ReachabilityMatrix {
    type Output = ReachabilityMatrix;

    fn mul(self, rhs: ReachabilityMatrix) -> ReachabilityMatrix {
        &self * &rhs
    }
}

/// Element-wise OR: the union of both edge sets.
impl Add for &ReachabilityMatrix {
    type Output = ReachabilityMatrix;

    fn add(self, rhs: &ReachabilityMatrix) -> ReachabilityMatrix {
        assert_eq!(self.size, rhs.size, "matrix sizes differ");
        let rows = self
            .rows
            .iter()
            .zip(&rhs.rows)
            .map(|(a, b)| a.iter().zip(b).map(|(&x, &y)| x || y).collect())
            .collect();
        ReachabilityMatrix { size: self.size, rows }
    }
}

impl Add for ReachabilityMatrix {
    type Output = ReachabilityMatrix;

    fn add(self, rhs: ReachabilityMatrix) -> ReachabilityMatrix {
        &self + &rhs
    }
}

/// Rows of 0s and 1s.
impl fmt::Display for ReachabilityMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in &self.rows {
            let line: String = row.iter().map(|&b| if b { '1' } else { '0' }).collect();
            writeln!(f, "{line}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::VecDeque;

    fn graph(size: usize, edges: &[(usize, usize)]) -> ReachabilityMatrix {
        let mut m = ReachabilityMatrix::new(size);
        for &(i, j) in edges {
            m.connect(i, j);
        }
        m
    }

    fn bfs_distance(m: &ReachabilityMatrix, a: usize, b: usize) -> Option<usize> {
        let mut dist = vec![None; m.size()];
        dist[a] = Some(0);
        let mut queue = VecDeque::from([a]);
        while let Some(u) = queue.pop_front() {
            for v in 0..m.size() {
                if m.get(u, v) && dist[v].is_none() {
                    dist[v] = Some(dist[u].unwrap() + 1);
                    queue.push_back(v);
                }
            }
        }
        dist[b]
    }

    /// It runs from a to b, and each step is a real move along an edge.
    fn is_path(m: &ReachabilityMatrix, path: &[usize], a: usize, b: usize) -> bool {
        path.first() == Some(&a)
            && path.last() == Some(&b)
            && path.windows(2).all(|w| w[0] != w[1] && m.get(w[0], w[1]))
    }

    #[test]
    fn new_is_the_identity() {
        let m = ReachabilityMatrix::new(3);
        assert_eq!(m.to_string(), "100\n010\n001\n");
        assert_eq!(&m * &m, m);
    }

    #[test]
    fn powers_mean_within_d_steps() {
        // 0 - 1 - 2
        let a = graph(3, &[(0, 1), (1, 2)]);
        assert_eq!(a.to_string(), "110\n111\n011\n");
        assert_eq!((&a * &a).to_string(), "111\n111\n111\n");
        assert_eq!(a.pow(0), ReachabilityMatrix::new(3));
        assert_eq!(a.pow(1), a);
    }

    #[test]
    #[should_panic]
    fn diagonal_cannot_be_cleared() {
        ReachabilityMatrix::new(2).set(1, 1, false);
    }

    #[test]
    fn sum_is_union_of_edges() {
        let sum = &graph(3, &[(0, 1)]) + &graph(3, &[(1, 2)]);
        assert_eq!(sum, graph(3, &[(0, 1), (1, 2)]));
    }

    #[test]
    fn powers_become_idempotent() {
        // Directed cycle 0 -> 1 -> 2 -> 3 -> 0.
        let mut c = ReachabilityMatrix::new(4);
        for i in 0..4 {
            c.set(i, (i + 1) % 4, true);
        }
        assert!(c.pow(2).get(0, 2) && !c.pow(2).get(0, 3));
        assert!((0..4).all(|j| c.pow(3).get(0, j)));
        assert_eq!(c.pow(3), c.pow(4));
    }

    #[test]
    fn shortest_path_on_a_line() {
        let a = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert_eq!(a.shortest_path(0, 4), Some(vec![0, 1, 2, 3, 4]));
        assert_eq!(a.shortest_path(3, 1), Some(vec![3, 2, 1]));
        assert_eq!(a.shortest_path(2, 2), Some(vec![2]));
    }

    #[test]
    fn even_cycle_and_isolated_vertex() {
        // A 4-cycle is bipartite: exact-step powers would alternate forever.
        let a = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        assert_eq!(a.shortest_path(0, 2).map(|p| p.len()), Some(3));
        assert_eq!(a.shortest_path(0, 4), None); // vertex 4 is isolated
    }

    #[test]
    fn respects_edge_direction() {
        let mut a = ReachabilityMatrix::new(3);
        a.set(0, 1, true);
        a.set(1, 2, true);
        assert_eq!(a.shortest_path(0, 2), Some(vec![0, 1, 2]));
        assert_eq!(a.shortest_path(2, 0), None);
    }

    #[test]
    fn takes_the_shortcut() {
        // Long way 0-1-2-3-4-5, shortcut 0-6-5.
        let a = graph(7, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 6), (6, 5)]);
        assert_eq!(a.shortest_path(0, 5), Some(vec![0, 6, 5]));
    }

    #[test]
    fn matches_bfs_on_pseudo_random_graphs() {
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for _ in 0..30 {
            let n = 2 + (next() % 15) as usize;
            let mut m = ReachabilityMatrix::new(n);
            for i in 0..n {
                for j in 0..n {
                    if i != j && next() % 5 == 0 {
                        m.set(i, j, true); // directed
                    }
                }
            }
            for a in 0..n {
                for b in 0..n {
                    let path = m.shortest_path(a, b);
                    assert_eq!(path.as_ref().map(|p| p.len() - 1), bfs_distance(&m, a, b), "{a} -> {b}\n{m}");
                    if let Some(p) = path {
                        assert!(is_path(&m, &p, a, b), "bad path {p:?}\n{m}");
                    }
                }
            }
        }
    }
}
