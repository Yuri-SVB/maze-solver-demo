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
    /// included), or `None` if `b` is unreachable. Uses repeated squaring.
    pub fn shortest_path(&self, a: usize, b: usize) -> Option<Vec<usize>> {
        self.shortest_path_with(a, b, Method::Squaring)
    }

    pub fn shortest_path_with(&self, a: usize, b: usize, method: Method) -> Option<Vec<usize>> {
        match method {
            Method::Linear => self.shortest_path_linear(a, b),
            Method::Squaring => self.shortest_path_squaring(a, b),
        }
    }

    /// Shortest path by repeated squaring: O(log d) matrix products instead of
    /// d, and O(log d) stored matrices instead of d.
    ///
    /// 1. Square: `S_t = A^(2^t)`, so `S_0 = A`, `S_1 = A^2`, `S_2 = A^4`, ...
    ///    Stop at the first `S_T` with `[a][b]` set: the distance `d` then
    ///    satisfies `2^(T-1) < d <= 2^T`. If instead `S_t^2 == S_t`, the powers
    ///    are idempotent (`S_t` is the full closure) with `[a][b]` false, so
    ///    `b` is unreachable.
    /// 2. Bisect (binary lifting): find the largest `s < 2^T` with
    ///    `A^s[a][b]` false. Going through `t = T-1 .. 0`, try adding `2^t`
    ///    to `s`, and keep it if `[a][b]` is still false. Only row `a` of
    ///    `A^s` is needed, so each try is a row-vector product, O(n^2). Then
    ///    `d = s + 1`.
    /// 3. Rebuild the path by splitting at midpoints. For a path `u -> v` of
    ///    length `k`, let `h` be the largest power of two below `k`. Any `m`
    ///    with `A^h[u][m]` and `A^(k-h)[m][v]` is exactly `h` steps from `u`
    ///    on some shortest path (otherwise `u -> v` would be shorter than
    ///    `k`). Recurse on both halves.
    pub fn shortest_path_squaring(&self, a: usize, b: usize) -> Option<Vec<usize>> {
        assert!(a < self.size && b < self.size, "vertex out of range");
        if a == b {
            return Some(vec![a]);
        }
        let mut squares = vec![self.clone()];
        loop {
            let last = squares.last().expect("squares starts non-empty");
            if last.get(a, b) {
                break;
            }
            let next = last * last;
            if &next == last {
                return None;
            }
            squares.push(next);
        }

        let top = squares.len() - 1;
        let mut row = unit(self.size, a); // row a of A^s
        let mut s = 0;
        for t in (0..top).rev() {
            let candidate = squares[t].row_times(&row);
            if !candidate[b] {
                row = candidate;
                s += 1 << t;
            }
        }

        let mut path = vec![a];
        extend_path(&squares, a, b, s + 1, &mut path);
        Some(path)
    }

    /// `x * self` for a row vector `x`: the vertices reachable within the
    /// matrix's step count from anything in `x`.
    fn row_times(&self, x: &[bool]) -> Vec<bool> {
        let mut out = vec![false; self.size];
        for (k, _) in x.iter().enumerate().filter(|(_, &on)| on) {
            for (o, &r) in out.iter_mut().zip(&self.rows[k]) {
                *o |= r;
            }
        }
        out
    }

    /// `self * x` for a column vector `x`: the vertices that reach something
    /// in `x` within the matrix's step count.
    fn times_column(&self, x: &[bool]) -> Vec<bool> {
        self.rows.iter().map(|row| row.iter().zip(x).any(|(&r, &c)| r && c)).collect()
    }

    /// Shortest path by computing the powers `A^0, A^1, A^2, ...` one at a
    /// time, keeping only row `a` of each: d row-vector products and d
    /// stored rows of n entries.
    ///
    /// Row `a` of `A^s` is the set of vertices within `s` steps of `a`. It is
    /// all the walk back ever reads, and the next one follows from it alone:
    /// row `a` of `A^(s+1)` = (row `a` of `A^s`) * `A`.
    ///
    /// 1. Compute the rows until the first one with `b` set; that `s` is the
    ///    shortest distance. Rows only ever gain entries. If one equals the
    ///    previous row with `b` still unset, every later row is the same
    ///    (each is the previous times `A`), so `b` is unreachable. This
    ///    happens within `size` steps.
    /// 2. Walk backwards. If the current vertex is at distance `s + 1`, its
    ///    predecessor is any `v` with `A[v][current]` set and `v` in row `s`:
    ///    `v` is at most `s` steps from `a`, so exactly `s`. The current
    ///    vertex itself never qualifies, because it isn't in row `s`.
    pub fn shortest_path_linear(&self, a: usize, b: usize) -> Option<Vec<usize>> {
        assert!(a < self.size && b < self.size, "vertex out of range");
        let mut rows = vec![unit(self.size, a)]; // rows[s] = row a of A^s
        loop {
            let last = rows.last().expect("rows starts non-empty");
            if last[b] {
                break;
            }
            let next = self.row_times(last);
            if &next == last {
                return None;
            }
            rows.push(next);
        }

        let distance = rows.len() - 1;
        let mut path = vec![b];
        let mut current = b;
        for s in (0..distance).rev() {
            current = (0..self.size)
                .find(|&v| rows[s][v] && self.get(v, current))
                .expect("a vertex at distance s + 1 has a neighbour at distance s");
            path.push(current);
        }
        path.reverse();
        Some(path)
    }
}

/// How `shortest_path` searches for the distance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// `A^0, A^1, A^2, ...` one power at a time, keeping only row `a`.
    Linear,
    /// `A, A^2, A^4, ...` by squaring, then bisection.
    Squaring,
}

/// The unit row vector for vertex `i`: row `i` of `A^0`.
fn unit(size: usize, i: usize) -> Vec<bool> {
    let mut v = vec![false; size];
    v[i] = true;
    v
}

/// Append the vertices after `u` on a shortest path `u -> v` of length `k`,
/// using `squares[t] = A^(2^t)`. Needs `k <= 2^(squares.len() - 1)`.
fn extend_path(squares: &[ReachabilityMatrix], u: usize, v: usize, k: usize, path: &mut Vec<usize>) {
    if k == 0 {
        return;
    }
    if k == 1 {
        path.push(v);
        return;
    }
    let j = (usize::BITS - 1 - (k - 1).leading_zeros()) as usize; // h = 2^j, the largest power of two < k
    let h = 1 << j;
    let rest = k - h; // 1..=h
    let reaches_v: Vec<bool> = if rest.is_power_of_two() {
        let r = &squares[rest.trailing_zeros() as usize];
        (0..r.size()).map(|m| r.get(m, v)).collect()
    } else {
        // Column v of A^rest: multiply the unit column through the squares for each bit of rest.
        let mut col = unit(squares[0].size(), v);
        for (t, square) in squares.iter().enumerate().take(j + 1) {
            if rest & (1 << t) != 0 {
                col = square.times_column(&col);
            }
        }
        col
    };
    let m = (0..squares[j].size())
        .find(|&m| squares[j].get(u, m) && reaches_v[m])
        .expect("a shortest path has a vertex h steps in");
    extend_path(squares, u, m, h, path);
    extend_path(squares, m, v, rest, path);
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

    /// Solve with both methods: they must agree on the distance and both
    /// return valid paths. Returns the squaring result.
    fn both(m: &ReachabilityMatrix, a: usize, b: usize) -> Option<Vec<usize>> {
        let linear = m.shortest_path_with(a, b, Method::Linear);
        let squaring = m.shortest_path_with(a, b, Method::Squaring);
        assert_eq!(linear.as_ref().map(Vec::len), squaring.as_ref().map(Vec::len), "{a} -> {b}\n{m}");
        for path in [&linear, &squaring].into_iter().flatten() {
            assert!(is_path(m, path, a, b), "bad path {path:?}\n{m}");
        }
        squaring
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
        assert_eq!(both(&a, 0, 4), Some(vec![0, 1, 2, 3, 4]));
        assert_eq!(both(&a, 3, 1), Some(vec![3, 2, 1]));
        assert_eq!(both(&a, 2, 2), Some(vec![2]));
    }

    #[test]
    fn even_cycle_and_isolated_vertex() {
        // A 4-cycle is bipartite: exact-step powers would alternate forever.
        let a = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        assert_eq!(both(&a, 0, 2).map(|p| p.len()), Some(3));
        assert_eq!(both(&a, 0, 4), None); // vertex 4 is isolated
    }

    #[test]
    fn respects_edge_direction() {
        let mut a = ReachabilityMatrix::new(3);
        a.set(0, 1, true);
        a.set(1, 2, true);
        assert_eq!(both(&a, 0, 2), Some(vec![0, 1, 2]));
        assert_eq!(both(&a, 2, 0), None);
    }

    #[test]
    fn takes_the_shortcut() {
        // Long way 0-1-2-3-4-5, shortcut 0-6-5.
        let a = graph(7, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 6), (6, 5)]);
        assert_eq!(both(&a, 0, 5), Some(vec![0, 6, 5]));
    }

    #[test]
    fn long_line_every_pair() {
        // Distances up to 39 exercise every bit pattern in the bisection.
        let n = 40;
        let edges: Vec<(usize, usize)> = (0..n - 1).map(|i| (i, i + 1)).collect();
        let a = graph(n, &edges);
        for i in 0..n {
            for j in 0..n {
                let expected: Vec<usize> = if i <= j { (i..=j).collect() } else { (j..=i).rev().collect() };
                assert_eq!(both(&a, i, j), Some(expected));
            }
        }
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
                    let path = both(&m, a, b);
                    assert_eq!(path.as_ref().map(|p| p.len() - 1), bfs_distance(&m, a, b), "{a} -> {b}\n{m}");
                    if let Some(p) = path {
                        assert!(is_path(&m, &p, a, b), "bad path {p:?}\n{m}");
                    }
                }
            }
        }
    }
}
