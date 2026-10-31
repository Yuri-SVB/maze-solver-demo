# Maze Solver

Finds shortest paths in random lattice mazes using boolean reachability-matrix
powers.

```sh
cd maze-solver
cargo run --release                         # random maze, random seed
cargo run --release -- --seed 2             # reproduce a maze
cargo run --release -- --width 20 --height 15 --density 0.55
cargo run --release -- --seed 2 --method linear   # compare the two searches
cargo test
```

## Pieces

| File            | What it holds                                                      |
|-----------------|--------------------------------------------------------------------|
| `src/matrix.rs` | `ReachabilityMatrix`: `Vec<Vec<bool>>` + size (diagonal always true), boolean `*` and `+`, `pow`, `shortest_path` (squaring or linear) |
| `src/lattice.rs`| `Lattice`: w x h tiles with left/right/up/down connections, `reachability()`, `render()` |
| `src/rng.rs`    | SplitMix64, so there are no dependencies and seeds reproduce mazes |
| `src/main.rs`   | CLI demo                                                           |

## How the solver works

Over the boolean semiring, `+` is OR and `*` is AND, so

```
(A * B)[i][j] = OR over k of (A[i][k] AND B[k][j])
```

**The diagonal is always true.** Every tile can "stay put", so `A[i][i]`
is true by construction: `new(n)` is the identity, and `set` refuses to
clear a diagonal entry. With that,

```
A^d[i][j] = true  iff  j is reachable from i in at most d steps
```

for every `d`, including `A^0 = I` and `A^1 = A`. This matters for mazes.
Without the diagonal, `A^d` would mean "exactly d steps". A grid is
bipartite, so exact-step powers alternate between two patterns forever
and never become idempotent. "At most d steps" only ever gains entries.

There are two searches, chosen with `Method` (`--method` in the CLI).

### Squaring (default)

1. **Square.** Compute `S_0 = A, S_1 = A^2, S_2 = A^4, ...`, with
   `S_(t+1) = S_t * S_t`, and keep each one.
   - Stop at the first `S_T` with `[a][b]` set. The distance `d` then
     satisfies `2^(T-1) < d <= 2^T`.
   - If `S_t * S_t == S_t` while `[a][b]` is unset, `S_t` is the full
     closure (idempotent), so `b` is unreachable.
2. **Bisect (binary lifting).** Find the largest `s < 2^T` with `A^s[a][b]`
   false. For `t = T-1 down to 0`, try `s + 2^t` and keep it if `[a][b]` is
   still false. Only row `a` of `A^s` is ever needed, so each try is a
   row-vector times matrix product. Then `d = s + 1`.
3. **Rebuild the path by midpoints.** For `u -> v` of length `k`, let `h`
   be the largest power of two below `k`. Pick any `m` with `A^h[u][m]` and
   `A^(k-h)[m][v]`. Since `u -> v` can't be shorter than `k`, `m` is exactly
   `h` steps from `u` on a shortest path. Recurse on `u -> m` and `m -> v`.
   `A^h` is a stored square. Column `v` of `A^(k-h)` is the unit column
   pushed through the squares for each bit of `k - h`.

### Linear

Row `a` of `A^s` is the set of tiles within `s` steps of `a`. That row is
all the walk back reads, and the next row follows from it alone:
`row_a(A^(s+1)) = row_a(A^s) * A`. So the linear search never forms a full
power. It keeps one row of `n` entries per step.

1. Compute `row_a(A^0) = e_a, row_a(A^1), row_a(A^2), ...` and keep each.
   - The first `d` whose row has `b` set is the shortest distance.
   - If a row equals the previous one while `b` is still unset, every later
     row is the same, so `b` is unreachable. This happens within `n` steps.
2. Walk backwards from `b`. At step `s`, the predecessor of the current tile
   is any `v` with `A[v][current]` set and `v` in row `s`. Because the
   current tile is exactly `s + 1` steps from `a`, `v` is exactly `s` steps
   away, so the result is a shortest path. The self-loop never gets picked:
   the current tile isn't in row `s`.

Each step is a row-vector product, which ORs together the rows of `A` for
the tiles reached so far. In effect this is breadth-first search written
as linear algebra.

## Limits

The matrices are `n x n` with `n = width * height`. Measured in a release
build, on a serpentine maze whose only path visits every tile (the worst
case for path length):

| maze  | path | linear (row a only) | squaring |
|-------|------|---------------------|----------|
| 20x20 | 399  | 1.5 ms              | 7.7 ms   |
| 40x40 | 1599 | 66 ms               | 0.38 s   |
| 50x50 | 2499 | 0.29 s              | 1.7 s    |

On random 20x20 mazes linear takes 30-190 µs and squaring 1-6 ms. Linear
stores `d` rows (`d * n` bools, ~6 MB at 50x50). Squaring stores about
`log2(d)` full matrices and does dense `n^3` products. The CLI allows
2500 tiles for linear and 1600 for squaring.
