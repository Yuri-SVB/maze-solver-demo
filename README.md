# Maze Solver

Finds shortest paths in random lattice mazes using boolean reachability-matrix
powers.

```sh
cd maze-solver
cargo run --release                         # random maze, random seed
cargo run --release -- --seed 2             # reproduce a maze
cargo run --release -- --width 20 --height 15 --density 0.55
cargo test
```

## Pieces

| File            | What it holds                                                      |
|-----------------|--------------------------------------------------------------------|
| `src/matrix.rs` | `ReachabilityMatrix`: `Vec<Vec<bool>>` + size (diagonal always true), boolean `*` and `+`, `pow`, `shortest_path` |
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

1. Compute `A^0 = I, A^1, A^2, ...` and keep each one.
   - The first `d` with `A^d[a][b]` set is the shortest distance.
   - If `A^d == A^(d-1)` while `[a][b]` is still unset, the power is
     idempotent and `b` is unreachable. This happens within `n` steps.
2. Walk backwards from `b`. At step `s`, the predecessor of the current tile
   is any `v` with `A[v][current]` and `A^s[a][v]` set. Because the current
   tile is exactly `s + 1` steps from `a`, `v` is exactly `s` steps away, so
   the result is a shortest path. The self-loop never gets picked: the
   current tile isn't within `s` steps.

The product skips the left operand's false entries, and the solver
multiplies the sparse `A` on the left (`A * A^(d-1)`, which equals
`A^(d-1) * A` because both are powers of `A`). That makes each step
O(n² + edges·n) instead of O(n³).

## Limits

The solver keeps every intermediate power, one `n x n` matrix per step of
the path, with `n = width * height`. The CLI caps `n` at 400. The worst case
there, a 399-step serpentine path through a 20x20 grid, takes about 0.3 s
and 70 MB in a release build. Only row `a` of each stored power is used
when walking back, so storing just that row would cut memory by a factor
of `n`.
