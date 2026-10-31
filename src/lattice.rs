//! A `width x height` lattice of square tiles. Each tile may be connected to
//! its left, right, up and down neighbours.
//!
//! Tiles are numbered row by row: tile `(x, y)` has index `y * width + x`,
//! with `(0, 0)` at the top left and `y` growing downwards.

use crate::matrix::{Method, ReachabilityMatrix};
use crate::rng::Rng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Direction {
    Left,
    Right,
    Up,
    Down,
}

impl Direction {
    pub const ALL: [Direction; 4] = [Direction::Left, Direction::Right, Direction::Up, Direction::Down];

    fn delta(self) -> (isize, isize) {
        match self {
            Direction::Left => (-1, 0),
            Direction::Right => (1, 0),
            Direction::Up => (0, -1),
            Direction::Down => (0, 1),
        }
    }
}

/// Connections are symmetric, so each one is stored once: on the tile to
/// its left (`right`) or the tile above it (`down`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lattice {
    width: usize,
    height: usize,
    right: Vec<bool>,
    down: Vec<bool>,
}

impl Lattice {
    /// A lattice with every connection closed.
    pub fn new(width: usize, height: usize) -> Self {
        assert!(width > 0 && height > 0, "lattice must be at least 1x1");
        let n = width * height;
        Lattice { width, height, right: vec![false; n], down: vec![false; n] }
    }

    /// Opens each possible connection independently with probability
    /// `density`. Near 0.5 there's a fair chance the start and exit aren't
    /// connected at all.
    pub fn random(width: usize, height: usize, density: f64, rng: &mut Rng) -> Self {
        let mut lattice = Self::new(width, height);
        for y in 0..height {
            for x in 0..width {
                if x + 1 < width && rng.chance(density) {
                    lattice.set_connected(x, y, Direction::Right, true);
                }
                if y + 1 < height && rng.chance(density) {
                    lattice.set_connected(x, y, Direction::Down, true);
                }
            }
        }
        lattice
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn len(&self) -> usize {
        self.width * self.height
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn index(&self, x: usize, y: usize) -> usize {
        assert!(x < self.width && y < self.height, "tile ({x}, {y}) is outside the lattice");
        y * self.width + x
    }

    pub fn coords(&self, index: usize) -> (usize, usize) {
        assert!(index < self.len(), "tile {index} is outside the lattice");
        (index % self.width, index / self.width)
    }

    /// The tile next to `(x, y)` in `dir`, if it's inside the lattice.
    pub fn neighbor(&self, x: usize, y: usize, dir: Direction) -> Option<(usize, usize)> {
        let (dx, dy) = dir.delta();
        let nx = x.checked_add_signed(dx).filter(|&v| v < self.width)?;
        let ny = y.checked_add_signed(dy).filter(|&v| v < self.height)?;
        Some((nx, ny))
    }

    pub fn is_connected(&self, x: usize, y: usize, dir: Direction) -> bool {
        let Some((nx, ny)) = self.neighbor(x, y, dir) else {
            return false;
        };
        match dir {
            Direction::Right => self.right[self.index(x, y)],
            Direction::Left => self.right[self.index(nx, ny)],
            Direction::Down => self.down[self.index(x, y)],
            Direction::Up => self.down[self.index(nx, ny)],
        }
    }

    /// Open or close the connection between `(x, y)` and its neighbour in
    /// `dir`. Both tiles see the change. Panics at the edge of the lattice.
    pub fn set_connected(&mut self, x: usize, y: usize, dir: Direction, open: bool) {
        let (nx, ny) = self.neighbor(x, y, dir).expect("no neighbour in that direction");
        let slot = match dir {
            Direction::Right => &mut self.right[y * self.width + x],
            Direction::Left => &mut self.right[ny * self.width + nx],
            Direction::Down => &mut self.down[y * self.width + x],
            Direction::Up => &mut self.down[ny * self.width + nx],
        };
        *slot = open;
    }

    /// One vertex per tile, with an undirected edge for each open connection
    /// (and, as always for a reachability matrix, each tile reaching itself).
    pub fn reachability(&self) -> ReachabilityMatrix {
        let mut m = ReachabilityMatrix::new(self.len());
        for y in 0..self.height {
            for x in 0..self.width {
                for dir in [Direction::Right, Direction::Down] {
                    if self.is_connected(x, y, dir) {
                        let (nx, ny) = self.neighbor(x, y, dir).expect("connected implies inside");
                        m.connect(self.index(x, y), self.index(nx, ny));
                    }
                }
            }
        }
        m
    }

    /// Shortest path between two tiles, as tile coordinates.
    pub fn shortest_path(&self, from: (usize, usize), to: (usize, usize)) -> Option<Vec<(usize, usize)>> {
        self.shortest_path_with(from, to, Method::Squaring)
    }

    pub fn shortest_path_with(
        &self,
        from: (usize, usize),
        to: (usize, usize),
        method: Method,
    ) -> Option<Vec<(usize, usize)>> {
        let (a, b) = (self.index(from.0, from.1), self.index(to.0, to.1));
        let path = self.reachability().shortest_path_with(a, b, method)?;
        Some(path.into_iter().map(|i| self.coords(i)).collect())
    }

    /// Draw the maze with `#` walls, marking `start` with S, `exit` with E and
    /// the tiles and openings along `path` with `.`. Each canvas cell is
    /// printed two characters wide so the maze looks square.
    pub fn render(&self, start: (usize, usize), exit: (usize, usize), path: &[(usize, usize)]) -> String {
        let (cols, rows) = (2 * self.width + 1, 2 * self.height + 1);
        let mut canvas = vec![vec!['#'; cols]; rows];
        for y in 0..self.height {
            for x in 0..self.width {
                canvas[2 * y + 1][2 * x + 1] = ' ';
                if self.is_connected(x, y, Direction::Right) {
                    canvas[2 * y + 1][2 * x + 2] = ' ';
                }
                if self.is_connected(x, y, Direction::Down) {
                    canvas[2 * y + 2][2 * x + 1] = ' ';
                }
            }
        }
        for (i, &(x, y)) in path.iter().enumerate() {
            canvas[2 * y + 1][2 * x + 1] = '.';
            if let Some(&(px, py)) = i.checked_sub(1).and_then(|p| path.get(p)) {
                // The opening between consecutive tiles sits halfway between them.
                canvas[y + py + 1][x + px + 1] = '.';
            }
        }
        canvas[2 * start.1 + 1][2 * start.0 + 1] = 'S';
        canvas[2 * exit.1 + 1][2 * exit.0 + 1] = 'E';

        let mut out = String::new();
        for row in canvas {
            for c in row {
                out.push_str(match c {
                    '#' => "██",
                    '.' => "··",
                    'S' => "S ",
                    'E' => "E ",
                    _ => "  ",
                });
            }
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 3x2 lattice shaped like a U:
    /// (0,0) (1,0)-(2,0)
    ///   |           |
    /// (0,1)-(1,1)-(2,1)
    fn u_shape() -> Lattice {
        let mut l = Lattice::new(3, 2);
        l.set_connected(0, 0, Direction::Down, true);
        l.set_connected(0, 1, Direction::Right, true);
        l.set_connected(1, 1, Direction::Right, true);
        l.set_connected(2, 1, Direction::Up, true);
        l.set_connected(2, 0, Direction::Left, true);
        l
    }

    #[test]
    fn index_and_coords_round_trip() {
        let l = Lattice::new(4, 3);
        assert_eq!(l.index(3, 2), 11);
        assert_eq!(l.coords(11), (3, 2));
        assert!((0..12).all(|i| { let (x, y) = l.coords(i); l.index(x, y) == i }));
    }

    #[test]
    fn connections_are_symmetric() {
        let l = u_shape();
        assert!(l.is_connected(0, 0, Direction::Down) && l.is_connected(0, 1, Direction::Up));
        assert!(l.is_connected(1, 0, Direction::Right) && l.is_connected(2, 0, Direction::Left));
        assert!(!l.is_connected(0, 0, Direction::Right) && !l.is_connected(1, 0, Direction::Left));
        assert!(!l.is_connected(0, 0, Direction::Up)); // edge of the lattice
    }

    #[test]
    #[should_panic]
    fn cannot_connect_past_the_edge() {
        Lattice::new(2, 2).set_connected(1, 0, Direction::Right, true);
    }

    #[test]
    fn reachability_has_one_vertex_per_tile() {
        let m = u_shape().reachability();
        assert_eq!(m.size(), 6);
        assert!((0..6).all(|i| m.get(i, i)));
        let edges = (0..6).flat_map(|i| (0..6).map(move |j| (i, j))).filter(|&(i, j)| i != j && m.get(i, j)).count();
        assert_eq!(edges, 2 * 5); // 5 undirected connections
        assert!(m.get(0, 3) && m.get(3, 0)); // (0,0) <-> (0,1)
        assert!(!m.get(0, 1)); // (0,0) and (1,0) are walled off
    }

    #[test]
    fn path_goes_around_the_wall() {
        let l = u_shape();
        assert_eq!(
            l.shortest_path((0, 0), (1, 0)),
            Some(vec![(0, 0), (0, 1), (1, 1), (2, 1), (2, 0), (1, 0)])
        );
    }

    #[test]
    fn no_path_through_closed_walls() {
        let mut l = u_shape();
        l.set_connected(1, 1, Direction::Right, false);
        assert_eq!(l.shortest_path((0, 0), (2, 0)), None);
    }

    #[test]
    fn random_density_extremes() {
        let mut rng = Rng::new(1);
        let closed = Lattice::random(5, 4, 0.0, &mut rng);
        assert_eq!(closed.shortest_path((0, 0), (4, 3)), None);
        let open = Lattice::random(5, 4, 1.0, &mut rng);
        assert_eq!(open.shortest_path((0, 0), (4, 3)).map(|p| p.len()), Some(8)); // Manhattan + 1
    }

    #[test]
    fn render_marks_start_exit_and_path() {
        let l = u_shape();
        let path = l.shortest_path((0, 0), (1, 0)).unwrap();
        let picture = l.render((0, 0), (1, 0), &path);
        let expected = [
            "██████████████",
            "██S ██E ····██",
            "██··██████··██",
            "██··········██",
            "██████████████",
        ]
        .join("\n")
            + "\n";
        assert_eq!(picture, expected);
    }
}
