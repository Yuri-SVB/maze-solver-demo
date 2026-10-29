//! Maze solver engine: lattice mazes turned into boolean reachability matrices,
//! solved by repeated "or/and" matrix multiplication.

pub mod lattice;
pub mod matrix;
pub mod rng;

pub use lattice::{Direction, Lattice};
pub use matrix::ReachabilityMatrix;
pub use rng::Rng;
