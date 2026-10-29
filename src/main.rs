//! Demo: generate a random maze, pick a start and an exit, and print the
//! shortest path if there is one.
//!
//!     cargo run --release -- [--width W] [--height H] [--density P] [--seed S]

use maze_solver::{Lattice, Rng};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

const USAGE: &str = "\
usage: maze-solver [--width W] [--height H] [--density P] [--seed S]

  --width W     tiles across        (default 16)
  --height H    tiles down          (default 10)
  --density P   chance each wall between neighbours is open, 0..1 (default 0.6)
  --seed S      random seed, to reproduce a maze (default: from the clock)

The solver keeps one n x n matrix per step of the path, with n = W * H,
so W * H is capped at 400 to bound memory.";

struct Options {
    width: usize,
    height: usize,
    density: f64,
    seed: u64,
}

fn parse_args() -> Result<Options, String> {
    let clock = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_nanos() as u64).unwrap_or(0);
    let mut opts = Options { width: 16, height: 10, density: 0.6, seed: clock };
    let mut args = std::env::args().skip(1);
    while let Some(flag) = args.next() {
        if flag == "-h" || flag == "--help" {
            return Err(String::new());
        }
        let value = args.next().ok_or_else(|| format!("{flag} needs a value"))?;
        let bad = || format!("invalid value for {flag}: {value}");
        match flag.as_str() {
            "--width" => opts.width = value.parse().map_err(|_| bad())?,
            "--height" => opts.height = value.parse().map_err(|_| bad())?,
            "--density" => opts.density = value.parse().map_err(|_| bad())?,
            "--seed" => opts.seed = value.parse().map_err(|_| bad())?,
            _ => return Err(format!("unknown option {flag}")),
        }
    }
    if opts.width == 0 || opts.height == 0 || opts.width * opts.height < 2 {
        return Err("the maze needs at least 2 tiles".into());
    }
    if opts.width * opts.height > 400 {
        return Err("width * height must be at most 400".into());
    }
    if !(0.0..=1.0).contains(&opts.density) {
        return Err("density must be between 0 and 1".into());
    }
    Ok(opts)
}

fn main() -> ExitCode {
    let opts = match parse_args() {
        Ok(o) => o,
        Err(msg) => {
            if !msg.is_empty() {
                eprintln!("error: {msg}\n");
            }
            eprintln!("{USAGE}");
            return if msg.is_empty() { ExitCode::SUCCESS } else { ExitCode::from(2) };
        }
    };

    let mut rng = Rng::new(opts.seed);
    let maze = Lattice::random(opts.width, opts.height, opts.density, &mut rng);
    let start = maze.coords(rng.below(maze.len()));
    let exit = loop {
        let tile = maze.coords(rng.below(maze.len()));
        if tile != start {
            break tile;
        }
    };

    println!(
        "{}x{} maze, density {}, seed {}  (rerun with --seed {} to reproduce)",
        opts.width, opts.height, opts.density, opts.seed, opts.seed
    );
    println!("start S = {start:?}, exit E = {exit:?}\n");

    match maze.shortest_path(start, exit) {
        Some(path) => {
            print!("{}", maze.render(start, exit, &path));
            let steps: Vec<String> = path.iter().map(|(x, y)| format!("({x},{y})")).collect();
            println!("\nShortest path: {} steps", path.len() - 1);
            println!("{}", steps.join(" -> "));
        }
        None => {
            print!("{}", maze.render(start, exit, &[]));
            println!("\nNo path from S to E.");
        }
    }
    ExitCode::SUCCESS
}
