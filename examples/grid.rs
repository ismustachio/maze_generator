extern crate maze_generator;
use maze_generator::Grid;

fn main() {
    let g = Grid::new(14, 14);
    g.to_png("./images/grid.png".to_string(), 10, 1);
    println!("Grid:\n{}", g);
}
