use maze_generator::{Grid, RecursiveBacktracker};

fn main() {
    let g = Grid::new(20, 20);
    let bg = RecursiveBacktracker::on(&g, None);
    bg.to_png("./images/RecursiveBacktracker.png".to_string(), 10, 1);
    println!("RecursiveBacktracker:\n{}", bg);
}
