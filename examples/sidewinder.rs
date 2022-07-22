use maze_generator::{Grid, Sidewinder};

fn main() {
    let g = Grid::new(10, 10);
    let sg = Sidewinder::on(&g);
    sg.to_png("./images/sidewinder.png".to_string(), 10, 1);
    println!("Sidewinder:\n{}", sg);
}
