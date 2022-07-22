use maze_generator::{Grid, Wilsons};

fn main() {
    let g = Grid::new(10, 10);
    let wg = Wilsons::on(&g);
    wg.to_png("./images/wilsons.png".to_string(), 10, 1);
    println!("Wilsons:\n{}", wg);
}
