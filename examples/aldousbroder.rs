use maze_generator::{AldousBroder, Grid};

fn main() {
    let g = Grid::new(10, 10);
    let ag = AldousBroder::on(&g);
    ag.to_png("./images/aldousbroder.png".to_string(), 10, 1);
    println!("AldousBroder:\n{}", ag);
}
