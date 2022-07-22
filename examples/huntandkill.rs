use maze_generator::{Grid, HuntAndKill};

fn main() {
    let g = Grid::new(20, 20);
    let bg = HuntAndKill::on(&g);
    bg.to_png("./images/huntandkill.png".to_string(), 10, 1);
    println!("HuntAndKill:\n{}", bg);
}
