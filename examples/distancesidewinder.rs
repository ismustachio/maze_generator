use maze_generator::{DistanceTable, Grid, Sidewinder};

fn main() {
    let g = Grid::new(20, 20);
    let mut sg = Sidewinder::on(&g);
    let dt = DistanceTable::new(&sg, sg[(0, 0)]);
    let max_tup = dt.max_table_entry();
    sg.set_cell_bodies(&dt);
    println!("Sidewinder:\n{}", sg);
    sg.to_png("./images/distancesidewinder.png".to_string(), 10, max_tup.1);
}
