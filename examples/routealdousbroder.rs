use maze_generator::{AldousBroder, DistanceTable, Grid};

fn main() {
    let g = Grid::new(14, 14);
    let mut ag = AldousBroder::on(&g);
    let dt = DistanceTable::new(&ag, ag[(0, 0)]);
    let max_tup = dt.max_table_entry();
    if let Some(route_dt) = dt.route_from_start(&ag[(7, 7)], &ag) {
        ag.set_cell_bodies(&route_dt);
        println!("AldousBroder:\n{}", ag);
        ag.to_png("./images/routealdousebroder.png".to_string(), 10, max_tup.1);
    } else {
        println!("Fatal error no routes available");
    }
}
