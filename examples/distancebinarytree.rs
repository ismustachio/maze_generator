use maze_generator::{BinaryTree, DistanceTable, Grid};

fn main() {
    let g = Grid::new(10, 10);
    let mut bg = BinaryTree::on(&g);
    let dt = DistanceTable::new(&bg, bg[(0, 0)]);
    let max_tup = dt.max_table_entry();
    bg.set_cell_bodies(&dt);
    println!("Binary Tree:\n{}", bg);
    bg.to_png("./images/distancebinarytree.png".to_string(), 10, max_tup.1);
}
