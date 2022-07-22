use maze_generator::{BinaryTree, Grid};

fn main() {
    let g = Grid::new(20, 20);
    let bg = BinaryTree::on(&g);
    bg.to_png("./images/binarytree.png".to_string(), 10, 1);
    println!("Binary Tree:\n{}", bg);
}
