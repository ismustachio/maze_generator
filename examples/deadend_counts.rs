use maze_generator::{AldousBroder, BinaryTree, Grid, HuntAndKill, Sidewinder, Wilsons};
use std::collections::HashMap;

fn main() {
    let algorithms: Vec<fn(&Grid) -> Grid> = [
        BinaryTree::on,
        Sidewinder::on,
        AldousBroder::on,
        Wilsons::on,
        HuntAndKill::on,
    ]
    .to_vec();
    let algorithms_name: Vec<&str> = [
        "BinaryTree",
        "Sidewinder",
        "AldousBroder",
        "Wilsons",
        "HuntAndKill",
    ]
    .to_vec();
    let tries = 100;
    let size = 20;
    let mut averages = HashMap::new();
    for (index, algorithm) in algorithms.iter().enumerate() {
        println!("Running {}", algorithms_name[index]);
        let mut deadend_counts = vec![];
        for _ in 0..tries {
            let grid = Grid::new(size, size);
            let agrid = algorithm(&grid);
            deadend_counts.push(agrid.deadends().len());
        }
        let len = deadend_counts.len();
        if let Some(sum) = deadend_counts.into_iter().reduce(|a, b| a + b) {
            averages.insert(algorithms_name[index], sum / len);
        }
    }
    let total_size = size * size;
    println!(
        "Average dead ends per {}x{} maze {} cells",
        size, size, total_size
    );
    for algorithm in algorithms_name {
        if let Some(average) = averages.get(algorithm) {
            let percentage = average * 100 / total_size;
            println!(
                "{}    {}/{}  {}%",
                algorithm, average, total_size, percentage
            );
        }
    }
}
