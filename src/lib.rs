use image::{ImageBuffer, Pixel, RgbImage};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::ops::Index;
use std::ops::IndexMut;

#[derive(Copy, Clone, PartialEq, Debug, Hash, Eq)]
pub struct Position {
    pub row: usize,
    pub col: usize,
}

#[derive(Copy, Clone, PartialEq, Debug, Hash, Eq)]
pub struct Cell {
    pub pos: Position,
    pub has_north: bool,
    pub has_south: bool,
    pub has_west: bool,
    pub has_east: bool,
}

pub struct DistanceTable {
    pub root: Cell,
    pub table: HashMap<Cell, usize>,
}

pub struct Grid {
    pub rows: usize,
    pub cols: usize,
    pub cells: Vec<Cell>,
    pub cell_links: HashMap<Cell, HashSet<Cell>>,
    pub cell_bodies: HashMap<Cell, String>,
}

pub struct Mask {
    pub rows: usize,
    pub cols: usize,
    pub bits: Vec<bool>,
}

impl DistanceTable {
    pub fn new(grid: &Grid, root: Cell) -> DistanceTable {
        let mut dt = DistanceTable {
            root,
            table: HashMap::from([(root, 0)]),
        };
        let mut frontier = vec![root];
        let mut new_frontier: Vec<Cell> = vec![];
        while !frontier.is_empty() {
            new_frontier.clear();
            for cell in &frontier {
                if grid.cell_links.get(cell).is_none() {
                    continue;
                }
                let links = grid.cell_links.get(cell).expect("should have cell links");
                for neighbor in links.iter() {
                    if dt.table.get(neighbor).is_none() {
                        let cell_distance =
                            dt.table.get(cell).expect("Should have a cell distance");
                        dt.table.insert(*neighbor, cell_distance + 1);
                        new_frontier.push(*neighbor);
                    }
                }
            }
            frontier.clear();
            frontier = new_frontier.clone();
        }
        dt
    }

    pub fn max_table_entry(&self) -> (Cell, usize) {
        let mut max_distance = 0;
        let mut max_cell = self.root;
        for _ in 0..2 {
            for (cell, distance) in &self.table {
                if distance > &max_distance {
                    max_cell = *cell;
                    max_distance = *distance;
                }
            }
        }
        (max_cell, max_distance)
    }

    pub fn route_from_start(&self, start_cell: &Cell, grid: &Grid) -> Option<DistanceTable> {
        let mut current = start_cell;
        let mut route = DistanceTable {
            root: self.root,
            table: HashMap::from([(self.root, 0)]),
        };
        if let Some(distance) = self.table.get(start_cell) {
            route.table.insert(*start_cell, *distance);
        } else {
            return None;
        }

        let mut neighbor_distance: usize = 0;
        let mut current_distance: usize = 0;
        while *current != self.root {
            let links = grid
                .cell_links
                .get(current)
                .expect("should have a set of cell links");
            if links.is_empty() {
                return None;
            }
            for link in links.iter() {
                if let Some(distance) = self.table.get(link) {
                    neighbor_distance = *distance;
                }
                if let Some(distance) = self.table.get(current) {
                    current_distance = *distance;
                }
                if neighbor_distance < current_distance {
                    route.table.insert(*link, neighbor_distance);
                    current = link;
                }
            }
        }
        Some(route)
    }
}

// BinaryTree is for each cell in the grid, randomly
// decide wether to link north or east.
pub struct BinaryTree;
impl BinaryTree {
    pub fn on(grid: &Grid) -> Grid {
        let mut out = Grid::new(grid.rows, grid.cols);
        let mut neighbors = vec![];
        let mut rng = rand::thread_rng();
        for cell in grid.cells.iter() {
            if cell.has_north {
                neighbors.push(grid.get_north(cell).expect("Should have a north"));
            }
            if cell.has_east {
                neighbors.push(grid.get_east(cell).expect("Should have an east"));
            }
            if !neighbors.is_empty() {
                let index = rng.gen_range(0..neighbors.len());
                out.link_cells(cell, &neighbors[index], true);
                neighbors.clear();
            }
        }
        out
    }
}

// Sidewinder randomly chooses at each step to
// save the current run of cells and either link
// the current cell to it's east neighbor, or
// pick a random cell from the current run of cells
// and link the current cell to it's north.
pub struct Sidewinder;
impl Sidewinder {
    pub fn on(grid: &Grid) -> Grid {
        let mut out = Grid::new(grid.rows, grid.cols);
        let mut run = vec![];
        let mut rng = rand::thread_rng();
        // TODO: will changing the uniform distribution have an effect on the shape?...
        let die = Uniform::from(0..2);
        for y in 0..grid.rows {
            run.clear();
            for x in 0..grid.cols {
                let cell = grid[(y, x)];
                run.push(cell);
                if !cell.has_east || (cell.has_north && die.sample(&mut rng) % 2 == 0) {
                    if !run.is_empty() {
                        let neighbor = run[rng.gen_range(0..run.len())];
                        if neighbor.has_north {
                            let north =
                                grid.get_north(&neighbor).expect("cell should have a north");
                            out.link_cells(&neighbor, &north, true);
                        }
                        run.clear();
                    }
                } else {
                    let east = grid.get_east(&cell).expect("cell should have a east");
                    out.link_cells(&cell, &east, true);
                }
            }
        }
        out
    }
}

// Start anywhere in the grid and choose a random neighbor.
// Move to that neighbor and if it hasn't previously been
// visited, link it to the prior cell. Repeat this process
// until every cell has been visited.
// developed by: Prof David Aldous, Distinguished Scientist at Google
// Andrei Broder
pub struct AldousBroder;
impl AldousBroder {
    pub fn on(grid: &Grid) -> Grid {
        let mut out = Grid::new(grid.rows, grid.cols);
        let mut unvisited = grid.size() - 1;
        // start a random cell on grid
        let mut cell = grid.random_cell();
        let mut rng = rand::thread_rng();
        while unvisited > 0 {
            let neighbors = grid.neighbors(&cell);
            if !neighbors.is_empty() {
                // select a random neighbor and link if not already
                let index = rng.gen_range(0..neighbors.len());
                let neighbor = neighbors[index];
                if let Some(links) = out.cell_links.get(&neighbor) {
                    if links.is_empty() {
                        out.link_cells(&cell, &neighbor, true);
                        unvisited -= 1;
                    }
                }
                cell = neighbor;
            }
        }
        out
    }
}

//
//
// Developed by David Bruce Wilson, a principle researcher at
// Microsoft
pub struct Wilsons;
impl Wilsons {
    pub fn on(grid: &Grid) -> Grid {
        let mut out = Grid::new(grid.rows, grid.cols);
        let mut unvisited: HashSet<Cell> = grid.cells.iter().cloned().collect();
        let mut rng = rand::thread_rng();
        unvisited.remove(&grid.cells[0]);

        while !unvisited.is_empty() {
            // all we want here is reference to the set cells
            let unvisited_vecs: Vec<&Cell> = unvisited.iter().collect();
            // select a random neighbor and link if not already
            let index = rng.gen_range(0..unvisited.len());
            let mut cell = *unvisited
                .get(unvisited_vecs[index])
                .expect("should have a cell");
            let mut path = vec![cell];

            while unvisited.contains(&cell) {
                let neighbors = grid.neighbors(&cell);
                if neighbors.is_empty() {
                    break;
                }
                let mut position: i32 = -1;
                let index = rng.gen_range(0..neighbors.len());
                cell = neighbors[index];
                for (index, c) in path.iter().enumerate() {
                    if *c == cell {
                        position = index as i32;
                        break;
                    }
                }
                if position >= 0 {
                    path = path[0..(position as usize)].to_vec();
                } else {
                    path.push(cell)
                }
            }
            if path.len() >= 2 {
                for index in 0..=(path.len() - 2) {
                    out.link_cells(&path[index], &path[index + 1], true);
                    unvisited.remove(&path[index]);
                }
            }
        }
        out
    }
}

pub struct HuntAndKill;
impl HuntAndKill {
    pub fn on(grid: &Grid) -> Grid {
        let mut out = Grid::new(grid.rows, grid.cols);
        let mut current: Option<Cell> = Some(grid.random_cell());
        let mut rng = rand::thread_rng();
        while current != None {
            let mut unvisited_neighbors = vec![];
            for neighbor in &out.neighbors(&current.expect("should have a cell")) {
                if let Some(links) = out.cell_links.get(neighbor) {
                    if links.is_empty() {
                        unvisited_neighbors.push(*neighbor);
                    }
                }
            }
            if !unvisited_neighbors.is_empty() {
                let index = rng.gen_range(0..unvisited_neighbors.len());
                let neighbor = unvisited_neighbors[index];
                out.link_cells(&current.expect("should have a cell"), &neighbor, true);
                current = Some(neighbor);
            } else {
                current = None;
                for cell in &grid.cells {
                    let mut visited_neighbors = vec![];
                    for neighbor in &grid.neighbors(cell) {
                        if let Some(links) = out.cell_links.get(neighbor) {
                            if !links.is_empty() {
                                visited_neighbors.push(*neighbor);
                            }
                        }
                    }

                    if let Some(links) = out.cell_links.get(cell) {
                        if links.is_empty() && !visited_neighbors.is_empty() {
                            let index = rng.gen_range(0..visited_neighbors.len());
                            current = Some(*cell);
                            let neighbor = visited_neighbors[index];
                            out.link_cells(&current.expect("should have a cell"), &neighbor, true);
                        }
                    }
                }
            }
        }
        out
    }
}

pub struct RecursiveBacktracker;
impl RecursiveBacktracker {
    pub fn on(grid: &Grid, start: Option<Cell>) -> Grid {
        let mut out = Grid::new(grid.rows, grid.cols);
        let mut rng = rand::thread_rng();
        let mut stack = vec![];
        if let Some(s) = start {
            stack.push(s);
        } else {
            stack.push(out.random_cell());
        }

        while !stack.is_empty() {
            let current = stack[stack.len() - 1];
            let mut neighbors = vec![];
            for neighbor in out.neighbors(&current) {
                if let Some(links) = out.cell_links.get(&neighbor) {
                    if links.is_empty() {
                        neighbors.push(neighbor);
                    }
                }
            }
            if neighbors.is_empty() {
                stack.pop();
            } else {
                let neighbor = neighbors[rng.gen_range(0..neighbors.len())];
                out.link_cells(&current, &neighbor, true);
                stack.push(neighbor);
            }
        }
        out
    }
}

impl fmt::Display for Grid {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut padding = 1;
        if !self.cell_bodies.is_empty() {
            let mut values = vec![];
            for v in self.cell_bodies.values() {
                values.push(v.trim().parse::<usize>().expect("should be a valid int"));
            }
            values.sort();
            padding = values[values.len() - 1].to_string().len();
        }
        let mut head = String::from("+");
        let ticks = "-".repeat(padding);
        let spaces = " ".repeat(padding);
        head.push_str(&ticks);
        head.push('-');
        head.push_str(&ticks);
        head.push('-');
        let mut east_boundary: String;
        let mut south_boundary: String;
        let mut south_boundary_ticks = ticks.clone();
        south_boundary_ticks.push('-');
        south_boundary_ticks.push('-');
        south_boundary_ticks.push_str(&ticks);
        let mut south_boundary_spaces = spaces.clone();
        south_boundary_spaces.push(' ');
        south_boundary_spaces.push(' ');
        south_boundary_spaces.push_str(&spaces);

        let mut body = head.repeat(self.cols);
        body.push('+');
        body.push('\n');
        for y in 0..self.rows {
            let mut top = String::from("|");
            let mut bottom = String::from("+");
            for x in 0..self.cols {
                let cell = self[(y, x)];
                if cell.has_east
                    && self.is_linked(
                        &cell,
                        &self.get_east(&cell).expect("should have a east cell"),
                    )
                {
                    east_boundary = " ".to_string();
                } else {
                    east_boundary = "|".to_string();
                }
                if cell.has_south
                    && self.is_linked(
                        &cell,
                        &self.get_south(&cell).expect("should have a south cell"),
                    )
                {
                    south_boundary = south_boundary_spaces.clone();
                } else {
                    south_boundary = south_boundary_ticks.clone();
                }
                let content = self.contents_of(&cell, padding);
                top.push_str(&content.to_owned());
                top.push_str(&east_boundary.to_owned());
                bottom.push_str(&south_boundary.to_owned());
                bottom.push('+');
            }
            body.push_str(&top.to_owned());
            body.push('\n');
            body.push_str(&bottom.to_owned());
            body.push('\n');
        }
        writeln!(f, "{}", body)
    }
}

impl Index<(usize, usize)> for Grid {
    type Output = Cell;
    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        if row > self.rows || col > self.cols {
            panic!("Index out of bounds");
        }
        &self.cells[self.cols * row + col]
    }
}

impl IndexMut<(usize, usize)> for Grid {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Cell {
        if row > self.rows || col > self.cols {
            panic!("Index out of bounds");
        }
        &mut self.cells[self.cols * row + col]
    }
}

impl Grid {
    pub fn new(rows: usize, cols: usize) -> Grid {
        let mut cells = vec![];
        let cell_bodies = HashMap::new();
        let mut cell_links = HashMap::new();
        for r in 0..rows {
            for c in 0..cols {
                let cell = Cell {
                    pos: Position { row: r, col: c },
                    has_north: r != 0,
                    has_south: (r + 1) != rows,
                    has_west: c != 0,
                    has_east: (c + 1) != cols,
                };
                cells.push(cell);
                cell_links.insert(cell, HashSet::new());
            }
        }
        Grid {
            rows,
            cols,
            cells,
            cell_links,
            cell_bodies,
        }
    }

    pub fn size(&self) -> usize {
        self.rows * self.cols
    }

    pub fn get_north(&self, cell: &Cell) -> Option<Cell> {
        if cell.has_north {
            Some(self[(cell.pos.row - 1, cell.pos.col)])
        } else {
            None
        }
    }

    pub fn get_south(&self, cell: &Cell) -> Option<Cell> {
        if cell.has_south {
            Some(self[(cell.pos.row + 1, cell.pos.col)])
        } else {
            None
        }
    }

    pub fn get_west(&self, cell: &Cell) -> Option<Cell> {
        if cell.has_west {
            Some(self[(cell.pos.row, cell.pos.col - 1)])
        } else {
            None
        }
    }

    pub fn get_east(&self, cell: &Cell) -> Option<Cell> {
        if cell.has_east {
            Some(self[(cell.pos.row, cell.pos.col + 1)])
        } else {
            None
        }
    }

    fn contents_of(&self, cell: &Cell, padding: usize) -> String {
        if let Some(body) = self.cell_bodies.get(cell) {
            body.to_string()
        } else {
            // treat nothing as one digit
            let mut body = " ".repeat(padding);
            if (padding - 1) == 0 {
                body.push(' ');
            }
            for _ in 0..padding - 1 {
                body.push(' ');
            }
            body.push(' ');
            body.push_str(&" ".repeat(padding));
            body
        }
    }

    fn draw_background(&self, img: &mut RgbImage, cell_size: usize, max_distance: usize) {
        for cell in &self.cells {
            let x1 = cell.pos.col * cell_size;
            let y1 = cell.pos.row * cell_size;
            let mut x2 = (cell.pos.col + 1) * cell_size;
            let background = self.background_color(cell, max_distance);
            let color = image::Rgb::from_slice(&background);
            if x2 == img.width() as usize {
                x2 -= 1;
            }
            let mut y2 = (cell.pos.row + 1) * cell_size;
            if y2 == img.height() as usize {
                y2 -= 1;
            }

            for y in y1..=y2 {
                for x in x1..=x2 {
                    img.put_pixel(y as u32, x as u32, *color);
                }
            }
        }
    }

    fn draw_cells(&self, img: &mut RgbImage, cell_size: usize) {
        let background: [u8; 3] = [0, 0, 0];
        let color = image::Rgb::from_slice(&background);
        for cell in &self.cells {
            let x1 = cell.pos.col * cell_size;
            let y1 = cell.pos.row * cell_size;
            let mut x2 = (cell.pos.col + 1) * cell_size;
            if x2 == img.width() as usize {
                x2 -= 1;
            }
            let mut y2 = (cell.pos.row + 1) * cell_size;
            if y2 == (img.height() as usize) {
                y2 -= 1;
            }

            if !cell.has_north {
                for x in x1..=x2 {
                    img.put_pixel(x as u32, y1 as u32, *color);
                }
            }

            if !cell.has_west {
                for y in y1..=y2 {
                    img.put_pixel(x1 as u32, y as u32, *color);
                }
            }

            if let Some(east) = self.get_east(cell) {
                if !self.is_linked(cell, &east) {
                    for y in y1..=y2 {
                        img.put_pixel(x2 as u32, y as u32, *color);
                    }
                }
            }

            if let Some(south) = self.get_south(cell) {
                if !self.is_linked(cell, &south) {
                    for x in x1..=x2 {
                        img.put_pixel(x as u32, y2 as u32, *color);
                    }
                }
            }
        }
    }

    pub fn to_png(&self, path: String, cell_size: usize, max_distance: usize) {
        let mut img: RgbImage = ImageBuffer::new(
            (self.cols * cell_size) as u32,
            (self.rows * cell_size) as u32,
        );
        self.draw_background(&mut img, cell_size, max_distance);
        self.draw_cells(&mut img, cell_size);
        img.save(path).expect("Failed to save png image");
    }

    fn background_color(&self, cell: &Cell, max_distance: usize) -> [u8; 3] {
        if let Some(distance) = self.cell_bodies.get(cell) {
            let distance = distance.trim().parse().unwrap_or(0);
            let intensity: u8 = ((max_distance - distance) as f64 / max_distance as f64) as u8;
            let dark: u8 = ((255 * intensity) as f64).round() as u8;
            let bright: u8 = (128.0 + (127 * intensity) as f64).round() as u8;
            [dark, bright, dark]
        } else {
            [0, 0, 0]
        }
    }

    fn random_cell(&self) -> Cell {
        let mut rng = rand::thread_rng();
        let col = rng.gen_range(0..self.cols);
        let row = rng.gen_range(0..self.rows);
        self[(row, col)]
    }

    fn is_linked(&self, cell1: &Cell, cell2: &Cell) -> bool {
        let mut linked = false;
        if let Some(links) = self.cell_links.get(cell1) {
            linked = links.get(cell2).is_some();
        }
        linked
    }

    fn link_cells(&mut self, cell1: &Cell, cell2: &Cell, bidi: bool) {
        let mut new_links = HashSet::new();
        if let Some(links) = self.cell_links.get(cell1) {
            if !links.contains(cell2) {
                new_links = links.clone();
                new_links.insert(*cell2);
                self.cell_links.insert(*cell1, new_links);
            }
        } else {
            new_links.insert(*cell2);
            self.cell_links.insert(*cell1, new_links);
        }

        if bidi {
            self.link_cells(cell2, cell1, false);
        }
    }

    fn unlink_cells(&mut self, cell1: &Cell, cell2: &Cell, bidi: bool) {
        let links = self.cell_links.get(cell1).unwrap();
        if links.contains(cell2) {
            let mut new_links: HashSet<Cell> = links.clone();
            new_links.remove(cell2);
            self.cell_links.insert(*cell1, new_links);
        }
        if bidi {
            self.link_cells(cell2, cell1, false);
        }
    }

    fn neighbors(&self, cell: &Cell) -> Vec<Cell> {
        let mut neighbors = vec![];
        if cell.has_north {
            neighbors.push(self.get_north(cell).expect("should have a north"));
        }
        if cell.has_south {
            neighbors.push(self.get_south(cell).expect("should have a south"));
        }
        if cell.has_east {
            neighbors.push(self.get_east(cell).expect("should have a east"));
        }
        if cell.has_west {
            neighbors.push(self.get_west(cell).expect("should have a west"));
        }
        neighbors
    }

    pub fn deadends(&self) -> Vec<Cell> {
        let mut list = vec![];
        for cell in &self.cells {
            if let Some(links) = self.cell_links.get(cell) {
                if links.len() == 1 {
                    list.push(*cell);
                }
            }
        }
        list
    }

    pub fn set_cell_bodies(&mut self, distance_table: &DistanceTable) {
        let max = distance_table.max_table_entry();
        let padding = max.1.to_string().len();
        for (k, v) in distance_table.table.iter() {
            let mut body = " ".repeat(padding);
            let vv = v;
            if vv.to_string().trim().len() < padding {
                for _ in 0..padding - 1 {
                    body.push(' ');
                }
            }
            body.push_str(v.to_string().trim());
            body.push_str(&" ".repeat(padding));
            self.cell_bodies.insert(*k, body.to_string());
        }
    }
}

impl Mask {
    pub fn new(rows: usize, cols: usize) -> Mask {
        let bits = (0..(rows * cols)).map(|_| true).collect();
        Mask { rows, cols, bits }
    }

    pub fn count(&self) -> usize {
        let mut c: usize = 0;
        for row in 0..self.rows {
            for col in 0..self.cols {
                if self[(row, col)] {
                    c += 1;
                }
            }
        }
        c
    }

    pub fn random_location(&self) -> (usize, usize) {
        let mut rng = rand::thread_rng();
        loop {
            let row = rng.gen_range(0..self.rows);
            let col = rng.gen_range(0..self.cols);
            if self[(row, col)] {
                return (row, col);
            }
        }
    }
}

impl Index<(usize, usize)> for Mask {
    type Output = bool;
    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        if row > self.rows || col > self.cols {
            panic!("Index out of bounds");
        }
        &self.bits[self.cols * row + col]
    }
}

impl IndexMut<(usize, usize)> for Mask {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut bool {
        if row > self.rows || col > self.cols {
            panic!("Index out of bounds");
        }
        &mut self.bits[self.cols * row + col]
    }
}

#[test]
fn new_grid() {
    let wanted_cells = vec![
        Cell {
            pos: Position { row: 0, col: 0 },
            has_north: false,
            has_south: true,
            has_west: false,
            has_east: true,
        },
        Cell {
            pos: Position { row: 0, col: 1 },
            has_north: false,
            has_south: true,
            has_west: true,
            has_east: false,
        },
        Cell {
            pos: Position { row: 1, col: 0 },
            has_north: true,
            has_south: false,
            has_west: false,
            has_east: true,
        },
        Cell {
            pos: Position { row: 1, col: 1 },
            has_north: true,
            has_south: false,
            has_west: true,
            has_east: false,
        },
    ];
    let g = Grid::new(2, 2);
    println!("{}", g);
    assert!(g.rows == 3);
    assert!(g.cols == 2);
    assert!(g.cells[0] == wanted_cells[0]);
    assert!(g.cells[1] == wanted_cells[1]);
    assert!(g.cells[2] == wanted_cells[2]);
    assert!(g.cells[3] == wanted_cells[3]);
}
