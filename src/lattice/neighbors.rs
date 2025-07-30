pub trait Vector3Ext {
    fn scale(&self, scalar: f64) -> [f64; 3];
    fn norm(&self) -> f64;
    fn add(&self, other: &[f64; 3]) -> [f64; 3];
    fn sub(&self, other: &[f64; 3]) -> [f64; 3];
}

impl Vector3Ext for [f64; 3] {
    fn scale(&self, scalar: f64) -> [f64; 3] {
        [self[0] * scalar, self[1] * scalar, self[2] * scalar]
    }

    fn norm(&self) -> f64 {
        (self[0].powi(2) + self[1].powi(2) + self[2].powi(2)).sqrt()
    }

    fn add(&self, other: &[f64; 3]) -> [f64; 3] {
        [self[0] + other[0], self[1] + other[1], self[2] + other[2]]
    }

    fn sub(&self, other: &[f64; 3]) -> [f64; 3] {
        [self[0] - other[0], self[1] - other[1], self[2] - other[2]]
    }
}
use itertools::iproduct;

pub struct Atoms {
    pub cell: [[f64; 3]; 3],
    pub positions: Vec<[f64; 3]>,
    pub pbc: [bool; 3],
    pub tolerance: f64,
}

#[derive(Debug, Clone)]
pub struct Neighbor {
    pub from: usize,
    pub to: usize,
    pub offset: [isize; 3],
}

#[derive(Debug, Clone)]
pub struct Distance {
    neighbor: Neighbor,
    pub distance: f64,
}

impl Atoms {
    pub fn calc_distance_from_to(&self, from: usize, to: usize, max_n: isize) -> Vec<Distance> {
        let mut result = vec![];

        let x_i_range = if self.pbc[0] { -max_n..max_n + 1 } else { 0..1 };
        let y_i_range = if self.pbc[1] { -max_n..max_n + 1 } else { 0..1 };
        let z_i_range = if self.pbc[2] { -max_n..max_n + 1 } else { 0..1 };

        for (x_i, y_i, z_i) in iproduct!(x_i_range, y_i_range, z_i_range) {
            if x_i == 0 && y_i == 0 && z_i == 0 && to == from {
                continue;
            }
            let to_position = self.positions[to]
                .add(&self.cell[0].scale(x_i as f64))
                .add(&self.cell[1].scale(y_i as f64))
                .add(&self.cell[2].scale(z_i as f64));

            let distance = Distance {
                neighbor: Neighbor {
                    from,
                    to,
                    offset: [x_i, y_i, z_i],
                },
                distance: self.positions[from].sub(&to_position).norm(), // distance: self.positions[from].sub(&to_position).norm()
            };
            result.push(distance);
        }
        result
    }

    pub fn calc_distance_from(&self, from: usize, max_n: isize) -> Vec<Distance> {
        let mut result = vec![];
        for to in 0..self.positions.len() {
            result.extend(self.calc_distance_from_to(from, to, max_n));
        }
        result
    }

    fn get_neighbor_from(&self, mut neighbors: Vec<Distance>, order: usize) -> Vec<Neighbor> {
        let mut result = vec![];

        neighbors.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap());

        let mut current_order = 0;
        let mut last_distance = 0.;

        for neighbor in neighbors {
            if (neighbor.distance - last_distance).abs() > self.tolerance {
                current_order += 1;
                last_distance = neighbor.distance;
            }

            if current_order == order {
                result.push(neighbor.neighbor);
            } else if current_order > order {
                break;
            }
        }
        result
    }

    pub fn find_neighbors_from(&self, index: usize, order: usize) -> Vec<Neighbor> {
        let neighbors = self.calc_distance_from(index, order as isize);
        // println!("{}", neighbors.len());

        self.get_neighbor_from(neighbors, order)
    }

    pub fn find_neighbors_from_to(
        &self,
        index: usize,
        to_index: usize,
        order: usize,
    ) -> Vec<Neighbor> {
        let neighbors = self.calc_distance_from_to(index, to_index, order as isize);
        // println!("{}", neighbors.len());

        self.get_neighbor_from(neighbors, order)
    }

    pub fn find_neighbors_all(&self, order: usize) -> Vec<Neighbor> {
        let mut neighbors = vec![];

        for i in 0..self.positions.len() {
            neighbors.extend(self.calc_distance_from(i, order as isize));
        }
        self.get_neighbor_from(neighbors, order)
    }
}

#[cfg(test)]
mod tests {
    use crate::lattice::neighbors::Atoms;

    #[test]
    fn test_find_neighbors_simple() {
        let atoms = Atoms {
            cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[0., 0., 0.0]],
            pbc: [true, true, false],
            tolerance: 0.000001,
        };

        let neighbors = atoms.find_neighbors_from(0, 2);

        // println!("{neighbors:?}");
        assert_eq!(neighbors.len(), 4)
    }

    #[test]
    fn test_find_neighbors_simple2() {
        let atoms = Atoms {
            cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
            pbc: [true, true, false],
            tolerance: 0.0001,
        };

        let neighbors = atoms.find_neighbors_from(0, 1);

        // println!("{neighbors:?}");
        assert_eq!(neighbors.len(), 4);
        let neighbors = atoms.find_neighbors_from(0, 2);

        println!("{neighbors:?}");
        assert_eq!(neighbors.len(), 4);
    }

    #[test]
    fn test_find_neighbors_all_simple() {
        let atoms = Atoms {
            cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[0., 0., 0.0]],
            pbc: [true, true, false],
            tolerance: 0.000001,
        };

        let neighbors = atoms.find_neighbors_all(1);

        // println!("{neighbors:?}");
        assert_eq!(neighbors.len(), 4)
    }

    #[test]
    fn test_find_neighbors_all_simple2() {
        let atoms = Atoms {
            cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
            pbc: [true, true, false],
            tolerance: 0.0001,
        };
        let neighbors = atoms.find_neighbors_all(1);

        // println!("{neighbors:?}");
        assert_eq!(neighbors.len(), 16)
    }

    #[test]
    fn test_find_neighbors_from_to_simple2() {
        let atoms = Atoms {
            cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
            pbc: [true, true, false],
            tolerance: 0.0001,
        };
        let neighbors = atoms.find_neighbors_from_to(0, 1, 1);

        // println!("{neighbors:?}");
        assert_eq!(neighbors.len(), 2)
    }
}
