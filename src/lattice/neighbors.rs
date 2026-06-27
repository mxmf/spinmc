pub trait Vector3Ext {
    fn scale(&self, scalar: f64) -> [f64; 3];
    fn dot(&self, other: &[f64; 3]) -> f64;
    fn norm(&self) -> f64;
    fn add(&self, other: &[f64; 3]) -> [f64; 3];
    fn sub(&self, other: &[f64; 3]) -> [f64; 3];
}

impl Vector3Ext for [f64; 3] {
    fn scale(&self, scalar: f64) -> [f64; 3] {
        [self[0] * scalar, self[1] * scalar, self[2] * scalar]
    }

    fn dot(&self, other: &[f64; 3]) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
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
    pub neighbor: Neighbor,
    pub distance: f64,
}

impl Atoms {
    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn reciprocal_cell(&self) -> Option<[[f64; 3]; 3]> {
        let volume = self.cell[0].dot(&Self::cross(self.cell[1], self.cell[2]));
        if volume.abs() <= f64::EPSILON {
            return None;
        }

        Some([
            Self::cross(self.cell[1], self.cell[2]).scale(1.0 / volume),
            Self::cross(self.cell[2], self.cell[0]).scale(1.0 / volume),
            Self::cross(self.cell[0], self.cell[1]).scale(1.0 / volume),
        ])
    }

    fn offset_bounds_in_radius(&self, delta: [f64; 3], radius: f64) -> [[isize; 2]; 3] {
        let Some(reciprocal) = self.reciprocal_cell() else {
            return [[0, 0]; 3];
        };
        let mut bounds = [[0isize, 0isize]; 3];

        for axis in 0..3 {
            if !self.pbc[axis] {
                continue;
            }

            let center = reciprocal[axis].dot(&delta);
            let width = reciprocal[axis].norm() * radius;
            bounds[axis] = [
                (center - width).floor() as isize - 1,
                (center + width).ceil() as isize + 1,
            ];
        }

        bounds
    }

    fn cell_offset(&self, offset: [isize; 3]) -> [f64; 3] {
        self.cell[0]
            .scale(offset[0] as f64)
            .add(&self.cell[1].scale(offset[1] as f64))
            .add(&self.cell[2].scale(offset[2] as f64))
    }

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

    pub fn calc_distance_range_from_to(
        &self,
        from: usize,
        to: usize,
        min_distance: f64,
        max_distance: f64,
    ) -> Vec<Distance> {
        let delta = self.positions[from].sub(&self.positions[to]);
        let bounds = self.offset_bounds_in_radius(delta, max_distance);
        let mut result = vec![];

        for (x_i, y_i, z_i) in iproduct!(
            bounds[0][0]..=bounds[0][1],
            bounds[1][0]..=bounds[1][1],
            bounds[2][0]..=bounds[2][1]
        ) {
            let offset = [x_i, y_i, z_i];
            if offset == [0, 0, 0] && from == to {
                continue;
            }

            let distance = delta.sub(&self.cell_offset(offset)).norm();
            if distance >= min_distance - self.tolerance
                && distance <= max_distance + self.tolerance
            {
                result.push(Distance {
                    neighbor: Neighbor { from, to, offset },
                    distance,
                });
            }
        }

        result
    }

    pub fn calc_distance_range_from(
        &self,
        from: usize,
        min_distance: f64,
        max_distance: f64,
    ) -> Vec<Distance> {
        let mut result = vec![];
        for to in 0..self.positions.len() {
            result.extend(self.calc_distance_range_from_to(from, to, min_distance, max_distance));
        }
        result
    }

    pub fn calc_distance_range_all(&self, min_distance: f64, max_distance: f64) -> Vec<Distance> {
        let mut result = vec![];
        for from in 0..self.positions.len() {
            result.extend(self.calc_distance_range_from(from, min_distance, max_distance));
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

    fn initial_search_radius(&self) -> f64 {
        self.cell.iter().map(|v| v.norm()).fold(1e-6, f64::max)
    }

    fn find_neighbors_by_radius<F>(&self, order: usize, distances: F) -> Vec<Neighbor>
    where
        F: Fn(f64) -> Vec<Distance>,
    {
        if order == 0 {
            return vec![];
        }

        let mut radius = self.initial_search_radius();
        for _ in 0..32 {
            let result = self.get_neighbor_from(distances(radius), order);
            if !result.is_empty() {
                return result;
            }
            radius *= 2.0;
        }
        vec![]
    }

    pub fn find_neighbors_from(&self, index: usize, order: usize) -> Vec<Neighbor> {
        self.find_neighbors_by_radius(order, |radius| {
            self.calc_distance_range_from(index, 0.0, radius)
        })
    }

    pub fn find_neighbors_from_to(
        &self,
        index: usize,
        to_index: usize,
        order: usize,
    ) -> Vec<Neighbor> {
        self.find_neighbors_by_radius(order, |radius| {
            self.calc_distance_range_from_to(index, to_index, 0.0, radius)
        })
    }

    pub fn find_neighbors_all(&self, order: usize) -> Vec<Neighbor> {
        self.find_neighbors_by_radius(order, |radius| self.calc_distance_range_all(0.0, radius))
    }
}

#[cfg(test)]
#[path = "neighbors_tests.rs"]
mod tests;
