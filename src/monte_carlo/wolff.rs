use crate::spin::{SpinState, SpinVector};
use std::collections::VecDeque;

use super::MonteCarlo;
pub struct Wolff<R: rand::Rng> {
    pub rng: R,
    pub beta: f64,
}
impl<S: SpinState, R: rand::Rng> MonteCarlo<S, R> for Wolff<R> {
    fn step(&mut self, grid: &mut crate::lattice::Grid<S, R>) -> usize {
        let init_spin_index = self.rng.random_range(0..grid.size);
        let axis = match grid.spins[init_spin_index].spinvector() {
            SpinVector::Ising(_) => grid.spins[init_spin_index].direction(),
            _ => grid.spins[init_spin_index]
                .spinvector()
                .random(&mut self.rng),
        };
        let mut visited = vec![false; grid.size];
        let mut cluster = Vec::new();
        let mut queue = VecDeque::new();
        visited[init_spin_index] = true;
        queue.push_back(init_spin_index);

        while let Some(site) = queue.pop_front() {
            cluster.push(site);
            for (neighbor, j) in grid.calc_inputs[site]
                .exchange_neighbor_index
                .iter()
                .zip(grid.calc_inputs[site].exchanges.iter())
            {
                if visited[*neighbor] {
                    continue;
                }

                let neighbor_spin = &grid.spins[*neighbor];

                if let SpinVector::Ising(_) = grid.spins[init_spin_index].spinvector() {
                    if neighbor_spin.dot(&grid.spins[site]) < 0. {
                        continue;
                    }
                };

                let p = 1.0
                    - (-2.0
                        * self.beta
                        * j
                        * (grid.spins[site].spinvector().dot(&axis))
                        * (neighbor_spin.spinvector().dot(&axis)))
                    .exp();

                if self.rng.random::<f64>() < p {
                    visited[*neighbor] = true;
                    queue.push_back(*neighbor);
                }
            }
        }

        for index in &cluster {
            grid.spins[*index].flip(&axis);
        }

        cluster.len()
    }
}
