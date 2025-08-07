use crate::{calculators::HamiltonianConfig, spin::SpinState};
use std::collections::VecDeque;

use super::MonteCarlo;
pub struct Wolff<R: rand::Rng> {
    pub rng: R,
    pub beta: f64,
    pub ham_config: HamiltonianConfig,
}
impl<S: SpinState, R: rand::Rng> MonteCarlo<S, R> for Wolff<R> {
    fn step(&mut self, grid: &mut crate::lattice::Grid<S, R>) -> usize {
        let init_spin_index = self.rng.random_range(0..grid.size);

        let axis = -grid.spins[init_spin_index].perturb(&mut self.rng, 1.0);

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

                if !neighbor_spin.is_aligned(&grid.spins[site]) {
                    continue;
                }

                let p = grid.spins[site].wolff_probability(
                    neighbor_spin,
                    &axis,
                    self.beta,
                    *j,
                    grid.calc_inputs[site].magnitude,
                    grid.calc_inputs[*neighbor].magnitude,
                );
                if self.rng.random::<f64>() < p {
                    visited[*neighbor] = true;
                    queue.push_back(*neighbor);
                }
            }
        }

        if self.ham_config.anisotropy_enable {
            let mut delta_e = 0.;
            for index in &cluster {
                let flip_spin = grid.spins[*index].flip(&axis);

                delta_e += flip_spin.ion_anisotropy_energy_diff(
                    &grid.calc_inputs[*index],
                    &grid.hamiltonian,
                    &grid.spins[*index],
                )
            }

            if delta_e >= 0.0 && self.rng.random::<f64>() > (-self.beta * delta_e).exp() {
                return 0;
            }
        }

        for index in &cluster {
            grid.spins[*index] = grid.spins[*index].flip(&axis);
        }

        // if self.ham_config.anisotropy_enable {
        //     if delta_e >= 0.0 && self.rng.random::<f64>() > (-self.beta * delta_e).exp() {
        //         for index in &cluster {
        //             grid.spins[*index] = grid.spins[*index].flip(&axis);
        //         }
        //     }
        // }

        cluster.len()
    }
}
