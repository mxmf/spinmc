use crate::spin::SpinState;

use super::MonteCarlo;

pub struct Metropolis<R: rand::Rng> {
    pub rng: R,
    pub beta: f64,
}

impl<S: SpinState, R: rand::Rng> MonteCarlo<S, R> for Metropolis<R> {
    fn step(&mut self, grid: &mut crate::lattice::Grid<S, R>) {
        for i in 0..grid.size {
            let proposed_spin =
                grid.spins[i].propose_perturbation(&mut self.rng, grid.calc_inputs[i].magnitude);
            let delta_e = proposed_spin.energy_diff(
                &grid.calc_inputs[i],
                &grid.hamiltonian,
                &grid.spins,
                &grid.spins[i],
            );
            let spin = &mut grid.spins[i];
            if delta_e < 0.0 || self.rng.random::<f64>() < (-self.beta * delta_e).exp() {
                *spin = proposed_spin;
            }
        }
    }
}
