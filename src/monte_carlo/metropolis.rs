use crate::spin::SpinState;
use rand::RngExt;

use super::MonteCarlo;

pub struct Metropolis<R: rand::Rng> {
    pub rng: R,
    pub beta: f64,
}

impl<S: SpinState, R: rand::Rng> MonteCarlo<S, R> for Metropolis<R> {
    fn step(&mut self, grid: &mut crate::lattice::Grid<S, R>) -> usize {
        for i in 0..grid.size {
            let proposed_spin = grid.spins[i].perturb(&mut self.rng, grid.calc_inputs[i].magnitude);
            let delta_e = proposed_spin.energy_diff(
                &grid.calc_inputs[i],
                &grid.hamiltonian,
                &grid.spins,
                &grid.spins[i],
            );
            let spin = &mut grid.spins[i];
            if accepts_metropolis_move(delta_e, self.beta, &mut self.rng) {
                *spin = proposed_spin;
            }
        }
        grid.size
    }
}

fn accepts_metropolis_move<R: rand::Rng>(delta_e: f64, beta: f64, rng: &mut R) -> bool {
    if delta_e < 0.0 {
        true
    } else if beta.is_infinite() {
        false
    } else {
        rng.random::<f64>() < (-beta * delta_e).exp()
    }
}

#[cfg(test)]
#[path = "metropolis_tests.rs"]
mod tests;
