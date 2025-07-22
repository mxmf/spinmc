mod metropolis;
mod stats;
use crate::lattice::Grid;
use crate::spin::SpinState;

pub use metropolis::Metropolis;
pub use stats::Stats;

pub trait MonteCarlo<S: SpinState, R: rand::Rng> {
    fn step(&mut self, grid: &mut Grid<S, R>);
}
