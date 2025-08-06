mod metropolis;
mod stats;
mod wolff;
use crate::lattice::Grid;
use crate::spin::SpinState;

pub use metropolis::Metropolis;
pub use stats::{StatResult, Stats, StatsConfig};
pub use wolff::Wolff;

pub trait MonteCarlo<S: SpinState, R: rand::Rng> {
    fn step(&mut self, grid: &mut Grid<S, R>) -> usize;
}

pub enum AnyMC<R: rand::Rng> {
    Metropolis(Metropolis<R>),
    Wolff(Wolff<R>),
}

impl<S: SpinState, R: rand::Rng> MonteCarlo<S, R> for AnyMC<R> {
    fn step(&mut self, grid: &mut crate::lattice::Grid<S, R>) -> usize {
        match self {
            AnyMC::Metropolis(mc) => mc.step(grid),
            AnyMC::Wolff(mc) => mc.step(grid),
        }
    }
}
