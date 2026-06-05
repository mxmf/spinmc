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

impl<R: rand::Rng> AnyMC<R> {
    pub fn beta(&self) -> f64 {
        match self {
            AnyMC::Metropolis(m) => m.beta,
            AnyMC::Wolff(w) => w.beta,
        }
    }
    pub fn set_beta(&mut self, beta: f64) {
        match self {
            AnyMC::Metropolis(m) => m.beta = beta,
            AnyMC::Wolff(w) => w.beta = beta,
        }
    }
}
