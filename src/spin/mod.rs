mod heisenberg;
mod ising;
mod xy;

use crate::calculators::{CalcInput, Hamiltonian};
use std::ops::{Add, Div, Mul};
use std::ops::{Neg, Sub};
use std::{iter::Sum, ops::AddAssign};

pub use heisenberg::HeisenbergSpin;
pub use ising::IsingSpin;
pub use xy::XYSpin;

#[cfg(not(feature = "snapshots"))]
mod private {
    pub trait H5Type {}
    impl<T> H5Type for T {}
}
#[cfg(not(feature = "snapshots"))]
use private::H5Type;

#[cfg(feature = "snapshots")]
use hdf5_metno::H5Type;

pub trait SpinState:
    Default
    + Clone
    + Copy
    + Send
    + Sync
    + H5Type
    + 'static
    + Add
    + AddAssign
    + for<'a> AddAssign<&'a Self>
    + Neg<Output = Self>
    + Sub
    + Div<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Sum
{
    fn zero() -> Self;
    fn along_x(magnitude: f64) -> Self;
    fn along_y(magnitude: f64) -> Self;
    fn along_z(magnitude: f64) -> Self;
    fn random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self;

    fn perturb<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self;

    fn dot(&self, other: &Self) -> f64;
    fn norm(&self) -> f64;
    fn norm_sqr(&self) -> f64;

    fn energy(&self, calc_input: &CalcInput<Self>, ham: &Hamiltonian, spins: &[Self]) -> f64 {
        ham.compute(self, calc_input, spins)
    }

    fn energy_diff(
        &self,
        calc_input: &CalcInput<Self>,
        ham: &Hamiltonian,
        spins: &[Self],
        old_spin: &Self,
    ) -> f64 {
        ham.compute(self, calc_input, spins) - ham.compute(old_spin, calc_input, spins)
    }

    fn is_aligned(&self, axis: &Self) -> bool;

    fn flip(&mut self, axis: &Self);
}
