pub mod heisenberg;

pub mod ising;
pub mod xy;

use std::iter::Sum;
use std::ops::Add;

pub use ising::IsingSpin;

use crate::calculators::{CalcInput, Hamiltonian};

#[derive(Clone, Debug)]
pub enum SpinVector {
    Ising(f64),
    XY(f64, f64),
    Heisenberg(f64, f64, f64),
}

impl SpinVector {
    pub fn zero(&self) -> Self {
        match self {
            SpinVector::Ising(_) => SpinVector::Ising(0.0),
            SpinVector::XY(_, _) => SpinVector::XY(0.0, 0.0),
            SpinVector::Heisenberg(_, _, _) => SpinVector::Heisenberg(0.0, 0.0, 0.0),
        }
    }

    pub fn norm(&self) -> f64 {
        match self {
            SpinVector::Ising(s) => s.abs(),
            SpinVector::XY(x, y) => (x.powi(2) + y.powi(2)).sqrt(),
            SpinVector::Heisenberg(x, y, z) => (x.powi(2) + y.powi(2) + z.powi(2)).sqrt(),
        }
    }
}

impl Add for SpinVector {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        match (self, other) {
            (SpinVector::Ising(a), SpinVector::Ising(b)) => SpinVector::Ising(a + b),
            (SpinVector::XY(ax, ay), SpinVector::XY(bx, by)) => SpinVector::XY(ax + bx, ay + by),
            (SpinVector::Heisenberg(ax, ay, az), SpinVector::Heisenberg(bx, by, bz)) => {
                SpinVector::Heisenberg(ax + bx, ay + by, az + bz)
            }
            _ => panic!("Cannot add different spin types"),
        }
    }
}

impl Sum for SpinVector {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = None;

        for item in iter {
            result = match result {
                Some(res) => Some(res + item),
                None => Some(item),
            };
        }

        result.unwrap_or(SpinVector::Heisenberg(0.0, 0.0, 0.0))
    }
}

pub trait SpinState: Clone + Send + Sync + 'static {
    fn new_x(magnitude: f64) -> Self;
    fn new_y(magnitude: f64) -> Self;
    fn new_z(magnitude: f64) -> Self;
    fn new_random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self;

    fn magnitude(&self) -> f64;

    fn direction(&self) -> SpinVector;

    fn spinvector(&self) -> SpinVector;

    fn random<R: rand::Rng>(&self, rng: &mut R) -> Self;

    fn propose_constrained_perturbation<R: rand::Rng>(&self, rng: &mut R) -> Self;

    fn dot(&self, other: &Self) -> f64;

    fn energy(&self, calc_input: &CalcInput, ham: &Hamiltonian, spins: &[Self]) -> f64 {
        ham.compute(self, calc_input, spins)
    }
    fn energy_diff(
        &self,
        calc_input: &CalcInput,
        ham: &Hamiltonian,
        spins: &[Self],
        old_spin: &Self,
    ) -> f64 {
        ham.compute(self, calc_input, spins) - ham.compute(old_spin, calc_input, spins)
    }
}
