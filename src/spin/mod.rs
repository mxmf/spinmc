pub mod heisenberg;
pub mod ising;
pub mod xy;

use std::ops::{Add, Div, Mul};
use std::{iter::Sum, ops::AddAssign};

use crate::calculators::{CalcInput, Hamiltonian};
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

#[derive(Clone, Debug)]
pub enum SpinVector {
    Ising(f64),
    XY(f64, f64),
    Heisenberg(f64, f64, f64),
}

impl AddAssign<&SpinVector> for SpinVector {
    fn add_assign(&mut self, rhs: &Self) {
        match (self, rhs) {
            (SpinVector::Ising(a), SpinVector::Ising(b)) => {
                *a += b;
            }
            (SpinVector::XY(ax, ay), SpinVector::XY(bx, by)) => {
                *ax += bx;
                *ay += by;
            }
            (SpinVector::Heisenberg(ax, ay, az), SpinVector::Heisenberg(bx, by, bz)) => {
                *ax += bx;
                *ay += by;
                *az += bz;
            }
            _ => panic!("Cannot add-assign different spin types"),
        }
    }
}

impl SpinVector {
    pub fn zero(kind: &SpinVector) -> Self {
        match kind {
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
    pub fn norm_sqr(&self) -> f64 {
        match self {
            SpinVector::Ising(s) => s.powi(2),
            SpinVector::XY(x, y) => x.powi(2) + y.powi(2),
            SpinVector::Heisenberg(x, y, z) => x.powi(2) + y.powi(2) + z.powi(2),
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

impl Mul<f64> for SpinVector {
    type Output = SpinVector;

    fn mul(self, rhs: f64) -> SpinVector {
        match self {
            SpinVector::Ising(s) => SpinVector::Ising(s * rhs),
            SpinVector::XY(x, y) => SpinVector::XY(x * rhs, y * rhs),
            SpinVector::Heisenberg(x, y, z) => SpinVector::Heisenberg(x * rhs, y * rhs, z * rhs),
        }
    }
}

impl Div<f64> for &SpinVector {
    type Output = SpinVector;

    fn div(self, rhs: f64) -> SpinVector {
        if rhs == 0.0 {
            panic!("Divide by zero in SpinVector");
        }
        match self {
            SpinVector::Ising(s) => SpinVector::Ising(s / rhs),
            SpinVector::XY(x, y) => SpinVector::XY(x / rhs, y / rhs),
            SpinVector::Heisenberg(x, y, z) => SpinVector::Heisenberg(x / rhs, y / rhs, z / rhs),
        }
    }
}

impl Div<f64> for SpinVector {
    type Output = SpinVector;

    fn div(self, rhs: f64) -> SpinVector {
        match self {
            SpinVector::Ising(s) => SpinVector::Ising(s / rhs),
            SpinVector::XY(x, y) => SpinVector::XY(x / rhs, y / rhs),
            SpinVector::Heisenberg(x, y, z) => SpinVector::Heisenberg(x / rhs, y / rhs, z / rhs),
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

        result.unwrap()
    }
}

pub trait SpinState: Default + Clone + Send + Sync + H5Type + 'static {
    fn new_x(magnitude: f64) -> Self;
    fn new_y(magnitude: f64) -> Self;
    fn new_z(magnitude: f64) -> Self;
    fn new_random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self;

    fn zero() -> SpinVector;

    fn magnitude(&self) -> f64;

    fn direction(&self) -> SpinVector;

    fn spinvector(&self) -> SpinVector;

    fn random<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self;

    fn propose_perturbation<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self;

    fn dot(&self, other: &Self) -> f64;

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
}
