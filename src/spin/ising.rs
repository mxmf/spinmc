use std::{
    iter::Sum,
    ops::{Add, AddAssign, Div, Mul, Neg, Sub},
};

use crate::spin::SpinState;

#[derive(Default, Debug, Clone, Copy)]
#[repr(C)]
#[cfg_attr(feature = "snapshots", derive(hdf5_metno::H5Type))]
pub struct IsingSpin {
    state: f64,
}

impl SpinState for IsingSpin {
    fn zero() -> Self {
        Self { state: 0. }
    }
    fn along_x(_magnitude: f64) -> Self {
        panic!("IsingSpin does not support creating spins along the x-axis")
    }

    fn along_y(_magnitude: f64) -> Self {
        panic!("IsingSpin does not support creating spins along the y-axis")
    }
    fn along_z(magnitude: f64) -> Self {
        Self { state: magnitude }
    }
    fn random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self {
        let value = if rng.random_bool(0.5) {
            magnitude
        } else {
            -magnitude
        };
        Self { state: value }
    }

    fn perturb<R: rand::Rng>(&self, _rng: &mut R, _magnitude: f64) -> Self {
        Self { state: -self.state }
    }

    fn dot(&self, other: &Self) -> f64 {
        self.state * other.state
    }

    fn norm(&self) -> f64 {
        self.state.abs()
    }

    fn norm_sqr(&self) -> f64 {
        self.state * self.state
    }

    fn energy_diff(
        &self,
        calc_input: &crate::calculators::CalcInput<IsingSpin>,
        ham: &crate::calculators::Hamiltonian,
        spins: &[Self],
        _old_spin: &Self,
    ) -> f64 {
        2. * self.local_energy(calc_input, ham, spins)
    }

    fn is_aligned(&self, axis: &Self) -> bool {
        self.state.signum() == axis.state.signum()
    }

    fn wolff_probability(
        &self,
        _other: &Self,
        _axis: &Self,
        beta: f64,
        j: f64,
        self_magnitude: f64,
        other_magnitude: f64,
    ) -> f64 {
        1.0 - (-2.0 * beta * j * self_magnitude * other_magnitude).exp()
    }

    fn flip(&mut self, _axis: &Self) {
        self.state = -self.state
    }

    fn to_array(&self) -> [f64; 3] {
        [0., 0., self.state]
    }
}

// +=
impl AddAssign<IsingSpin> for IsingSpin {
    fn add_assign(&mut self, rhs: IsingSpin) {
        self.state += rhs.state
    }
}
impl AddAssign<&IsingSpin> for IsingSpin {
    fn add_assign(&mut self, rhs: &IsingSpin) {
        self.state += rhs.state;
    }
}

impl Sum for IsingSpin {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl Add for IsingSpin {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            state: self.state + other.state,
        }
    }
}
impl Neg for IsingSpin {
    type Output = Self;
    fn neg(self) -> Self {
        Self { state: -self.state }
    }
}

impl Sub for IsingSpin {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            state: self.state - other.state,
        }
    }
}

impl Mul<f64> for IsingSpin {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            state: self.state * rhs,
        }
    }
}

impl Div<f64> for IsingSpin {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self {
            state: self.state / rhs,
        }
    }
}

impl Div<f64> for &IsingSpin {
    type Output = IsingSpin;
    fn div(self, rhs: f64) -> Self::Output {
        IsingSpin {
            state: self.state / rhs,
        }
    }
}
