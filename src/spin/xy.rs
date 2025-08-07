use rand_distr::{Distribution, UnitCircle};
use std::{
    iter::Sum,
    ops::{Add, AddAssign, Div, Mul, Neg, Sub},
};

use crate::spin::SpinState;

#[derive(Default, Debug, Clone, Copy)]
#[repr(C)]
#[cfg_attr(feature = "snapshots", derive(hdf5_metno::H5Type))]
pub struct XYSpin {
    x: f64,
    y: f64,
}

impl SpinState for XYSpin {
    fn zero() -> Self {
        Self { x: 0., y: 0. }
    }
    fn along_x(magnitude: f64) -> Self {
        Self {
            x: magnitude,
            y: 0.,
        }
    }

    fn along_y(magnitude: f64) -> Self {
        Self {
            x: 0.,
            y: magnitude,
        }
    }
    fn along_z(_magnitude: f64) -> Self {
        panic!("XYSpin does not support creating spins along the z-axis")
    }
    fn random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self {
        let [x, y]: [f64; 2] = UnitCircle.sample(rng);
        Self {
            x: x * magnitude,
            y: y * magnitude,
        }
    }

    fn perturb<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self {
        Self::random(rng, magnitude)
    }

    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    fn norm(&self) -> f64 {
        self.norm_sqr().sqrt()
    }

    fn norm_sqr(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    fn is_aligned(&self, _axis: &Self) -> bool {
        true
    }

    fn flip(&mut self, axis: &Self) -> Self {
        *self - *axis * 2. * (self.dot(axis))
    }

    fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, 0.]
    }
}

// +=
impl AddAssign<XYSpin> for XYSpin {
    fn add_assign(&mut self, rhs: XYSpin) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}
impl AddAssign<&XYSpin> for XYSpin {
    fn add_assign(&mut self, rhs: &XYSpin) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Sum for XYSpin {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl Add for XYSpin {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}
impl Neg for XYSpin {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl Sub for XYSpin {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<f64> for XYSpin {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Div<f64> for XYSpin {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Div<f64> for &XYSpin {
    type Output = XYSpin;
    fn div(self, rhs: f64) -> Self::Output {
        XYSpin {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}
