use rand_distr::{Distribution, UnitSphere};
use std::{
    iter::Sum,
    ops::{Add, AddAssign, Div, Mul, Neg, Sub},
};

use crate::spin::SpinState;

#[derive(Default, Debug, Clone, Copy)]
#[repr(C)]
#[cfg_attr(feature = "snapshots", derive(hdf5_metno::H5Type))]
pub struct HeisenbergSpin {
    x: f64,
    y: f64,
    z: f64,
}

impl SpinState for HeisenbergSpin {
    fn zero() -> Self {
        Self {
            x: 0.,
            y: 0.,
            z: 0.,
        }
    }
    fn along_x(magnitude: f64) -> anyhow::Result<Self> {
        Ok(Self {
            x: magnitude,
            y: 0.,
            z: 0.,
        })
    }

    fn along_y(magnitude: f64) -> anyhow::Result<Self> {
        Ok(Self {
            x: 0.,
            y: magnitude,
            z: 0.,
        })
    }
    fn along_z(magnitude: f64) -> anyhow::Result<Self> {
        Ok(Self {
            x: 0.,
            y: 0.,
            z: magnitude,
        })
    }
    fn random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self {
        let [x, y, z]: [f64; 3] = UnitSphere.sample(rng);
        Self {
            x: x * magnitude,
            y: y * magnitude,
            z: z * magnitude,
        }
    }

    fn perturb<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self {
        Self::random(rng, magnitude)
    }

    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn norm(&self) -> f64 {
        self.norm_sqr().sqrt()
    }

    fn norm_sqr(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn is_aligned(&self, _axis: &Self) -> bool {
        true
    }

    fn flip(&mut self, axis: &Self) -> Self {
        *self - *axis * 2. * (self.dot(axis))
    }

    fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

// +=
impl AddAssign<HeisenbergSpin> for HeisenbergSpin {
    fn add_assign(&mut self, rhs: HeisenbergSpin) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z
    }
}
impl AddAssign<&HeisenbergSpin> for HeisenbergSpin {
    fn add_assign(&mut self, rhs: &HeisenbergSpin) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z
    }
}

impl Sum for HeisenbergSpin {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl Add for HeisenbergSpin {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}
impl Neg for HeisenbergSpin {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Sub for HeisenbergSpin {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for HeisenbergSpin {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f64> for HeisenbergSpin {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Div<f64> for &HeisenbergSpin {
    type Output = HeisenbergSpin;
    fn div(self, rhs: f64) -> Self::Output {
        HeisenbergSpin {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}
