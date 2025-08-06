use crate::spin::{SpinState, SpinVector};
use rand_distr::{Distribution, UnitSphere};

#[derive(Default, Debug, Clone)]
#[repr(C)]
#[cfg_attr(feature = "snapshots", derive(hdf5_metno::H5Type))]
pub struct HeisenbergSpin {
    state: [f64; 3],
}

impl SpinState for HeisenbergSpin {
    fn zero() -> SpinVector {
        SpinVector::Heisenberg(0., 0., 0.)
    }
    fn new_x(magnitude: f64) -> Self {
        Self {
            state: [magnitude, 0., 0.],
        }
    }
    fn new_y(magnitude: f64) -> Self {
        Self {
            state: [0., magnitude, 0.],
        }
    }
    fn new_z(magnitude: f64) -> Self {
        Self {
            state: [0., 0., magnitude],
        }
    }
    fn new_random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self {
        let unit: [f64; 3] = UnitSphere.sample(rng);
        Self {
            state: [
                unit[0] * magnitude,
                unit[1] * magnitude,
                unit[2] * magnitude,
            ],
        }
    }
    fn magnitude(&self) -> f64 {
        (self.state[0] * self.state[0]
            + self.state[1] * self.state[1]
            + self.state[2] * self.state[2])
            .sqrt()
    }

    fn direction(&self) -> SpinVector {
        SpinVector::Heisenberg(
            self.state[0] / self.magnitude(),
            self.state[1] / self.magnitude(),
            self.state[2] / self.magnitude(),
        )
    }

    fn spinvector(&self) -> SpinVector {
        SpinVector::Heisenberg(self.state[0], self.state[1], self.state[2])
    }

    fn random<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self {
        Self::new_random(rng, magnitude)
    }

    fn propose_perturbation<R: rand::Rng>(&self, rng: &mut R, magnitude: f64) -> Self {
        self.random(rng, magnitude)
    }

    fn dot(&self, other: &Self) -> f64 {
        self.state[0] * other.state[0] + self.state[1] * other.state[1]
    }

    fn energy_diff(
        &self,
        calc_input: &crate::calculators::CalcInput<HeisenbergSpin>,
        ham: &crate::calculators::Hamiltonian,
        spins: &[Self],
        old_spin: &Self,
    ) -> f64 {
        self.energy(calc_input, ham, spins) - old_spin.energy(calc_input, ham, spins)
    }
    fn flip(&mut self, axis: &SpinVector) {
        let det_vec = (axis.clone() * 2. * (self.spinvector().dot(&axis))).to_vec();
        self.state = [
            self.state[0] - det_vec[0],
            self.state[1] - det_vec[1],
            self.state[2] - det_vec[2],
        ]
    }
}
