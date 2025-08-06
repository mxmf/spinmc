use crate::spin::{SpinState, SpinVector};

#[derive(Default, Debug, Clone)]
#[repr(C)]
#[cfg_attr(feature = "snapshots", derive(hdf5_metno::H5Type))]
pub struct IsingSpin {
    state: f64,
}

impl SpinState for IsingSpin {
    fn zero() -> SpinVector {
        SpinVector::Ising(0.)
    }
    fn new_x(magnitude: f64) -> Self {
        Self { state: magnitude }
    }

    fn new_y(magnitude: f64) -> Self {
        Self { state: magnitude }
    }
    fn new_z(magnitude: f64) -> Self {
        Self { state: magnitude }
    }
    fn new_random<R: rand::Rng>(rng: &mut R, magnitude: f64) -> Self {
        let sign = if rng.random_bool(0.5) { 1.0 } else { -1.0 };
        let value = sign * magnitude;
        Self { state: value }
    }

    fn magnitude(&self) -> f64 {
        self.state.abs()
    }

    fn direction(&self) -> SpinVector {
        SpinVector::Ising(self.state.signum())
    }
    fn spinvector(&self) -> SpinVector {
        SpinVector::Ising(self.state)
    }

    fn random<R: rand::Rng>(&self, rng: &mut R, _magnitude: f64) -> Self {
        let sign = if rng.random_bool(0.5) { 1.0 } else { -1.0 };
        Self {
            state: self.magnitude() * sign,
        }
    }

    fn propose_perturbation<R: rand::Rng>(&self, _rng: &mut R, _magnitude: f64) -> Self {
        Self { state: -self.state }
    }

    fn dot(&self, other: &Self) -> f64 {
        self.state * other.state
    }

    fn energy_diff(
        &self,
        calc_input: &crate::calculators::CalcInput<IsingSpin>,
        ham: &crate::calculators::Hamiltonian,
        spins: &[Self],
        _old_spin: &Self,
    ) -> f64 {
        2. * self.energy(calc_input, ham, spins)
    }

    fn flip(&mut self, _axis: &SpinVector) {
        self.state = -self.state
    }
}
