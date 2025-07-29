use crate::spin::SpinState;

#[derive(Clone, Debug, Default)]
pub struct CalcInput {
    pub exchange_neighbors: Option<Vec<(usize, f64)>>,
    pub dm_neighbors: Option<Vec<(usize, [f64; 3], f64)>>,
    pub magnetic_field: Option<[f64; 3]>,
    pub easy_axis: Option<[f64; 3]>,
}

fn exchange_energy<S: SpinState>(spin: &S, calc_input: &CalcInput, spins: &[S]) -> f64 {
    if let Some(list) = calc_input.exchange_neighbors.as_ref() {
        list.iter().map(|(n, j)| -j * spin.dot(&spins[*n])).sum()
    } else {
        0.0
    }
}

fn zeeman_energy<S: SpinState>(_: &S, _: &CalcInput) -> f64 {
    unimplemented!();
}

fn anisotropy_energy<S: SpinState>(_: &S, _: &CalcInput) -> f64 {
    unimplemented!();
}

fn dm_energy<S: SpinState>(_: &S, _: &CalcInput, _: &[S]) -> f64 {
    unimplemented!();
}

#[derive(Clone, Copy, Debug)]
pub struct HamiltonianConfig {
    pub exchange_enable: bool,
    pub anisotropy_enable: bool,
    pub zeeman_enable: bool,
    pub dm_enable: bool,
}

#[derive(Clone, Debug)]
pub struct Hamiltonian {
    config: HamiltonianConfig,
}
impl Hamiltonian {
    pub fn new(config: HamiltonianConfig) -> Self {
        Self { config }
    }

    pub fn compute<S: SpinState>(&self, spin: &S, calc_input: &CalcInput, spins: &[S]) -> f64 {
        let mut result = 0.0;
        if self.config.exchange_enable {
            result += exchange_energy(spin, calc_input, spins);
        }

        if self.config.zeeman_enable {
            result += zeeman_energy(spin, calc_input);
        }

        if self.config.anisotropy_enable {
            result += anisotropy_energy(spin, calc_input);
        }
        if self.config.dm_enable {
            result += dm_energy(spin, calc_input, spins)
        }
        result
    }
}
