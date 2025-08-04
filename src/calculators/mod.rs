use crate::spin::SpinState;
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub struct CalcInput<S: SpinState> {
    pub magnitude: f64,
    pub exchange_neighbors: Option<Vec<(*const S, f64)>>,
    pub dm_neighbors: Option<Vec<(usize, [f64; 3], f64)>>,
    pub magnetic_field: Option<[f64; 3]>,
    pub easy_axis: Option<[f64; 3]>,
}

impl<S: SpinState> Default for CalcInput<S> {
    fn default() -> Self {
        CalcInput {
            magnitude: 0.0,
            exchange_neighbors: None,
            dm_neighbors: None,
            magnetic_field: None,
            easy_axis: None,
        }
    }
}

impl<S: SpinState> CalcInput<S> {
    pub fn validate_exchange_neighbor(&self) -> bool {
        let pairs = &self.exchange_neighbors;
        if let Some(vec) = pairs {
            let mut seen = HashSet::new();
            for (index, _) in vec {
                if !seen.insert(index) {
                    panic!(
                        "Duplicate neighbor indices found in your exchange coupling configuration. Please ensure all neighbor indices are unique."
                    );
                }
            }
        }
        true
    }
}

fn exchange_energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
    if let Some(list) = calc_input.exchange_neighbors.as_ref() {
        list.iter()
            .map(|(n, j)| {
                unsafe {
                    let neighbor = &*(*n); // *n 是 *const S，解引用为 &S
                    -j * spin.dot(neighbor)
                }
            })
            .sum()
    } else {
        0.0
    }
}

fn zeeman_energy<S: SpinState>(_: &S, _: &CalcInput<S>) -> f64 {
    unimplemented!();
}

fn anisotropy_energy<S: SpinState>(_: &S, _: &CalcInput<S>) -> f64 {
    unimplemented!();
}

fn dm_energy<S: SpinState>(_: &S, _: &CalcInput<S>, _: &[S]) -> f64 {
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

    pub fn compute<S: SpinState>(&self, spin: &S, calc_input: &CalcInput<S>, spins: &[S]) -> f64 {
        let mut result = 0.0;
        if self.config.exchange_enable {
            result += exchange_energy(spin, calc_input);
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
