use crate::{config::Config, spin::SpinState};
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub struct CalcInput<S: SpinState> {
    pub magnitude: f64,
    pub exchange_neighbors: Option<Vec<(*const S, f64)>>,
    pub exchanges: Vec<f64>,
    pub exchange_neighbor_index: Vec<usize>,
    pub dm_neighbors: Option<Vec<(usize, [f64; 3], f64)>>,
    pub magnetic_field: Option<[f64; 3]>,
    pub easy_axis: Option<[f64; 3]>,
    pub anisotropy: (f64, [f64; 3]),
}

impl<S: SpinState> Default for CalcInput<S> {
    fn default() -> Self {
        CalcInput {
            magnitude: 0.0,
            exchange_neighbor_index: vec![],
            exchanges: vec![],
            exchange_neighbors: None,
            dm_neighbors: None,
            magnetic_field: None,
            easy_axis: None,
            anisotropy: (0., [0., 0., 1.]),
        }
    }
}

impl<S: SpinState> CalcInput<S> {
    pub fn validate_exchange_neighbor(&self) {
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

            if vec.len() != self.exchange_neighbor_index.len() || vec.len() != self.exchanges.len()
            {
                panic!()
            }
        }
    }
}

/// Compute total exchange energy, with 1/2 factor to avoid double counting.
/// Should only be used in total energy evaluation.
fn exchange_energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
    if let Some(list) = calc_input.exchange_neighbors.as_ref() {
        let energy: f64 = list
            .iter()
            .map(|(n, j)| {
                unsafe {
                    let neighbor = &*(*n); // *n 是 *const S，解引用为 &S
                    -j * spin.dot(neighbor)
                }
            })
            .sum();
        energy / 2.
        // Each exchange interaction between sites i and j is counted twice
        // (once from i → j and once from j → i), so we divide by 2 to avoid double-counting.
    } else {
        0.0
    }
}

/// Compute local exchange energy for one site. No 1/2 factor.
/// Used in energy_diff or local site updates.
fn local_exchange_energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
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

fn anisotropy_energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
    let (strength, axis) = calc_input.anisotropy;

    let spin_array = spin.to_array();

    let dot = spin_array[0] * axis[0] + spin_array[1] * axis[1] + spin_array[2] * axis[2];

    strength * dot * dot
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
    pub config: HamiltonianConfig,
}
impl Hamiltonian {
    pub fn new(config: &Config) -> Self {
        let exchange_enable = !config.exchange_params.is_empty();
        let anisotropy_enable = !config.anisotropy_params.is_empty();

        let ham_config = HamiltonianConfig {
            exchange_enable,
            anisotropy_enable,
            zeeman_enable: false,
            dm_enable: false,
        };
        Self { config: ham_config }
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
    pub fn local_compute<S: SpinState>(
        &self,
        spin: &S,
        calc_input: &CalcInput<S>,
        spins: &[S],
    ) -> f64 {
        let mut result = 0.0;
        if self.config.exchange_enable {
            result += local_exchange_energy(spin, calc_input);
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

    pub fn compute_anisotropy<S: SpinState>(&self, spin: &S, calc_input: &CalcInput<S>) -> f64 {
        anisotropy_energy(spin, calc_input)
    }
}
