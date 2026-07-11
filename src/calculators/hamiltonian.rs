use crate::{
    calculators::{input::CalcInput, onsite, pairs},
    config::Config,
    spin::SpinState,
};

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
        let exchange_enable = !config.parsed_exchange.is_empty();
        let anisotropy_enable = !config.parsed_anisotropy.is_empty();

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
            result += pairs::energy(spin, calc_input);
        }
        result += onsite::energy(
            spin,
            calc_input,
            self.config.anisotropy_enable,
            self.config.zeeman_enable,
        );
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
            result += pairs::local_energy(spin, calc_input);
        }
        result += onsite::energy(
            spin,
            calc_input,
            self.config.anisotropy_enable,
            self.config.zeeman_enable,
        );
        if self.config.dm_enable {
            result += dm_energy(spin, calc_input, spins)
        }
        result
    }

    pub fn compute_anisotropy<S: SpinState>(&self, spin: &S, calc_input: &CalcInput<S>) -> f64 {
        onsite::anisotropy_energy(spin, calc_input)
    }
}

fn dm_energy<S: SpinState>(_: &S, _: &CalcInput<S>, _: &[S]) -> f64 {
    unimplemented!();
}
