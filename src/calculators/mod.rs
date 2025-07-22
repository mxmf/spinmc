use crate::spin::SpinState;

#[derive(Clone, Debug, Default)]
pub struct CalcInput {
    pub exchange_neighbors: Option<Vec<(usize, f64)>>,
    pub dm_neighbors: Option<Vec<(usize, [f64; 3], f64)>>,
    pub magnetic_field: Option<[f64; 3]>,
    pub easy_axis: Option<[f64; 3]>,
}

#[derive(Clone, Debug)]
pub struct Exchange;
impl Exchange {
    pub fn energy<S: SpinState>(spin: &S, calc_input: &CalcInput, spins: &[S]) -> f64 {
        if let Some(list) = calc_input.exchange_neighbors.as_ref() {
            list.iter().map(|(n, j)| -j * spin.dot(&spins[*n])).sum()
        } else {
            0.0
        }
    }
}

#[derive(Clone, Debug)]
pub struct Anisotropy;
impl Anisotropy {
    fn energy<S: SpinState>(_: &S, _: &CalcInput) -> f64 {
        unimplemented!();
    }
}

#[derive(Clone, Debug)]
pub struct Zeeman;
impl Zeeman {
    fn energy<S: SpinState>(_: &S, _: &CalcInput) -> f64 {
        unimplemented!();
    }
}

#[derive(Clone, Debug)]
pub struct DM;
impl DM {
    fn energy<S: SpinState>(_: &S, _: &CalcInput) -> f64 {
        unimplemented!();
    }
}

#[derive(Clone, Copy, Debug)]
pub struct HamiltonianConfig {
    pub exchange_enable: bool,
    pub anisotropy_enable: bool,
    pub zeeman_enable: bool,
    pub dm_enable: bool,
}

#[derive(Clone, Debug)]
pub enum EnergyModelType {
    Exchange(Exchange),
    Anisotropy(Anisotropy),
    Zeeman(Zeeman),
    DM(DM),
}

impl EnergyModelType {
    pub fn energy<S: SpinState>(&self, spin: &S, input: &CalcInput, spins: &[S]) -> f64 {
        match self {
            Self::Exchange(_) => Exchange::energy(spin, input, spins),
            Self::Anisotropy(_) => Anisotropy::energy(spin, input),
            Self::Zeeman(_) => Zeeman::energy(spin, input),
            Self::DM(_) => DM::energy(spin, input),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Hamiltonian {
    models: Vec<EnergyModelType>,
}
impl Hamiltonian {
    pub fn new(config: HamiltonianConfig) -> Self {
        let mut models = Vec::new();
        if config.exchange_enable {
            models.push(EnergyModelType::Exchange(Exchange));
        }
        if config.anisotropy_enable {
            models.push(EnergyModelType::Anisotropy(Anisotropy));
        }
        if config.zeeman_enable {
            models.push(EnergyModelType::Zeeman(Zeeman));
        }
        if config.dm_enable {
            models.push(EnergyModelType::DM(DM));
        }
        Self { models }
    }

    pub fn compute<S: SpinState>(&self, spin: &S, calc_input: &CalcInput, spins: &[S]) -> f64 {
        self.models
            .iter()
            .map(|model| model.energy(spin, calc_input, spins))
            .sum()
    }
}
