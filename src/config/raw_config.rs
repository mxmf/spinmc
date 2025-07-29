use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
pub struct Grid {
    pub dim: [usize; 3],
    pub sublattices: usize,
    pub spin_magnitudes: Vec<f64>,
    pub pbc: [bool; 3],
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Exchange {
    pub from_sub: usize,
    pub to_sub: usize,

    pub offsets: Option<Vec<[isize; 3]>>,
    pub neighbor_order: Option<usize>,

    pub strength: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct DM {
    pub from_sub: usize,
    pub to_sub: usize,

    pub offsets: Option<Vec<[isize; 3]>>,
    pub neighbor_order: Option<usize>,

    pub strength: [f64; 3],
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Zeeman {
    pub saxis: [f64; 3],
    pub strength: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Anisotropy {
    pub saxis: Vec<[f64; 3]>,
    pub strength: f64,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "snake_case")]
pub enum InitialState {
    Random,
    X,
    Y,
    Z,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "snake_case")]
pub enum Model {
    Ising,
    Xy,
    Heisenberg,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Simulation {
    pub initial_state: InitialState,
    pub model: Model,
    pub n_equil: usize,
    pub n_steps: usize,
    pub temperatures: Option<Vec<f64>>,
    pub temp_start: Option<f64>,
    pub temp_end: Option<f64>,
    pub temp_step: Option<f64>,
    pub num_threads: usize,

    // default kb = 8.617333262145×10^−5 eV/K
    pub kb: Option<f64>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct Output {
    pub outfile: String,
    pub energy: Option<bool>,
    pub heat_capacity: Option<bool>,
    pub magnetization: Option<bool>,
    pub susceptibility: Option<bool>,
    pub magnetization_abs: Option<bool>,
    pub susceptibility_abs: Option<bool>,
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct RawConfig {
    pub grid: Grid,
    pub simulation: Simulation,
    pub output: Output,
    pub exchange: Option<Vec<Exchange>>,
    pub dm: Option<Vec<DM>>,
    pub zeeman: Option<Zeeman>,
    pub anisotropy: Option<Anisotropy>,
}

impl RawConfig {
    pub fn load_from_file(path: &str) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let mut config: RawConfig = toml::from_str(&content)?;
        config.validate()?;
        Ok(config)
    }
}
impl RawConfig {
    pub fn validate(&mut self) -> anyhow::Result<()> {
        self.validate_grid()?;
        self.validate_output()?;
        Ok(())
    }

    fn validate_grid(&mut self) -> anyhow::Result<()> {
        if self.grid.sublattices != self.grid.spin_magnitudes.len() {
            anyhow::bail!(
                "spin_magnitude length ({}) does not match sublattices ({})",
                self.grid.spin_magnitudes.len(),
                self.grid.sublattices
            );
        }

        Ok(())
    }
    fn validate_output(&mut self) -> anyhow::Result<()> {
        match self.output {
            Output {
                energy: None,
                heat_capacity: None,
                magnetization: None,
                susceptibility: None,
                magnetization_abs: None,
                susceptibility_abs: None,
                ..
            } => {
                panic!("No output fields specified: Please enable at least one observable.")
            }

            Output {
                energy: Some(false),
                heat_capacity: Some(false),
                magnetization: Some(false),
                susceptibility: Some(false),
                magnetization_abs: Some(false),
                susceptibility_abs: Some(false),
                ..
            } => {
                panic!("No output fields specified: Please enable at least one observable.")
            }

            _ => {}
        }
        Ok(())
    }
}
