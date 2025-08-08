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
    pub from_sub: Option<usize>,
    pub to_sub: Option<usize>,

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

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Anisotropy {
    pub saxis: Vec<[f64; 3]>,
    pub strength: Vec<f64>,
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

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
    Metropolis,
    Wolff,
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
    pub algorithm: Option<Algorithm>,
    // default kb = 8.617333262145×10^−5 eV/K
    pub kb: Option<f64>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Structure {
    pub cell: [[f64; 3]; 3],
    pub positions: Vec<[f64; 3]>,
    pub tolerance: Option<f64>,
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
    pub group_magnetization: Option<bool>,
    pub group_susceptibility: Option<bool>,
    pub group: Option<Vec<Vec<usize>>>,
    pub stats_interval: Option<usize>,
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
    pub structure: Option<Structure>,

    #[cfg(feature = "snapshots")]
    pub snapshots: Option<crate::snapshots::Snapshots>,
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
        self.validate_anistropy()?;
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
                anyhow::bail!("No output fields specified: Please enable at least one observable.");
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
                anyhow::bail!("No output fields specified: Please enable at least one observable.");
            }

            _ => {}
        }

        if (self.output.group_magnetization == Some(true)
            || self.output.group_susceptibility == Some(true))
            && self.output.group.is_none()
        {
            anyhow::bail!("group_magnetization/group_susceptibility set but group is None");
        }

        if let Some(group) = &self.output.group {
            for i in group {
                for j in i {
                    if *j >= self.grid.sublattices {
                        anyhow::bail!("group_magnetization member{j} not in sublattices, check it");
                    }
                }
            }
        };

        Ok(())
    }

    fn validate_anistropy(&self) -> anyhow::Result<()> {
        if let Some(anisotropy) = &self.anisotropy {
            if anisotropy.saxis.len() != anisotropy.strength.len() {
                anyhow::bail!(
                    "anisotropy saxis length ({}) does not match anisotropy strength length ({})",
                    anisotropy.saxis.len(),
                    anisotropy.strength.len(),
                );
            }
            if anisotropy.strength.len() != self.grid.sublattices {
                anyhow::bail!(
                    "anisotropy strength length ({}) does not match sublattices ({})",
                    anisotropy.strength.len(),
                    self.grid.sublattices
                );
            }
        }
        Ok(())
    }
}
