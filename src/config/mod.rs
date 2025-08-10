use serde::{Deserialize, Serialize};
use std::fmt;

mod anisotropy;
mod exchange;
mod grid;
mod output;
mod simulation;
mod structure;

pub use anisotropy::{Anisotropy, ParsedAnisotropy};
pub use exchange::{Exchange, ParsedExchange};
pub use grid::Grid;
pub use output::Output;
pub use simulation::Simulation;
pub use structure::Structure;

#[cfg(feature = "snapshots")]
mod snapshots;
#[cfg(feature = "snapshots")]
pub use snapshots::Snapshots;
#[cfg(feature = "snapshots")]
pub use snapshots::save_snapshots_to_hdf5;

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct Config {
    pub grid: Grid,
    pub simulation: Simulation,
    pub structure: Option<Structure>,
    pub output: Output,
    #[cfg(feature = "snapshots")]
    pub snapshots: Option<Snapshots>,

    pub exchange: Vec<Exchange>,
    #[serde(skip)]
    pub parsed_exchange: Vec<ParsedExchange>,

    pub anisotropy: Option<Anisotropy>,
    #[serde(skip)]
    pub parsed_anisotropy: Vec<ParsedAnisotropy>,
}

impl Config {
    pub fn new(content: &str) -> anyhow::Result<Self> {
        let mut config: Config = toml::from_str(content)?;
        for exchange in &config.exchange {
            let exchange_params =
                exchange.parse(&config.structure, config.grid.periodic_boundary)?;
            config.parsed_exchange.extend(exchange_params);
        }
        if let Some(anisotropy) = &config.anisotropy {
            config.parsed_anisotropy = anisotropy.parse()?;
        }
        config.validate()?;
        Ok(config)
    }

    fn validate(&mut self) -> anyhow::Result<()> {
        self.grid.validate()?;
        self.simulation.validate()?;
        for exchange in &self.exchange {
            exchange.validate(self.grid.sublattices)?;
        }
        if let Some(stru) = &self.structure {
            stru.validate(self.grid.sublattices)?;
        }
        self.output.validate()?;
        #[cfg(feature = "snapshots")]
        if let Some(snap) = &self.snapshots {
            snap.validate()?;
        }
        Ok(())
    }
}

impl fmt::Display for Config {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{}", self.grid)?;
        if let Some(stru) = &self.structure {
            writeln!(f, "{stru}")?;
        }
        writeln!(f, "{}", self.simulation)?;
        writeln!(f, "{}", self.output)?;

        if !&self.parsed_exchange.is_empty() {
            writeln!(f, "\nExchange Parameters:")?;
            writeln!(
                f,
                "  {:<4} | {:<3} | {:>3} {:>3} {:>3}  | {:>12}",
                "from", "to", "x", "y", "z", "strength (eV)"
            )?;
            for exchange in &self.parsed_exchange {
                writeln!(f, "{exchange}")?;
            }
        }

        if !&self.parsed_anisotropy.is_empty() {
            writeln!(f, "\nAnisotropy Parameters:")?;
            writeln!(
                f,
                "  {:<6} | {:>10}     | {:>12}",
                " ion ", "saxis", "strength (eV)"
            )?;
            for (i, anisotropy) in self.parsed_anisotropy.iter().enumerate() {
                write!(f, "  ion{i:<4}| {anisotropy}")?;
            }
        }
        #[cfg(feature = "snapshots")]
        if let Some(snapshots) = &self.snapshots {
            writeln!(f, "{snapshots}")?;
        }

        Ok(())
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
#[serde(rename_all = "snake_case")]
pub enum Algorithm {
    Metropolis,
    Wolff,
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
