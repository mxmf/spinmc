use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub struct Output {
    #[serde(default = "default_outfile")]
    pub savefile: String,
    #[serde(default = "default_false")]
    pub energy: bool,
    #[serde(default = "default_false")]
    pub heat_capacity: bool,
    #[serde(default = "default_false")]
    pub magnetization: bool,
    #[serde(default = "default_false")]
    pub susceptibility: bool,
    #[serde(default = "default_false")]
    pub magnetization_abs: bool,
    #[serde(default = "default_false")]
    pub susceptibility_abs: bool,
    #[serde(default = "default_false")]
    pub group_magnetization: bool,
    #[serde(default = "default_false")]
    pub group_susceptibility: bool,
    #[serde(default)]
    pub group: Vec<Vec<usize>>,
    #[serde(default = "default_stats_interval")]
    pub stats_interval: usize,
}

fn default_false() -> bool {
    false
}

fn default_outfile() -> String {
    "result.txt".to_string()
}

fn default_stats_interval() -> usize {
    1
}

impl Output {
    pub fn validate(&self) -> anyhow::Result<()> {
        if let (false, false, false, false, false, false, false, false) = (
            self.energy,
            self.heat_capacity,
            self.magnetization,
            self.susceptibility,
            self.magnetization_abs,
            self.susceptibility_abs,
            self.group_magnetization,
            self.group_susceptibility,
        ) {
            anyhow::bail!("No output fields specified: Please enable at least one observable.")
        }
        Ok(())
    }
}

impl fmt::Display for Output {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\nOutput:")?;
        writeln!(f, "  Output File: {}", self.savefile)?;
        writeln!(f, "  stats_interval: {}", self.stats_interval)?;

        writeln!(f, "  Energy [E = <H> / N ]: {}", self.energy)?;
        writeln!(
            f,
            "  Heat Capacity [ C = (⟨E²⟩ - ⟨E⟩²) / (N kB T²) ]: {}",
            self.heat_capacity
        )?;
        writeln!(
            f,
            "  Magnetization [ M = ⟨Σ s_i⟩ / N ] : {}",
            self.magnetization
        )?;
        writeln!(
            f,
            "  Susceptibility [ χ = (⟨M²⟩ - ⟨M⟩²) / (N kB T) ]: {}",
            self.susceptibility
        )?;
        writeln!(
            f,
            "  Magnetization_abs [M = ⟨|Σ s_i|⟩ / N]: {}",
            self.magnetization_abs
        )?;
        writeln!(
            f,
            "  susceptibility_abs [  χ(|M|) = (⟨|M|²⟩ - ⟨|M|⟩²) / (N kB T) ]: {}",
            self.susceptibility_abs
        )?;
        writeln!(f, "  Group Magnetization: {}", self.group_magnetization)?;
        writeln!(f, "  Group Susceptibility: {}", self.group_susceptibility)?;
        if self.group_magnetization || self.group_susceptibility {
            writeln!(f, "  Groups:")?;
            for (i, group) in self.group.iter().enumerate() {
                writeln!(f, "    Group {i}: {group:?}")?;
            }
        }

        Ok(())
    }
}
