use crate::config::Config;
use crate::lattice::Grid;
use crate::spin::{SpinState, SpinVector};
use std::fmt;

#[derive(Clone, Debug)]
pub struct StatsConfig {
    pub energy: bool,
    pub heat_capacity: bool,
    pub magnetization: bool,
    pub susceptibility: bool,
    pub magnetization_abs: bool,
    pub susceptibility_abs: bool,
}

impl fmt::Display for StatsConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "# T (K)")?;
        if self.energy {
            write!(f, "\tEnergy (eV)")?;
        }
        if self.heat_capacity {
            write!(f, "\t$C$ (eV/K)")?;
        }
        if self.magnetization {
            write!(f, "\tM ($\\mu_B$)")?;
        }
        if self.susceptibility {
            write!(f, "\t$\\chi$ ")?;
        }
        if self.magnetization_abs {
            write!(f, "\t|M| ($\\mu_B$)")?;
        }
        if self.susceptibility_abs {
            write!(f, "\t$|\\chi|$ ")?;
        }
        Ok(())
    }
}
#[derive(Debug, Default)]
pub struct StatResult {
    pub t: f64,
    pub energy: Option<f64>,
    pub specific_heat: Option<f64>,
    pub magnetization: Option<f64>,      // |<M>| / N
    pub susceptibility: Option<f64>,     // ( < M^2 > - <M>^2)/(N * k_B * T)
    pub magnetization_abs: Option<f64>,  // < |M| >/ N
    pub susceptibility_abs: Option<f64>, // ( < |M|^2 > - <M>^2)/(N * k_B * T)
}

impl fmt::Display for StatResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.4}", self.t)?;
        if let Some(e) = self.energy {
            write!(f, "\t{e:.8}")?;
        }
        if let Some(c) = self.specific_heat {
            write!(f, "\t{c:.8}")?;
        }

        if let Some(m) = self.magnetization {
            write!(f, "\t{m:.8}")?;
        }
        if let Some(chi) = self.susceptibility {
            write!(f, "\t{chi:.8}")?;
        }
        if let Some(m_abs) = self.magnetization_abs {
            write!(f, "\t{m_abs:.8}")?;
        }
        if let Some(ch_absi) = self.susceptibility_abs {
            write!(f, "\t{ch_absi:.8}")?;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct Stats {
    pub energy_sum: f64,
    pub energy2_sum: f64,
    pub m_sum: SpinVector, // ∑ M
    pub m_2_sum: f64,      // ∑ M^2
    pub m_abs_sum: f64,    // ∑ |M|
    pub steps: usize,
    pub size: f64,
    pub kb: f64,
    pub t: f64,
    pub stats_config: StatsConfig,
}

impl Stats {
    pub fn new<S: SpinState>(config: &Config, t: f64) -> Self {
        let stats_config = StatsConfig {
            energy: config.energy,
            heat_capacity: config.heat_capacity,
            magnetization: config.magnetization,
            susceptibility: config.susceptibility,
            magnetization_abs: config.magnetization_abs,
            susceptibility_abs: config.susceptibility_abs,
        };
        Self {
            energy_sum: 0.,
            energy2_sum: 0.,
            m_sum: S::zero(),
            m_2_sum: 0.,
            m_abs_sum: 0.,
            steps: 0,
            kb: config.kb,
            t,

            size: (config.dim[0] * config.dim[1] * config.dim[2] * config.sublattices) as f64,
            stats_config,
        }
    }

    pub fn record<S: SpinState, R: rand::Rng>(&mut self, grid: &Grid<S, R>) {
        if self.stats_config.energy {
            let energy = grid.total_energy();
            self.energy_sum += energy;

            if self.stats_config.heat_capacity {
                self.energy2_sum += energy * energy;
            }
        }

        if self.stats_config.magnetization
            || self.stats_config.susceptibility
            || self.stats_config.magnetization_abs
            || self.stats_config.susceptibility_abs
        {
            let spin_vec = grid.total_spin_vector();

            if self.stats_config.magnetization || self.stats_config.susceptibility {
                self.m_sum += &spin_vec;
            }

            if self.stats_config.magnetization_abs || self.stats_config.susceptibility_abs {
                self.m_abs_sum += spin_vec.norm();
            }
            if self.stats_config.susceptibility || self.stats_config.susceptibility_abs {
                self.m_2_sum += spin_vec.norm_sqr();
            }
        }

        self.steps += 1;
    }

    pub fn result(&self) -> StatResult {
        let energy = if self.stats_config.energy {
            Some(self.energy_sum / self.steps as f64 / self.size)
        } else {
            None
        };

        let specific_heat = if self.stats_config.heat_capacity {
            let e_avg = self.energy_sum / self.steps as f64;
            let e2_avg = self.energy2_sum / self.steps as f64;
            Some((e2_avg - e_avg * e_avg) / (self.kb * self.t * self.t) / self.size)
        } else {
            None
        };

        let magnetization = if self.stats_config.magnetization {
            Some((&self.m_sum / self.steps as f64).norm() / self.size)
        } else {
            None
        };

        let susceptibility = if self.stats_config.susceptibility {
            let m_avg = &self.m_sum / self.steps as f64; //<M>
            // let m_avg = self.m_norm_sum / self.steps as f64; // < |M| >
            let m2_avg = self.m_2_sum / self.steps as f64; // < |M|^2>
            Some((m2_avg - m_avg.norm_sqr()) / (self.kb * self.t) / self.size)
        } else {
            None
        };

        let magnetization_abs = if self.stats_config.magnetization_abs {
            Some(self.m_abs_sum / self.steps as f64 / self.size)
        } else {
            None
        };

        let susceptibility_abs = if self.stats_config.susceptibility_abs {
            let m_abs_avg = self.m_abs_sum / self.steps as f64; //<M>
            let m2_avg = self.m_2_sum / self.steps as f64; // < |M|^2>
            Some((m2_avg - m_abs_avg * m_abs_avg) / (self.kb * self.t) / self.size)
        } else {
            None
        };

        StatResult {
            t: self.t,
            energy,
            specific_heat,
            magnetization,
            susceptibility,
            magnetization_abs,
            susceptibility_abs,
        }
    }
}
