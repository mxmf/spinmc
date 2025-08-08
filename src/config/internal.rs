use super::raw_config::Algorithm;
use super::raw_config::Model;
use crate::lattice::Atoms;
use std::fmt;

pub use super::raw_config::InitialState;
use super::raw_config::RawConfig;

#[derive(Debug, Clone)]
pub struct ExchangeParams {
    pub from_sub: usize,
    pub to_sub: usize,
    pub offsets: Vec<[isize; 3]>,
    pub strength: f64,
}
#[derive(Debug, Clone)]
pub struct AnisotropyParams {
    pub saxis: [f64; 3],
    pub strength: f64,
}

#[derive(Debug, Clone)]
pub struct Config {
    // grid
    pub dim: [usize; 3],
    pub sublattices: usize,
    pub spin_magnitudes: Vec<f64>,
    pub pbc: [bool; 3],

    //rules
    pub exchange_params: Vec<ExchangeParams>,
    pub anisotropy_params: Vec<AnisotropyParams>,

    // simulation
    pub initial_state: InitialState,
    pub model: Model,
    pub n_equil: usize,
    pub n_steps: usize,
    pub temperatures: Vec<f64>,
    pub num_threads: usize,
    pub kb: f64,
    pub algorithm: Algorithm,

    // output
    pub outfile: String,
    pub energy: bool,
    pub heat_capacity: bool,
    pub magnetization: bool,
    pub susceptibility: bool,
    pub magnetization_abs: bool,
    pub susceptibility_abs: bool,
    pub group_magnetization: bool,
    pub group_susceptibility: bool,
    pub group: Vec<Vec<usize>>,

    //snapshot
    #[cfg(feature = "snapshots")]
    pub snapshot_enable: bool,
    #[cfg(feature = "snapshots")]
    pub snapshot_params: crate::snapshots::ParsedSnapshots,
}

impl Config {
    pub fn new(path: &str) -> anyhow::Result<Self> {
        let raw_config = RawConfig::load_from_file(path)?;
        let temperatures = Self::resolve_temperatures(&raw_config)?;
        let exchange_params = Self::resolve_exchange(&raw_config)?;
        let anisotropy_params = Self::resolve_anisotropy(&raw_config)?;
        let kb = raw_config.simulation.kb.unwrap_or(0.00008617333262145178);

        let energy = raw_config.output.energy.unwrap_or_default();
        let heat_capacity = raw_config.output.heat_capacity.unwrap_or_default();
        let magnetization = raw_config.output.magnetization.unwrap_or_default();
        let susceptibility = raw_config.output.susceptibility.unwrap_or_default();
        let magnetization_abs = raw_config.output.magnetization_abs.unwrap_or_default();
        let susceptibility_abs = raw_config.output.susceptibility_abs.unwrap_or_default();
        let group_magnetization = raw_config.output.group_magnetization.unwrap_or_default();
        let group_susceptibility = raw_config.output.group_susceptibility.unwrap_or_default();
        let group = raw_config.output.group.unwrap_or_default();

        #[cfg(feature = "snapshots")]
        let (snapshots_enable, snapshots_params) = if let Some(snapshot) = raw_config.snapshots {
            (
                true,
                snapshot.parse(raw_config.simulation.n_equil, raw_config.simulation.n_steps),
            )
        } else {
            (false, crate::snapshots::ParsedSnapshots::default())
        };

        let algorithm = raw_config.simulation.algorithm.unwrap_or(Algorithm::Wolff);

        Ok(Self {
            dim: raw_config.grid.dim,
            sublattices: raw_config.grid.sublattices,
            spin_magnitudes: raw_config.grid.spin_magnitudes,
            pbc: raw_config.grid.pbc,
            exchange_params,
            anisotropy_params,
            initial_state: raw_config.simulation.initial_state,
            model: raw_config.simulation.model,
            n_equil: raw_config.simulation.n_equil,
            n_steps: raw_config.simulation.n_steps,
            temperatures,
            num_threads: raw_config.simulation.num_threads,
            algorithm,
            kb,
            outfile: raw_config.output.outfile,
            energy,
            heat_capacity,
            magnetization,
            susceptibility,
            magnetization_abs,
            susceptibility_abs,
            group_magnetization,
            group_susceptibility,
            group,
            #[cfg(feature = "snapshots")]
            snapshot_enable: snapshots_enable,
            #[cfg(feature = "snapshots")]
            snapshot_params: snapshots_params,
        })
    }

    fn resolve_temperatures(raw_config: &RawConfig) -> anyhow::Result<Vec<f64>> {
        let sim = &raw_config.simulation;
        if sim.temperatures.is_some()
            && (sim.temp_start.is_some() || sim.temp_end.is_some() || sim.temp_step.is_some())
        {
            anyhow::bail!(
                "Cannot specify both `temperatures` and `temp_start`/`temp_end`/`temp_step` — please choose one."
            );
        }

        if let Some(temperatures) = &sim.temperatures {
            return Ok(temperatures.clone());
        }

        match (sim.temp_start, sim.temp_end, sim.temp_step) {
            (Some(start), Some(end), Some(step)) => {
                let mut vec = vec![];
                let mut t = start;
                while t <= end + 1e-8 {
                    vec.push(t);
                    t += step;
                }
                Ok(vec)
            }
            (None, None, None) => anyhow::bail!(
                "You must provide either `temperatures` or all of `temp_start`, `temp_end`, `temp_step`."
            ),
            _ => anyhow::bail!(
                "`temp_start`, `temp_end`, and `temp_step` must all be specified together."
            ),
        }
    }

    fn resolve_exchange(raw_config: &RawConfig) -> anyhow::Result<Vec<ExchangeParams>> {
        let mut exchange_params = vec![];
        if let Some(exchange_raws) = &raw_config.exchange {
            for exchange_raw in exchange_raws {
                let (from_sub, to_sub, offsets, neighbor_order, strength, structure) = (
                    &exchange_raw.from_sub,
                    &exchange_raw.to_sub,
                    &exchange_raw.offsets,
                    &exchange_raw.neighbor_order,
                    &exchange_raw.strength,
                    &raw_config.structure,
                );

                match (from_sub, to_sub, offsets, neighbor_order, structure) {
                    (_, _, Some(_), Some(_), _) => anyhow::bail!(
                        "Invalid configuration: do not set both `offsets` and `neighbor_order`; only one should be specified.",
                    ),

                    (_, _, None, None, _) => anyhow::bail!(
                        "Missing configuration: you must specify either `offsets` or `neighbor_order`.",
                    ),

                    (_, _, _, Some(_), None) => anyhow::bail!(
                        "Incomplete configuration: when using `neighbor_order`, `structure` must be set.",
                    ),

                    (from_sub, to_sub, Some(offsets), None, _) => {
                        if let (Some(from_sub), Some(to_sub)) = (from_sub, to_sub) {
                            exchange_params.push(ExchangeParams {
                                from_sub: *from_sub,
                                to_sub: *to_sub,
                                offsets: offsets.clone(),
                                strength: *strength,
                            });
                        } else {
                            anyhow::bail!(
                                "Incomplete configuration: when using `offsets`, both `from_sub` and `to_sub` must be specified."
                            )
                        }
                    }

                    (from_sub, to_sub, None, Some(neighbor_order), Some(structure)) => {
                        let atoms = Atoms {
                            cell: structure.cell,
                            positions: structure.positions.clone(),
                            pbc: raw_config.grid.pbc,
                            tolerance: structure.tolerance.unwrap_or(0.0001),
                        };

                        let neighbors = match (from_sub, to_sub) {
                            (Some(from), Some(to)) => {
                                atoms.find_neighbors_from_to(*from, *to, *neighbor_order)
                            }
                            (Some(from), None) => atoms.find_neighbors_from(*from, *neighbor_order),
                            (None, None) => atoms.find_neighbors_all(*neighbor_order),
                            (None, Some(_)) => anyhow::bail!(
                                "Invalid configuration: `from_sub` must be specified when using `neighbor_order`."
                            ),
                        };

                        for neighbor in neighbors {
                            exchange_params.push(ExchangeParams {
                                from_sub: neighbor.from,
                                to_sub: neighbor.to,
                                offsets: vec![neighbor.offset],
                                strength: *strength,
                            });
                        }
                    }
                }
            }
        }
        Ok(exchange_params)
    }

    fn resolve_anisotropy(raw_config: &RawConfig) -> anyhow::Result<Vec<AnisotropyParams>> {
        let mut result = vec![];

        if let Some(anisotropy) = &raw_config.anisotropy {
            for (saxis, strength) in anisotropy.saxis.iter().zip(anisotropy.strength.iter()) {
                let saxis_norm =
                    (saxis[0] * saxis[0] + saxis[1] * saxis[1] + saxis[2] * saxis[2]).sqrt();
                if saxis_norm == 0.0 {
                    anyhow::bail!("Anisotropy direction vector {:?} has zero length", saxis);
                }

                let ani = AnisotropyParams {
                    saxis: [
                        saxis[0] / saxis_norm,
                        saxis[1] / saxis_norm,
                        saxis[2] / saxis_norm,
                    ],
                    strength: *strength,
                };
                result.push(ani);
            }
        }
        Ok(result)
    }
}

impl fmt::Display for Config {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut hamiltonian = String::new();
        if !self.exchange_params.is_empty() {
            hamiltonian += "-∑⟨i,j⟩ Jᵢⱼ Sᵢ · Sⱼ";
        }
        // TODO more hamiltonian term

        writeln!(f, "\n=== Simulation Configuration ===")?;
        writeln!(f, "Hamiltonian: H = {hamiltonian}")?;
        writeln!(f, "\nGrid:")?;
        writeln!(f, "  Dimension: {:?}, PBC: {:?}", self.dim, self.pbc)?;
        writeln!(f, "  Sublattices: {}", self.sublattices)?;
        writeln!(f, "  Spin Magnitudes: {:?}", self.spin_magnitudes)?;

        writeln!(f, "\nExchange Parameters:")?;
        writeln!(
            f,
            "{:<4} | {:<3} | {:>3} {:>3} {:>3}  | {:>12}",
            "from", "to", "x", "y", "z", "strength (eV)"
        )?;
        for params in &self.exchange_params {
            let (from_sub, to_sub, strength) = (params.from_sub, params.to_sub, params.strength);

            for offset in &params.offsets {
                writeln!(
                    f,
                    "{from_sub:<4} | {to_sub:<3} | {:>3} {:>3} {:>3}  | {strength:>8.12}",
                    offset[0], offset[1], offset[2]
                )?;
            }
        }

        if !&self.anisotropy_params.is_empty() {
            writeln!(f, "\nAnisotropy Parameters:")?;
            writeln!(
                f,
                "{:<6} | {:>10}     | {:>12}",
                "ion", "saxis", "strength (eV)"
            )?;
            for (i, ani) in self.anisotropy_params.iter().enumerate() {
                let [x, y, z] = ani.saxis;
                writeln!(
                    f,
                    "ion{i:<4}| {x:>4} {y:>4} {z:>4} | {:>8.12}",
                    ani.strength
                )?;
            }
        }
        // for ions in

        writeln!(f, "\nSimulation Parameters:")?;
        writeln!(f, "  Initial State: {:?}", self.initial_state)?;
        writeln!(f, "  Model: {:?}", self.model)?;
        writeln!(f, "  Equilibration Steps: {}", self.n_equil)?;
        writeln!(f, "  Simulation Steps: {}", self.n_steps)?;
        writeln!(f, "  Boltzmann Constant (kB): {} (eV/K)", self.kb)?;
        writeln!(f, "  Algorithm: {:?}", self.algorithm)?;
        writeln!(f, "  Threads: {}", self.num_threads)?;

        write!(f, "Temperatures (K):\n  ")?;
        for t in &self.temperatures {
            write!(f, "{t:.4}   ")?;
        }
        writeln!(f)?;

        writeln!(f, "\nOutput:")?;
        writeln!(f, "  Output File: {}", self.outfile)?;

        writeln!(
            f,
            "  Energy [E = <H> / N = < {} >/N]: {}",
            hamiltonian, self.energy
        )?;
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

        #[cfg(feature = "snapshots")]
        if self.snapshot_enable {
            writeln!(f, "{}", self.snapshot_params)?;
        } else {
            writeln!(f, "\nSnapshots: Disble")?;
        }

        Ok(())
    }
}
