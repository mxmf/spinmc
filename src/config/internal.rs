use super::raw_config::Model;
use crate::lattice::Atoms;

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
pub struct Config {
    // grid
    pub dim: [usize; 3],
    pub sublattices: usize,
    pub spin_magnitudes: Vec<f64>,
    pub pbc: [bool; 3],

    //rules
    pub exchange_params: Vec<ExchangeParams>,

    // simulation
    pub initial_state: InitialState,
    pub model: Model,
    pub n_equil: usize,
    pub n_steps: usize,
    pub temperatures: Vec<f64>,
    pub num_threads: usize,
    pub kb: f64,

    //output
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
}

impl Config {
    pub fn new(path: &str) -> anyhow::Result<Self> {
        let raw_config = RawConfig::load_from_file(path)?;
        let temperatures = Self::resolve_temperatures(&raw_config)?;
        let exchange_params = Self::resolve_exchange(&raw_config)?;
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

        Ok(Self {
            dim: raw_config.grid.dim,
            sublattices: raw_config.grid.sublattices,
            spin_magnitudes: raw_config.grid.spin_magnitudes,
            pbc: raw_config.grid.pbc,
            exchange_params,
            initial_state: raw_config.simulation.initial_state,
            model: raw_config.simulation.model,
            n_equil: raw_config.simulation.n_equil,
            n_steps: raw_config.simulation.n_steps,
            temperatures,
            num_threads: raw_config.simulation.num_threads,
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
        })
    }

    fn resolve_temperatures(raw_config: &RawConfig) -> anyhow::Result<Vec<f64>> {
        let sim = &raw_config.simulation;
        if sim.temperatures.is_some()
            && (sim.temp_start.is_some() || sim.temp_end.is_some() || sim.temp_step.is_some())
        {
            return Err(anyhow::anyhow!(
                "Cannot specify both `temperatures` and `temp_start`/`temp_end`/`temp_step` — please choose one."
            ));
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
            (None, None, None) => Err(anyhow::anyhow!(
                "You must provide either `temperatures` or all of `temp_start`, `temp_end`, `temp_step`."
            )),
            _ => Err(anyhow::anyhow!(
                "`temp_start`, `temp_end`, and `temp_step` must all be specified together."
            )),
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
                    (_, _, Some(_), Some(_), _) => {
                        panic!(
                            "Invalid configuration: do not set both `offsets` and `neighbor_order`; only one should be specified."
                        )
                    }
                    (_, _, None, None, _) => {
                        panic!(
                            "Missing configuration: you must specify either `offsets` or `neighbor_order`."
                        )
                    }
                    (_, _, _, Some(_), None) => {
                        panic!(
                            "Incomplete configuration: when using `neighbor_order`, `structure` must be set."
                        )
                    }
                    (from_sub, to_sub, Some(offsets), None, _) => {
                        if let (Some(from_sub), Some(to_sub)) = (from_sub, to_sub) {
                            exchange_params.push(ExchangeParams {
                                from_sub: *from_sub,
                                to_sub: *to_sub,
                                offsets: offsets.clone(),
                                strength: *strength,
                            });
                        } else {
                            panic!(
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
                            (None, Some(_)) => {
                                panic!(
                                    "Invalid configuration: `from_sub` must be specified when using `neighbor_order`."
                                );
                            }
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
}
