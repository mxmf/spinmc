use super::raw_config::Model;

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
}

impl Config {
    pub fn new(path: &str) -> anyhow::Result<Self> {
        let raw_config = RawConfig::load_from_file(path)?;
        let temperatures = Self::resolve_temperatures(&raw_config)?;
        let exchange_params = Self::resolve_exchange(&raw_config)?;
        let kb = raw_config.simulation.kb.unwrap_or(0.00008617333262145178);

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
                let (from_sub, to_sub, offsets, neighbor_order, strength) = (
                    &exchange_raw.from_sub,
                    &exchange_raw.to_sub,
                    &exchange_raw.offsets,
                    &exchange_raw.neighbor_order,
                    &exchange_raw.strength,
                );

                match (offsets, neighbor_order) {
                    (Some(offsets), None) => {
                        exchange_params.push(ExchangeParams {
                            from_sub: *from_sub,
                            to_sub: *to_sub,
                            offsets: offsets.clone(),
                            strength: *strength,
                        });
                    }
                    (None, Some(_neighbor_order)) => {
                        unimplemented!("neighbor to offsets")
                    }

                    (Some(_), Some(_)) => {
                        return Err(anyhow::anyhow!(
                            "Cannot specify both `offsets` and `neighbor_order` — they are mutually exclusive."
                        ));
                    }
                    (None, None) => {
                        return Err(anyhow::anyhow!(
                            "You must specify either `offsets` or `neighbor_order` in exchange rule."
                        ));
                    }
                }
            }
        }
        Ok(exchange_params)
    }
}
