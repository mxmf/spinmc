use anyhow::Result;
use clap::Parser;
use rand_core::SeedableRng;
use rand_pcg::Pcg64Mcg;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use tracing::info;
use tracing_subscriber::FmtSubscriber;

use mc_curie::{
    config::{self, Config},
    lattice::Grid,
    monte_carlo::{Metropolis, MonteCarlo, StatResult, Stats, StatsConfig},
    spin::IsingSpin,
};

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Path to config file (TOML format)
    #[arg(short, long, default_value = "config.toml")]
    config: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let run_config = Config::new(&args.config)?;

    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");
    info!("{run_config}");

    let stats_config = StatsConfig {
        energy: run_config.energy,
        heat_capacity: run_config.heat_capacity,
        magnetization: run_config.magnetization,
        susceptibility: run_config.susceptibility,
        magnetization_abs: run_config.magnetization_abs,
        susceptibility_abs: run_config.susceptibility_abs,
        group_magnetization: run_config.group_magnetization,
        group_susceptibility: run_config.group_susceptibility,
        group_num: run_config.group.len(),
    };

    let results = run_simulation(&run_config, &stats_config)?;

    let file = File::create(&run_config.outfile)?;

    let mut writer = BufWriter::new(file);

    writeln!(writer, "{stats_config}",)?;

    for result in results.iter() {
        writeln!(writer, "{result}",)?;
    }

    Ok(())
}

type Type = anyhow::Result<Vec<StatResult>>;

fn run_simulation(run_config: &Config, stats_config: &StatsConfig) -> Type {
    ThreadPoolBuilder::new()
        .num_threads(run_config.num_threads)
        .build_global()
        .unwrap();

    info!("Start run simulations");
    let results: Vec<StatResult> = run_config
        .temperatures
        .par_iter()
        .map(|t| {
            // TODO add more  rng method
            let rng = Pcg64Mcg::from_rng(&mut rand::rng());

            let (mut grid, mut stats) = match run_config.model {
                config::Model::Ising => {
                    let stats = Stats::new::<IsingSpin>(run_config, *t, stats_config.clone());
                    let grid = Grid::<IsingSpin, _>::new(run_config.clone(), rng.clone());
                    (grid, stats)
                }
                _ => {
                    unimplemented!("xy and Heisenberg model")
                }
            };

            let beta = 1. / (run_config.kb * t);
            let mut mc = Metropolis { rng, beta };

            for _ in 0..run_config.n_equil {
                mc.step(&mut grid);
            }

            for _ in 0..run_config.n_steps {
                mc.step(&mut grid);
                stats.record(&grid);
            }
            info!("Simulation on temperature {t} K fininshed");
            stats.result()
        })
        .collect();

    Ok(results)
}
