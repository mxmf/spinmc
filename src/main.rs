use anyhow::Result;
use clap::Parser;
use mc_curie::spin::SpinState;
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
    spin::{IsingSpin, XYSpin},
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

    let results = run_parallel_simulations(&run_config, &stats_config)?;

    let file = File::create(&run_config.outfile)?;

    let mut writer = BufWriter::new(file);

    writeln!(writer, "{stats_config}",)?;

    for result in results.iter() {
        writeln!(writer, "{result}",)?;
    }

    Ok(())
}

type Type = anyhow::Result<Vec<StatResult>>;

fn run_parallel_simulations(run_config: &Config, stats_config: &StatsConfig) -> Type {
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

            match run_config.model {
                config::Model::Ising => {
                    let stats = Stats::new::<IsingSpin>(run_config, *t, stats_config.clone());
                    let mut grid = Grid::<IsingSpin, _>::new(run_config.clone(), rng.clone());
                    run_single_simulate::<IsingSpin, _>(&mut grid, stats, run_config, *t, rng)
                }
                config::Model::Xy => {
                    let stats = Stats::new::<XYSpin>(run_config, *t, stats_config.clone());
                    let mut grid = Grid::<XYSpin, _>::new(run_config.clone(), rng.clone());
                    run_single_simulate::<XYSpin, _>(&mut grid, stats, run_config, *t, rng)
                }
                _ => {
                    unimplemented!("xy and Heisenberg model")
                }
            }
        })
        .collect();

    Ok(results)
}

fn run_single_simulate<S: SpinState, R: rand::Rng>(
    grid: &mut Grid<S, R>,
    mut stats: Stats,
    run_config: &Config,
    t: f64,
    rng: R,
) -> StatResult {
    let beta = 1. / (run_config.kb * t);
    let mut mc = Metropolis { rng, beta };

    #[cfg(feature = "snapshots")]
    let (mut equil_snapshots, mut steps_snapshots) = (vec![], vec![]);

    for _step in 0..run_config.n_equil {
        mc.step(grid);
        #[cfg(feature = "snapshots")]
        {
            if run_config.snapshot_enable
                && run_config.snapshot_params.snapshot_equil_interval > 0
                && _step % run_config.snapshot_params.snapshot_equil_interval == 0
            {
                equil_snapshots.push(grid.spins_to_array());
            }
        }
    }

    for _step in 0..run_config.n_steps {
        mc.step(grid);
        stats.record(grid);

        #[cfg(feature = "snapshots")]
        if run_config.snapshot_enable
            && run_config.snapshot_params.snapshot_equil_interval > 0
            && _step % run_config.snapshot_params.snapshot_equil_interval == 0
        {
            steps_snapshots.push(grid.spins_to_array());
        }
    }
    info!("Simulation at temperature {t:.4} K fininshed");

    #[cfg(feature = "snapshots")]
    if run_config.snapshot_enable {
        let snapshot_dir = &run_config.snapshot_params.snapshot_dir;
        std::fs::create_dir_all(snapshot_dir).unwrap();
        let file_name = format!("{snapshot_dir}/T_{t:.4}.h5");
        match mc_curie::snapshots::save_snapshots_to_hdf5(
            &file_name,
            &equil_snapshots,
            &steps_snapshots,
        ) {
            Ok(_) => info!("Saved snapshots to file {file_name} successfully"),
            Err(e) => {
                info!("Failed to save snapshots to file {file_name} because {e}")
            }
        };
    };

    stats.result()
}
