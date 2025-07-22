use anyhow::Result;
use clap::Parser;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

use MC_Curie::{
    config::{self, Config},
    lattice::Grid,
    monte_carlo::{Metropolis, MonteCarlo, Stats},
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

    let results = run_simulation(&run_config)?;

    // todo Output struct and output writer
    let file = File::create(&run_config.outfile)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "# t: 温度, E: 平均能量(eV), Cv: 比热容(eV/K), M: 平均磁化强度, χ: 磁化率(eV^-1)"
    )?;
    for (t, e, cv, m, chi) in results.iter() {
        writeln!(writer, "{t:.4}\t{e:.6}\t{cv:.6}\t{m:.6}\t{chi:.6}",)?;
    }

    Ok(())
}

// todo output struct
type Type = anyhow::Result<Vec<(f64, f64, f64, f64, f64)>>;

fn run_simulation(run_config: &Config) -> Type {
    ThreadPoolBuilder::new()
        .num_threads(run_config.num_threads)
        .build_global()
        .unwrap();

    let results = run_config
        .temperatures
        .par_iter()
        .map(|t| {
            let rng = rand::rng();
            let mut grid = match run_config.model {
                config::Model::Ising => Grid::<IsingSpin, _>::new(run_config.clone(), rng.clone()),
                _ => {
                    unimplemented!("xy and Heisenberg model")
                }
            };

            let beta = 1. / (run_config.kb * t);
            let mut mc = Metropolis { rng, beta };

            for _ in 0..run_config.n_equil {
                mc.step(&mut grid);
            }

            let mut stats = Stats::new(grid.size);
            for _ in 0..run_config.n_steps {
                mc.step(&mut grid);
                stats.record(&grid);
            }

            (
                *t,
                stats.mean_energy(),
                stats.specific_heat(beta),
                stats.mean_magnetization(),
                stats.susceptibility(beta),
            )
        })
        .collect();
    Ok(results)
}
