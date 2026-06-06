use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar, ProgressStyle};
use rand::Rng;
use rand_core::SeedableRng;
use rand_pcg::Pcg64Mcg;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use tracing::info;

use crate::{
    config::{self, Algorithm, Config},
    lattice::Grid,
    monte_carlo::{AnyMC, Metropolis, MonteCarlo, StatResult, Stats, StatsConfig, Wolff},
    spin::{HeisenbergSpin, IsingSpin, SpinState, XYSpin},
};

/// Allows sharing &mut [Stats<S>] across threads via raw pointer.
/// Each thread writes to disjoint elements (guaranteed by temp_to_replica being a permutation).
struct StatsRef<S: SpinState>(*mut Stats<S>);
unsafe impl<S: SpinState + Send> Send for StatsRef<S> {}
impl<S: SpinState> Clone for StatsRef<S> {
    fn clone(&self) -> Self {
        StatsRef(self.0)
    }
}

pub fn run(content: &str) -> anyhow::Result<()> {
    let run_config = Config::new(content)?;
    info!("{run_config}");

    let stats_config = StatsConfig {
        energy: run_config.output.energy,
        heat_capacity: run_config.output.heat_capacity,
        magnetization: run_config.output.magnetization,
        susceptibility: run_config.output.susceptibility,
        magnetization_abs: run_config.output.magnetization_abs,
        susceptibility_abs: run_config.output.susceptibility_abs,
        group_magnetization: run_config.output.group_magnetization,
        group_susceptibility: run_config.output.group_susceptibility,
        group_magnetization_abs: run_config.output.group_magnetization_abs,
        group_susceptibility_abs: run_config.output.group_susceptibility_abs,
        group_num: run_config.output.group.len(),
    };

    ThreadPoolBuilder::new()
        .num_threads(run_config.simulation.num_threads)
        .build_global()
        .unwrap();

    let results = match run_config.simulation.model {
        config::Model::Ising => run_simulations::<IsingSpin>(&run_config, &stats_config),
        config::Model::Xy => run_simulations::<XYSpin>(&run_config, &stats_config),
        config::Model::Heisenberg => run_simulations::<HeisenbergSpin>(&run_config, &stats_config),
    }?;

    let file = File::create(&run_config.output.savefile)?;

    let mut writer = BufWriter::new(&file);

    writeln!(writer, "{stats_config}")?;

    for result in results.iter() {
        writeln!(writer, "{result}")?;
    }

    info!(
        "Simulation completed. Results saved to file: {}",
        run_config.output.savefile
    );
    Ok(())
}

struct Systems<S: SpinState> {
    stats: Vec<Stats<S>>,
    grids: Vec<Grid<S, Pcg64Mcg>>,
    algos: Vec<AnyMC<Pcg64Mcg>>,
}

fn build_systems<S: SpinState>(
    config: &Config,
    stats_config: &StatsConfig,
) -> anyhow::Result<Systems<S>> {
    let kb = config.simulation.boltzmann_constant;
    let mut stats = Vec::new();
    let mut grids = Vec::new();
    let mut algos = Vec::new();

    for &t in &config.simulation.temperatures {
        let beta = 1.0 / (kb * t);
        let rng = Pcg64Mcg::from_rng(&mut rand::rng());
        let grid = Grid::<S, Pcg64Mcg>::new(config, rng.clone())?;
        let mc = match config.simulation.algorithm {
            Algorithm::Wolff => AnyMC::Wolff(Wolff {
                rng,
                beta,
                ham_config: grid.hamiltonian.config,
            }),
            Algorithm::Metropolis => AnyMC::Metropolis(Metropolis { rng, beta }),
        };
        stats.push(Stats::<S>::new(config, t, stats_config.clone()));
        grids.push(grid);
        algos.push(mc);
    }
    Ok(Systems {
        stats,
        grids,
        algos,
    })
}

// Non-PT: into_par_iter().map() — each thread owns its data
fn run_independent<S: SpinState, R: rand::Rng + Clone + SeedableRng + Send>(
    config: &Config,
    stats: Vec<Stats<S>>,
    algos: Vec<AnyMC<R>>,
    grids: Vec<Grid<S, R>>,
) -> Vec<StatResult> {
    let equil_steps = config.simulation.equilibration_steps;
    let meas_steps = config.simulation.measurement_steps;
    let stats_interval = config.output.stats_interval;
    let total_steps = equil_steps + meas_steps;
    let num_threads = rayon::current_num_threads();
    let progress_interval = (total_steps / 100).max(1) as u64;

    let multi = MultiProgress::new();
    let pb = multi.add(ProgressBar::new(stats.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} temperatures ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    let sub_pbs: Vec<ProgressBar> = (0..num_threads)
        .map(|i| {
            let sp = multi.insert(i + 1, ProgressBar::new(total_steps as u64));
            sp.set_style(
                ProgressStyle::default_bar()
                    .template("  {msg} [{bar:20.cyan/dim}] {pos}/{len}")
                    .unwrap()
                    .progress_chars("█░"),
            );
            sp
        })
        .collect();
    let sub_counter = AtomicUsize::new(0);

    stats
        .into_par_iter()
        .zip(algos.into_par_iter())
        .zip(grids.into_par_iter())
        .progress_with(pb)
        .enumerate()
        .map(|(idx, ((mut stat, mut mc), mut grid))| {
            let bar_id = sub_counter.fetch_add(1, Ordering::Relaxed) % num_threads;
            let sub_pb = &sub_pbs[bar_id];
            sub_pb.reset();
            sub_pb.set_message(format!("T={:.4}", config.simulation.temperatures[idx]));
            #[cfg(not(feature = "snapshots"))]
            let _ = idx;
            #[cfg(feature = "snapshots")]
            let t = config.simulation.temperatures[idx];
            #[cfg(feature = "snapshots")]
            let (mut equil_snapshots, mut measure_snapshots) = (vec![], vec![]);

            for _step in 0..equil_steps {
                mc.step(&mut grid);
                if _step as u64 % progress_interval == 0 {
                    sub_pb.set_position(_step as u64);
                }

                #[cfg(feature = "snapshots")]
                if let Some(snapshots) = &config.snapshots
                    && snapshots.equilibration_interval > 0
                    && _step % snapshots.equilibration_interval == 0
                {
                    equil_snapshots.push(grid.spins_to_array());
                }
            }
            sub_pb.set_position(equil_steps as u64);
            for step in 0..meas_steps {
                mc.step(&mut grid);
                if step % stats_interval == 0 {
                    stat.record(&grid);
                }
                if step as u64 % progress_interval == 0 {
                    sub_pb.set_position(equil_steps as u64 + step as u64);
                }

                #[cfg(feature = "snapshots")]
                if let Some(snapshots) = &config.snapshots
                    && snapshots.measurement_interval > 0
                    && step % snapshots.measurement_interval == 0
                {
                    measure_snapshots.push(grid.spins_to_array());
                }
            }

            #[cfg(feature = "snapshots")]
            if let Some(snapshots) = &config.snapshots {
                let snapshot_dir = &snapshots.save_directory;
                std::fs::create_dir_all(snapshot_dir).unwrap();
                let file_name = format!("{snapshot_dir}/T_{t:.4}.h5");
                match config::save_snapshots_to_hdf5(
                    &file_name,
                    &equil_snapshots,
                    &measure_snapshots,
                ) {
                    Ok(_) => info!("Saved snapshots to file {file_name} successfully"),
                    Err(e) => {
                        info!("Failed to save snapshots to file {file_name} because {e}")
                    }
                };
            };
            sub_pb.set_position(total_steps as u64);
            sub_pb.finish_with_message(format!("T={:.4} ✓", config.simulation.temperatures[idx]));

            stat.result()
        })
        .collect()
}

// PT: batched MC (pt_interval steps per fork-join) + parallel measurement + even-odd swap
fn run_pt<S: SpinState, R: rand::Rng + Clone + SeedableRng + Send + Sync>(
    config: &Config,
    stats: &mut [Stats<S>],
    algos: &mut [AnyMC<R>],
    grids: &mut [Grid<S, R>],
) {
    let n_temps = config.simulation.temperatures.len();
    let pt_interval = config.simulation.pt_interval;
    let equil_steps = config.simulation.equilibration_steps;
    let meas_steps = config.simulation.measurement_steps;
    let stats_interval = config.output.stats_interval;
    let total_steps = equil_steps + meas_steps;

    let mut temp_to_replica: Vec<usize> = (0..n_temps).collect();

    let pb = ProgressBar::new(total_steps as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} sweeps ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let mut sweep = 0usize;
    while sweep < total_steps {
        let batch_end = (sweep + pt_interval).min(total_steps);
        let batch_nsteps = batch_end - sweep;
        let start = sweep;

        // Batch: pt_interval MC steps + measurements, ONE fork-join, cache stays hot
        let sr = StatsRef(stats.as_mut_ptr());
        let ttr_ptr = temp_to_replica.as_ptr() as usize;
        grids
            .par_iter_mut()
            .zip(algos.par_iter_mut())
            .enumerate()
            .for_each_with((sr, ttr_ptr), |(sr, ttrp), (r, (grid, mc))| {
                let stats = unsafe { std::slice::from_raw_parts_mut(sr.0, n_temps) };
                let ttr = unsafe { std::slice::from_raw_parts(*ttrp as *const usize, n_temps) };

                for offset in 0..batch_nsteps {
                    let s = start + offset;
                    mc.step(grid);
                    if s >= equil_steps {
                        let do_meas = stats_interval == 0 || s % stats_interval == 0;
                        if do_meas {
                            for t in 0..n_temps {
                                if ttr[t] == r {
                                    stats[t].record(grid);
                                    break;
                                }
                            }
                        }
                    }
                }
            });

        // PT swap between batches (serial)
        for s in start..batch_end {
            if s > 0 && s % pt_interval == 0 {
                let mut rng = Pcg64Mcg::from_rng(&mut rand::rng());
                let swap_start = (s / pt_interval) % 2;
                for t in (swap_start..n_temps - 1).step_by(2) {
                    let i = temp_to_replica[t];
                    let j = temp_to_replica[t + 1];
                    let e_i = grids[i].total_energy();
                    let e_j = grids[j].total_energy();
                    let beta_i = algos[i].beta();
                    let beta_j = algos[j].beta();
                    let delta = (e_i - e_j) * (beta_i - beta_j);
                    if delta >= rng.random::<f64>().ln() {
                        algos[i].set_beta(beta_j);
                        algos[j].set_beta(beta_i);
                        temp_to_replica.swap(t, t + 1);
                    }
                }
            }
        }

        sweep = batch_end;
        pb.set_position(sweep as u64);
    }
    pb.finish_with_message("done");
}

// Dispatch
fn run_simulations<S: SpinState>(
    config: &Config,
    stats_config: &StatsConfig,
) -> anyhow::Result<Vec<StatResult>> {
    let Systems {
        mut stats,
        mut grids,
        mut algos,
    } = build_systems::<S>(config, stats_config)?;

    if config.simulation.pt_interval > 0 {
        info!(
            "PT enabled, swap every {} sweeps",
            config.simulation.pt_interval
        );
        run_pt(config, &mut stats, &mut algos, &mut grids);
        Ok(stats.into_iter().map(|s| s.result()).collect())
    } else {
        info!("Non-PT mode, independent temperatures");
        Ok(run_independent(config, stats, algos, grids))
    }
}
