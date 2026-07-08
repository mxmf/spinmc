use anyhow::Context;
use indicatif::{
    MultiProgress, ParallelProgressIterator, ProgressBar, ProgressDrawTarget, ProgressStyle,
};
use rand::RngExt;
use rand_core::SeedableRng;
use rand_pcg::Pcg64Mcg;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, IsTerminal, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use tracing::info;

use crate::{
    config::{self, Algorithm, Config},
    lattice::Grid,
    monte_carlo::{AnyMC, Metropolis, MonteCarlo, StatResult, Stats, StatsConfig, Wolff},
    spin::{HeisenbergSpin, IsingSpin, SpinState, XYSpin},
};

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

    let pool = build_thread_pool(run_config.simulation.num_threads)?;

    let results = pool.install(|| match run_config.simulation.model {
        config::Model::Ising => run_simulations::<IsingSpin>(&run_config, &stats_config),
        config::Model::Xy => run_simulations::<XYSpin>(&run_config, &stats_config),
        config::Model::Heisenberg => run_simulations::<HeisenbergSpin>(&run_config, &stats_config),
    })?;

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

fn build_thread_pool(num_threads: usize) -> anyhow::Result<rayon::ThreadPool> {
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .context("Failed to build Rayon thread pool")
}

struct Systems<S: SpinState> {
    stats: Vec<Stats<S>>,
    grids: Vec<Grid<S, Pcg64Mcg>>,
    algos: Vec<AnyMC<Pcg64Mcg>>,
}

#[derive(Clone, Copy)]
struct ProgressConfig {
    use_bars: bool,
    log_interval: usize,
}

fn progress_config(config: &Config, total_steps: usize) -> ProgressConfig {
    progress_config_with_terminal(config, total_steps, std::io::stderr().is_terminal())
}

fn progress_config_with_terminal(
    config: &Config,
    total_steps: usize,
    stderr_is_terminal: bool,
) -> ProgressConfig {
    let use_bars = config.output.progress_bar && stderr_is_terminal;
    let log_interval = if config.output.progress_log_interval > 0 {
        config.output.progress_log_interval
    } else if use_bars {
        0
    } else {
        auto_progress_log_interval(total_steps)
    };

    ProgressConfig {
        use_bars,
        log_interval,
    }
}

fn auto_progress_log_interval(total_steps: usize) -> usize {
    (total_steps / 20).max(1)
}

fn should_log_progress(completed: usize, total: usize, interval: usize) -> bool {
    interval > 0 && (completed == total || completed.is_multiple_of(interval))
}

fn maybe_hide_progress_bar(pb: &ProgressBar, progress: ProgressConfig) {
    if !progress.use_bars {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }
}

fn progress_style(template: &str, chars: &str) -> ProgressStyle {
    match ProgressStyle::default_bar().template(template) {
        Ok(style) => style.progress_chars(chars),
        Err(_) => ProgressStyle::default_bar().progress_chars(chars),
    }
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
) -> anyhow::Result<Vec<StatResult>> {
    let equil_steps = config.simulation.equilibration_steps;
    let meas_steps = config.simulation.measurement_steps;
    let stats_interval = config.output.stats_interval;
    let total_steps = equil_steps + meas_steps;
    let num_threads = rayon::current_num_threads();
    let progress = progress_config(config, total_steps);
    let temp_count = stats.len();

    let multi = MultiProgress::new();
    let pb = multi.add(ProgressBar::new(stats.len() as u64));
    maybe_hide_progress_bar(&pb, progress);
    pb.set_style(progress_style(
        "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} temperatures",
        "#>-",
    ));
    let sub_pbs: Vec<ProgressBar> = (0..num_threads)
        .map(|i| {
            let sp = multi.insert(i + 1, ProgressBar::new(total_steps as u64));
            maybe_hide_progress_bar(&sp, progress);
            sp.set_style(progress_style(
                "  {msg} [{elapsed_precise}/{eta}] [{bar:20.cyan/dim}] {pos}/{len}",
                "█░",
            ));
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
        .map(|(idx, ((mut stat, mut mc), mut grid))| -> anyhow::Result<StatResult> {
            let bar_id = sub_counter.fetch_add(1, Ordering::Relaxed) % num_threads;
            let sub_pb = &sub_pbs[bar_id];
            sub_pb.reset();
            sub_pb.set_message(format!("T={:.4}", config.simulation.temperatures[idx]));
            #[cfg(feature = "snapshots")]
            let t = config.simulation.temperatures[idx];
            #[cfg(feature = "snapshots")]
            let (mut equil_snapshots, mut measure_snapshots) = (vec![], vec![]);

            for _step in 0..equil_steps {
                mc.step(&mut grid);
                let completed = _step + 1;
                sub_pb.set_position(completed as u64);
                if should_log_progress(completed, total_steps, progress.log_interval) {
                    info!(
                        "Progress: temperature={}/{}, T={:.4}, phase=equilibration, sweep={}/{}, {:.1}%",
                        idx + 1,
                        temp_count,
                        config.simulation.temperatures[idx],
                        completed,
                        total_steps,
                        completed as f64 * 100.0 / total_steps as f64
                    );
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
                let completed = equil_steps + step + 1;
                sub_pb.set_position(completed as u64);
                if should_log_progress(completed, total_steps, progress.log_interval) {
                    info!(
                        "Progress: temperature={}/{}, T={:.4}, phase=measurement, sweep={}/{}, {:.1}%",
                        idx + 1,
                        temp_count,
                        config.simulation.temperatures[idx],
                        completed,
                        total_steps,
                        completed as f64 * 100.0 / total_steps as f64
                    );
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
                std::fs::create_dir_all(snapshot_dir).with_context(|| {
                    format!("Failed to create snapshot directory: {snapshot_dir}")
                })?;
                let file_name = format!("{snapshot_dir}/T_{t:.4}.npz");
                match config::save_snapshots_to_npz(
                    &file_name,
                    &equil_snapshots,
                    &measure_snapshots,
                    snapshots.compression_level,
                ) {
                    Ok(_) => info!("Saved snapshots to file {file_name} successfully"),
                    Err(e) => {
                        info!("Failed to save snapshots to file {file_name} because {e}")
                    }
                };
            };
            sub_pb.set_position(total_steps as u64);
            sub_pb.finish_with_message(format!("T={:.4} ✓", config.simulation.temperatures[idx]));

            Ok(stat.result())
        })
        .collect()
}

// PT: batched MC (pt_interval steps per fork-join) + parallel measurement + even-odd swap
// Each parallel thread writes to stats[r] directly (indexed by replica).
// When a PT swap accepts, stats entries are swapped alongside betas so that
// each Stats object follows the temperature it represents.
fn run_pt<S: SpinState, R: rand::Rng + Clone + SeedableRng + Send + Sync>(
    config: &Config,
    stats: &mut [Stats<S>],
    algos: &mut [AnyMC<R>],
    grids: &mut [Grid<S, R>],
) -> anyhow::Result<Vec<StatResult>> {
    let n_temps = config.simulation.temperatures.len();
    let pt_interval = config.simulation.pt_interval;
    let equil_steps = config.simulation.equilibration_steps;
    let meas_steps = config.simulation.measurement_steps;
    let stats_interval = config.output.stats_interval;
    let total_steps = equil_steps + meas_steps;
    let progress = progress_config(config, total_steps);

    // temp_to_replica[t] = replica index currently simulating temperature t
    let mut temp_to_replica: Vec<usize> = (0..n_temps).collect();

    #[cfg(feature = "snapshots")]
    let (mut equil_snapshots, mut measure_snapshots) = {
        let mut equil = Vec::with_capacity(n_temps);
        let mut meas = Vec::with_capacity(n_temps);
        for _ in 0..n_temps {
            equil.push(Vec::new());
            meas.push(Vec::new());
        }
        (equil, meas)
    };

    let pb = ProgressBar::new(total_steps as u64);
    maybe_hide_progress_bar(&pb, progress);
    pb.set_style(progress_style(
        "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} sweeps ({eta})",
        "#>-",
    ));

    let mut sweep = 0usize;
    while sweep < total_steps {
        let batch_end = (sweep + pt_interval).min(total_steps);
        let batch_nsteps = batch_end - sweep;
        let start = sweep;

        // Parallel batch: each replica writes to its own stats[r].
        // Stats are indexed by replica (same as grids/algos), so
        // par_iter_mut gives each thread exclusive access to its own stat.
        grids
            .par_iter_mut()
            .zip(algos.par_iter_mut())
            .zip(stats.par_iter_mut())
            .for_each(|((grid, mc), stat)| {
                for offset in 0..batch_nsteps {
                    let s = start + offset;
                    mc.step(grid);
                    if s >= equil_steps {
                        let do_meas = stats_interval == 0 || s.is_multiple_of(stats_interval);
                        if do_meas {
                            stat.record(grid);
                        }
                    }
                }
            });

        // PT swap between batches (serial).
        if batch_end < total_steps && batch_end.is_multiple_of(pt_interval) {
            let mut rng = Pcg64Mcg::from_rng(&mut rand::rng());
            let swap_start = (batch_end / pt_interval) % 2;
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
                    // Swap stats alongside temperatures so each Stats
                    // object tracks the same temperature throughout.
                    stats.swap(i, j);
                    temp_to_replica.swap(t, t + 1);
                }
            }
        }

        // Collect snapshots after batch + swap, grouped by temperature.
        // Note: when pt_interval > snapshot_interval, intermediate snapshots
        // within a batch are not captured.
        #[cfg(feature = "snapshots")]
        if let Some(snaps) = &config.snapshots {
            let is_equil = batch_end <= equil_steps;
            let interval = if is_equil {
                snaps.equilibration_interval
            } else {
                snaps.measurement_interval
            };
            if interval > 0 && batch_end.is_multiple_of(interval) {
                for t in 0..n_temps {
                    let r = temp_to_replica[t];
                    let snap = grids[r].spins_to_array();
                    if is_equil {
                        equil_snapshots[t].push(snap);
                    } else {
                        measure_snapshots[t].push(snap);
                    }
                }
            }
        }

        sweep = batch_end;
        pb.set_position(sweep as u64);
        if should_log_progress(sweep, total_steps, progress.log_interval) {
            info!(
                "Progress: parallel_tempering sweep={}/{}, {:.1}%",
                sweep,
                total_steps,
                sweep as f64 * 100.0 / total_steps as f64
            );
        }
    }
    pb.finish_with_message("done");

    // Save snapshots by temperature (not by replica).
    #[cfg(feature = "snapshots")]
    if let Some(snaps) = &config.snapshots {
        let snapshot_dir = &snaps.save_directory;
        std::fs::create_dir_all(snapshot_dir)
            .with_context(|| format!("Failed to create snapshot directory: {snapshot_dir}"))?;
        for t in 0..n_temps {
            let file_name = format!(
                "{snapshot_dir}/T_{:.4}.npz",
                config.simulation.temperatures[t]
            );
            match config::save_snapshots_to_npz(
                &file_name,
                &equil_snapshots[t],
                &measure_snapshots[t],
                snaps.compression_level,
            ) {
                Ok(_) => info!("Saved snapshots to file {file_name} successfully"),
                Err(e) => info!("Failed to save snapshots to file {file_name} because {e}"),
            };
        }
    }

    // Return results in temperature order.
    // After all swaps, stats[temp_to_replica[t]] holds data for temperature t.
    Ok((0..n_temps)
        .map(|t| {
            let r = temp_to_replica[t];
            stats[r].result()
        })
        .collect())
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
        run_pt(config, &mut stats, &mut algos, &mut grids)
    } else {
        info!("Non-PT mode, independent temperatures");
        run_independent(config, stats, algos, grids)
    }
}

#[cfg(test)]
#[path = "runner_tests.rs"]
mod tests;
