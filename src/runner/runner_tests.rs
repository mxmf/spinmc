use super::*;

fn energy_stats_config(group_num: usize) -> StatsConfig {
    StatsConfig {
        energy: true,
        heat_capacity: false,
        magnetization: false,
        susceptibility: false,
        magnetization_abs: false,
        susceptibility_abs: false,
        group_magnetization: false,
        group_susceptibility: false,
        group_magnetization_abs: false,
        group_susceptibility_abs: false,
        group_num,
    }
}

fn unique_temp_file(prefix: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "{prefix}_{}_{}.txt",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ))
}

fn non_comment_lines(content: &str) -> Vec<&str> {
    content
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .collect()
}

fn progress_config_toml(extra_output: &str) -> String {
    format!(
        r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 10
measurement_steps = 10
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
group = [[0]]
{extra_output}
"#
    )
}

#[test]
fn auto_progress_log_interval_normal() {
    assert_eq!(auto_progress_log_interval(100), 5);
}

#[test]
fn auto_progress_log_interval_less_than_20() {
    assert_eq!(auto_progress_log_interval(10), 1);
}

#[test]
fn auto_progress_log_interval_zero() {
    assert_eq!(auto_progress_log_interval(0), 1);
}

#[test]
fn should_log_progress_at_multiple() {
    assert!(should_log_progress(10, 100, 5));
}

#[test]
fn should_log_progress_at_completion() {
    assert!(should_log_progress(100, 100, 5));
}

#[test]
fn should_log_progress_not_at_non_multiple() {
    assert!(!should_log_progress(11, 100, 5));
}

#[test]
fn should_log_progress_interval_zero_never() {
    assert!(!should_log_progress(10, 100, 0));
}

#[test]
fn progress_config_uses_explicit_log_interval() {
    let toml = progress_config_toml(
        r#"
progress_bar = false
progress_log_interval = 7
"#,
    );
    let config = Config::new(&toml).unwrap();
    let progress = progress_config_with_terminal(&config, 20, true);
    assert!(!progress.use_bars);
    assert_eq!(progress.log_interval, 7);
}

#[test]
fn progress_config_tty_uses_progress_bar_without_log_interval() {
    let toml = progress_config_toml("");
    let config = Config::new(&toml).unwrap();
    let progress = progress_config_with_terminal(&config, 100, true);
    assert!(progress.use_bars);
    assert_eq!(progress.log_interval, 0);
}

#[test]
fn progress_config_non_tty_uses_auto_log_interval() {
    let toml = progress_config_toml("");
    let config = Config::new(&toml).unwrap();
    let progress = progress_config_with_terminal(&config, 100, false);
    assert!(!progress.use_bars);
    assert_eq!(progress.log_interval, 5);
}

#[test]
fn build_systems_correct_number_of_replicas() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 10
measurement_steps = 10
temperatures = [1.0, 2.0, 3.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let stats_config = energy_stats_config(1);
    let sys = build_systems::<IsingSpin>(&config, &stats_config).unwrap();
    assert_eq!(sys.stats.len(), 3);
    assert_eq!(sys.grids.len(), 3);
    assert_eq!(sys.algos.len(), 3);
}

#[test]
fn build_systems_beta_calculation() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 10
measurement_steps = 10
temperatures = [2.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let stats_config = energy_stats_config(1);
    let sys = build_systems::<IsingSpin>(&config, &stats_config).unwrap();
    let beta = sys.algos[0].beta();
    let expected = 1.0 / (8.617333262145e-5 * 2.0);
    assert!((beta - expected).abs() < 1e-10);
}

#[test]
fn build_systems_wolff_algorithm() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 10
measurement_steps = 10
temperatures = [1.0]
num_threads = 1
algorithm = "wolff"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let stats_config = energy_stats_config(1);
    let sys = build_systems::<IsingSpin>(&config, &stats_config).unwrap();
    assert_eq!(sys.algos.len(), 1);
    assert!(matches!(sys.algos[0], AnyMC::Wolff(_)));
}

#[test]
fn run_simulations_parallel_tempering_returns_temperature_order() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 2
measurement_steps = 3
temperatures = [1.0, 2.0]
num_threads = 1
pt_interval = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
stats_interval = 1
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let stats_config = energy_stats_config(1);
    let results = run_simulations::<IsingSpin>(&config, &stats_config).unwrap();

    assert_eq!(results.len(), 2);
    assert_eq!(results[0].t, 1.0);
    assert_eq!(results[1].t, 2.0);
    assert!(results.iter().all(|result| result.energy.is_some()));
    assert!(
        results
            .iter()
            .all(|result| result.energy.unwrap().is_finite())
    );
}

#[test]
fn run_end_to_end_ising() {
    let savefile = unique_temp_file("spinmc_test_result");
    let savefile_str = savefile.to_str().unwrap();
    let toml = format!(
        r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 5
measurement_steps = 5
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
savefile = "{savefile_str}"
group = [[0]]
"#
    );
    run(&toml).unwrap();
    // Verify output file was created and is not empty
    let content = std::fs::read_to_string(&savefile).unwrap();
    assert!(content.contains("#T(K)"));
    assert!(content.contains("Energy"));
    let result_lines = non_comment_lines(&content);
    assert_eq!(result_lines.len(), 1);
    assert!(result_lines[0].starts_with("1.000000"));
    // Clean up
    let _ = std::fs::remove_file(savefile);
}

#[test]
fn run_invalid_toml_errors() {
    assert!(run("invalid toml {{{").is_err());
}
