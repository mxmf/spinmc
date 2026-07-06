use super::*;

#[test]
fn config_new_minimal_ising() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0, 2.0]
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
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    assert_eq!(config.grid.dimensions, [2, 2, 1]);
    assert_eq!(config.parsed_exchange.len(), 1);
    assert!(config.parsed_anisotropy.is_empty());
}

#[test]
fn config_new_with_anisotropy() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
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
offsets = [[1, 0, 0]]
strength = 1.0

[anisotropy]
axis = [[0.0, 0.0, 1.0]]
strength = [2.0]

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    assert_eq!(config.parsed_anisotropy.len(), 1);
    assert_eq!(config.parsed_anisotropy[0].strength, 2.0);
}

#[test]
fn config_new_invalid_toml() {
    assert!(Config::new("invalid {{{ toml").is_err());
}

#[test]
fn config_new_missing_model() {
    let toml = r#"
[simulation]
initial_state = "random"
equilibration_steps = 100
measurement_steps = 1000
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
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    assert!(Config::new(toml).is_err());
}

#[test]
fn config_new_with_temperature_range() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
num_threads = 1
algorithm = "metropolis"

[[simulation.temperature_range]]
start = 1.0
end = 3.0
step = 1.0

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    assert_eq!(config.simulation.temperatures.len(), 3);
}

#[test]
fn config_new_cell_structure() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 2
spin_magnitudes = [1.0, 1.0]
periodic_boundary = [true, true, true]

[structure]
cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    assert!(config.structure.is_some());
}

#[test]
fn display_config() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
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
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let s = format!("{config}");
    assert!(s.contains("Grid"));
    assert!(s.contains("Simulation"));
}

#[test]
fn display_config_with_structure_exchange_and_anisotropy() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 2
spin_magnitudes = [1.0, 1.0]
periodic_boundary = [true, true, true]

[structure]
cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]

[[exchange]]
from_sublattice = 0
to_sublattice = 1
offsets = [[0, 0, 0]]
strength = 1.0

[anisotropy]
axis = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]
strength = [2.0, 3.0]

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let s = format!("{config}");
    assert!(s.contains("Structure"));
    assert!(s.contains("Exchange Parameters"));
    assert!(s.contains("Anisotropy Parameters"));
    assert!(s.contains("ion0"));
    assert!(s.contains("ion1"));
}
