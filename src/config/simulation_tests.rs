use super::*;

fn simulation_with_temperature_range(start: f64, end: f64, step: f64) -> Simulation {
    Simulation {
        initial_state: InitialState::Random,
        model: Model::Ising,
        equilibration_steps: 100,
        measurement_steps: 1000,
        temperatures: vec![],
        temperature_range: vec![TemperatureRange { start, end, step }],
        num_threads: 1,
        pt_interval: 0,
        algorithm: Algorithm::Metropolis,
        boltzmann_constant: 1.0,
    }
}

#[test]
fn validate_temperatures_ok() {
    let mut s = Simulation {
        initial_state: InitialState::Random,
        model: Model::Ising,
        equilibration_steps: 100,
        measurement_steps: 1000,
        temperatures: vec![1.0, 2.0, 3.0],
        temperature_range: vec![],
        num_threads: 1,
        pt_interval: 0,
        algorithm: Algorithm::Metropolis,
        boltzmann_constant: 1.0,
    };
    assert!(s.validate().is_ok());
}

#[test]
fn validate_temperature_range_expands() {
    let mut s = simulation_with_temperature_range(1.0, 3.0, 1.0);
    assert!(s.validate().is_ok());
    assert_eq!(s.temperatures.len(), 3);
    assert!((s.temperatures[0] - 1.0).abs() < 1e-10);
    assert!((s.temperatures[1] - 2.0).abs() < 1e-10);
    assert!((s.temperatures[2] - 3.0).abs() < 1e-10);
}

#[test]
fn validate_neither_specified_errors() {
    let mut s = Simulation {
        initial_state: InitialState::Random,
        model: Model::Ising,
        equilibration_steps: 100,
        measurement_steps: 1000,
        temperatures: vec![],
        temperature_range: vec![],
        num_threads: 1,
        pt_interval: 0,
        algorithm: Algorithm::Metropolis,
        boltzmann_constant: 1.0,
    };
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("Either 'temperatures' or 'temperature_range'"));
}

#[test]
fn validate_both_specified_errors() {
    let mut s = Simulation {
        initial_state: InitialState::Random,
        model: Model::Ising,
        equilibration_steps: 100,
        measurement_steps: 1000,
        temperatures: vec![1.0],
        temperature_range: vec![TemperatureRange {
            start: 1.0,
            end: 3.0,
            step: 1.0,
        }],
        num_threads: 1,
        pt_interval: 0,
        algorithm: Algorithm::Metropolis,
        boltzmann_constant: 1.0,
    };
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("Only one of 'temperatures' or 'temperature_range'"));
}

#[test]
fn validate_temperature_range_zero_step_errors() {
    let mut s = simulation_with_temperature_range(1.0, 3.0, 0.0);
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("step"));
    assert!(err.contains("positive"));
}

#[test]
fn validate_temperature_range_negative_step_errors() {
    let mut s = simulation_with_temperature_range(1.0, 3.0, -1.0);
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("step"));
    assert!(err.contains("positive"));
}

#[test]
fn validate_temperature_range_start_greater_than_end_errors() {
    let mut s = simulation_with_temperature_range(3.0, 1.0, 1.0);
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("start"));
    assert!(err.contains("end"));
}

#[test]
fn display() {
    let mut s = Simulation {
        initial_state: InitialState::Z,
        model: Model::Heisenberg,
        equilibration_steps: 100,
        measurement_steps: 1000,
        temperatures: vec![1.0, 2.0],
        temperature_range: vec![],
        num_threads: 4,
        pt_interval: 0,
        algorithm: Algorithm::Wolff,
        boltzmann_constant: 1.0,
    };
    s.validate().unwrap();
    let output = format!("{s}");
    assert!(output.contains("Simulation"));
    assert!(output.contains("Heisenberg"));
    assert!(output.contains("Wolff"));
}

#[test]
fn display_parallel_tempering_enabled() {
    let s = Simulation {
        initial_state: InitialState::Z,
        model: Model::Ising,
        equilibration_steps: 10,
        measurement_steps: 20,
        temperatures: vec![1.0],
        temperature_range: vec![],
        num_threads: 2,
        pt_interval: 5,
        algorithm: Algorithm::Metropolis,
        boltzmann_constant: 1.0,
    };
    let output = format!("{s}");
    assert!(output.contains("enabled"));
    assert!(output.contains("swap every 5"));
}
