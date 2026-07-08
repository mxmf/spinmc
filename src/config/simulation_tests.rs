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

fn simulation_with_temperatures(temperatures: Vec<f64>) -> Simulation {
    Simulation {
        initial_state: InitialState::Random,
        model: Model::Ising,
        equilibration_steps: 100,
        measurement_steps: 1000,
        temperatures,
        temperature_range: vec![],
        num_threads: 1,
        pt_interval: 0,
        algorithm: Algorithm::Metropolis,
        boltzmann_constant: 1.0,
    }
}

#[test]
fn validate_temperatures_ok() {
    let mut s = simulation_with_temperatures(vec![1.0, 2.0, 3.0]);
    assert!(s.validate().is_ok());
}

#[test]
fn validate_zero_temperature_ok() {
    let mut s = simulation_with_temperatures(vec![0.0, 1.0]);
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
fn validate_temperature_range_non_finite_errors() {
    for (start, end, step, expected) in [
        (f64::NAN, 3.0, 1.0, "start"),
        (1.0, f64::INFINITY, 1.0, "end"),
        (1.0, 3.0, f64::NAN, "step"),
    ] {
        let mut s = simulation_with_temperature_range(start, end, step);
        let err = s.validate().unwrap_err().to_string();
        assert!(err.contains(expected));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_temperature_range_zero_start_expands() {
    let mut s = simulation_with_temperature_range(0.0, 2.0, 1.0);
    assert!(s.validate().is_ok());
    assert_eq!(s.temperatures, vec![0.0, 1.0, 2.0]);
}

#[test]
fn validate_temperature_range_negative_start_errors() {
    let mut s = simulation_with_temperature_range(-1.0, 3.0, 1.0);
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("start"));
    assert!(err.contains("non-negative"));
}

#[test]
fn validate_negative_temperature_errors() {
    let mut s = simulation_with_temperatures(vec![1.0, -1.0]);
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("temperatures[1]"));
    assert!(err.contains("non-negative"));
}

#[test]
fn validate_non_finite_temperature_errors() {
    for temperature in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut s = simulation_with_temperatures(vec![temperature]);
        let err = s.validate().unwrap_err().to_string();
        assert!(err.contains("temperatures[0]"));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_zero_measurement_steps_errors() {
    let mut s = simulation_with_temperatures(vec![1.0]);
    s.measurement_steps = 0;
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("measurement_steps"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_zero_num_threads_errors() {
    let mut s = simulation_with_temperatures(vec![1.0]);
    s.num_threads = 0;
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("num_threads"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_non_positive_boltzmann_constant_errors() {
    let mut s = simulation_with_temperatures(vec![1.0]);
    s.boltzmann_constant = 0.0;
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("boltzmann_constant"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_non_finite_boltzmann_constant_errors() {
    for boltzmann_constant in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut s = simulation_with_temperatures(vec![1.0]);
        s.boltzmann_constant = boltzmann_constant;
        let err = s.validate().unwrap_err().to_string();
        assert!(err.contains("boltzmann_constant"));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_step_count_overflow_errors() {
    let mut s = simulation_with_temperatures(vec![1.0]);
    s.equilibration_steps = usize::MAX;
    s.measurement_steps = 1;
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("equilibration_steps"));
    assert!(err.contains("measurement_steps"));
    assert!(err.contains("too large"));
}

#[test]
fn validate_parallel_tempering_requires_multiple_temperatures() {
    let mut s = simulation_with_temperatures(vec![1.0]);
    s.pt_interval = 10;
    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("parallel tempering"));
    assert!(err.contains("at least two temperatures"));
}

#[test]
fn validate_parallel_tempering_allows_multiple_temperatures() {
    let mut s = simulation_with_temperatures(vec![1.0, 2.0]);
    s.pt_interval = 10;
    assert!(s.validate().is_ok());
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
