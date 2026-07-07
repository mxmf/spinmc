use super::*;

fn output_with_energy() -> Output {
    Output {
        savefile: "result.txt".into(),
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
        group: vec![],
        stats_interval: 1,
        progress_bar: true,
        progress_log_interval: 0,
    }
}

#[test]
fn validate_all_false_errors() {
    let mut o = output_with_energy();
    o.energy = false;
    let err = o.validate(1).unwrap_err().to_string();
    assert!(err.contains("No output fields specified"));
}

#[test]
fn validate_one_enough() {
    let o = output_with_energy();
    assert!(o.validate(1).is_ok());
}

#[test]
fn validate_empty_savefile_errors() {
    let mut o = output_with_energy();
    o.savefile = "  ".into();
    let err = o.validate(1).unwrap_err().to_string();
    assert!(err.contains("savefile"));
    assert!(err.contains("empty"));
}

#[test]
fn validate_zero_stats_interval_errors() {
    let mut o = output_with_energy();
    o.stats_interval = 0;
    let err = o.validate(1).unwrap_err().to_string();
    assert!(err.contains("stats_interval"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_group_observable_requires_groups() {
    let mut o = output_with_energy();
    o.group_magnetization = true;
    let err = o.validate(1).unwrap_err().to_string();
    assert!(err.contains("group output"));
    assert!(err.contains("at least one group"));
}

#[test]
fn validate_empty_group_errors() {
    let mut o = output_with_energy();
    o.group = vec![vec![]];
    let err = o.validate(1).unwrap_err().to_string();
    assert!(err.contains("group[0]"));
    assert!(err.contains("empty"));
}

#[test]
fn validate_group_sublattice_out_of_range_errors() {
    let mut o = output_with_energy();
    o.group = vec![vec![0, 2]];
    let err = o.validate(2).unwrap_err().to_string();
    assert!(err.contains("group[0]"));
    assert!(err.contains("2"));
    assert!(err.contains("out of range"));
}

#[test]
fn validate_group_output_ok() {
    let mut o = output_with_energy();
    o.group_magnetization = true;
    o.group = vec![vec![0], vec![1]];
    assert!(o.validate(2).is_ok());
}

#[test]
fn display() {
    let o = Output {
        savefile: "out.txt".into(),
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
        group: vec![],
        stats_interval: 1,
        progress_bar: true,
        progress_log_interval: 0,
    };
    let s = format!("{o}");
    assert!(s.contains("out.txt"));
    assert!(s.contains("Energy"));
}

#[test]
fn display_all_observables_and_groups() {
    let o = Output {
        savefile: "out.txt".into(),
        energy: true,
        heat_capacity: true,
        magnetization: true,
        susceptibility: true,
        magnetization_abs: true,
        susceptibility_abs: true,
        group_magnetization: true,
        group_susceptibility: true,
        group_magnetization_abs: true,
        group_susceptibility_abs: true,
        group: vec![vec![0], vec![1, 2]],
        stats_interval: 3,
        progress_bar: false,
        progress_log_interval: 7,
    };
    let s = format!("{o}");
    for expected in [
        "Heat Capacity",
        "Magnetization",
        "Susceptibility",
        "Magnetization_abs",
        "susceptibility_abs",
        "Group Magnetization",
        "Group Susceptibility",
        "Group |Magnetization|",
        "Group |Susceptibility|",
        "Group 0",
        "Group 1",
    ] {
        assert!(s.contains(expected), "missing {expected} in {s}");
    }
}
