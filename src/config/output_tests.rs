use super::*;

#[test]
fn validate_all_false_errors() {
    let o = Output {
        savefile: Default::default(),
        energy: false,
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
    let err = o.validate().unwrap_err().to_string();
    assert!(err.contains("No output fields specified"));
}

#[test]
fn validate_one_enough() {
    let o = Output {
        savefile: Default::default(),
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
    assert!(o.validate().is_ok());
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
