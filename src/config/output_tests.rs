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
