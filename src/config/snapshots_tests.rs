use super::*;

#[test]
fn validate_trivially_ok() {
    let s = Snapshots {
        equilibration_interval: 0,
        measurement_interval: 0,
        compression_level: 0,
        save_directory: "snapshots".into(),
    };
    assert!(s.validate().is_ok());
}

#[test]
fn display() {
    let s = Snapshots {
        equilibration_interval: 0,
        measurement_interval: 0,
        compression_level: 0,
        save_directory: "snapshots".into(),
    };
    let out = format!("{s}");
    assert!(out.contains("Snapshots"));
}
