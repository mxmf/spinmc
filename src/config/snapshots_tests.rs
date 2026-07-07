use super::*;
use ndarray::{Array5, ArrayD};
use ndarray_npy::NpzReader;
use std::fs::File;

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
fn validate_accepts_supported_compression_levels() {
    for compression_level in 0..=9 {
        let s = Snapshots {
            equilibration_interval: 0,
            measurement_interval: 0,
            compression_level,
            save_directory: "snapshots".into(),
        };

        assert!(
            s.validate().is_ok(),
            "compression level {compression_level} should be supported"
        );
    }
}

#[test]
fn validate_rejects_compression_level_above_nine() {
    let s = Snapshots {
        equilibration_interval: 0,
        measurement_interval: 0,
        compression_level: 10,
        save_directory: "snapshots".into(),
    };

    let err = s.validate().unwrap_err().to_string();
    assert!(err.contains("compression_level"));
}

#[test]
fn display() {
    let s = Snapshots {
        equilibration_interval: 5,
        measurement_interval: 7,
        compression_level: 9,
        save_directory: "custom_snapshots".into(),
    };
    let out = format!("{s}");
    assert!(out.contains("Snapshots"));
    assert!(out.contains("Equilibration Interval: 5 steps"));
    assert!(out.contains("Sampling Interval: 7 steps"));
    assert!(out.contains("Compression Level: 9"));
    assert!(out.contains("Snapshot Directory: custom_snapshots"));
}

#[test]
fn deserialize_uses_default_save_directory() {
    let s: Snapshots = toml::from_str(
        r#"
equilibration_interval = 2
measurement_interval = 3
"#,
    )
    .unwrap();

    assert_eq!(s.equilibration_interval, 2);
    assert_eq!(s.measurement_interval, 3);
    assert_eq!(s.compression_level, 0);
    assert_eq!(s.save_directory, "snapshots");
}

// --- NPZ round-trip ---

fn unique_npz_path(prefix: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "{prefix}_{}_{}.npz",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ))
}

fn dummy_snap(shape: [usize; 5], offset: f64) -> Array5<f64> {
    Array5::from_shape_fn(shape, |idx| {
        offset
            + idx.0 as f64
            + idx.1 as f64 * 10.0
            + idx.2 as f64 * 100.0
            + idx.3 as f64 * 1000.0
            + idx.4 as f64 * 10000.0
    })
}

fn read_npz_arrays(filename: &str) -> (ArrayD<f64>, ArrayD<f64>) {
    let mut reader = NpzReader::new(File::open(filename).unwrap()).unwrap();
    let equil = reader.by_name("equil").unwrap();
    let steps = reader.by_name("steps").unwrap();
    (equil, steps)
}

#[test]
fn npz_round_trip_all_supported_compression_levels() {
    let equil = vec![
        dummy_snap([1, 2, 2, 2, 3], 0.0),
        dummy_snap([1, 2, 2, 2, 3], 100_000.0),
    ];
    let steps = vec![dummy_snap([1, 2, 2, 2, 3], 200_000.0)];

    for compression_level in 0..=9 {
        let file = unique_npz_path(&format!(
            "test_spinmc_npz_round_trip_level_{compression_level}"
        ));
        let filename = file.to_str().unwrap();

        save_snapshots_to_npz(filename, &equil, &steps, compression_level).unwrap();

        let (read_equil, read_steps) = read_npz_arrays(filename);

        assert_eq!(read_equil.shape(), &[2, 1, 2, 2, 2, 3]);
        assert_eq!(read_steps.shape(), &[1, 1, 2, 2, 2, 3]);
        assert_eq!(read_equil[[0, 0, 1, 1, 1, 2]], equil[0][[0, 1, 1, 1, 2]]);
        assert_eq!(read_equil[[1, 0, 1, 1, 1, 2]], equil[1][[0, 1, 1, 1, 2]]);
        assert_eq!(read_steps[[0, 0, 1, 1, 1, 2]], steps[0][[0, 1, 1, 1, 2]]);

        std::fs::remove_file(filename).unwrap();
    }
}

#[test]
fn npz_empty_snapshots() {
    let equil: Vec<Array5<f64>> = vec![];
    let steps: Vec<Array5<f64>> = vec![];

    let file = unique_npz_path("test_spinmc_npz_empty");
    let filename = file.to_str().unwrap();

    save_snapshots_to_npz(filename, &equil, &steps, 0).unwrap();

    let (read_equil, read_steps) = read_npz_arrays(filename);

    assert_eq!(read_equil.shape(), &[0, 0, 0, 0, 0, 3]);
    assert_eq!(read_steps.shape(), &[0, 0, 0, 0, 0, 3]);

    std::fs::remove_file(filename).unwrap();
}

#[test]
fn npz_empty_side_uses_non_empty_side_shape() {
    let snap = dummy_snap([2, 3, 2, 1, 3], 10.0);
    let file = unique_npz_path("test_spinmc_npz_empty_side");
    let filename = file.to_str().unwrap();

    save_snapshots_to_npz(filename, std::slice::from_ref(&snap), &[], 6).unwrap();

    let (read_equil, read_steps) = read_npz_arrays(filename);

    assert_eq!(read_equil.shape(), &[1, 2, 3, 2, 1, 3]);
    assert_eq!(read_steps.shape(), &[0, 2, 3, 2, 1, 3]);

    std::fs::remove_file(filename).unwrap();
}

#[test]
fn npz_named_arrays() {
    let snap = dummy_snap([1, 2, 2, 1, 3], 0.0);
    let file = unique_npz_path("test_spinmc_npz_names");
    let filename = file.to_str().unwrap();

    save_snapshots_to_npz(filename, &[snap], &[], 8).unwrap();

    let mut reader = NpzReader::new(File::open(filename).unwrap()).unwrap();
    let names = reader.names().unwrap();

    assert!(names.contains(&"equil".to_string()));
    assert!(names.contains(&"steps".to_string()));

    std::fs::remove_file(filename).unwrap();
}
