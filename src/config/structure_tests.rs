use super::*;

fn poscar_path() -> String {
    // Relative to the workspace root (where tests are run from)
    "examples/heisenberg_2d_cri3_poscar/POSCAR".to_string()
}

#[test]
fn parse_from_poscar_file() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let s = sc.parse().unwrap();
    assert_eq!(s.positions.len(), 8); // CrI3 has 8 atoms
    assert_eq!(s.cell.len(), 3);
}

#[test]
fn parse_from_poscar_with_magnetic_indices() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: Some(0.01),
        magnetic_indices: Some(vec![0, 1]), // only Cr atoms
    };
    let parsed = sc.parse().unwrap();
    assert_eq!(parsed.positions.len(), 2);
}

#[test]
fn parse_from_nonexistent_file_errors() {
    let sc = StructureConf {
        file: Some("nonexistent_file.vasp".into()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("failed to open POSCAR"));
}

#[test]
fn parse_magnetic_indices_out_of_range_errors() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: Some(vec![100]),
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("out of range"));
}

#[test]
fn parse_magnetic_indices_duplicate_errors() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: Some(vec![0, 0]),
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("duplicate"));
}

#[test]
fn validate_file_mode_without_magnetic_indices_positions_mismatch() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    // CrI3 has 8 atoms, sublattices=5 -> mismatch
    let err = sc.validate(5).unwrap_err().to_string();
    assert!(err.contains("number of atoms in file"));
}

#[test]
fn validate_file_mode_ok() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    // CrI3 has 8 atoms
    assert!(sc.validate(8).is_ok());
}

#[test]
fn parse_cell_and_positions_ok() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        tolerance: Some(0.01),
        magnetic_indices: None,
    };
    let s = sc.parse().unwrap();
    assert_eq!(s.cell.len(), 3);
    assert_eq!(s.positions.len(), 2);
    assert_eq!(s.tolerance, Some(0.01));
}

#[test]
fn parse_cell_and_positions_without_tolerance() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let s = sc.parse().unwrap();
    assert_eq!(s.tolerance, None);
}

#[test]
fn parse_missing_both_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("structure must be provided"));
}

#[test]
fn parse_file_and_cell_both_errors() {
    let sc = StructureConf {
        file: Some("test.poscar".into()),
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("both `file` and `cell`"));
}

#[test]
fn parse_file_and_positions_both_errors() {
    let sc = StructureConf {
        file: Some("test.poscar".into()),
        format: None,
        cell: None,
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("both `file` and `positions`"));
}

#[test]
fn parse_cell_without_positions_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("missing `positions`"));
}

#[test]
fn parse_positions_without_cell_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: None,
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("missing `cell`"));
}

#[test]
fn parse_magnetic_indices_with_cell_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: Some(vec![0]),
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("`magnetic_indices` is only valid"));
}

#[test]
fn parse_format_with_cell_errors() {
    let sc = StructureConf {
        file: None,
        format: Some("vasp".into()),
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("`format` is only valid"));
}

#[test]
fn validate_cell_mode_positions_mismatch() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    // positions=1, sublattices=2 → mismatch
    let err = sc.validate(2).unwrap_err().to_string();
    assert!(err.contains("number of provided positions"));
}

#[test]
fn validate_cell_mode_ok() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        tolerance: None,
        magnetic_indices: None,
    };
    assert!(sc.validate(2).is_ok());
}

#[test]
fn validate_magnetic_indices_length_mismatch() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: Some(vec![0, 1]),
    };
    // magnetic_indices=2, sublattices=3 → mismatch
    let err = sc.validate(3).unwrap_err().to_string();
    assert!(err.contains("number of magnetic_indices"));
}

#[test]
fn display() {
    let sc = StructureConf {
        file: Some("test.vasp".into()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let s = format!("{sc}");
    assert!(s.contains("test.vasp"));
}
