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
        magnetic_indices: Some(vec![6, 7]), // only Cr atoms
    };
    let parsed = sc.parse().unwrap();
    assert_eq!(parsed.positions.len(), 2);
    let full_structure = parsed.full_structure.as_ref().unwrap();
    assert_eq!(full_structure.atoms.len(), 8);
    assert_eq!(full_structure.atoms[0].element, "I");
    assert_eq!(full_structure.atoms[7].element, "Cr");
    assert_eq!(parsed.positions[0], full_structure.atoms[6].position);
}

#[test]
fn parse_from_poscar_without_magnetic_indices_keeps_full_structure() {
    let sc = StructureConf {
        file: Some(poscar_path()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let s = sc.parse().unwrap();
    let full_structure = s.full_structure.as_ref().unwrap();
    assert_eq!(full_structure.cell, s.cell);
    assert_eq!(full_structure.atoms.len(), s.positions.len());
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
    assert!(err.contains("number of atoms"));
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
fn parse_cell_and_positions_uses_placeholder_elements_in_full_structure() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let s = sc.parse().unwrap();
    let full_structure = s.full_structure.as_ref().unwrap();
    assert_eq!(full_structure.cell, s.cell);
    assert_eq!(full_structure.atoms.len(), 2);
    assert!(full_structure.atoms.iter().all(|atom| atom.element == "X"));
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
fn parse_magnetic_indices_with_cell_filters_positions() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        tolerance: None,
        magnetic_indices: Some(vec![0]),
    };
    let s = sc.parse().unwrap();
    assert_eq!(s.positions.len(), 1);
    assert_eq!(s.full_structure.as_ref().unwrap().atoms.len(), 2);
}

#[test]
fn parse_magnetic_indices_with_cell_out_of_range_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: Some(vec![5]),
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("out of range"));
}

#[test]
fn parse_magnetic_indices_with_cell_duplicate_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        tolerance: None,
        magnetic_indices: Some(vec![0, 0]),
    };
    let err = sc.parse().unwrap_err().to_string();
    assert!(err.contains("duplicate"));
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
fn validate_cell_mode_with_magnetic_indices_allows_extra_positions() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]),
        tolerance: None,
        magnetic_indices: Some(vec![0, 2]),
    };
    assert!(sc.validate(2).is_ok());
}

#[test]
fn validate_empty_file_path_errors() {
    let sc = StructureConf {
        file: Some("  ".into()),
        format: None,
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("file path"));
    assert!(err.contains("empty"));
}

#[test]
fn validate_empty_format_errors() {
    let sc = StructureConf {
        file: Some("POSCAR".into()),
        format: Some(" ".into()),
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("format"));
    assert!(err.contains("empty"));
}

#[test]
fn validate_unsupported_format_errors() {
    let sc = StructureConf {
        file: Some("structure.xyz".into()),
        format: Some("xyz".into()),
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("unsupported"));
    assert!(err.contains("xyz"));
}

#[test]
fn validate_format_with_cell_errors() {
    let sc = StructureConf {
        file: None,
        format: Some("vasp".into()),
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("`format` is only valid"));
}

#[test]
fn validate_non_finite_tolerance_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: Some(f64::NAN),
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("tolerance"));
    assert!(err.contains("finite"));
}

#[test]
fn validate_negative_tolerance_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: Some(-0.1),
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("tolerance"));
    assert!(err.contains("non-negative"));
}

#[test]
fn validate_non_finite_cell_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, f64::NAN, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("cell[0][1]"));
    assert!(err.contains("finite"));
}

#[test]
fn validate_zero_volume_cell_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("cell vectors"));
    assert!(err.contains("non-zero volume"));
}

#[test]
fn validate_empty_positions_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("positions"));
    assert!(err.contains("at least one"));
}

#[test]
fn validate_non_finite_positions_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, f64::INFINITY, 0.0]]),
        tolerance: None,
        magnetic_indices: None,
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("positions[0][1]"));
    assert!(err.contains("finite"));
}

#[test]
fn validate_magnetic_indices_duplicate_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        tolerance: None,
        magnetic_indices: Some(vec![0, 0]),
    };
    let err = sc.validate(2).unwrap_err().to_string();
    assert!(err.contains("duplicate"));
}

#[test]
fn validate_magnetic_indices_out_of_range_errors() {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0]]),
        tolerance: None,
        magnetic_indices: Some(vec![1]),
    };
    let err = sc.validate(1).unwrap_err().to_string();
    assert!(err.contains("out of range"));
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

#[test]
fn display_file_mode_includes_format_and_parse_error() {
    let sc = StructureConf {
        file: Some("missing.poscar".into()),
        format: Some("POSCAR".into()),
        cell: None,
        positions: None,
        tolerance: None,
        magnetic_indices: None,
    };

    let s = format!("{sc}");

    assert!(s.contains("Mode: from file"));
    assert!(s.contains("File: missing.poscar"));
    assert!(s.contains("Format: POSCAR"));
    assert!(s.contains("[Failed to parse:"));
}
