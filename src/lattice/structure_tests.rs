use super::{is_vasp_format, load_from_file, parse_poscar_from_str};
use std::fs;

fn parse(content: &str) -> super::Structure {
    parse_poscar_from_str(content.trim(), "test").expect("failed to parse POSCAR")
}

// ─── is_vasp_format ────────────────────────────────────────────────

#[test]
fn vasp_format_explicit() {
    assert!(is_vasp_format(&Some("poscar".into()), "anyfile"));
    assert!(is_vasp_format(&Some("POSCAR".into()), "anyfile"));
    assert!(is_vasp_format(&Some("vasp".into()), "anyfile"));
    assert!(is_vasp_format(&Some("VASP".into()), "anyfile"));
    assert!(is_vasp_format(&Some("contcar".into()), "anyfile"));
    assert!(is_vasp_format(&Some("CONTCAR".into()), "anyfile"));
}

#[test]
fn not_vasp_format_explicit() {
    assert!(!is_vasp_format(&Some("xyz".into()), "anyfile"));
    assert!(!is_vasp_format(&Some("cif".into()), "anyfile"));
    assert!(!is_vasp_format(&Some("pdb".into()), "anyfile"));
}

#[test]
fn vasp_format_detect_vasp_extension() {
    assert!(is_vasp_format(&None, "structure.vasp"));
    assert!(is_vasp_format(&None, "STRUCTURE.VASP"));
    assert!(is_vasp_format(&None, "path/to/structure.VASP"));
}

#[test]
fn vasp_format_detect_poscar_contcar() {
    assert!(is_vasp_format(&None, "POSCAR"));
    assert!(is_vasp_format(&None, "CONTCAR"));
    assert!(is_vasp_format(&None, "path/to/POSCAR"));
    assert!(is_vasp_format(&None, "path/to/CONTCAR"));
}

#[test]
fn vasp_format_not_detected() {
    assert!(!is_vasp_format(&None, "structure.xyz"));
    assert!(!is_vasp_format(&None, "structure.cif"));
    assert!(!is_vasp_format(&None, "structure.pdb"));
    assert!(!is_vasp_format(&None, "data.txt"));
}

#[test]
fn explicit_format_takes_precedence() {
    // format explicitly set to xyz even though filename looks like POSCAR
    assert!(!is_vasp_format(&Some("xyz".into()), "POSCAR"));
    assert!(!is_vasp_format(&Some("cif".into()), "structure.vasp"));
}

#[test]
fn load_from_file_with_explicit_xyz_format() {
    let source = std::fs::read_to_string("examples/heisenberg_2d_cri3_poscar/POSCAR").unwrap();
    let source = parse_poscar_from_str(&source, "examples/heisenberg_2d_cri3_poscar/POSCAR")
        .expect("failed to parse source POSCAR");
    let path = std::env::temp_dir().join(format!("spinmc-cri3-poscar-{}.xyz", std::process::id()));
    fs::write(&path, extended_xyz_from_structure(&source)).unwrap();

    let structure = load_from_file(path.to_str().unwrap(), Some("XYZ".into())).unwrap();

    assert_positions_close(&structure.positions, &source.positions);
    assert_cell_close(structure.cell, source.cell);
    assert!(structure.frame.is_some());

    fs::remove_file(path).unwrap();
}

#[test]
fn load_from_file_infers_xyz_format_from_extension() {
    let source = std::fs::read_to_string("examples/heisenberg_2d_cri3_poscar/POSCAR").unwrap();
    let source = parse_poscar_from_str(&source, "examples/heisenberg_2d_cri3_poscar/POSCAR")
        .expect("failed to parse source POSCAR");
    let path = std::env::temp_dir().join(format!(
        "spinmc-cri3-poscar-infer-{}.xyz",
        std::process::id()
    ));
    fs::write(&path, extended_xyz_from_structure(&source)).unwrap();

    let structure = load_from_file(path.to_str().unwrap(), None).unwrap();

    assert_positions_close(&structure.positions, &source.positions);
    assert_cell_close(structure.cell, source.cell);
    assert!(structure.frame.is_some());

    fs::remove_file(path).unwrap();
}

fn extended_xyz_from_structure(structure: &super::Structure) -> String {
    let cell = structure.cell;
    let mut xyz = format!(
        "{}\nLattice=\"{} {} {} {} {} {} {} {} {}\" Properties=species:S:1:pos:R:3\n",
        structure.positions.len(),
        cell[0][0],
        cell[1][0],
        cell[2][0],
        cell[0][1],
        cell[1][1],
        cell[2][1],
        cell[0][2],
        cell[1][2],
        cell[2][2],
    );
    for position in &structure.positions {
        xyz.push_str(&format!(
            "X {:.16} {:.16} {:.16}\n",
            position[0], position[1], position[2]
        ));
    }
    xyz
}

fn assert_positions_close(left: &[[f64; 3]], right: &[[f64; 3]]) {
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        assert_vector_close(*left, *right);
    }
}

fn assert_cell_close(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) {
    for (left, right) in left.iter().zip(right.iter()) {
        assert_vector_close(*left, *right);
    }
}

fn assert_vector_close(left: [f64; 3], right: [f64; 3]) {
    for (left, right) in left.iter().zip(right) {
        assert!((left - right).abs() < 1e-10, "{left} != {right}");
    }
}

// ─── parse_poscar_from_str: happy paths ────────────────────────────

#[test]
fn basic_direct() {
    let s = parse(
        "\
cubic BN
1.0
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
    );
    assert_eq!(s.cell, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    assert_eq!(s.positions.len(), 1);
    assert!((s.positions[0][0] - 0.5).abs() < 1e-10);
    assert!((s.positions[0][1] - 0.5).abs() < 1e-10);
    assert!((s.positions[0][2] - 0.5).abs() < 1e-10);
}

#[test]
fn direct_to_cartesian() {
    let s = parse(
        "\
test
1.0
  2 0 0
  0 2 0
  0 0 2
  1
direct
  0.25 0.25 0.25
",
    );
    assert!((s.positions[0][0] - 0.5).abs() < 1e-10);
    assert!((s.positions[0][1] - 0.5).abs() < 1e-10);
    assert!((s.positions[0][2] - 0.5).abs() < 1e-10);
}

#[test]
fn cartesian_coords() {
    let s = parse(
        "\
test
2.0
  1 0 0
  0 1 0
  0 0 1
  1
Cartesian
  1.0 2.0 3.0
",
    );
    assert!((s.positions[0][0] - 2.0).abs() < 1e-10);
    assert!((s.positions[0][1] - 4.0).abs() < 1e-10);
    assert!((s.positions[0][2] - 6.0).abs() < 1e-10);
}

#[test]
fn no_element_names() {
    let s = parse(
        "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  2
Direct
  0.0 0.0 0.0
  0.5 0.5 0.5
",
    );
    assert_eq!(s.positions.len(), 2);
    assert!((s.positions[0][0] - 0.0).abs() < 1e-10);
    assert!((s.positions[1][0] - 0.5).abs() < 1e-10);
}

#[test]
fn with_element_names() {
    let s = parse(
        "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  B N
  1 1
Direct
  0.0 0.0 0.0
  0.25 0.25 0.25
",
    );
    assert_eq!(s.positions.len(), 2);
    assert!((s.positions[0][0] - 0.0).abs() < 1e-10);
    assert!((s.positions[1][0] - 0.25).abs() < 1e-10);
}

#[test]
fn scale_factor() {
    let s = parse(
        "\
test
3.57
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.0 0.0 0.0
",
    );
    assert!((s.cell[0][0] - 3.57).abs() < 1e-10);
}

#[test]
fn three_scale_factors() {
    let s = parse(
        "\
test
2.0 3.0 4.0
  1 0 0
  0 1 0
  0 0 1
  1
Cartesian
  1.0 1.0 1.0
",
    );
    // lattice: each row component * respective scale
    assert!((s.cell[0][0] - 2.0).abs() < 1e-10);
    assert!((s.cell[1][1] - 3.0).abs() < 1e-10);
    assert!((s.cell[2][2] - 4.0).abs() < 1e-10);
    // Cartesian pos: x*2, y*3, z*4
    assert!((s.positions[0][0] - 2.0).abs() < 1e-10);
    assert!((s.positions[0][1] - 3.0).abs() < 1e-10);
    assert!((s.positions[0][2] - 4.0).abs() < 1e-10);
}

#[test]
fn selective_dynamics_with_direct() {
    let s = parse(
        "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  1
Selective dynamics
Direct
  0.5 0.5 0.5 T T F
",
    );
    assert_eq!(s.positions.len(), 1);
}

#[test]
fn selective_dynamics_with_cartesian() {
    let s = parse(
        "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  1
s
c
  1.0 2.0 3.0 F F F
",
    );
    assert_eq!(s.positions.len(), 1);
}

#[test]
fn coord_mode_direct_variants() {
    for mode in &["Direct", "direct", "D", "d", "Xyz", ""] {
        let content = format!(
            "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  1
{mode}
  0.5 0.5 0.5
"
        );
        let s = parse(&content);
        assert!(
            (s.positions[0][0] - 0.5).abs() < 1e-10,
            "failed for mode={mode:?}"
        );
    }
}

#[test]
fn coord_mode_cartesian_variants() {
    for mode in &["Cartesian", "cartesian", "C", "c", "K", "k"] {
        let content = format!(
            "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  1
{mode}
  2.0 3.0 4.0
"
        );
        let s = parse(&content);
        assert!(
            (s.positions[0][0] - 2.0).abs() < 1e-10,
            "failed for mode={mode:?}"
        );
    }
}

#[test]
fn multiple_atoms_with_velocities_tail() {
    let s = parse(
        "\
test
1.0
  1 0 0
  0 1 0
  0 0 1
  2
Direct
  0.0 0.0 0.0
  0.5 0.5 0.5

  0.1 0.2 0.3
  0.4 0.5 0.6
",
    );
    assert_eq!(s.positions.len(), 2);
}

// ─── parse_poscar_from_str: error paths ────────────────────────────

#[test]
fn error_too_few_lines() {
    let err = parse_poscar_from_str("too\nshort", "test").unwrap_err();
    assert!(err.to_string().contains("too short"));
}

#[test]
fn error_invalid_scale_not_a_number() {
    let err = parse_poscar_from_str(
        "\
comment
abc
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("scale"));
}

#[test]
fn error_scale_too_many_numbers() {
    let err = parse_poscar_from_str(
        "\
comment
1.0 2.0 3.0 4.0
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("scale"));
}

#[test]
fn error_negative_scale_not_supported() {
    let err = parse_poscar_from_str(
        "\
comment
-1.0
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("must be positive"));
}

#[test]
fn error_zero_scale_not_supported() {
    let err = parse_poscar_from_str(
        "\
comment
1.0 0.0 1.0
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("must be positive"));
}

#[test]
fn error_lattice_vector_too_few_columns() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("lattice vector"));
}

#[test]
fn error_lattice_vector_non_numeric() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0 a
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("lattice vector"));
}

#[test]
fn error_coordinate_non_numeric() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 abc 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("coordinate"));
}

#[test]
fn error_coordinate_too_few_columns() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0 0
  0 1 0
  0 0 1
  1
Direct
  0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("invalid coordinate line"));
}

#[test]
fn error_insufficient_coordinate_lines() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0 0
  0 1 0
  0 0 1
  2
Direct
  0.5 0.5 0.5
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("insufficient"));
}

#[test]
fn error_missing_coordinate_data() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0 0
  0 1 0
  0 0 1
  1
",
        "test",
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("missing coordinate data")
            || err.to_string().contains("too short")
    );
}

#[test]
fn error_missing_coordinate_data_after_element_counts() {
    let err = parse_poscar_from_str(
        "\
comment
1.0
  1 0 0
  0 1 0
  0 0 1
H
1
",
        "test",
    )
    .unwrap_err();
    assert!(err.to_string().contains("missing coordinate data"));
}
