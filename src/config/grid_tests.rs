use super::*;

#[test]
fn validate_ok() {
    let g = Grid {
        dimensions: [2, 2, 1],
        sublattices: 2,
        spin_magnitudes: vec![1.0, 1.0],
        periodic_boundary: [true, true, true],
    };
    assert!(g.validate().is_ok());
}

#[test]
fn validate_mismatched_length() {
    let g = Grid {
        dimensions: [2, 2, 1],
        sublattices: 3,
        spin_magnitudes: vec![1.0, 1.0],
        periodic_boundary: [true, true, true],
    };
    let err = g.validate().unwrap_err().to_string();
    assert!(err.contains("spin_magnitude length"));
    assert!(err.contains("sublattices"));
}

#[test]
fn display() {
    let g = Grid {
        dimensions: [2, 2, 1],
        sublattices: 1,
        spin_magnitudes: vec![1.0],
        periodic_boundary: [true, true, false],
    };
    let s = format!("{g}");
    assert!(s.contains("Grid"));
    assert!(s.contains("2, 2, 1"));
}
