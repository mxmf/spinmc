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
fn validate_zero_dimension_errors() {
    let g = Grid {
        dimensions: [2, 0, 1],
        sublattices: 1,
        spin_magnitudes: vec![1.0],
        periodic_boundary: [true, true, true],
    };
    let err = g.validate().unwrap_err().to_string();
    assert!(err.contains("dimensions[1]"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_zero_sublattices_errors() {
    let g = Grid {
        dimensions: [2, 2, 1],
        sublattices: 0,
        spin_magnitudes: vec![],
        periodic_boundary: [true, true, true],
    };
    let err = g.validate().unwrap_err().to_string();
    assert!(err.contains("sublattices"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_non_positive_spin_magnitude_errors() {
    let g = Grid {
        dimensions: [2, 2, 1],
        sublattices: 2,
        spin_magnitudes: vec![1.0, 0.0],
        periodic_boundary: [true, true, true],
    };
    let err = g.validate().unwrap_err().to_string();
    assert!(err.contains("spin_magnitudes[1]"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_non_finite_spin_magnitude_errors() {
    for magnitude in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let g = Grid {
            dimensions: [2, 2, 1],
            sublattices: 1,
            spin_magnitudes: vec![magnitude],
            periodic_boundary: [true, true, true],
        };
        let err = g.validate().unwrap_err().to_string();
        assert!(err.contains("spin_magnitudes[0]"));
        assert!(err.contains("finite"));
    }
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
