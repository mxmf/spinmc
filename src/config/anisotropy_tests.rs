use super::*;

#[test]
fn validate_ok() {
    let a = Anisotropy {
        axis: vec![[0.0, 0.0, 1.0]],
        strength: vec![1.0],
    };
    assert!(a.validate(1).is_ok());
}

#[test]
fn validate_axis_strength_length_mismatch() {
    let a = Anisotropy {
        axis: vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        strength: vec![1.0],
    };
    let err = a.validate(2).unwrap_err().to_string();
    assert!(err.contains("axis and strength arrays"));
}

#[test]
fn validate_strength_length_not_equal_sublattices() {
    let a = Anisotropy {
        axis: vec![[0.0, 0.0, 1.0]],
        strength: vec![1.0],
    };
    let err = a.validate(2).unwrap_err().to_string();
    assert!(err.contains("strength arrays"));
    assert!(err.contains("2"));
}

#[test]
fn validate_zero_length_axis_errors() {
    let a = Anisotropy {
        axis: vec![[0.0, 0.0, 0.0]],
        strength: vec![1.0],
    };
    let err = a.validate(1).unwrap_err().to_string();
    assert!(err.contains("axis[0]"));
    assert!(err.contains("non-zero length"));
}

#[test]
fn validate_non_finite_axis_errors() {
    for component in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let a = Anisotropy {
            axis: vec![[0.0, component, 1.0]],
            strength: vec![1.0],
        };
        let err = a.validate(1).unwrap_err().to_string();
        assert!(err.contains("axis[0][1]"));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_non_finite_strength_errors() {
    for strength in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let a = Anisotropy {
            axis: vec![[0.0, 0.0, 1.0]],
            strength: vec![strength],
        };
        let err = a.validate(1).unwrap_err().to_string();
        assert!(err.contains("strength[0]"));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_zero_strength_ok() {
    let a = Anisotropy {
        axis: vec![[0.0, 0.0, 1.0]],
        strength: vec![0.0],
    };
    assert!(a.validate(1).is_ok());
}

#[test]
fn parse_normalizes_axis() {
    let a = Anisotropy {
        axis: vec![[2.0, 0.0, 0.0]],
        strength: vec![3.0],
    };
    let result = a.parse().unwrap();
    assert_eq!(result.len(), 1);
    assert!((result[0].axis[0] - 1.0).abs() < 1e-10);
    assert!((result[0].axis[1] - 0.0).abs() < 1e-10);
    assert!((result[0].axis[2] - 0.0).abs() < 1e-10);
    assert_eq!(result[0].strength, 3.0);
}

#[test]
fn parse_zero_length_axis_errors() {
    let a = Anisotropy {
        axis: vec![[0.0, 0.0, 0.0]],
        strength: vec![1.0],
    };
    let err = a.parse().unwrap_err().to_string();
    assert!(err.contains("zero length"));
}

#[test]
fn parse_multiple() {
    let a = Anisotropy {
        axis: vec![[0.0, 1.0, 0.0], [0.0, 0.0, 3.0]],
        strength: vec![2.0, 5.0],
    };
    let result = a.parse().unwrap();
    assert_eq!(result.len(), 2);
    // First: axis normalized (0,1,0), strength=2
    assert!((result[0].axis[1] - 1.0).abs() < 1e-10);
    assert_eq!(result[0].strength, 2.0);
    // Second: axis normalized (0,0,1), strength=5
    assert!((result[1].axis[2] - 1.0).abs() < 1e-10);
    assert_eq!(result[1].strength, 5.0);
}

#[test]
fn display_parsed_anisotropy() {
    let p = ParsedAnisotropy {
        axis: [0.0, 0.0, 1.0],
        strength: 2.5,
    };
    let s = format!("{p}");
    assert!(s.contains("2.5"));
}
