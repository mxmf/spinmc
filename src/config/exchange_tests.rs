use super::*;
use crate::config::structure::StructureConf;
use crate::lattice::Structure;

fn make_square_structure() -> Structure {
    let sc = StructureConf {
        file: None,
        format: None,
        cell: Some([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]]),
        positions: Some(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        tolerance: Some(0.0001),
        magnetic_indices: None,
    };
    sc.parse().unwrap()
}

#[test]
fn parse_neighbor_order_with_structure() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.5),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
    assert_eq!(result[0].j_scalar, Some(1.5));
}

#[test]
fn parse_neighbor_order_from_only() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: None,
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(2.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
}

#[test]
fn parse_neighbor_order_all() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
}

#[test]
fn parse_neighbor_order_to_without_from_errors() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e
        .parse(&Some(structure), [true, true, false])
        .unwrap_err()
        .to_string();
    assert!(err.contains("from_sublattice"));
}

#[test]
fn parse_distance_range_with_structure() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.5, 1.5]),
        j_scalar: Some(3.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
    assert_eq!(result[0].j_scalar, Some(3.0));
}

#[test]
fn parse_distance_range_from_only() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: None,
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.5, 1.5]),
        j_scalar: Some(3.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
    assert!(result.iter().all(|exchange| exchange.from_sub == 0));
}

#[test]
fn parse_distance_range_all_sublattices() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.5, 1.5]),
        j_scalar: Some(3.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
    assert!(result.iter().any(|exchange| exchange.from_sub == 0));
    assert!(result.iter().any(|exchange| exchange.from_sub == 1));
}

#[test]
fn parse_distance_range_to_without_from_errors() {
    let structure = make_square_structure();
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.5, 1.5]),
        j_scalar: Some(3.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e
        .parse(&Some(structure), [true, true, false])
        .unwrap_err()
        .to_string();
    assert!(err.contains("from_sublattice"));
}

#[test]
fn parse_without_selection_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: None,
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.parse(&None, [true, true, true]).unwrap_err().to_string();
    assert!(err.contains("must specify one"));
}

#[test]
fn parse_with_multiple_selection_modes_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.parse(&None, [true, true, true]).unwrap_err().to_string();
    assert!(err.contains("only one"));
}

#[test]
fn validate_offsets_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    assert!(e.validate(2).is_ok());
}

#[test]
fn validate_neighbor_order_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    assert!(e.validate(2).is_ok());
}

#[test]
fn validate_distance_range_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.0, 5.0]),
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    assert!(e.validate(2).is_ok());
}

#[test]
fn validate_none_specified_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: None,
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("exactly one"));
}

#[test]
fn validate_multiple_specified_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("exactly one"));
}

#[test]
fn validate_non_finite_j_scalar_errors() {
    for strength in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let e = Exchange {
            from_sublattice: Some(0),
            to_sublattice: Some(0),
            offsets: Some(vec![[1, 0, 0]]),
            neighbor_order: None,
            distance_range: None,
            j_scalar: Some(strength),
            j_diagonal: None,
            j_tensor: None,
            dm: None,
        };
        let err = e.validate(1).unwrap_err().to_string();
        assert!(err.contains("j_scalar"));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_zero_strength_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(0.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    assert!(e.validate(1).is_ok());
}

#[test]
fn validate_dm_only_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: None,
        j_diagonal: None,
        j_tensor: None,
        dm: Some([0.0, 0.0, 1.0]),
    };
    assert!(e.validate(1).is_ok());
}

#[test]
fn validate_dm_with_j_scalar_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: Some([0.0, 0.0, 1.0]),
    };
    assert!(e.validate(1).is_ok());
}

#[test]
fn validate_dm_with_j_tensor_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: None,
        j_diagonal: None,
        j_tensor: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        dm: Some([0.0, 0.0, 1.0]),
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("j_tensor"));
    assert!(err.contains("dm"));
}

#[test]
fn validate_non_finite_dm_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: None,
        j_diagonal: None,
        j_tensor: None,
        dm: Some([0.0, f64::NAN, 1.0]),
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("dm[1]"));
    assert!(err.contains("finite"));
}

#[test]
fn validate_empty_offsets_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: Some(vec![]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("offsets"));
    assert!(err.contains("at least one"));
}

#[test]
fn validate_offsets_without_from_to_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("from_sublattice"));
    assert!(err.contains("to_sublattice"));
}

#[test]
fn validate_neighbor_order_zero_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: None,
        neighbor_order: Some(0),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("neighbor_order"));
    assert!(err.contains("greater than zero"));
}

#[test]
fn validate_neighbor_order_to_without_from_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: Some(0),
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("from_sublattice"));
    assert!(err.contains("neighbor_order"));
}

#[test]
fn validate_from_sublattice_out_of_range() {
    let e = Exchange {
        from_sublattice: Some(2),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(2).unwrap_err().to_string();
    assert!(err.contains("from_sublattice"));
    assert!(err.contains("out of range"));
}

#[test]
fn validate_to_sublattice_out_of_range() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(2),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(2).unwrap_err().to_string();
    assert!(err.contains("to_sublattice"));
    assert!(err.contains("out of range"));
}

#[test]
fn validate_distance_range_negative_min() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([-1.0, 5.0]),
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("minimum"));
    assert!(err.contains("non-negative"));
}

#[test]
fn validate_distance_range_non_finite_errors() {
    for (distance_range, expected) in [
        ([f64::NAN, 5.0], "minimum"),
        ([0.0, f64::INFINITY], "maximum"),
    ] {
        let e = Exchange {
            from_sublattice: Some(0),
            to_sublattice: Some(0),
            offsets: None,
            neighbor_order: None,
            distance_range: Some(distance_range),
            j_scalar: Some(1.0),
            j_diagonal: None,
            j_tensor: None,
            dm: None,
        };
        let err = e.validate(1).unwrap_err().to_string();
        assert!(err.contains(expected));
        assert!(err.contains("finite"));
    }
}

#[test]
fn validate_distance_range_max_less_than_min() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([5.0, 3.0]),
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("maximum"));
    assert!(err.contains("minimum"));
}

#[test]
fn validate_distance_range_to_without_from_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: Some(0),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.0, 5.0]),
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("from_sublattice"));
    assert!(err.contains("distance_range"));
}

#[test]
fn parse_offsets_without_from_to_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.parse(&None, [true, true, true]).unwrap_err().to_string();
    assert!(err.contains("from_sublattice"));
    assert!(err.contains("to_sublattice"));
}

#[test]
fn parse_offsets_ok() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: Some(vec![[1, 0, 0], [0, 1, 0]]),
        neighbor_order: None,
        distance_range: None,
        j_scalar: Some(2.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let result = e.parse(&None, [true, true, true]).unwrap();
    assert_eq!(result.len(), 2);
    assert_eq!(result[0].offset, [1, 0, 0]);
    assert_eq!(result[0].j_scalar, Some(2.0));
    assert_eq!(result[1].offset, [0, 1, 0]);
}

#[test]
fn parse_neighbor_order_without_structure_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: Some(1),
        distance_range: None,
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.parse(&None, [true, true, true]).unwrap_err().to_string();
    assert!(err.contains("structure"));
}

#[test]
fn parse_distance_range_without_structure_errors() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(1),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([0.0, 5.0]),
        j_scalar: Some(1.0),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let err = e.parse(&None, [true, true, true]).unwrap_err().to_string();
    assert!(err.contains("structure"));
}

#[test]
fn display_parsed_exchange() {
    let p = ParsedExchange {
        from_sub: 0,
        to_sub: 1,
        offset: [1, 0, -1],
        j_scalar: Some(0.5),
        j_diagonal: None,
        j_tensor: None,
        dm: None,
    };
    let s = format!("{p}");
    assert!(s.contains("0.5"));
    assert!(s.contains("1"));
}
