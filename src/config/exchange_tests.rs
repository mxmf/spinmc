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
        strength: 1.5,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
    assert_eq!(result[0].strength, 1.5);
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
        strength: 2.0,
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 3.0,
    };
    let result = e.parse(&Some(structure), [true, true, false]).unwrap();
    assert!(!result.is_empty());
    assert_eq!(result[0].strength, 3.0);
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
        strength: 3.0,
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
        strength: 3.0,
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
        strength: 3.0,
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 1.0,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("exactly one"));
}

#[test]
fn validate_from_sublattice_out_of_range() {
    let e = Exchange {
        from_sublattice: Some(2),
        to_sublattice: Some(0),
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 1.0,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("minimum"));
    assert!(err.contains("non-negative"));
}

#[test]
fn validate_distance_range_max_less_than_min() {
    let e = Exchange {
        from_sublattice: Some(0),
        to_sublattice: Some(0),
        offsets: None,
        neighbor_order: None,
        distance_range: Some([5.0, 3.0]),
        strength: 1.0,
    };
    let err = e.validate(1).unwrap_err().to_string();
    assert!(err.contains("maximum"));
    assert!(err.contains("minimum"));
}

#[test]
fn parse_offsets_without_from_to_errors() {
    let e = Exchange {
        from_sublattice: None,
        to_sublattice: None,
        offsets: Some(vec![[1, 0, 0]]),
        neighbor_order: None,
        distance_range: None,
        strength: 1.0,
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
        strength: 2.0,
    };
    let result = e.parse(&None, [true, true, true]).unwrap();
    assert_eq!(result.len(), 2);
    assert_eq!(result[0].offset, [1, 0, 0]);
    assert_eq!(result[0].strength, 2.0);
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
        strength: 1.0,
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
        strength: 1.0,
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
        strength: 0.5,
    };
    let s = format!("{p}");
    assert!(s.contains("0.5"));
    assert!(s.contains("1"));
}
