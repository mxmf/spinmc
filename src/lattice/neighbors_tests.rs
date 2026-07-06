use super::Atoms;
use super::Vector3Ext;

fn direct_to_cart(cell: [[f64; 3]; 3], direct: [f64; 3]) -> [f64; 3] {
    [
        direct[0] * cell[0][0] + direct[1] * cell[1][0] + direct[2] * cell[2][0],
        direct[0] * cell[0][1] + direct[1] * cell[1][1] + direct[2] * cell[2][1],
        direct[0] * cell[0][2] + direct[1] * cell[1][2] + direct[2] * cell[2][2],
    ]
}

// --- Vector3Ext ---

#[test]
fn vector3ext_scale() {
    let v = [1.0, 2.0, 3.0];
    assert_eq!(v.scale(2.0), [2.0, 4.0, 6.0]);
    assert_eq!(v.scale(0.0), [0.0, 0.0, 0.0]);
}

#[test]
fn vector3ext_dot() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 5.0, 6.0];
    assert_eq!(a.dot(&b), 32.0);
    assert_eq!(a.dot(&[0.0, 1.0, 0.0]), 2.0);
}

#[test]
fn vector3ext_norm() {
    assert!(([3.0, 4.0, 0.0].norm() - 5.0).abs() < 1e-10);
    assert_eq!([0.0, 0.0, 0.0].norm(), 0.0);
}

#[test]
fn vector3ext_add() {
    assert_eq!([1.0, 2.0, 3.0].add(&[4.0, 5.0, 6.0]), [5.0, 7.0, 9.0]);
}

#[test]
fn vector3ext_sub() {
    assert_eq!([5.0, 7.0, 9.0].sub(&[4.0, 5.0, 6.0]), [1.0, 2.0, 3.0]);
}

// --- calc_distance variants ---

fn two_site_unit_atoms() -> Atoms {
    Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0], [0.5, 0., 0.]],
        pbc: [true, true, true],
        tolerance: 0.0001,
    }
}

#[test]
fn test_calc_distance_from_to_max_n() {
    let atoms = two_site_unit_atoms();
    let distances = atoms.calc_distance_from_to(0, 1, 1);
    assert_eq!(distances.len(), 27);
    let nearest: Vec<_> = distances
        .iter()
        .filter(|distance| (distance.distance - 0.5).abs() < 1e-10)
        .map(|distance| distance.neighbor.offset)
        .collect();
    assert_eq!(nearest, vec![[-1, 0, 0], [0, 0, 0]]);
}

#[test]
fn test_calc_distance_from_max_n() {
    let atoms = two_site_unit_atoms();
    let distances = atoms.calc_distance_from(0, 1);
    assert_eq!(distances.len(), 53);
    assert!(distances.iter().all(|distance| distance.neighbor.from == 0));
    assert_eq!(
        distances
            .iter()
            .filter(|distance| distance.neighbor.to == 0)
            .count(),
        26
    );
    assert_eq!(
        distances
            .iter()
            .filter(|distance| distance.neighbor.to == 1)
            .count(),
        27
    );
}

#[test]
fn test_calc_distance_from_to_self_excluded() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, true],
        tolerance: 0.0001,
    };
    let distances = atoms.calc_distance_from_to(0, 0, 2);
    // Should not include the zero-distance self pair
    assert_eq!(distances.len(), 124);
    assert!(distances.iter().all(|d| d.distance > 0.0));
    assert!(
        distances
            .iter()
            .all(|distance| distance.neighbor.offset != [0, 0, 0])
    );
}

// --- cell_offset ---

#[test]
fn test_cell_offset() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.5, 2.0, 0.0], [0.0, 0.25, 3.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, true],
        tolerance: 0.0001,
    };
    let offset = atoms.cell_offset([1, 2, 3]);
    assert_eq!(offset, [2.0, 4.75, 9.0]);
}

// --- reciprocal_cell zero volume ---

#[test]
fn test_reciprocal_cell_zero_volume() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, true],
        tolerance: 0.0001,
    };
    // Cell vectors (1,0,0) and (1,0,0) are linearly dependent → volume=0
    let bounds = atoms.offset_bounds_in_radius([0.0, 0.0, 0.0], 5.0);
    assert_eq!(bounds, [[0, 0]; 3]);
}

// --- find_neighbors variants ---

#[test]
fn test_get_neighbor_from() {
    let atoms = two_site_unit_atoms();
    let distances = atoms.calc_distance_from(0, 1);
    let neighbors = atoms.get_neighbor_from(distances, 1);
    let offsets: Vec<_> = neighbors.iter().map(|neighbor| neighbor.offset).collect();
    assert_eq!(offsets, vec![[-1, 0, 0], [0, 0, 0]]);
    assert!(neighbors.iter().all(|neighbor| neighbor.from == 0));
    assert!(neighbors.iter().all(|neighbor| neighbor.to == 1));
}

#[test]
fn test_find_neighbors_by_radius_found() {
    let atoms = two_site_unit_atoms();
    // find_neighbors_from uses find_neighbors_by_radius internally
    let neighbors = atoms.find_neighbors_from(0, 1);
    let offsets: Vec<_> = neighbors.iter().map(|neighbor| neighbor.offset).collect();
    assert_eq!(offsets, vec![[-1, 0, 0], [0, 0, 0]]);
}

#[test]
fn test_find_neighbors_order_zero_is_empty() {
    let atoms = two_site_unit_atoms();

    assert!(atoms.find_neighbors_from(0, 0).is_empty());
    assert!(atoms.find_neighbors_from_to(0, 1, 0).is_empty());
}

// --- offset_bounds_in_radius edge cases ---

#[test]
fn test_offset_bounds_in_radius_zero_radius() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, true],
        tolerance: 0.0001,
    };
    let bounds = atoms.offset_bounds_in_radius([0.0, 0.0, 0.0], 0.0);
    assert_eq!(bounds, [[-1, 1], [-1, 1], [-1, 1]]);
}

#[test]
fn test_offset_bounds_in_radius_orthogonal_cell() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, true],
        tolerance: 0.000001,
    };

    let bounds = atoms.offset_bounds_in_radius([0.0, 0.0, 0.0], 1.0);

    assert_eq!(bounds, [[-2, 2], [-2, 2], [-2, 2]]);
}

#[test]
fn test_offset_bounds_in_radius_respects_non_periodic_axis() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let bounds = atoms.offset_bounds_in_radius([0.0, 0.0, 0.0], 2.0);

    assert_eq!(bounds[2], [0, 0]);
}

#[test]
fn test_offset_bounds_in_radius_covers_skewed_cell_offset() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.49, 0.01, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let bounds = atoms.offset_bounds_in_radius([0.0, 0.0, 0.0], 0.04);

    assert!(bounds[0][0] <= 1 && 1 <= bounds[0][1]);
    assert!(bounds[1][0] <= -2 && -2 <= bounds[1][1]);
}

#[test]
fn test_calc_distance_range_from_to_filters_by_distance() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let mut distances = atoms.calc_distance_range_from_to(0, 0, 0.9, 1.1);
    distances.sort_by_key(|distance| distance.neighbor.offset);
    let offsets: Vec<[isize; 3]> = distances
        .iter()
        .map(|distance| distance.neighbor.offset)
        .collect();

    assert_eq!(offsets, vec![[-1, 0, 0], [0, -1, 0], [0, 1, 0], [1, 0, 0]]);
}

#[test]
fn test_calc_distance_range_from_to_skips_self_origin() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, true],
        tolerance: 0.000001,
    };

    let distances = atoms.calc_distance_range_from_to(0, 0, 0.0, 0.1);

    assert!(distances.is_empty());
}

#[test]
fn test_calc_distance_range_from_to_respects_open_boundary() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0], [0., 0., 0.9]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let distances = atoms.calc_distance_range_from_to(0, 1, 0.05, 0.2);

    assert!(distances.is_empty());
}

#[test]
fn test_calc_distance_range_from_collects_all_targets() {
    let atoms = Atoms {
        cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
        pbc: [true, true, false],
        tolerance: 0.0001,
    };

    let distances = atoms.calc_distance_range_from(0, 0.9, 1.1);

    assert_eq!(distances.len(), 4);
    assert!(distances.iter().all(|distance| distance.neighbor.from == 0));
}

#[test]
fn test_calc_distance_range_all_collects_all_sources() {
    let atoms = Atoms {
        cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
        pbc: [true, true, false],
        tolerance: 0.0001,
    };

    let distances = atoms.calc_distance_range_all(0.9, 1.1);

    assert_eq!(distances.len(), 16);
}

#[test]
fn test_find_neighbors_simple() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let neighbors = atoms.find_neighbors_from(0, 2);

    assert_eq!(neighbors.len(), 4)
}

#[test]
fn test_find_neighbors_simple2() {
    let atoms = Atoms {
        cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
        pbc: [true, true, false],
        tolerance: 0.0001,
    };

    let neighbors = atoms.find_neighbors_from(0, 1);

    assert_eq!(neighbors.len(), 4);
    let neighbors = atoms.find_neighbors_from(0, 2);

    assert_eq!(neighbors.len(), 4);
}

#[test]
fn test_find_neighbors_all_simple() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let neighbors = atoms.find_neighbors_all(1);

    assert_eq!(neighbors.len(), 4)
}

#[test]
fn test_find_neighbors_all_simple2() {
    let atoms = Atoms {
        cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
        pbc: [true, true, false],
        tolerance: 0.0001,
    };
    let neighbors = atoms.find_neighbors_all(1);

    assert_eq!(neighbors.len(), 16)
}

#[test]
fn test_find_neighbors_from_to_simple2() {
    let atoms = Atoms {
        cell: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
        pbc: [true, true, false],
        tolerance: 0.0001,
    };
    let neighbors = atoms.find_neighbors_from_to(0, 1, 1);

    assert_eq!(neighbors.len(), 2)
}

#[test]
fn test_find_neighbors_from_handles_skewed_cell_large_offset() {
    let atoms = Atoms {
        cell: [[1.0, 0.0, 0.0], [0.49, 0.01, 0.0], [0.0, 0.0, 1.0]],
        positions: vec![[0., 0., 0.0]],
        pbc: [true, true, false],
        tolerance: 0.000001,
    };

    let neighbors = atoms.find_neighbors_from(0, 1);
    let offsets: Vec<[isize; 3]> = neighbors.iter().map(|neighbor| neighbor.offset).collect();

    assert!(offsets.contains(&[1, -2, 0]));
    assert!(offsets.contains(&[-1, 2, 0]));
}

#[test]
fn test_find_neighbors_cri3_single_layer_cr_cr_orders_1_to_3() {
    let cell = [
        [7.0037999153, 0.0000000000, 0.0000000000],
        [-3.5018999577, 6.0654686497, 0.0000000000],
        [0.0000000000, 0.0000000000, 23.1228008270],
    ];
    let positions = vec![
        direct_to_cart(cell, [0.000000000, 0.000000000, 0.067589996]),
        direct_to_cart(cell, [0.333330027, 0.666670062, 0.067469997]),
    ];
    let atoms = Atoms {
        cell,
        positions,
        pbc: [true, true, false],
        tolerance: 0.0001,
    };

    for (order, expected) in [
        (1, vec![[-1, -1, 0], [0, -1, 0], [0, 0, 0]]),
        (2, vec![[-1, -2, 0], [-1, 0, 0], [1, 0, 0]]),
        (
            3,
            vec![
                [-2, -2, 0],
                [-2, -1, 0],
                [0, -2, 0],
                [0, 1, 0],
                [1, -1, 0],
                [1, 1, 0],
            ],
        ),
    ] {
        let mut offsets: Vec<[isize; 3]> = atoms
            .find_neighbors_from_to(0, 1, order)
            .iter()
            .map(|neighbor| neighbor.offset)
            .collect();
        offsets.sort();

        assert_eq!(offsets, expected);
    }
}
