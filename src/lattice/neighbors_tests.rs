use super::Atoms;

fn direct_to_cart(cell: [[f64; 3]; 3], direct: [f64; 3]) -> [f64; 3] {
    [
        direct[0] * cell[0][0] + direct[1] * cell[1][0] + direct[2] * cell[2][0],
        direct[0] * cell[0][1] + direct[1] * cell[1][1] + direct[2] * cell[2][1],
        direct[0] * cell[0][2] + direct[1] * cell[1][2] + direct[2] * cell[2][2],
    ]
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
