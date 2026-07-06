use super::*;
use crate::calculators::{CalcInput, Hamiltonian, HamiltonianConfig};
use crate::spin::{IsingSpin, SpinState};
use rand::SeedableRng;
use rand::rngs::SmallRng;

// --- coord_to_index ---

#[test]
fn coord_to_index_basic() {
    // dim = [2, 3, 4], base = 2*3*4 = 24 per sublattice
    // sublattice 1, coord (0,1,2): 1*24 + 0*12 + 1*4 + 2 = 30
    let idx = coord_to_index([0, 1, 2], 1, [2, 3, 4]);
    assert_eq!(idx, 30);
}

#[test]
fn coord_to_index_origin() {
    let idx = coord_to_index([0, 0, 0], 0, [3, 3, 3]);
    assert_eq!(idx, 0);
}

#[test]
fn coord_to_index_last_site() {
    // dim [2,2,2], base=8, sublattice 1, coord (1,1,1): 1*8 + 1*4 + 1*2 + 1 = 15
    let idx = coord_to_index([1, 1, 1], 1, [2, 2, 2]);
    assert_eq!(idx, 15);
}

// --- safe_coord_to_index ---

#[test]
fn safe_coord_to_index_in_bounds() {
    let idx = safe_coord_to_index([0, 0, 0], 0, [3, 3, 3], 1, [true, true, true]);
    assert_eq!(idx, Some(0));
}

#[test]
fn safe_coord_to_index_pbc_wrap_negative() {
    // x=-1, wrap to x=2 (dim=3)
    let idx = safe_coord_to_index([-1, 0, 0], 0, [3, 3, 3], 1, [true, true, true]);
    // coord = (2,0,0), sublattice=0, base=27
    // 0*27 + 2*9 + 0*3 + 0 = 18
    assert_eq!(idx, Some(18));
}

#[test]
fn safe_coord_to_index_pbc_wrap_positive() {
    let idx = safe_coord_to_index([3, 0, 0], 0, [3, 3, 3], 1, [true, true, true]);
    // wraps to (0,0,0)
    assert_eq!(idx, Some(0));
}

#[test]
fn safe_coord_to_index_non_periodic_out_of_bounds() {
    let idx = safe_coord_to_index([-1, 0, 0], 0, [3, 3, 3], 1, [false, true, true]);
    assert_eq!(idx, None);
}

#[test]
fn safe_coord_to_index_non_periodic_positive_oob() {
    let idx = safe_coord_to_index([3, 0, 0], 0, [3, 3, 3], 1, [false, true, true]);
    assert_eq!(idx, None);
}

#[test]
fn safe_coord_to_index_sublattice_oob() {
    // sublattice 1 >= sublattices 1 → None
    let idx = safe_coord_to_index([0, 0, 0], 1, [3, 3, 3], 1, [true, true, true]);
    assert_eq!(idx, None);
}

#[test]
fn safe_coord_to_index_partial_pbc() {
    // only x periodic; y=-1 out of bounds, not periodic → None
    let idx = safe_coord_to_index([0, -1, 0], 0, [3, 3, 3], 2, [true, false, true]);
    assert_eq!(idx, None);
}

// --- Grid helper ---

fn make_ising_grid_2x1_with_exchange(j: f64) -> Grid<IsingSpin, SmallRng> {
    let dim = [2usize, 1, 1];
    let size = 2;
    let spins = vec![
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    // SAFETY: raw pointers obtained here stay valid because `spins` Vec's heap allocation
    // doesn't move when the Vec is moved into the Grid below.
    let p0: *const IsingSpin = &spins[0];
    let p1: *const IsingSpin = &spins[1];

    let calc_inputs = vec![
        CalcInput {
            magnitude: 1.0,
            exchange_neighbors: Some(vec![(p1, j)]),
            exchanges: vec![j],
            exchange_neighbor_index: vec![1],
            ..Default::default()
        },
        CalcInput {
            magnitude: 1.0,
            exchange_neighbors: Some(vec![(p0, j)]),
            exchanges: vec![j],
            exchange_neighbor_index: vec![0],
            ..Default::default()
        },
    ];

    let hamiltonian = Hamiltonian {
        config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        },
    };
    let rng = SmallRng::seed_from_u64(0);

    Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs,
        rng,
        hamiltonian,
        group_index: vec![vec![0, 1]],
    }
}

// --- total_energy ---

#[test]
fn total_energy_ferromagnetic_2_spins() {
    let grid = make_ising_grid_2x1_with_exchange(1.0);
    // Each spin sees neighbor: exchange_energy = -j * (+1)*(+1) / 2 = -0.5
    // Total = 2 * (-0.5) = -1.0
    let e = grid.total_energy();
    assert!((e + 1.0).abs() < 1e-10, "got {e}");
}

#[test]
fn total_energy_antiferromagnetic_2_spins() {
    let dim = [2usize, 1, 1];
    let size = 2;
    let spins = vec![
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(-1.0).unwrap(),
    ];
    let p0: *const IsingSpin = &spins[0];
    let p1: *const IsingSpin = &spins[1];

    let calc_inputs = vec![
        CalcInput {
            magnitude: 1.0,
            exchange_neighbors: Some(vec![(p1, 1.0)]),
            exchanges: vec![1.0],
            exchange_neighbor_index: vec![1],
            ..Default::default()
        },
        CalcInput {
            magnitude: 1.0,
            exchange_neighbors: Some(vec![(p0, 1.0)]),
            exchanges: vec![1.0],
            exchange_neighbor_index: vec![0],
            ..Default::default()
        },
    ];
    let grid = Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs,
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: true,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![vec![0, 1]],
    };
    let e = grid.total_energy();
    // Spin 0 (+1) sees neighbor (-1): -1 * (-1) / 2 = 0.5
    // Spin 1 (-1) sees neighbor (+1): -1 * (-1) / 2 = 0.5  (dot=-1, -j*dot = 1)
    // Wait: exchange_energy per site: sum(-j*dot)/2.
    // Spin 0: -1 * (+1)*(-1) = -1 * -1 = 1. 1/2 = 0.5
    // Spin 1: -1 * (-1)*(+1) = -1 * -1 = 1. 1/2 = 0.5
    // Total = 1.0
    assert!((e - 1.0).abs() < 1e-10, "got {e}");
}

// --- total_spin_vector ---

#[test]
fn total_spin_vector_all_up() {
    let dim = [2usize, 1, 1];
    let size = 2;
    let spins = vec![
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(2.0).unwrap(),
    ];
    let grid = Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![],
    };
    let total = grid.total_spin_vector();
    assert_eq!(total.to_array(), [0.0, 0.0, 3.0]);
}

#[test]
fn total_spin_vector_with_negative() {
    let dim = [2usize, 1, 1];
    let size = 2;
    let spins = vec![
        IsingSpin::along_z(3.0).unwrap(),
        IsingSpin::along_z(-1.0).unwrap(),
    ];
    let grid = Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![],
    };
    let total = grid.total_spin_vector();
    assert_eq!(total.to_array(), [0.0, 0.0, 2.0]);
}

// --- partial_spin_vector ---

#[test]
fn partial_spin_vector_group() {
    let dim = [2usize, 1, 1];
    let size = 2;
    let spins = vec![
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(3.0).unwrap(),
    ];
    let grid = Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![vec![0], vec![1], vec![0, 1]],
    };
    assert_eq!(grid.partial_spin_vector(0).to_array(), [0.0, 0.0, 1.0]);
    assert_eq!(grid.partial_spin_vector(1).to_array(), [0.0, 0.0, 3.0]);
    assert_eq!(grid.partial_spin_vector(2).to_array(), [0.0, 0.0, 4.0]);
}

// --- get_spin_by_coord ---

#[test]
fn get_spin_by_coord_in_bounds() {
    let dim = [2usize, 2, 1];
    let size = 4;
    let spins = vec![
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(2.0).unwrap(),
        IsingSpin::along_z(3.0).unwrap(),
        IsingSpin::along_z(4.0).unwrap(),
    ];
    let grid = Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![],
    };
    let s = grid.get_spin_by_coord(0, 0, 1, 0);
    // index = 0*4 + 0*2 + 1*1 + 0 = 1 → spin value 2.0
    assert_eq!(s.unwrap().to_array(), [0.0, 0.0, 2.0]);
}

#[test]
fn get_spin_by_coord_out_of_bounds() {
    let dim = [2usize, 2, 1];
    let size = 4;
    let spins = vec![
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(2.0).unwrap(),
        IsingSpin::along_z(3.0).unwrap(),
        IsingSpin::along_z(4.0).unwrap(),
    ];
    let grid = Grid {
        spins,
        size,
        dim,
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![],
    };
    assert!(grid.get_spin_by_coord(0, 3, 0, 0).is_none());
}

#[test]
fn get_spin_by_coord_rejects_each_out_of_bounds_axis() {
    let dim = [2usize, 2, 2];
    let spins = (0..8)
        .map(|i| IsingSpin::along_z(i as f64).unwrap())
        .collect();
    let grid = Grid {
        spins,
        size: 8,
        dim,
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![],
    };

    assert!(grid.get_spin_by_coord(0, -1, 0, 0).is_none());
    assert!(grid.get_spin_by_coord(0, 0, -1, 0).is_none());
    assert!(grid.get_spin_by_coord(0, 0, 0, -1).is_none());
    assert!(grid.get_spin_by_coord(0, 2, 0, 0).is_none());
    assert!(grid.get_spin_by_coord(0, 0, 2, 0).is_none());
    assert!(grid.get_spin_by_coord(0, 0, 0, 2).is_none());
}

#[test]
fn get_spin_by_coord_rejects_out_of_bounds_sublattice() {
    let grid = Grid {
        spins: vec![IsingSpin::along_z(1.0).unwrap()],
        size: 1,
        dim: [1, 1, 1],
        num_sublattices: 1,
        calc_inputs: vec![],
        rng: SmallRng::seed_from_u64(0),
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
        group_index: vec![],
    };

    assert!(grid.get_spin_by_coord(1, 0, 0, 0).is_none());
}

// --- rng is accessible ---

#[test]
fn grid_rng_accessible() {
    let rng = SmallRng::seed_from_u64(42);
    let grid = Grid {
        spins: vec![IsingSpin::along_z(1.0).unwrap()],
        size: 1,
        dim: [1, 1, 1],
        num_sublattices: 1,
        calc_inputs: vec![],
        rng,
        group_index: vec![],
        hamiltonian: Hamiltonian {
            config: HamiltonianConfig {
                exchange_enable: false,
                anisotropy_enable: false,
                zeeman_enable: false,
                dm_enable: false,
            },
        },
    };
    // Just verify rng is still accessible after moving into Grid
    let _ = grid.rng;
}
