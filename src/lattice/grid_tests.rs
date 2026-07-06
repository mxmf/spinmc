use super::*;
use crate::calculators::{CalcInput, Hamiltonian, HamiltonianConfig};
use crate::config::Config;
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

// --- Grid::new() with Config ---

fn minimal_config(toml_extra: &str) -> Config {
    let base = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
group = [[0]]
"#;
    let full = format!("{base}\n{toml_extra}");
    Config::new(&full).unwrap()
}

#[test]
fn grid_new_basic() {
    let config = minimal_config("");
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.size, 4);
    assert_eq!(grid.dim, [2, 2, 1]);
    assert_eq!(grid.num_sublattices, 1);
    assert_eq!(grid.spins.len(), 4);
}

#[test]
fn grid_new_exchange_neighbors_built() {
    let config = minimal_config("");
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    // Every site has exactly one +x neighbor; x=1 wraps to x=0 because PBC is enabled.
    for ci in &grid.calc_inputs {
        let neighbors = ci.exchange_neighbors.as_ref().unwrap();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(ci.exchange_neighbor_index.len(), 1);
        assert_eq!(ci.exchanges, vec![1.0]);
    }
}

#[test]
fn grid_new_can_compute_energy() {
    let config = minimal_config("");
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    let e = grid.total_energy();
    // Four sites each contribute -J*s_i*s_j/2 = -0.5 for the +x bond.
    assert!((e + 2.0).abs() < 1e-10, "got {e}");
}

#[test]
fn grid_new_group_index() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 2
spin_magnitudes = [1.0, 1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 1
offsets = [[0, 0, 0]]
strength = 1.0

[output]
energy = true
group = [[0, 1]]
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.num_sublattices, 2);
    assert_eq!(grid.size, 8);
    // group_index len = number of groups
    assert_eq!(grid.group_index.len(), 1);
    // One group containing both sublattices → 8 indices
    assert_eq!(grid.group_index[0].len(), 8);
}

#[test]
fn grid_new_initial_state_random() {
    let toml = r#"
[simulation]
initial_state = "random"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [10, 1, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    // Random Ising spins should all be ±1
    for s in &grid.spins {
        assert!((s.to_array()[2].abs() - 1.0).abs() < 1e-10);
    }
}

#[test]
fn grid_new_initial_state_xy() {
    let toml = r#"
[simulation]
initial_state = "x"
model = "xy"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 1, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<crate::spin::XYSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.spins.len(), 2);
    for spin in &grid.spins {
        assert_eq!(spin.to_array(), [1.0, 0.0, 0.0]);
    }
}

#[test]
fn grid_new_initial_state_heisenberg() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "heisenberg"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 1, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<crate::spin::HeisenbergSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.spins.len(), 2);
    for spin in &grid.spins {
        assert_eq!(spin.to_array(), [0.0, 0.0, 1.0]);
    }
}

#[test]
fn grid_new_initial_state_heisenberg_y() {
    let toml = r#"
[simulation]
initial_state = "y"
model = "heisenberg"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 1, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<crate::spin::HeisenbergSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.spins.len(), 2);
    for spin in &grid.spins {
        assert_eq!(spin.to_array(), [0.0, 1.0, 0.0]);
    }
}

#[test]
fn grid_new_non_periodic_boundary_excludes_oob_neighbors() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 1, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [false, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.calc_inputs[0].exchange_neighbor_index, vec![1]);
    assert!(grid.calc_inputs[1].exchange_neighbor_index.is_empty());
}

#[test]
fn grid_new_with_anisotropy() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[anisotropy]
axis = [[0.0, 0.0, 1.0]]
strength = [2.0]

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    // All calc_inputs should have anisotropy set
    for ci in &grid.calc_inputs {
        assert!((ci.anisotropy.0 - 2.0).abs() < 1e-10);
        assert_eq!(ci.anisotropy.1, [0.0, 0.0, 1.0]);
    }
}

#[test]
fn grid_new_multi_sublattice() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [1, 1, 1]
sublattices = 3
spin_magnitudes = [1.0, 2.0, 3.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[0, 0, 0]]
strength = 0.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.size, 3);
    assert_eq!(grid.num_sublattices, 3);
    // Magnitudes should match
    assert!((grid.calc_inputs[0].magnitude - 1.0).abs() < 1e-10);
    assert!((grid.calc_inputs[1].magnitude - 2.0).abs() < 1e-10);
    assert!((grid.calc_inputs[2].magnitude - 3.0).abs() < 1e-10);
    assert_eq!(grid.spins[0].to_array(), [0.0, 0.0, 1.0]);
    assert_eq!(grid.spins[1].to_array(), [0.0, 0.0, 2.0]);
    assert_eq!(grid.spins[2].to_array(), [0.0, 0.0, 3.0]);
}

#[test]
fn grid_new_partial_spin_vector_works() {
    let config = minimal_config("");
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    // With 1 group of 1 sublattice, partial_spin_vector should be accessible
    let v = grid.partial_spin_vector(0);
    // All spins are +1, 4 sites → total = 4.0 in z
    assert!((v.to_array()[2] - 4.0).abs() < 1e-10);
}

#[test]
fn grid_new_total_spin_vector_works() {
    let config = minimal_config("");
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    let total = grid.total_spin_vector();
    assert!((total.to_array()[2] - 4.0).abs() < 1e-10);
}

#[test]
fn grid_new_3d() {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "metropolis"

[grid]
dimensions = [3, 2, 2]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0]]
strength = 1.0

[output]
energy = true
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(0);
    let grid: Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();
    assert_eq!(grid.size, 12);
    assert_eq!(grid.dim, [3, 2, 2]);
}
