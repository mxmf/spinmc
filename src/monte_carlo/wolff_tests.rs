use super::*;
use crate::calculators::HamiltonianConfig;
use crate::config::Config;
use crate::lattice::Grid;
use crate::monte_carlo::Wolff;
use crate::spin::IsingSpin;
use rand::SeedableRng;
use rand::rngs::SmallRng;

fn make_ferro_grid() -> crate::lattice::Grid<IsingSpin, SmallRng> {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "wolff"

[grid]
dimensions = [4, 4, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(42);
    Grid::new(&config, rng).unwrap()
}

#[test]
fn wolff_high_beta_flips_entire_cluster() {
    let mut grid = make_ferro_grid();
    let initial_spins: Vec<f64> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    let mut wolff = Wolff {
        rng: SmallRng::seed_from_u64(99),
        beta: 1e6,
        ham_config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        },
    };
    let cluster_size = wolff.step(&mut grid);
    // At very high beta, Ising ferromagnet cluster should include all sites
    assert_eq!(cluster_size, grid.size);
    // After flipping, the spins should differ from initial
    let final_spins: Vec<f64> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    let flipped_all = initial_spins
        .iter()
        .zip(final_spins.iter())
        .all(|(a, b)| (a + b).abs() < 1e-10);
    assert!(
        flipped_all,
        "high-beta Wolff step should flip the whole cluster"
    );
}

#[test]
fn wolff_zero_beta_small_cluster() {
    let mut grid = make_ferro_grid();
    let mut wolff = Wolff {
        rng: SmallRng::seed_from_u64(99),
        beta: 0.0,
        ham_config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        },
    };
    // At beta=0, wolff_probability = 1 - exp(0) = 0, so no bonds activated
    // Cluster = just the seed spin, so cluster size should be 1
    let cluster_size = wolff.step(&mut grid);
    assert_eq!(cluster_size, 1);
}

#[test]
fn wolff_repeated_steps() {
    let mut grid = make_ferro_grid();
    let mut wolff = Wolff {
        rng: SmallRng::seed_from_u64(99),
        beta: 10.0,
        ham_config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        },
    };
    for _ in 0..10 {
        let size = wolff.step(&mut grid);
        assert!(size > 0, "cluster should never be empty without anisotropy");
    }
}

fn make_aniso_grid() -> crate::lattice::Grid<IsingSpin, SmallRng> {
    let toml = r#"
[simulation]
initial_state = "z"
model = "ising"
equilibration_steps = 100
measurement_steps = 1000
temperatures = [1.0]
num_threads = 1
algorithm = "wolff"

[grid]
dimensions = [2, 2, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[1, 0, 0], [0, 1, 0]]
strength = 1.0

[anisotropy]
axis = [[0.0, 0.0, 1.0]]
strength = [2.0]

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(42);
    Grid::new(&config, rng).unwrap()
}

#[test]
fn wolff_with_anisotropy_runs() {
    let mut grid = make_aniso_grid();
    let mut wolff = Wolff {
        rng: SmallRng::seed_from_u64(99),
        beta: 1.0,
        ham_config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: true,
            zeeman_enable: false,
            dm_enable: false,
        },
    };
    let size = wolff.step(&mut grid);
    assert!(size > 0);
    assert!(size <= grid.size);
}
