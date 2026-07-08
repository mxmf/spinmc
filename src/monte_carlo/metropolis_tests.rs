use super::*;
use crate::config::Config;
use crate::lattice::Grid;
use crate::monte_carlo::Metropolis;
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
algorithm = "metropolis"

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

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(42);
    Grid::new(&config, rng).unwrap()
}

#[test]
fn metropolis_step_zero_beta_accepts_all() {
    let mut grid = make_ferro_grid();
    let initial_spins: Vec<f64> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    let mut mc = Metropolis {
        rng: SmallRng::seed_from_u64(1),
        beta: 0.0,
    };
    let changed = mc.step(&mut grid);
    assert_eq!(changed, grid.size);
    // At beta=0, every move is accepted, all spins get flipped (Ising perturb flips sign)
    let final_spins: Vec<f64> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    let flipped_all = initial_spins
        .iter()
        .zip(final_spins.iter())
        .all(|(a, b)| (a + b).abs() < 1e-10);
    assert!(flipped_all, "beta=0 should flip all spins");
}

#[test]
fn metropolis_step_high_beta_preserves_ground_state() {
    let mut grid = make_ferro_grid();
    let initial_energy = grid.total_energy();
    let mut mc = Metropolis {
        rng: SmallRng::seed_from_u64(1),
        beta: 1e10,
    };
    mc.step(&mut grid);
    // Ground state should be preserved (very low T)
    assert!(
        (grid.total_energy() - initial_energy).abs() < 1e-10,
        "ground state should not change at low T"
    );
}

#[test]
fn metropolis_step_returns_grid_size() {
    let mut grid = make_ferro_grid();
    let mut mc = Metropolis {
        rng: SmallRng::seed_from_u64(1),
        beta: 1.0,
    };
    let changed = mc.step(&mut grid);
    assert_eq!(changed, grid.size);
}

#[test]
fn metropolis_single_spin_grid() {
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
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, true]

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[0, 0, 0]]
strength = 0.0

[output]
energy = true
group = [[0]]
"#;
    let config = Config::new(toml).unwrap();
    let rng = SmallRng::seed_from_u64(42);
    let mut grid: crate::lattice::Grid<IsingSpin, SmallRng> = Grid::new(&config, rng).unwrap();

    let mut mc = Metropolis {
        rng: SmallRng::seed_from_u64(1),
        beta: 1.0,
    };
    let changed = mc.step(&mut grid);
    assert_eq!(changed, 1);
    assert_eq!(grid.spins[0].to_array(), [0.0, 0.0, -1.0]);
}

#[test]
fn metropolis_zero_exchange_accepts_all_moves() {
    let mut grid = make_ferro_grid();
    for calc_input in &mut grid.calc_inputs {
        calc_input.exchanges.fill(0.0);
        if let Some(neighbors) = &mut calc_input.exchange_neighbors {
            for (_, j) in neighbors {
                *j = 0.0;
            }
        }
    }
    let spins_before: Vec<_> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    let mut mc = Metropolis {
        rng: SmallRng::seed_from_u64(1),
        beta: 1e10,
    };
    mc.step(&mut grid);
    let spins_after: Vec<_> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    for (before, after) in spins_before.iter().zip(spins_after.iter()) {
        assert!((before + after).abs() < 1e-10);
    }
}

#[test]
fn metropolis_zero_temperature_rejects_equal_energy_moves() {
    let mut grid = make_ferro_grid();
    for calc_input in &mut grid.calc_inputs {
        calc_input.exchanges.fill(0.0);
        if let Some(neighbors) = &mut calc_input.exchange_neighbors {
            for (_, j) in neighbors {
                *j = 0.0;
            }
        }
    }
    let spins_before: Vec<_> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    let mut mc = Metropolis {
        rng: SmallRng::seed_from_u64(1),
        beta: f64::INFINITY,
    };
    mc.step(&mut grid);
    let spins_after: Vec<_> = grid.spins.iter().map(|s| s.to_array()[2]).collect();
    assert_eq!(spins_after, spins_before);
}
