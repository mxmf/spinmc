use super::*;
use crate::calculators::HamiltonianConfig;
use crate::config::Config;
use crate::lattice::Grid;
use crate::spin::{IsingSpin, SpinState};
use rand::SeedableRng;
use rand::rngs::SmallRng;

fn single_spin_grid() -> Grid<IsingSpin, SmallRng> {
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
"#;
    let config = Config::new(toml).unwrap();
    Grid::new(&config, SmallRng::seed_from_u64(42)).unwrap()
}

#[test]
fn any_mc_beta_set_beta() {
    let mc = AnyMC::Metropolis(Metropolis {
        rng: SmallRng::seed_from_u64(0),
        beta: 2.0,
    });
    assert!((mc.beta() - 2.0).abs() < 1e-10);

    let mut mc = AnyMC::Wolff(Wolff {
        rng: SmallRng::seed_from_u64(0),
        beta: 3.0,
        ham_config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        },
    });
    assert!((mc.beta() - 3.0).abs() < 1e-10);
    mc.set_beta(5.0);
    assert!((mc.beta() - 5.0).abs() < 1e-10);
}

#[test]
fn any_mc_step_dispatches_metropolis() {
    let mut grid = single_spin_grid();
    let mut mc = AnyMC::Metropolis(Metropolis {
        rng: SmallRng::seed_from_u64(0),
        beta: 1.0,
    });

    assert_eq!(mc.step(&mut grid), 1);
    assert_eq!(grid.spins[0].to_array(), [0.0, 0.0, -1.0]);
}

#[test]
fn any_mc_step_dispatches_wolff() {
    let mut grid = single_spin_grid();
    let mut mc = AnyMC::Wolff(Wolff {
        rng: SmallRng::seed_from_u64(0),
        beta: 0.0,
        ham_config: HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        },
    });

    assert_eq!(mc.step(&mut grid), 1);
    assert_eq!(grid.spins[0].to_array(), [0.0, 0.0, -1.0]);
}
