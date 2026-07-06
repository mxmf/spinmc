use super::*;
use crate::spin::{HeisenbergSpin, IsingSpin, SpinState};

fn make_calc_input_with_exchange(spins: &[IsingSpin; 2], j: f64) -> CalcInput<IsingSpin> {
    CalcInput {
        magnitude: 1.0,
        exchange_neighbors: Some(vec![
            (&spins[0] as *const IsingSpin, j),
            (&spins[1] as *const IsingSpin, j),
        ]),
        exchanges: vec![j, j],
        exchange_neighbor_index: vec![0, 1],
        ..Default::default()
    }
}

fn make_calc_input_no_exchange() -> CalcInput<IsingSpin> {
    CalcInput::default()
}

// --- CalcInput ---

#[test]
fn calc_input_default() {
    let ci: CalcInput<IsingSpin> = Default::default();
    assert_eq!(ci.magnitude, 0.0);
    assert!(ci.exchange_neighbor_index.is_empty());
    assert!(ci.exchanges.is_empty());
    assert!(ci.exchange_neighbors.is_none());
    assert!(ci.dm_neighbors.is_none());
    assert!(ci.magnetic_field.is_none());
    assert!(ci.easy_axis.is_none());
    assert_eq!(ci.anisotropy, (0.0, [0.0, 0.0, 1.0]));
}

#[test]
fn validate_exchange_neighbor_none_ok() {
    let ci: CalcInput<IsingSpin> = Default::default();
    assert!(ci.validate_exchange_neighbor().is_ok());
}

#[test]
fn validate_exchange_neighbor_unique_ok() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = CalcInput {
        magnitude: 1.0,
        exchange_neighbors: Some(vec![
            (&spins[0] as *const IsingSpin, 1.0),
            (&spins[1] as *const IsingSpin, 1.0),
        ]),
        exchanges: vec![1.0, 1.0],
        exchange_neighbor_index: vec![0, 1],
        ..Default::default()
    };
    assert!(ci.validate_exchange_neighbor().is_ok());
}

#[test]
fn validate_exchange_neighbor_duplicate_error() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = CalcInput {
        magnitude: 1.0,
        exchange_neighbors: Some(vec![
            (&spins[0] as *const IsingSpin, 1.0),
            (&spins[0] as *const IsingSpin, 2.0),
        ]),
        exchanges: vec![1.0, 2.0],
        exchange_neighbor_index: vec![0, 1],
        ..Default::default()
    };
    assert!(ci.validate_exchange_neighbor().is_err());
}

#[test]
fn validate_exchange_neighbor_length_mismatch_error() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = CalcInput {
        magnitude: 1.0,
        exchange_neighbors: Some(vec![
            (&spins[0] as *const IsingSpin, 1.0),
            (&spins[1] as *const IsingSpin, 1.0),
        ]),
        exchanges: vec![1.0], // too short
        exchange_neighbor_index: vec![0, 1],
        ..Default::default()
    };
    assert!(ci.validate_exchange_neighbor().is_err());
}

// --- exchange_energy (with 1/2 factor) ---

#[test]
fn exchange_energy_ferromagnetic() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = make_calc_input_with_exchange(&spins, 2.0);
    let result = exchange_energy(&spins[0], &ci);
    // Two neighbors both parallel, each -j*dot = -2.0*1 = -2.0, sum = -4.0, /2 = -2.0
    assert!((result + 2.0).abs() < 1e-10);
}

#[test]
fn exchange_energy_antiferromagnetic() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(-1.0).unwrap(),
    ];
    let ci = make_calc_input_with_exchange(&spins, 2.0);
    let result = exchange_energy(&spins[0], &ci);
    // One neighbor parallel (-2*1=-2), one antiparallel (-2*(-1)=2), sum=0, /2=0
    assert!((result - 0.0).abs() < 1e-10);
}

#[test]
fn exchange_energy_no_neighbors_returns_zero() {
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = make_calc_input_no_exchange();
    assert_eq!(exchange_energy(&spin, &ci), 0.0);
}

// --- local_exchange_energy (no 1/2 factor) ---

#[test]
fn local_exchange_energy_ferromagnetic() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = make_calc_input_with_exchange(&spins, 2.0);
    let result = local_exchange_energy(&spins[0], &ci);
    // Two neighbors parallel: -4.0 total (no /2)
    assert!((result + 4.0).abs() < 1e-10);
}

#[test]
fn local_exchange_energy_no_neighbors_returns_zero() {
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = make_calc_input_no_exchange();
    assert_eq!(local_exchange_energy(&spin, &ci), 0.0);
}

#[test]
fn exchange_vs_local_differs_by_factor_two() {
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = make_calc_input_with_exchange(&spins, 2.0);
    let e_tot = exchange_energy(&spins[0], &ci);
    let e_loc = local_exchange_energy(&spins[0], &ci);
    assert!((2.0 * e_tot - e_loc).abs() < 1e-10);
}

// --- anisotropy_energy ---

#[test]
fn anisotropy_energy_parallel_to_axis() {
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = CalcInput {
        anisotropy: (3.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    // dot = 1, energy = -3 * 1^2 = -3
    assert!((anisotropy_energy(&spin, &ci) + 3.0).abs() < 1e-10);
}

#[test]
fn anisotropy_energy_perpendicular_to_axis() {
    let ci = CalcInput {
        anisotropy: (3.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    // Heisenberg along x: dot with z = 0, energy = 0
    let spin = HeisenbergSpin::along_x(1.0).unwrap();
    assert!((anisotropy_energy(&spin, &ci) + 0.0).abs() < 1e-10);
}

#[test]
fn anisotropy_energy_zero_strength() {
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = CalcInput {
        anisotropy: (0.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    assert_eq!(anisotropy_energy(&spin, &ci), 0.0);
}

#[test]
fn anisotropy_energy_heisenberg_parallel_to_axis() {
    let ci = CalcInput {
        anisotropy: (2.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    // Heisenberg along z, parallel to axis
    let spin = HeisenbergSpin::along_z(1.0).unwrap();
    assert!((anisotropy_energy(&spin, &ci) + 2.0).abs() < 1e-10);
}

#[test]
fn anisotropy_energy_heisenberg_perpendicular_to_axis() {
    let ci = CalcInput {
        anisotropy: (2.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    // Heisenberg along x, perpendicular to z-axis
    let spin = HeisenbergSpin::along_x(1.0).unwrap();
    assert!((anisotropy_energy(&spin, &ci) - 0.0).abs() < 1e-10);
}

// --- Hamiltonian ---

fn make_ham(exchange: bool, anisotropy: bool) -> Hamiltonian {
    Hamiltonian {
        config: HamiltonianConfig {
            exchange_enable: exchange,
            anisotropy_enable: anisotropy,
            zeeman_enable: false,
            dm_enable: false,
        },
    }
}

#[test]
fn hamiltonian_compute_all_disabled() {
    let ham = make_ham(false, false);
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = make_calc_input_with_exchange(
        &[
            IsingSpin::along_z(1.0).unwrap(),
            IsingSpin::along_z(1.0).unwrap(),
        ],
        2.0,
    );
    assert_eq!(ham.compute(&spin, &ci, &[]), 0.0);
}

#[test]
fn hamiltonian_compute_exchange_only() {
    let ham = make_ham(true, false);
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = make_calc_input_with_exchange(&spins, 2.0);
    let e_calc = ham.compute(&spins[0], &ci, &[]);
    let e_direct = exchange_energy(&spins[0], &ci);
    assert!((e_calc - e_direct).abs() < 1e-10);
}

#[test]
fn hamiltonian_compute_anisotropy_only() {
    let ham = make_ham(false, true);
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = CalcInput {
        anisotropy: (3.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    let e_calc = ham.compute(&spin, &ci, &[]);
    let e_direct = anisotropy_energy(&spin, &ci);
    assert!((e_calc - e_direct).abs() < 1e-10);
}

#[test]
fn hamiltonian_compute_exchange_plus_anisotropy() {
    let ham = make_ham(true, true);
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = CalcInput {
        magnitude: 1.0,
        exchange_neighbors: Some(vec![
            (&spins[0] as *const IsingSpin, 2.0),
            (&spins[1] as *const IsingSpin, 2.0),
        ]),
        exchanges: vec![2.0, 2.0],
        exchange_neighbor_index: vec![0, 1],
        anisotropy: (3.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    let e = ham.compute(&spins[0], &ci, &[]);
    let expected = exchange_energy(&spins[0], &ci) + anisotropy_energy(&spins[0], &ci);
    assert!((e - expected).abs() < 1e-10);
}

#[test]
fn hamiltonian_local_compute_uses_no_half_factor() {
    let ham = make_ham(true, false);
    let spins = [
        IsingSpin::along_z(1.0).unwrap(),
        IsingSpin::along_z(1.0).unwrap(),
    ];
    let ci = make_calc_input_with_exchange(&spins, 2.0);
    let e = ham.local_compute(&spins[0], &ci, &[]);
    assert!((e + 4.0).abs() < 1e-10);
}

#[test]
fn hamiltonian_compute_anisotropy_delegates() {
    let ham = make_ham(false, true);
    let spin = IsingSpin::along_z(1.0).unwrap();
    let ci = CalcInput {
        anisotropy: (5.0, [0.0, 0.0, 1.0]),
        ..Default::default()
    };
    let e = ham.compute_anisotropy(&spin, &ci);
    // -5 * (1)^2 = -5
    assert!((e + 5.0).abs() < 1e-10);
}
