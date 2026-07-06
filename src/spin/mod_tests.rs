use super::*;
use crate::calculators::{CalcInput, HamiltonianConfig};

fn exchange_calc_input<S: SpinState>(neighbor: &S) -> CalcInput<S> {
    CalcInput {
        magnitude: 1.0,
        exchange_neighbors: Some(vec![(neighbor as *const S, 1.0)]),
        exchanges: vec![1.0],
        exchange_neighbor_index: vec![0],
        dm_neighbors: None,
        magnetic_field: None,
        easy_axis: None,
        anisotropy: (0.0, [0.0, 0.0, 1.0]),
    }
}

fn anisotropy_calc_input<S: SpinState>() -> CalcInput<S> {
    CalcInput {
        magnitude: 1.0,
        exchange_neighbors: None,
        exchanges: vec![],
        exchange_neighbor_index: vec![],
        dm_neighbors: None,
        magnetic_field: None,
        easy_axis: None,
        anisotropy: (2.0, [0.0, 0.0, 1.0]),
    }
}

fn hamiltonian(exchange_enable: bool, anisotropy_enable: bool) -> Hamiltonian {
    Hamiltonian {
        config: HamiltonianConfig {
            exchange_enable,
            anisotropy_enable,
            zeeman_enable: false,
            dm_enable: false,
        },
    }
}

#[test]
fn default_wolff_probability_zero_beta() {
    let s = HeisenbergSpin::along_x(1.0).unwrap();
    let axis = HeisenbergSpin::along_z(1.0).unwrap();
    let p = s.wolff_probability(&s, &axis, 0.0, 1.0, 1.0, 1.0);
    assert!((p - 0.0).abs() < 1e-10);
}

#[test]
fn default_wolff_probability_large_beta_parallel() {
    let s = HeisenbergSpin::along_x(1.0).unwrap();
    let axis = HeisenbergSpin::along_x(1.0).unwrap();
    let p = s.wolff_probability(&s, &axis, 100.0, 1.0, 1.0, 1.0);
    assert!((p - 1.0).abs() < 1e-10);
}

#[test]
fn default_wolff_probability_orthogonal_axis() {
    let s = HeisenbergSpin::along_x(1.0).unwrap();
    let axis = HeisenbergSpin::along_z(1.0).unwrap();
    let p = s.wolff_probability(&s, &axis, 1.0, 1.0, 1.0, 1.0);
    assert!((p - 0.0).abs() < 1e-10);
}

#[test]
fn default_wolff_probability_zero_j() {
    let s = HeisenbergSpin::along_x(1.0).unwrap();
    let axis = HeisenbergSpin::along_x(1.0).unwrap();
    let p = s.wolff_probability(&s, &axis, 1.0, 0.0, 1.0, 1.0);
    assert!((p - 0.0).abs() < 1e-10);
}

#[test]
fn default_energy_delegates_to_hamiltonian_compute() {
    let spin = HeisenbergSpin::along_z(1.0).unwrap();
    let neighbor = HeisenbergSpin::along_z(1.0).unwrap();
    let ci = exchange_calc_input(&neighbor);
    let ham = hamiltonian(true, false);

    assert!((spin.energy(&ci, &ham, &[]) + 0.5).abs() < 1e-10);
}

#[test]
fn default_local_energy_delegates_to_hamiltonian_local_compute() {
    let spin = HeisenbergSpin::along_z(1.0).unwrap();
    let neighbor = HeisenbergSpin::along_z(1.0).unwrap();
    let ci = exchange_calc_input(&neighbor);
    let ham = hamiltonian(true, false);

    assert!((spin.local_energy(&ci, &ham, &[]) + 1.0).abs() < 1e-10);
}

#[test]
fn default_energy_diff_uses_local_energy_difference() {
    let new_spin = HeisenbergSpin::along_z(-1.0).unwrap();
    let old_spin = HeisenbergSpin::along_z(1.0).unwrap();
    let neighbor = HeisenbergSpin::along_z(1.0).unwrap();
    let ci = exchange_calc_input(&neighbor);
    let ham = hamiltonian(true, false);

    assert!((new_spin.energy_diff(&ci, &ham, &[], &old_spin) - 2.0).abs() < 1e-10);
}

#[test]
fn default_ion_anisotropy_energy_delegates() {
    let spin = HeisenbergSpin::along_z(1.0).unwrap();
    let ci = anisotropy_calc_input();
    let ham = hamiltonian(false, true);

    assert!((spin.ion_anisotropy_energy(&ci, &ham) + 2.0).abs() < 1e-10);
}

#[test]
fn default_ion_anisotropy_energy_diff_delegates() {
    let new_spin = HeisenbergSpin::along_z(1.0).unwrap();
    let old_spin = HeisenbergSpin::along_x(1.0).unwrap();
    let ci = anisotropy_calc_input();
    let ham = hamiltonian(false, true);

    assert!((new_spin.ion_anisotropy_energy_diff(&ci, &ham, &old_spin) + 2.0).abs() < 1e-10);
}
