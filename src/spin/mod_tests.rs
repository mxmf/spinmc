use super::*;

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
