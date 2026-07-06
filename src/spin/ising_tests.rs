use super::*;
use rand::SeedableRng;
use rand::rngs::SmallRng;

fn make_rng() -> SmallRng {
    SmallRng::seed_from_u64(42)
}

#[test]
fn zero_has_zero_state() {
    let s = IsingSpin::zero();
    assert_eq!(s.state, 0.0);
}

#[test]
fn along_x_returns_error() {
    let r = IsingSpin::along_x(1.0);
    assert!(r.is_err());
}

#[test]
fn along_y_returns_error() {
    let r = IsingSpin::along_y(1.0);
    assert!(r.is_err());
}

#[test]
fn along_z_positive() {
    let s = IsingSpin::along_z(2.5).unwrap();
    assert_eq!(s.state, 2.5);
}

#[test]
fn random_only_produces_plus_or_minus_magnitude() {
    let mag = 3.0;
    let mut rng = make_rng();
    for _ in 0..100 {
        let s = IsingSpin::random(&mut rng, mag);
        assert!(s.state == mag || s.state == -mag);
    }
}

#[test]
fn perturb_flips_sign() {
    let s = IsingSpin { state: 5.0 };
    let p = s.perturb(&mut make_rng(), 5.0);
    assert_eq!(p.state, -5.0);
}

#[test]
fn perturb_flips_negative_to_positive() {
    let s = IsingSpin { state: -3.0 };
    let p = s.perturb(&mut make_rng(), 3.0);
    assert_eq!(p.state, 3.0);
}

#[test]
fn dot_same_sign_equals_product() {
    let a = IsingSpin { state: 2.0 };
    let b = IsingSpin { state: 2.0 };
    assert_eq!(a.dot(&b), 4.0);
}

#[test]
fn dot_opposite_sign_negative_product() {
    let a = IsingSpin { state: 2.0 };
    let b = IsingSpin { state: -2.0 };
    assert_eq!(a.dot(&b), -4.0);
}

#[test]
fn norm_equals_absolute() {
    assert_eq!(IsingSpin { state: 3.0 }.norm(), 3.0);
    assert_eq!(IsingSpin { state: -3.0 }.norm(), 3.0);
    assert_eq!(IsingSpin { state: 0.0 }.norm(), 0.0);
}

#[test]
fn norm_sqr_equals_state_squared() {
    assert_eq!(IsingSpin { state: 3.0 }.norm_sqr(), 9.0);
    assert_eq!(IsingSpin { state: -3.0 }.norm_sqr(), 9.0);
}

#[test]
fn same_side_same_sign_true() {
    let a = IsingSpin { state: 1.0 };
    let b = IsingSpin { state: 2.0 };
    assert!(a.same_side(&b));
}

#[test]
fn same_side_opposite_sign_false() {
    let a = IsingSpin { state: 1.0 };
    let b = IsingSpin { state: -1.0 };
    assert!(!a.same_side(&b));
}

#[test]
fn same_side_zero_against_positive() {
    let a = IsingSpin { state: 0.0 };
    let b = IsingSpin { state: 1.0 };
    assert!(a.same_side(&b));
}

#[test]
fn same_side_zero_against_negative() {
    let a = IsingSpin { state: 0.0 };
    let b = IsingSpin { state: -1.0 };
    assert!(!a.same_side(&b));
}

#[test]
fn wolff_probability_zero_beta() {
    let s = IsingSpin { state: 1.0 };
    let p = s.wolff_probability(&s, &s, 0.0, 1.0, 1.0, 1.0);
    assert!((p - 0.0).abs() < 1e-10);
}

#[test]
fn wolff_probability_large_beta_approaches_one() {
    let s = IsingSpin { state: 1.0 };
    let p = s.wolff_probability(&s, &s, 100.0, 1.0, 1.0, 1.0);
    assert!((p - 1.0).abs() < 1e-10);
}

#[test]
fn wolff_probability_zero_j() {
    let s = IsingSpin { state: 1.0 };
    let p = s.wolff_probability(&s, &s, 1.0, 0.0, 1.0, 1.0);
    assert!((p - 0.0).abs() < 1e-10);
}

#[test]
fn flip_negates_state() {
    let s = IsingSpin { state: 3.0 };
    let axis = s;
    let flipped = s.flip(&axis);
    assert_eq!(flipped.state, -3.0);
}

#[test]
fn to_array_puts_state_in_z() {
    let s = IsingSpin { state: 2.0 };
    assert_eq!(s.to_array(), [0.0, 0.0, 2.0]);
}

#[test]
fn to_array_negative_state() {
    let s = IsingSpin { state: -1.5 };
    assert_eq!(s.to_array(), [0.0, 0.0, -1.5]);
}

#[test]
fn default_is_zero() {
    let s: IsingSpin = Default::default();
    assert_eq!(s.state, 0.0);
}

#[test]
fn neg() {
    let s = IsingSpin { state: 4.0 };
    assert_eq!((-s).state, -4.0);
}

#[test]
fn add() {
    let a = IsingSpin { state: 1.0 };
    let b = IsingSpin { state: 2.0 };
    assert_eq!((a + b).state, 3.0);
}

#[test]
fn add_assign() {
    let mut a = IsingSpin { state: 1.0 };
    a += IsingSpin { state: 2.0 };
    assert_eq!(a.state, 3.0);
}

#[test]
fn add_assign_ref() {
    let mut a = IsingSpin { state: 1.0 };
    let b = IsingSpin { state: 2.0 };
    a += &b;
    assert_eq!(a.state, 3.0);
}

#[test]
fn sub() {
    let a = IsingSpin { state: 5.0 };
    let b = IsingSpin { state: 2.0 };
    assert_eq!((a - b).state, 3.0);
}

#[test]
fn mul_f64() {
    let s = IsingSpin { state: 3.0 };
    assert_eq!((s * 2.0).state, 6.0);
}

#[test]
fn div_f64() {
    let s = IsingSpin { state: 6.0 };
    assert_eq!((s / 2.0).state, 3.0);
}

#[test]
fn div_f64_ref() {
    let s = IsingSpin { state: 6.0 };
    assert_eq!((&s / 2.0).state, 3.0);
}

#[test]
fn sum() {
    let spins = vec![
        IsingSpin { state: 1.0 },
        IsingSpin { state: 2.0 },
        IsingSpin { state: -3.0 },
    ];
    let total: IsingSpin = spins.into_iter().sum();
    assert_eq!(total.state, 0.0);
}

#[test]
fn sum_empty() {
    let spins: Vec<IsingSpin> = vec![];
    let total: IsingSpin = spins.into_iter().sum();
    assert_eq!(total.state, 0.0);
}

#[test]
fn div_f64_by_zero_produces_inf() {
    let s = IsingSpin { state: 1.0 };
    let result = s / 0.0;
    assert!(result.state.is_infinite());
}

#[test]
fn clone_and_copy() {
    let s = IsingSpin { state: 5.0 };
    let c = s;
    assert_eq!(c.state, 5.0);
    assert_eq!(s.state, 5.0);
}
