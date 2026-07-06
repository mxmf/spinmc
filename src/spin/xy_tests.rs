use super::*;
use rand::SeedableRng;
use rand::rngs::SmallRng;

fn make_rng() -> SmallRng {
    SmallRng::seed_from_u64(42)
}

#[test]
fn zero() {
    let s = XYSpin::zero();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 0.0);
}

#[test]
fn along_x() {
    let s = XYSpin::along_x(2.5).unwrap();
    assert_eq!(s.x, 2.5);
    assert_eq!(s.y, 0.0);
}

#[test]
fn along_y() {
    let s = XYSpin::along_y(3.0).unwrap();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 3.0);
}

#[test]
fn along_z_returns_error() {
    let r = XYSpin::along_z(1.0);
    assert!(r.is_err());
}

#[test]
fn random_on_unit_circle() {
    let mut rng = make_rng();
    let mag = 2.0;
    for _ in 0..50 {
        let s = XYSpin::random(&mut rng, mag);
        let ns = s.norm_sqr();
        assert!((ns - mag * mag).abs() < 1e-10);
    }
}

#[test]
fn perturb_changes_direction() {
    let mut rng = make_rng();
    let s = XYSpin { x: 1.0, y: 0.0 };
    let p = s.perturb(&mut rng, 1.0);
    assert!((p.norm_sqr() - 1.0).abs() < 1e-10);
}

#[test]
fn dot_parallel() {
    let a = XYSpin { x: 2.0, y: 0.0 };
    let b = XYSpin { x: 3.0, y: 0.0 };
    assert_eq!(a.dot(&b), 6.0);
}

#[test]
fn dot_orthogonal() {
    let a = XYSpin { x: 2.0, y: 0.0 };
    let b = XYSpin { x: 0.0, y: 3.0 };
    assert_eq!(a.dot(&b), 0.0);
}

#[test]
fn norm() {
    let s = XYSpin { x: 3.0, y: 4.0 };
    assert!((s.norm() - 5.0).abs() < 1e-10);
}

#[test]
fn norm_sqr() {
    let s = XYSpin { x: 3.0, y: 4.0 };
    assert_eq!(s.norm_sqr(), 25.0);
}

#[test]
fn same_side_always_true() {
    let a = XYSpin { x: 1.0, y: 0.0 };
    let b = XYSpin { x: -1.0, y: 0.0 };
    let c = XYSpin { x: 0.0, y: 1.0 };
    assert!(a.same_side(&b));
    assert!(a.same_side(&c));
}

#[test]
fn flip_reflection_perpendicular_axis() {
    let s = XYSpin { x: 1.0, y: 0.0 };
    let axis = XYSpin { x: 0.0, y: 1.0 };
    let old = s;
    let flipped = s.flip(&axis);
    assert!((flipped.x - 1.0).abs() < 1e-10);
    assert!((flipped.y - 0.0).abs() < 1e-10);
    let expected = old - axis * 2.0 * old.dot(&axis);
    assert!((flipped.x - expected.x).abs() < 1e-10);
    assert!((flipped.y - expected.y).abs() < 1e-10);
}

#[test]
fn flip_reflection_parallel_axis() {
    let s = XYSpin { x: 0.0, y: 1.0 };
    let axis = XYSpin { x: 0.0, y: 1.0 };
    let old = s;
    let flipped = s.flip(&axis);
    assert!((flipped.x).abs() < 1e-10);
    assert!((flipped.y + 1.0).abs() < 1e-10);
    let expected = old - axis * 2.0 * old.dot(&axis);
    assert!((flipped.x - expected.x).abs() < 1e-10);
    assert!((flipped.y - expected.y).abs() < 1e-10);
}

#[test]
fn flip_reflection_formula() {
    let mut s = XYSpin { x: 1.0, y: 1.0 };
    let norm = s.norm();
    s = XYSpin {
        x: s.x / norm,
        y: s.y / norm,
    };
    let axis = XYSpin { x: 1.0, y: 0.0 };
    let old = s;
    let flipped = s.flip(&axis);
    let expected = old - axis * 2.0 * old.dot(&axis);
    assert!((flipped.x - expected.x).abs() < 1e-10);
    assert!((flipped.y - expected.y).abs() < 1e-10);
}

#[test]
fn to_array_appends_zero_z() {
    let s = XYSpin { x: 1.0, y: 2.0 };
    assert_eq!(s.to_array(), [1.0, 2.0, 0.0]);
}

#[test]
fn default_is_zero() {
    let s: XYSpin = Default::default();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 0.0);
}

#[test]
fn neg() {
    let s = XYSpin { x: 1.0, y: -2.0 };
    let n = -s;
    assert_eq!(n.x, -1.0);
    assert_eq!(n.y, 2.0);
}

#[test]
fn add() {
    let a = XYSpin { x: 1.0, y: 2.0 };
    let b = XYSpin { x: 3.0, y: 4.0 };
    let c = a + b;
    assert_eq!(c.x, 4.0);
    assert_eq!(c.y, 6.0);
}

#[test]
fn add_assign_val() {
    let mut a = XYSpin { x: 1.0, y: 2.0 };
    a += XYSpin { x: 3.0, y: 4.0 };
    assert_eq!(a.x, 4.0);
    assert_eq!(a.y, 6.0);
}

#[test]
fn add_assign_ref() {
    let mut a = XYSpin { x: 1.0, y: 2.0 };
    let b = XYSpin { x: 3.0, y: 4.0 };
    a += &b;
    assert_eq!(a.x, 4.0);
    assert_eq!(a.y, 6.0);
}

#[test]
fn sub() {
    let a = XYSpin { x: 5.0, y: 3.0 };
    let b = XYSpin { x: 2.0, y: 1.0 };
    let c = a - b;
    assert_eq!(c.x, 3.0);
    assert_eq!(c.y, 2.0);
}

#[test]
fn mul_f64() {
    let s = XYSpin { x: 2.0, y: 3.0 };
    let m = s * 2.0;
    assert_eq!(m.x, 4.0);
    assert_eq!(m.y, 6.0);
}

#[test]
fn div_f64() {
    let s = XYSpin { x: 4.0, y: 6.0 };
    let d = s / 2.0;
    assert_eq!(d.x, 2.0);
    assert_eq!(d.y, 3.0);
}

#[test]
fn div_f64_ref() {
    let s = XYSpin { x: 4.0, y: 6.0 };
    let d = &s / 2.0;
    assert_eq!(d.x, 2.0);
    assert_eq!(d.y, 3.0);
}

#[test]
fn sum() {
    let spins = vec![
        XYSpin { x: 1.0, y: 0.0 },
        XYSpin { x: 0.0, y: 2.0 },
        XYSpin { x: 3.0, y: -1.0 },
    ];
    let total: XYSpin = spins.into_iter().sum();
    assert_eq!(total.x, 4.0);
    assert_eq!(total.y, 1.0);
}

#[test]
fn sum_empty() {
    let spins: Vec<XYSpin> = vec![];
    let total: XYSpin = spins.into_iter().sum();
    assert_eq!(total.x, 0.0);
    assert_eq!(total.y, 0.0);
}

#[test]
fn div_f64_by_zero_produces_inf() {
    let s = XYSpin { x: 1.0, y: 1.0 };
    let result = s / 0.0;
    assert!(result.x.is_infinite());
    assert!(result.y.is_infinite());
}

#[test]
fn clone_and_copy() {
    let s = XYSpin { x: 1.0, y: 2.0 };
    let c = s;
    assert_eq!(c.x, 1.0);
    assert_eq!(c.y, 2.0);
    assert_eq!(s.x, 1.0);
    assert_eq!(s.y, 2.0);
}
