use super::*;
use rand::SeedableRng;
use rand::rngs::SmallRng;

fn make_rng() -> SmallRng {
    SmallRng::seed_from_u64(42)
}

#[test]
fn zero() {
    let s = HeisenbergSpin::zero();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 0.0);
    assert_eq!(s.z, 0.0);
}

#[test]
fn along_x() {
    let s = HeisenbergSpin::along_x(2.5).unwrap();
    assert_eq!(s.x, 2.5);
    assert_eq!(s.y, 0.0);
    assert_eq!(s.z, 0.0);
}

#[test]
fn along_y() {
    let s = HeisenbergSpin::along_y(3.0).unwrap();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 3.0);
    assert_eq!(s.z, 0.0);
}

#[test]
fn along_z() {
    let s = HeisenbergSpin::along_z(4.0).unwrap();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 0.0);
    assert_eq!(s.z, 4.0);
}

#[test]
fn random_on_unit_sphere() {
    let mut rng = make_rng();
    let mag = 3.0;
    for _ in 0..50 {
        let s = HeisenbergSpin::random(&mut rng, mag);
        let ns = s.norm_sqr();
        assert!((ns - mag * mag).abs() < 1e-10);
    }
}

#[test]
fn perturb_changes_direction() {
    let mut rng = make_rng();
    let s = HeisenbergSpin {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let p = s.perturb(&mut rng, 1.0);
    assert!((p.norm_sqr() - 1.0).abs() < 1e-10);
}

#[test]
fn dot_parallel() {
    let a = HeisenbergSpin {
        x: 2.0,
        y: 0.0,
        z: 0.0,
    };
    let b = HeisenbergSpin {
        x: 3.0,
        y: 0.0,
        z: 0.0,
    };
    assert_eq!(a.dot(&b), 6.0);
}

#[test]
fn dot_orthogonal() {
    let a = HeisenbergSpin {
        x: 2.0,
        y: 0.0,
        z: 0.0,
    };
    let b = HeisenbergSpin {
        x: 0.0,
        y: 3.0,
        z: 0.0,
    };
    assert_eq!(a.dot(&b), 0.0);
}

#[test]
fn dot_three_component() {
    let a = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let b = HeisenbergSpin {
        x: 4.0,
        y: 5.0,
        z: 6.0,
    };
    assert_eq!(a.dot(&b), 32.0);
}

#[test]
fn norm() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 2.0,
    };
    assert!((s.norm() - 3.0).abs() < 1e-10);
}

#[test]
fn norm_sqr() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 2.0,
    };
    assert_eq!(s.norm_sqr(), 9.0);
}

#[test]
fn same_side_always_true() {
    let a = HeisenbergSpin {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let b = HeisenbergSpin {
        x: -1.0,
        y: 0.0,
        z: 0.0,
    };
    let c = HeisenbergSpin {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };
    assert!(a.same_side(&b));
    assert!(a.same_side(&c));
}

#[test]
fn flip_reflection() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let axis = HeisenbergSpin {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };
    let old = s;
    let flipped = s.flip(&axis);
    let expected = old - axis * 2.0 * old.dot(&axis);
    assert!((flipped.x - expected.x).abs() < 1e-10);
    assert!((flipped.y - expected.y).abs() < 1e-10);
    assert!((flipped.z - expected.z).abs() < 1e-10);
}

#[test]
fn flip_reflection_parallel_axis() {
    let s = HeisenbergSpin {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };
    let axis = HeisenbergSpin {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };
    let flipped = s.flip(&axis);
    assert!((flipped.x - 0.0).abs() < 1e-10);
    assert!((flipped.y - 0.0).abs() < 1e-10);
    assert!((flipped.z + 1.0).abs() < 1e-10);
}

#[test]
fn to_array() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    assert_eq!(s.to_array(), [1.0, 2.0, 3.0]);
}

#[test]
fn default_is_zero() {
    let s: HeisenbergSpin = Default::default();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 0.0);
    assert_eq!(s.z, 0.0);
}

#[test]
fn neg() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: -2.0,
        z: 3.0,
    };
    let n = -s;
    assert_eq!(n.x, -1.0);
    assert_eq!(n.y, 2.0);
    assert_eq!(n.z, -3.0);
}

#[test]
fn add() {
    let a = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let b = HeisenbergSpin {
        x: 4.0,
        y: 5.0,
        z: 6.0,
    };
    let c = a + b;
    assert_eq!(c.x, 5.0);
    assert_eq!(c.y, 7.0);
    assert_eq!(c.z, 9.0);
}

#[test]
fn add_assign_val() {
    let mut a = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    a += HeisenbergSpin {
        x: 4.0,
        y: 5.0,
        z: 6.0,
    };
    assert_eq!(a.x, 5.0);
    assert_eq!(a.y, 7.0);
    assert_eq!(a.z, 9.0);
}

#[test]
fn add_assign_ref() {
    let mut a = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let b = HeisenbergSpin {
        x: 4.0,
        y: 5.0,
        z: 6.0,
    };
    a += &b;
    assert_eq!(a.x, 5.0);
    assert_eq!(a.y, 7.0);
    assert_eq!(a.z, 9.0);
}

#[test]
fn sub() {
    let a = HeisenbergSpin {
        x: 5.0,
        y: 3.0,
        z: 1.0,
    };
    let b = HeisenbergSpin {
        x: 2.0,
        y: 1.0,
        z: 0.0,
    };
    let c = a - b;
    assert_eq!(c.x, 3.0);
    assert_eq!(c.y, 2.0);
    assert_eq!(c.z, 1.0);
}

#[test]
fn mul_f64() {
    let s = HeisenbergSpin {
        x: 2.0,
        y: 3.0,
        z: 4.0,
    };
    let m = s * 2.0;
    assert_eq!(m.x, 4.0);
    assert_eq!(m.y, 6.0);
    assert_eq!(m.z, 8.0);
}

#[test]
fn div_f64() {
    let s = HeisenbergSpin {
        x: 4.0,
        y: 6.0,
        z: 8.0,
    };
    let d = s / 2.0;
    assert_eq!(d.x, 2.0);
    assert_eq!(d.y, 3.0);
    assert_eq!(d.z, 4.0);
}

#[test]
fn div_f64_ref() {
    let s = HeisenbergSpin {
        x: 4.0,
        y: 6.0,
        z: 8.0,
    };
    let d = &s / 2.0;
    assert_eq!(d.x, 2.0);
    assert_eq!(d.y, 3.0);
    assert_eq!(d.z, 4.0);
}

#[test]
fn sum() {
    let spins = vec![
        HeisenbergSpin {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        },
        HeisenbergSpin {
            x: 0.0,
            y: 2.0,
            z: 0.0,
        },
        HeisenbergSpin {
            x: 0.0,
            y: 0.0,
            z: 3.0,
        },
    ];
    let total: HeisenbergSpin = spins.into_iter().sum();
    assert_eq!(total.x, 1.0);
    assert_eq!(total.y, 2.0);
    assert_eq!(total.z, 3.0);
}

#[test]
fn sum_empty() {
    let spins: Vec<HeisenbergSpin> = vec![];
    let total: HeisenbergSpin = spins.into_iter().sum();
    assert_eq!(total.x, 0.0);
    assert_eq!(total.y, 0.0);
    assert_eq!(total.z, 0.0);
}

#[test]
fn div_f64_by_zero_produces_inf() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };
    let result = s / 0.0;
    assert!(result.x.is_infinite());
    assert!(result.y.is_infinite());
    assert!(result.z.is_infinite());
}

#[test]
fn clone_and_copy() {
    let s = HeisenbergSpin {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let c = s;
    assert_eq!(c.x, 1.0);
    assert_eq!(c.y, 2.0);
    assert_eq!(c.z, 3.0);
    assert_eq!(s.x, 1.0);
    assert_eq!(s.y, 2.0);
    assert_eq!(s.z, 3.0);
}
