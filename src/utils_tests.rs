use super::*;

#[test]
fn large_positive_uses_scientific() {
    let result = fmt_fixed_width(2.5e6, 12);
    assert!(result.contains('e'));
}

#[test]
fn large_negative_uses_scientific() {
    let result = fmt_fixed_width(-3.0e7, 12);
    assert!(result.contains('e'));
}

#[test]
fn small_positive_uses_scientific() {
    let result = fmt_fixed_width(5e-5, 12);
    assert!(result.contains('e'));
}

#[test]
fn small_negative_uses_scientific() {
    let result = fmt_fixed_width(-8e-6, 12);
    assert!(result.contains('e'));
}

#[test]
fn normal_positive_fixed_point() {
    let result = fmt_fixed_width(3.125, 12);
    assert!(!result.contains('e'));
}

#[test]
fn normal_negative_fixed_point() {
    let result = fmt_fixed_width(-42.0, 12);
    assert!(!result.contains('e'));
}

#[test]
fn zero_fixed_point() {
    let result = fmt_fixed_width(0.0, 12);
    assert!(!result.contains('e'));
}

#[test]
fn boundary_1e6_uses_scientific() {
    let result = fmt_fixed_width(1e6, 12);
    assert!(result.contains('e'));
}

#[test]
fn boundary_neg_1e6_uses_scientific() {
    let result = fmt_fixed_width(-1e6, 12);
    assert!(result.contains('e'));
}

#[test]
fn boundary_1e4_neg_fixed_point() {
    let result = fmt_fixed_width(1e-4, 12);
    assert!(!result.contains('e'));
}

#[test]
fn boundary_neg_1e4_neg_fixed_point() {
    let result = fmt_fixed_width(-1e-4, 12);
    assert!(!result.contains('e'));
}

#[test]
fn just_below_1e6_fixed_point() {
    let result = fmt_fixed_width(999_999.0, 12);
    assert!(!result.contains('e'));
}

#[test]
fn just_below_1e4_neg_uses_scientific() {
    let result = fmt_fixed_width(9.9e-5, 12);
    assert!(result.contains('e'));
}

#[test]
fn width_is_respected() {
    let result = fmt_fixed_width(42.0, 15);
    assert_eq!(result.len(), 15);
}
