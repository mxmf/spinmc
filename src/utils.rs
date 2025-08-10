pub fn fmt_fixed_width(num: f64, width: usize) -> String {
    if num.abs() >= 1e6 || (num != 0.0 && num.abs() < 1e-4) {
        format!("{num:<width$.5e}")
    } else {
        format!("{num:<width$.6}")
    }
}
