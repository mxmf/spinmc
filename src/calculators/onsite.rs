use crate::{calculators::CalcInput, spin::SpinState};

pub(crate) fn energy<S: SpinState>(
    spin: &S,
    calc_input: &CalcInput<S>,
    anisotropy_enable: bool,
    zeeman_enable: bool,
) -> f64 {
    let mut result = 0.0;
    if anisotropy_enable {
        result += anisotropy_energy(spin, calc_input);
    }
    if zeeman_enable {
        result += zeeman_energy(spin, calc_input);
    }
    result
}

pub(crate) fn anisotropy_energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
    let (strength, axis) = calc_input.anisotropy;
    let spin_array = spin.to_array();
    let dot = spin_array[0] * axis[0] + spin_array[1] * axis[1] + spin_array[2] * axis[2];

    -strength * dot * dot
}

pub(crate) fn zeeman_energy<S: SpinState>(_: &S, _: &CalcInput<S>) -> f64 {
    unimplemented!();
}
