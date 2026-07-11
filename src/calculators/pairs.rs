use crate::{calculators::CalcInput, spin::SpinState};

/// Compute total pair energy for one site, with 1/2 factor to avoid
/// double-counting when summing over all sites.
pub(crate) fn energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
    local_energy(spin, calc_input) / 2.0
}

/// Compute local pair energy for one site. No 1/2 factor.
/// Used by local updates and energy differences.
pub(crate) fn local_energy<S: SpinState>(spin: &S, calc_input: &CalcInput<S>) -> f64 {
    scalar_exchange_energy(spin, calc_input)
}

pub(crate) fn scalar_exchange_energy<S: SpinState>(
    spin: &S,
    calc_input: &CalcInput<S>,
) -> f64 {
    let Some(neighbors) = calc_input.exchange_neighbors.as_ref() else {
        return 0.0;
    };

    neighbors
        .iter()
        .map(|(neighbor, j)| unsafe {
            let neighbor = &*(*neighbor);
            -j * spin.dot(neighbor)
        })
        .sum()
}
