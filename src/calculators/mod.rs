pub mod hamiltonian;
pub mod input;
pub mod onsite;
pub mod pairs;

pub use hamiltonian::{Hamiltonian, HamiltonianConfig};
pub use input::CalcInput;

#[cfg(test)]
pub(crate) use onsite::{anisotropy_energy, zeeman_energy};
#[cfg(test)]
pub(crate) use pairs::{energy as exchange_energy, local_energy as local_exchange_energy};

#[cfg(test)]
#[path = "calculators_tests.rs"]
mod tests;
