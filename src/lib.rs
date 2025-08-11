pub mod calculators;
pub mod config;
pub mod lattice;
pub mod monte_carlo;
pub mod runner;
pub mod spin;
pub mod utils;

#[cfg(feature = "python-extension")]
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "python-extension")]
use runner::run;
#[cfg(feature = "python-extension")]
use tracing_subscriber::FmtSubscriber;

#[cfg(feature = "python-extension")]
#[pyfunction]
fn run_from_py(content: &str) -> PyResult<()> {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap();
    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");
    run(content).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(())
}

#[cfg(feature = "python-extension")]
#[pymodule]
fn _spinmc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_from_py, m)?)?;
    Ok(())
}
