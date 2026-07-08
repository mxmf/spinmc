pub mod calculators;
pub mod config;
pub mod lattice;
pub mod monte_carlo;
pub mod runner;
pub mod spin;
pub mod utils;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[cfg(feature = "python-extension")]
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "python-extension")]
use runner::run;
#[cfg(feature = "python-extension")]
use tracing_subscriber::FmtSubscriber;

#[cfg(feature = "python-extension")]
#[pyfunction]
fn run_from_py(content: &str) -> PyResult<()> {
    let _ = ctrlc::set_handler(|| std::process::exit(2));
    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    let _ = tracing::subscriber::set_global_default(subscriber);
    run(content).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(())
}

#[cfg(feature = "python-extension")]
#[pymodule]
fn _spinmc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_from_py, m)?)?;
    Ok(())
}
