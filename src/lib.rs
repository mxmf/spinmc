pub mod calculators;
pub mod config;
pub mod lattice;
pub mod monte_carlo;
pub mod runner;
pub mod spin;
pub mod utils;
use pyo3::{exceptions::PyValueError, prelude::*};
use runner::run;
use tracing_subscriber::FmtSubscriber;

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

#[pymodule]
fn _spinmc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_from_py, m)?)?;
    Ok(())
}
