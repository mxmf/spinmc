use crate::spin::SpinState;
use std::collections::HashSet;
#[derive(Clone, Debug)]

pub struct CalcInput<S: SpinState> {
    pub magnitude: f64,
    pub exchange_neighbors: Option<Vec<(*const S, f64)>>,
    pub exchanges: Vec<f64>,
    pub exchange_neighbor_index: Vec<usize>,
    pub dm_neighbors: Option<Vec<(usize, [f64; 3], f64)>>,
    pub magnetic_field: Option<[f64; 3]>,
    pub easy_axis: Option<[f64; 3]>,
    pub anisotropy: (f64, [f64; 3]),
}
unsafe impl<S: SpinState> Send for CalcInput<S> {}
unsafe impl<S: SpinState> Sync for CalcInput<S> {}

impl<S: SpinState> Default for CalcInput<S> {
    fn default() -> Self {
        CalcInput {
            magnitude: 0.0,
            exchange_neighbor_index: vec![],
            exchanges: vec![],
            exchange_neighbors: None,
            dm_neighbors: None,
            magnetic_field: None,
            easy_axis: None,
            anisotropy: (0., [0., 0., 1.]),
        }
    }
}

impl<S: SpinState> CalcInput<S> {
    pub fn validate_exchange_neighbor(&self) -> anyhow::Result<()> {
        let pairs = &self.exchange_neighbors;
        if let Some(vec) = pairs {
            let mut seen = HashSet::new();
            for (index, _) in vec {
                if !seen.insert(index) {
                    anyhow::bail!(
                        "Duplicate neighbor indices found in your exchange coupling configuration. Please ensure all neighbor indices are unique."
                    );
                }
            }

            if vec.len() != self.exchange_neighbor_index.len() || vec.len() != self.exchanges.len()
            {
                anyhow::bail!(
                    "unexpected exchange data length: exchange_ptr={}, exchange_neighbors={}, exchanges={}",
                    vec.len(),
                    self.exchange_neighbor_index.len(),
                    self.exchanges.len()
                );
            }
        }
        Ok(())
    }
}
