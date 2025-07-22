use crate::lattice::Grid;
use crate::spin::SpinState;

#[derive(Default, Debug)]
pub struct Stats {
    pub energy_sum: f64,
    pub energy2_sum: f64,
    pub magnetization_sum: f64,
    pub magnetization2_sum: f64,
    pub steps: usize,
    pub size: f64,
}

impl Stats {
    pub fn new(size: usize) -> Self {
        Self {
            size: size as f64,
            ..Default::default()
        }
    }

    pub fn record<S: SpinState, R: rand::Rng>(&mut self, grid: &Grid<S, R>) {
        let e = grid.total_energy();
        let m = grid.total_spin_vector().norm();

        self.energy_sum += e;
        self.energy2_sum += e * e;
        self.magnetization_sum += m;
        self.magnetization2_sum += m * m;
        self.steps += 1;
    }

    pub fn mean_energy(&self) -> f64 {
        self.energy_sum / self.steps as f64 / self.size
    }

    pub fn specific_heat(&self, beta: f64) -> f64 {
        let e_avg = self.energy_sum / self.steps as f64;
        let e2_avg = self.energy2_sum / self.steps as f64;
        beta * beta * (e2_avg - e_avg * e_avg) / self.size
    }

    pub fn mean_magnetization(&self) -> f64 {
        self.magnetization_sum / self.steps as f64 / self.size
    }

    pub fn susceptibility(&self, beta: f64) -> f64 {
        let m_avg = self.magnetization_sum / self.steps as f64;
        let m2_avg = self.magnetization2_sum / self.steps as f64;
        beta * (m2_avg - m_avg * m_avg) / (self.size * self.size)
    }
}
