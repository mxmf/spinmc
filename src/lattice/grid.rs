use crate::calculators::{CalcInput, Hamiltonian, HamiltonianConfig};
use crate::config::{Config, InitialState};
use crate::spin::SpinState;
use itertools::{iproduct, zip_eq};

pub struct Grid<S: SpinState, R: rand::Rng> {
    pub size: usize,
    pub spins: Vec<S>,
    pub calc_inputs: Vec<CalcInput<S>>,
    pub rng: R,
    pub hamiltonian: Hamiltonian,
    pub group_index: Vec<Vec<usize>>,
}

impl<S: SpinState, R: rand::Rng> Grid<S, R> {
    pub fn new(config: Config, mut rng: R) -> Self {
        let dim = config.dim;
        let sublattices = config.sublattices;
        let mut spins = vec![];
        let total_sites = dim[0] * dim[1] * dim[2];
        let mut group_index = Vec::new();

        for sub_lattice_group in config.group {
            let mut indexs = Vec::new();

            for sub_lattice in sub_lattice_group {
                for (x, y, z) in iproduct!(0..dim[0], 0..dim[1], 0..dim[2]) {
                    indexs.push(coord_to_index(
                        [x as isize, y as isize, z as isize],
                        sub_lattice,
                        dim,
                        sublattices,
                    ))
                }
            }
            group_index.push(indexs);
        }

        let mut calc_inputs: Vec<CalcInput<S>> = vec![];
        for magnitude in &config.spin_magnitudes {
            let new_spin = match config.initial_state {
                InitialState::Random => S::new_z(*magnitude),
                InitialState::X => S::new_x(*magnitude),
                InitialState::Y => S::new_y(*magnitude),
                InitialState::Z => S::new_z(*magnitude),
            };
            spins.extend(std::iter::repeat_n(new_spin, total_sites));
            calc_inputs.extend(std::iter::repeat_n(CalcInput::default(), total_sites));
        }

        if let InitialState::Random = &config.initial_state {
            for spin in &mut spins {
                *spin = spin.random(&mut rng);
            }
        }

        let hamiltonian = Hamiltonian::new(HamiltonianConfig {
            exchange_enable: true,
            anisotropy_enable: false,
            zeeman_enable: false,
            dm_enable: false,
        });

        for (sublattice, x, y, z) in iproduct!(0..sublattices, 0..dim[0], 0..dim[1], 0..dim[2]) {
            let index = coord_to_index(
                [x as isize, y as isize, z as isize],
                sublattice,
                dim,
                sublattices,
            );

            let calc_input = &mut calc_inputs[index];

            let mut exchange_neighbors = vec![];
            for exchange_param in &config.exchange_params {
                for offset in &exchange_param.offsets {
                    if exchange_param.from_sub == sublattice {
                        let offset_coord = [
                            offset[0] + x as isize,
                            offset[1] + y as isize,
                            offset[2] + z as isize,
                        ];

                        let offset_index_opt = safe_coord_to_index(
                            offset_coord,
                            exchange_param.to_sub,
                            dim,
                            sublattices,
                            config.pbc,
                        );
                        if let Some(offset_index) = offset_index_opt {
                            exchange_neighbors
                                .push((&spins[offset_index] as *const S, exchange_param.strength));
                        }
                    }
                }
            }

            calc_input.exchange_neighbors = Some(exchange_neighbors);
            calc_input.validate_exchange_neighbor();
        }

        Self {
            rng,
            size: total_sites * sublattices,
            spins,
            calc_inputs,
            hamiltonian,
            group_index,
        }
    }

    pub fn total_energy(&self) -> f64 {
        zip_eq(self.spins.iter(), self.calc_inputs.iter())
            .map(|(spin, calc_input)| spin.energy(calc_input, &self.hamiltonian, &self.spins))
            .sum::<f64>()
            / 2.0
    }

    pub fn partial_spin_vector(&self, index: usize) -> crate::spin::SpinVector {
        self.group_index[index]
            .iter()
            .map(|i| self.spins[*i].spinvector())
            .sum()
    }

    pub fn total_spin_vector(&self) -> crate::spin::SpinVector {
        self.spins.iter().map(|spin| spin.spinvector()).sum()
    }
}

fn _index_to_coord(index: usize, dim: [usize; 3], sublattices: usize) -> ([isize; 3], usize) {
    let [lx, ly, _lz] = dim;
    let sub_lattice = index % sublattices;
    let linear_index = index / sublattices;
    let z = linear_index / (lx * ly);
    let y = (linear_index / lx) % ly;
    let x = linear_index % lx;

    ([x as isize, y as isize, z as isize], sub_lattice)
}

fn coord_to_index(
    coord: [isize; 3],
    sublattice: usize,
    dim: [usize; 3],
    sublattices: usize,
) -> usize {
    let [x, y, z] = coord;
    let (x, y, z) = (x as usize, y as usize, z as usize);

    let linear_index = x + y * dim[0] + z * dim[0] * dim[1];
    linear_index * sublattices + sublattice
}

fn safe_coord_to_index(
    mut coord: [isize; 3],
    sublattice: usize,
    dim: [usize; 3],
    sublattices: usize,
    pbc: [bool; 3],
) -> Option<usize> {
    for i in 0..3 {
        if coord[i] < 0 || coord[i] >= dim[i] as isize {
            if pbc[i] {
                coord[i] = (coord[i] + dim[i] as isize) % dim[i] as isize;
            } else {
                return None;
            }
        }
    }

    if sublattice >= sublattices {
        return None;
    }

    Some(coord_to_index(coord, sublattice, dim, sublattices))
}
