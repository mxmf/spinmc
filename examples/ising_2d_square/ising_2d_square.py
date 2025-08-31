from spinmc import Spinmc
import numpy as np

ising_2d_square = (
    Spinmc()
    .set_grid(
        dimensions=(50, 50, 1),
        sublattices=1,
        spin_magnitudes=[1.0],
        periodic_boundary=(True, True, False),
    )
    .add_exchange(
        strength=1,
        neighbor_order=1,
    )
    .set_structure(
        cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        positions=[[0, 0, 0]],
    )
    .set_simulation(
        initial_state="random",
        model="ising",
        equilibration_steps=100000,
        measurement_steps=100000,
        algorithm="metropolis",
        num_threads=10,
        temperatures=np.arange(0, 3.1, 0.1),
        boltzmann_constand=1,
    )
    .set_output()
).run()
