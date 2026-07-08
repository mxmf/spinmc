# SpinMC

SpinMC is an open-source Monte Carlo simulation package for classical spin models, designed for studying magnetic systems and phase transitions in statistical physics.

## ✨ Features

- **Multiple Spin Models**:
  - Ising model
  - XY model
  - Heisenberg model

- **Support Interactions**:
  - Exchange coupling (isotropic)
  - Automatic exchange neighbor generation by coordination shell or distance range
  - Single-ion anisotropy
  - [Planned] Anisotropic exchange coupling
  - [Planned] External magnetic field

- **Support Algorithms**:
  - Metropolis
  - Wolff cluster algorithm
  - [Planned] Parallel Tempering (Replica Exchange)

- **Simulation Capabilities**:
  - Energy, Heat capacity
  - Magnetization, Susceptibility
  - Absolute magnetization, absolute susceptibility
  - Group-wise (sublattice) magnetization & susceptibility
  - Group-wise (sublattice) Absolute magnetization & susceptibility
  - [planned] magnetic hysteresis loop
  - [planned] binder cumulant (u4)

## 📦 Installation

### Requirements

- Python 3.9+ (for pip installation)
- Rust 1.88+ (for source compilation)

### From PyPI

```bash
pip install spinmc
```

### From Pre-built Releases

1. Visit the [Releases page](../../releases)
2. Download the package for your OS
3. Unpack the archive and run the executable

### From source:

```bash
git clone https://github.com/mxmf/spinmc.git
cd spinmc
cargo build --release
```

## 🚀 Quick Start

1. Create a configuration file (e.g., `ising.toml`):

```toml
[grid]
dimensions = [50, 50, 1]
sublattices = 1
spin_magnitudes = [1.0]
periodic_boundary = [true, true, false]

[simulation]
initial_state = "random"
boltzmann_constant = 1
model = "ising"
equilibration_steps = 10000
measurement_steps = 100000
algorithm = "wolff"
num_threads = 10
temperature_range = [
  { start = 1, end = 3, step = 0.1 },
]

[output]
savefile = "result.txt"
energy = true
heat_capacity = true
magnetization = true
susceptibility = true
magnetization_abs = true
susceptibility_abs = true

[[exchange]]
from_sublattice = 0
to_sublattice = 0
offsets = [[0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]]
strength = 1.0
```

Exchange interactions can also be generated from the `[structure]` section instead of listing every offset manually. The structure can be written directly with `cell` and `positions`, or loaded from an external structure file:

```toml
[structure]
file = "POSCAR"
# format = "poscar"          # Optional; currently only POSCAR/CONTCAR and .vasp files are supported.
# magnetic_indices = [0, 2]  # Optional; select magnetic atoms when the file also contains non-magnetic atoms.
tolerance = 0.0001
```

Use `neighbor_order = 1` for the first coordination shell, or `distance_range = [min, max]` to include every pair whose distance falls within the given range in Å. Specify both `from_sublattice` and `to_sublattice` to target one pair, only `from_sublattice` to target all destinations from one source sublattice, or neither to apply the rule to all sublattice pairs.

See `examples/heisenberg_2d_cri3_poscar` for a CrI3 monolayer example loaded from a POSCAR file.

2. Run the simulation:

```bash
spinmc run -i ising.toml
```

3. Results will be saved in `result.txt` with the requested measurements.

4. If `spinmc` is installed via **Python**, you can plot the results with:

```bash
spinmc plot -i result.txt
```

## 📚 See more examples in the [examples folder](examples).
