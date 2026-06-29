# CrI3 monolayer from POSCAR

This example loads the CrI3 monolayer structure from `POSCAR` and uses
`magnetic_indices = [6, 7]` to select the two Cr atoms as magnetic sites.
The nearest-neighbor Cr-Cr exchange paths are generated automatically from
the structure with `neighbor_order = 1`.

Run it from this directory so the relative `POSCAR` path resolves correctly:

```bash
cd examples/heisenberg_2d_cri3_poscar
spinmc run -i config.toml
```
