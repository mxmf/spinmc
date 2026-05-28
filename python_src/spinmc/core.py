import inspect
from typing import Any, Literal, TypedDict
import numpy.typing as npt
import numpy as np


class StructureConfig(TypedDict, total=False):
    cell: list[list[float]]
    positions: list[list[float]]
    tolerance: float


class Spinmc:
    def __init__(self) -> None:
        self.params_dict: dict[str, Any] = {"exchange": []}

    def __save_dict(self, key: str) -> None:
        frame = inspect.currentframe()
        if frame is not None and frame.f_back is not None:
            args = frame.f_back.f_locals.copy()
            args.pop("self")
            filtered_args = {k: v for k, v in args.items() if v is not None}
            self.params_dict[key] = filtered_args

    def set_grid(
        self,
        dimensions: tuple[int, int, int],
        sublattices: int,
        spin_magnitudes: list[float],
        periodic_boundary: tuple[bool, bool, bool],
    ):
        self.__save_dict("grid")
        return self

    def set_structure(
        self,
        cell: list[list[float]],
        positions: list[list[float]],
        tolerance: float = 0.0001,
    ):
        self.__save_dict("structure")
        return self

    def set_simulation(
        self,
        initial_state: Literal["random", "x", "y", "z"],
        model: Literal["ising", "xy", "heisenberg"],
        equilibration_steps: int,
        measurement_steps: int,
        algorithm: Literal["wolff", "metropolis"],
        num_threads: int,
        temperatures: list[float] | npt.NDArray[np.floating],
        boltzmann_constand: float | None = None,
    ):
        self.__save_dict("simulation")
        self.params_dict["simulation"]["temperatures"] = [
            float(i) for i in self.params_dict["simulation"]["temperatures"]
        ]
        return self

    def add_exchange(
        self,
        strength: float,
        from_sublattice: int | None = None,
        to_sublattice: int | None = None,
        offsets: list[list[int]] | None = None,
        neighbor_order: int | None = None,
    ):
        self.__save_dict("exchange_tmp")
        self.params_dict["exchange"].append(self.params_dict.pop("exchange_tmp"))
        return self

    def set_anisotropy(
        self,
        axis: list[list[float]],
        strength: list[float],
    ):
        self.__save_dict("anisotropy")
        return self

    def set_output(
        self,
        outfile: str = "result.txt",
        stats_interval: int = 1,
        energy: bool = True,
        heat_capacit: bool = True,
        magnetization: bool = False,
        susceptibility: bool = False,
        magnetization_abs: bool = True,
        susceptibility_abs: bool = True,
        group_magnetization: bool = True,
        group_susceptibility: bool = True,
        group_magnetization_abs: bool = True,
        group_susceptibility_abs: bool = True,
        group: list[list[int]] = [],
    ):
        self.__save_dict("output")
        return self

    def set_snapshot(
        self,
        equilibration_interval: int,
        measurement_interval: int,
        compression_level: int,
        save_directory: str,
    ):
        self.__save_dict("snapshot")
        return self

    def run(self):
        from ._spinmc import run_from_py  # pyright: ignore[reportUnknownVariableType]
        import rtoml

        toml_str = rtoml.dumps(self.params_dict)

        run_from_py(toml_str)
