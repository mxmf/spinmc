from pathlib import Path
from typing import Annotated
import cappa
from dataclasses import dataclass
from cappa import Subcommands


@dataclass
class Run:
    input: Annotated[Path, cappa.Arg(short="-i", parse=Path)] = Path("./config.toml")

    def __call__(self):

        with open(self.input) as f:
            from ._spinmc import run_from_py  # ty:ignore[unresolved-import]  # pyright: ignore[reportUnknownVariableType]

            toml_str = f.read()
            run_from_py(toml_str)


@dataclass
class Plot:
    input: Annotated[Path, cappa.Arg(short="-i", parse=Path)] = Path("./result.txt")

    def __call__(self):

        import math

        import matplotlib.pyplot as plt
        import numpy as np
        import numpy.typing as npt
        from scipy.signal import find_peaks

        def find_tc(
            t: npt.NDArray[np.floating], y: npt.NDArray[np.floating], ylabel: str
        ):
            if ("chi" in ylabel) or ("$C$" in ylabel):
                peaks, _ = find_peaks(y)
                if len(peaks) == 0:
                    return t[np.argmax(y)]
                return t[peaks[np.argmax(y[peaks])]]
            else:
                dy = np.gradient(y, t)
                peaks, _ = find_peaks(np.abs(dy))
                if len(peaks) == 0:
                    return t[np.argmax(np.abs(dy))]
                return t[peaks[np.argmax(np.abs(dy[peaks]))]]

        data = np.loadtxt(self.input, unpack=True)
        t = data[0]
        ys = data[1:]

        with open(self.input, "r") as f_result:
            ylabels = f_result.readline().split("\t")[1:]

        assert len(ys) == len(ylabels)
        _fig, axs = plt.subplots(math.ceil(len(ys) / 2), 2, layout="constrained")
        axs = np.array(axs).flatten()

        for ax, y, ylabel in zip(axs, ys, ylabels):
            ax.plot(t, y, marker=".")
            ax.set_xlabel("T")
            ax.set_ylabel(ylabel)

            tc = find_tc(t, np.array(y), ylabel)

            ax.axvline(tc, color="red", linestyle="--", linewidth=1)
            ax.text(
                tc,
                0.9,
                rf"$T_c = {tc:.3f}$",
                color="red",
                transform=ax.get_xaxis_transform(),
            )
        plt.show()


@dataclass
class Spinmc:
    cmd: Subcommands[Run | Plot]


def app():
    from importlib.metadata import version as _ver

    cappa.invoke(Spinmc, version=_ver("spinmc"))
