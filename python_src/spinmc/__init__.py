from typing import Annotated
import typer
from pathlib import Path


app = typer.Typer(
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
    pretty_exceptions_enable=False,
)


@app.command()
def run(
    input_file: Annotated[Path, typer.Option("--input", "-i", exists=True)] = Path(
        "./config.toml"
    ),  # pyright: ignore[reportCallInDefaultInitializer]
):
    with open(input_file) as f:
        from ._spinmc import run_from_py  # pyright: ignore[reportUnknownVariableType]

        toml_str = f.read()
        run_from_py(toml_str)


@app.command()
def plot(
    input_file: Annotated[Path, typer.Option("--input", "-i", exists=True)] = Path(
        "./result.txt"
    ),  # pyright: ignore[reportCallInDefaultInitializer]
):
    import matplotlib.pyplot as plt
    import math
    import numpy as np
    from scipy.signal import find_peaks

    def find_tc(t, y, ylabel: str):
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

    data = np.loadtxt(input_file, unpack=True)
    t = data[0]
    ys = data[1:]

    with open(input_file, "r") as f_result:
        ylabels = f_result.readline().split("\t")[1:]

    assert len(ys) == len(ylabels)
    fig, axs = plt.subplots(math.ceil(len(ys) / 2), 2, layout="constrained")
    axs = np.array(axs).flatten()

    for ax, y, ylabel in zip(axs, ys, ylabels):
        ax.plot(t, y, marker=".")
        ax.set_xlabel("T")
        ax.set_ylabel(ylabel)

        tc = find_tc(t, y, ylabel)

        ax.axvline(tc, color="red", linestyle="--", linewidth=1)
        ax.text(
            tc,
            0.9,
            rf"$T_c = {tc:.3f}$",
            color="red",
            transform=ax.get_xaxis_transform(),
        )
    plt.show()
