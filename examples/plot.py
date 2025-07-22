import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


def estimate_tc_derivative(t, y):
    dy = np.gradient(y, t)
    peaks, _ = find_peaks(np.abs(dy))
    if len(peaks) == 0:
        return t[np.argmax(np.abs(dy))]
    return t[peaks[np.argmax(np.abs(dy[peaks]))]]


def estimate_tc_peak(t, y):
    peaks, _ = find_peaks(y)
    if len(peaks) == 0:
        return t[np.argmax(y)]
    return t[peaks[np.argmax(y[peaks])]]


if __name__ == "__main__":
    t, e, cv, m, chi = np.loadtxt("./result.txt", unpack=True)
    ys = (e, cv, m, chi)
    ylabels = ("E", r"$C_V$", "|M|", r"$\chi$")
    fig, axs = plt.subplots(2, 2, layout="constrained")
    axs = axs.flat

    for ax, y, ylabel in zip(axs, ys, ylabels):
        ax.plot(t, y, marker=".")
        ax.set_xlabel("T")
        ax.set_ylabel(ylabel)

        if ylabel in (r"$C_V$", r"$\chi$"):
            tc = estimate_tc_peak(t, y)
        else:
            tc = estimate_tc_derivative(t, y)

        ax.axvline(tc, color="red", linestyle="--", linewidth=1)
        ax.text(
            tc,
            0.9,
            rf"$T_c = {tc:.3f}$",
            color="red",
            transform=ax.get_xaxis_transform(),
        )

    plt.show()
