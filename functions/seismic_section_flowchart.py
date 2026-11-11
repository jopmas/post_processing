"""
Flowchart of seismic_section.py / seismic_section_parallel.py: from the Mandyoc
outputs to the synthetic seismic sections.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# Colours of the groups of steps
COLORS = {
    "input": "#dbe9f6",  # Mandyoc / tapioca files
    "model": "#dcefd9",  # model building
    "conv": "#fde5cc",  # convolutional section
    "fd": "#e6dcf0",  # tremelique
    "cache": "#f2f2f2",  # cache of the zoom
    "output": "#eeeeee",  # display and outputs
}


def _box(ax, x, y, w, h, title, text, color):
    """Box centred at (x, y) with a bold title and a smaller description."""
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc=COLORS[color], ec="0.35", lw=1))  # fmt: skip
    ax.text(x, y + h / 2 - 0.18, title, ha="center", va="top", fontsize=9.5,
            fontweight="bold")  # fmt: skip
    if text:
        ax.text(x, y + h / 2 - 0.5, text, ha="center", va="top", fontsize=8,
                linespacing=1.35)  # fmt: skip


def _arrow(ax, start, end, label=None, label_offset=(0.12, 0)):
    ax.annotate("", xy=end, xytext=start,
                arrowprops=dict(arrowstyle="-|>", color="0.3", lw=1.2,
                                shrinkA=0, shrinkB=0))  # fmt: skip
    if label:
        ax.text((start[0] + end[0]) / 2 + label_offset[0],
                (start[1] + end[1]) / 2 + label_offset[1], label, fontsize=8.5,
                style="italic", color="0.25", va="center")  # fmt: skip


def plot_flowchart(ax=None, parallel=False):
    """
    Plot the flowchart of the synthetic seismic section.

    Parameters
    ----------
    ax: `~matplotlib.axes.Axes`, optional
        Axis to draw on; a new figure is created if None.
    parallel: bool
        Flowchart of seismic_section_parallel.py (multithreaded tremelique and
        its own file of results) instead of seismic_section.py.

    Returns
    -------
    ax: `~matplotlib.axes.Axes`
    """
    prefix = "seismic_section_parallel" if parallel else "seismic_section"
    solver = "ParallelAcoustic (multithreaded)" if parallel else "tremelique Acoustic"
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 15), layout="constrained")
    ax.set_xlim(0, 12)
    ax.set_ylim(-0.9, 17.2)
    ax.axis("off")

    # Inputs and model building (trunk)
    _box(ax, 6, 16.3, 9.6, 1.3, "Mandyoc outputs (NetCDF from tapioca)",
         "lithology.nc (200 m), temperature.nc, strain.nc (1 km), surface.nc,\n"
         "interfaces.txt (rho0 of each layer), param.txt (thermal expansion α)", "input")  # fmt: skip
    _box(ax, 6, 14.35, 9.6, 1.3, "1. Load",
         "one time step (TIME), cropped around the window;\n"
         "depths referenced to sea level (sticky air removed)", "model")  # fmt: skip
    _box(ax, 6, 12.4, 9.6, 1.3, "2. Fine grid",
         "smoothed lithology indicators → largest wins (sharp contacts,\n"
         "no staircase); water above the seafloor (surface.nc)", "model")  # fmt: skip
    _box(ax, 6, 10.05, 9.6, 1.9, "3. Physical properties: density and Vp",
         "density rho0 (1 − α T); Vp per lithology − 0.4 m/s/K × T\n"
         "+ serpentinization: mantle below 400 °C, down to 4 km below the basement\n"
         "+ heterogeneity: von Kármán random media (laminated lower crust)\n"
         "+ faults: localized strain → damage zone + fabric along the band (brittle)",
         "model")  # fmt: skip
    for y0, y1 in [(15.65, 15.0), (13.7, 13.05), (11.75, 11.05)]:
        _arrow(ax, (6, y0), (6, y1))

    # Convolutional branch (left)
    _box(ax, 2.55, 6.7, 4.7, 1.8, "4. Convolutional section",
         "whole window (X0–X1)\n"
         "depth → TWT (true velocities)\n"
         "reflectivity (Z2 − Z1) / (Z2 + Z1)\n"
         "∗ Ricker 25 Hz, + noise, top mute", "conv")  # fmt: skip
    _arrow(ax, (4.2, 9.1), (2.55, 7.6))

    # tremelique branch (right)
    _box(ax, 8.75, 7.3, 6.1, 1.0, "5. Zoom model (XZ0–XZ1, DXZ)", "", "fd")
    ax.text(8.75, 7.08, "fast impedance-matched water, receivers on the seafloor",
            ha="center", va="center", fontsize=8)  # fmt: skip
    _arrow(ax, (7.8, 9.1), (8.75, 7.8))
    _box(ax, 8.75, 5.75, 6.1, 1.1, f"_output/{prefix}_zoom_<window>.nc",
         "exists, with the same parameters?", "cache")  # fmt: skip
    _arrow(ax, (8.75, 6.8), (8.75, 6.3))
    _box(ax, 7.6, 3.85, 3.8, 1.75, "Exploding reflector",
         f"{solver}\nhalf velocities, one simulation\n→ zero-offset data + noise",
         "fd")  # fmt: skip
    _arrow(ax, (7.6, 5.2), (7.6, 4.725), "no", (0.1, 0))
    _box(ax, 7.6, 1.95, 3.8, 1.55, "Reverse-time migration",
         "traces reversed in time,\nsmoothed model, Laplacian\n→ depth image; save .nc",
         "fd")  # fmt: skip
    _arrow(ax, (7.6, 2.975), (7.6, 2.725))
    _box(ax, 10.7, 3.85, 2.2, 1.2, "Load .nc", "skip the\nsimulations", "cache")
    _arrow(ax, (10.3, 5.2), (10.7, 4.45), "yes", (0.12, 0))

    # Display and outputs (trunk)
    _box(ax, 6, 0.0, 11.4, 1.2, "6. Display and outputs (_output folder of the scenario)",
         "depth → TWT with true velocities; AGC, grey scale; figures (model, convolutional,\n"
         "zoom sections, final figure) and the NetCDF of the zoom", "output")  # fmt: skip
    _arrow(ax, (2.55, 5.8), (2.55, 0.6))
    _arrow(ax, (7.6, 1.175), (7.6, 0.6))
    _arrow(ax, (10.7, 3.25), (10.7, 0.6))
    return ax
