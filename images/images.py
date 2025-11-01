# FILTR demo diagram generator
# This script creates illustrative diagrams/plots used in the presentation.
# It saves PNGs to /mnt/data/filtr_demo and also displays them inline.
#
# Rules honored for this environment:
# - Uses matplotlib only (no seaborn).
# - Each chart uses a single figure (no subplots).
# - Does not set explicit colors or styles.

import os
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

OUT_DIR = "./images"
os.makedirs(OUT_DIR, exist_ok=True)


def savefig(name):
    path = os.path.join(OUT_DIR, name)
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    print(f"Saved: {path}")
    return path


# 1) Row-by-row iteration diagram
def plot_row_iteration(n=10, m=16):
    plt.figure(figsize=(6, 4), facecolor="white")
    ax = plt.gca()
    ax.set_facecolor("white")

    # Draw grid
    ax.set_xlim(-0.5, m - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(np.arange(-0.5, m, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=0.5)
    ax.set_xticks([])
    ax.set_yticks([])

    # Arrows across each row
    for i in range(n):
        for j in range(m - 1):
            ax.arrow(
                j,
                i,
                0.9,
                0,
                head_width=0.18,
                head_length=0.18,
                length_includes_head=True,
                color="black",
            )
        if i < n - 1:
            ax.arrow(
                m - 1,
                i,
                0,
                0.9,
                head_width=0.18,
                head_length=0.18,
                length_includes_head=True,
                color="black",
            )

    plt.title("Row-by-row iteration", color="black")
    p = savefig("row_iteration.png")
    plt.show()
    return p


# 2) Antidiagonal-by-antidiagonal iteration diagram
def plot_antidiagonal_iteration(n=10, m=16):
    plt.figure(figsize=(6, 4), facecolor="white")
    ax = plt.gca()
    ax.set_facecolor("white")

    # Draw grid
    ax.set_xlim(-0.5, m - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(np.arange(-0.5, m, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=0.5)
    ax.set_xticks([])
    ax.set_yticks([])

    # Arrows across each anti-diagonal
    for s in range(n + m - 1):
        coords = [(i, s - i) for i in range(n) if 0 <= s - i < m]
        coords.sort(key=lambda t: (t[0], t[1]))
        for (i1, j1), (i2, j2) in zip(coords[:-1], coords[1:]):
            ax.arrow(
                j1,
                i1,
                (j2 - j1) * 0.9,
                (i2 - i1) * 0.9,
                head_width=0.18,
                head_length=0.18,
                length_includes_head=True,
                color="black",
            )

    plt.title("Antidiagonal-by-antidiagonal iteration", color="black")
    p = savefig("antidiagonal_iteration.png")
    plt.show()
    return p


# 3) Speedup plot for schedule change (synthetic data)
def plot_speedup():
    rng = np.random.default_rng(0)  # reproducible noise
    sizes = np.array([256, 512, 1024, 2048, 4096])

    # Base runtimes
    row_time = 1.9 * np.sqrt(sizes) + 10
    anti_time = 0.9 * np.sqrt(sizes) + 8

    # Add small Gaussian noise (e.g. ±5%)
    row_time += rng.normal(0, 0.005 * row_time, size=row_time.shape)
    anti_time += rng.normal(0, 0.005 * anti_time, size=anti_time.shape)

    speedup = row_time / anti_time

    fig = plt.figure(figsize=(6, 4), facecolor="white")
    plt.plot(sizes, speedup, marker="o", color="black")
    plt.xlabel("Problem size (N)")
    plt.ylabel("Speedup (row / antidiagonal)")
    plt.title("Speedup from schedule change")
    plt.xscale("log", base=2)
    plt.grid(True)

    p = savefig("speedup_plot.png")
    plt.show()
    return p


# 4) Banded alignment diagram
def plot_banded_alignment(n=20, m=30, bandwidth=6):
    # Matrix with outside-band masked
    fig = plt.figure(figsize=(7, 4))
    mask = np.ones((n, m))  # 1 = white (outside band), 0 = grey (band)
    for i in range(n):
        for j in range(m):
            if abs(j - i * (m - 1) / (n - 1)) <= bandwidth:
                mask[i, j] = 0  # band region

    # Use a colormap: 0 -> grey, 1 -> white
    cmap = ListedColormap(["#cccccc", "#ffffff"])

    plt.imshow(mask, aspect="equal", cmap=cmap, vmin=0, vmax=1)
    plt.xticks(np.arange(-0.5, m, 1), [])
    plt.yticks(np.arange(-0.5, n, 1), [])
    plt.grid(which="both")
    plt.title("Banded Needleman–Wunsch")
    p = savefig("banded_alignment.png")
    plt.show()
    return p


def plot_xdrop_diagram():
    # Use the provided x_drop_mat
    mat = np.array(
        [
            [0, -1, -10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-1, 3, 2, -10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-2, 2, 6, 5, -10, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [-10, 1, 5, 9, 8, -10, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, -10, 4, 8, 12, 11, 10, -10, 0, 0, 0, 0, 0, 0],
            [0, 0, -10, 7, 11, 10, 14, 13, 12, 11, -10, 0, 0, 0],
            [0, 0, 0, -10, 10, 9, 13, 12, 11, -10, 0, 0, 0, 0],
            [0, 0, 0, 0, 9, -10, 12, 11, 15, 14, 13, -10, 0, 0],
            [0, 0, 0, 0, -10, 0, 11, 15, 14, 13, -10, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -10, 14, 13, 17, 16, -10, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 13, 12, 16, 20, 19, -10, 0],
            [0, 0, 0, 0, 0, 0, 0, 12, -10, 15, 19, 23, 22, 21],
            [0, 0, 0, 0, 0, 0, 0, -10, 0, -10, 18, 22, 21, 20],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -10, 21, 20, -10],
        ]
    )
    # Create a mask: 1 for nonzero (grey), 0 for zero (white)
    mask = (mat != 0).astype(int)
    cmap = ListedColormap(["#ffffff", "#cccccc"])  # white, grey

    fig = plt.figure(figsize=(7, 4))
    plt.imshow(mask, aspect="equal", cmap=cmap, vmin=0, vmax=1)
    plt.xticks(np.arange(-0.5, mat.shape[1], 1), [])
    plt.yticks(np.arange(-0.5, mat.shape[0], 1), [])
    plt.grid(which="both")
    plt.title("X-drop=3 diagram, AAAACCTGAAAGG x AAAAACGTAAAAA")
    p = savefig("xdrop_diagram.png")
    plt.show()
    return p


# 6) Runtime comparison (synthetic)
# def plot_runtime_comparison():
#     methods = ["NW (full)", "NW (banded)", "NW (x-drop)"]
#     runtime_ms = [1200, 520, 380]
#     fig = plt.figure(figsize=(6, 4))
#     x = np.arange(len(methods))
#     plt.bar(x, runtime_ms)
#     plt.xticks(x, methods, rotation=15)
#     plt.ylabel("Cells Explored")
#     plt.title("Cells Explore comparison")
#     plt.grid(axis="y")
#     p = savefig("cells_explored_comparison.png")
#     plt.show()
#     return p


# 7) Cells explored comparison (synthetic)
def plot_cells_explored():
    methods = ["NW (full)", "NW (banded)", "NW (x-drop)"]
    cells_thousands = [16, 6.5, 4.1]  # synthetic
    fig = plt.figure(figsize=(6, 4))
    x = np.arange(len(methods))
    plt.bar(x, cells_thousands)
    plt.xticks(x, methods, rotation=15)
    plt.ylabel("Cells explored (thousands)")
    plt.title("Search-space size")
    plt.grid(axis="y")
    p = savefig("cells_explored.png")
    plt.show()
    return p


# 8) Performance summary table (synthetic)
# def plot_performance_table():
#     fig = plt.figure(figsize=(7, 2.8))
#     plt.axis("off")
#     col_labels = ["Method", "Runtime (ms)", "Cells (M)", "Speedup vs Full"]
#     rows = [
#         ["NW (full)", "1200", "16.0", "1.00×"],
#         ["NW (banded)", "520", "6.5", "2.31×"],
#         ["NW (x-drop)", "380", "4.1", "3.16×"],
#     ]
#     table = plt.table(cellText=rows, colLabels=col_labels, loc="center")
#     table.scale(1, 1.6)
#     plt.title("Performance summary (synthetic)")
#     p = savefig("performance_table.png")
#     plt.show()
#     return p


paths = []
# paths.append(plot_row_iteration())
# paths.append(plot_antidiagonal_iteration())
# paths.append(plot_speedup())
# paths.append(plot_banded_alignment())
paths.append(plot_xdrop_diagram())
# paths.append(plot_cells_explored())
# paths.append(plot_performance_table())

print("\nGenerated files:")
for p in paths:
    print(p)
