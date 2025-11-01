import random
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
import math
from ed import edit_distance
from banded_ed import banded_edit_distance
from x_drop import x_drop
from z_drop import z_drop_edit_distance
from wfa_ed import wfa_edit_distance


def random_dna(length):
    return ''.join(random.choice("ACGT") for _ in range(length))


def mutate_sequence(seq, divergence):
    """
    Introduce random substitutions to create divergence%.
    Only uses substitutions to maintain equal sequence length for square matrices.
    """
    seq = list(seq)
    n = len(seq)
    num_changes = int(n * divergence)
    
    # Only use substitutions to keep sequences the same length
    indices = random.sample(range(n), min(num_changes, n))
    for idx in indices:
        # Change to a different base
        current = seq[idx]
        new_base = random.choice([b for b in "ACGT" if b != current])
        seq[idx] = new_base
    
    return ''.join(seq)


def print_matrix(matrix):
    for row in matrix:
        print("\t".join(f"{val:5}" if val != float("-inf") else " -inf" for val in row))

def count_cells(matrix):
    # print_matrix(matrix)
    """Count cells that were actually computed (not pruned)."""
    count = 0
    for row in matrix:
        for val in row:
            if val != math.inf:
                count += 1
    return count



def compare_methods(seq_len=1000, divergence_levels=np.linspace(0, 0.8, 11)):
    results = {
        "divergence": [],
        "edit_distance": [],
        "banded": [],
        "xdrop": [],
        "zdrop": [],
        "wfa": [],
        "cells_edit": [],
        "cells_banded": [],
        "cells_xdrop": [],
        "cells_zdrop": [],
        "cells_wfa": [],
    }

    for div in divergence_levels:
        seq1 = random_dna(seq_len)
        seq2 = mutate_sequence(seq1, div)
        
        # Verify sequences are equal length for square matrix
        print(f"\n{'='*60}")
        print(f"Divergence: {div*100:.1f}%")
        print(f"Seq1 length: {len(seq1)}, Seq2 length: {len(seq2)}")
        assert len(seq1) == len(seq2), f"Sequences must be equal length! seq1={len(seq1)}, seq2={len(seq2)}"

        # Ground truth
        D_full, dist_true, cells_full = edit_distance(seq1, seq2)
        cells_full = len(seq1) * len(seq2)
        print(f"Matrix dimensions: {len(D_full)}x{len(D_full[0])} (should be {len(seq1)+1}x{len(seq2)+1})")
        print(f"Full DP: {dist_true}")
        print(f"Full DP cells: {cells_full}")

        # Banded
        D_band, dist_band, cells_band = banded_edit_distance(seq1, seq2, bandwidth=2)
        cells_band = count_cells(D_band)
        print(f"Banded DP: {dist_band}")
        print(f"Banded DP cells: {cells_band}")

        # X-drop (with edit distance scoring: match=0, mismatch=1, gap=1)
        # Use x_drop threshold of 10 (reasonable for edit distance)
        dist_x, dp_x = x_drop(seq1, seq2, x_drop=10)
        cells_x = count_cells(dp_x)
        print(f"X-drop DP: {dist_x}")
        print(f"X-drop DP cells: {cells_x}")


        # WFA
        W, dist_wfa, cells_wfa = wfa_edit_distance(seq1, seq2)
        print(f"WFA DP: {dist_wfa}")
        print(f"WFA DP cells: {cells_wfa}")

        # Store results
        results["divergence"].append(div * 100)
        results["edit_distance"].append(dist_true)
        results["banded"].append(dist_band)
        results["xdrop"].append(dist_x)
        results["wfa"].append(dist_wfa)

        results["cells_edit"].append(cells_full)
        results["cells_banded"].append(cells_band)
        results["cells_xdrop"].append(cells_x)
        results["cells_wfa"].append(cells_wfa)

        print(f"Divergence {div*100:.1f}% -> true={dist_true}, WFA={dist_wfa}, Xdrop={dist_x}, Banded={dist_band}")

    return results


def plot_results(results):
    div = np.array(results["divergence"])
    true_dist = np.array(results["edit_distance"])

    # Accuracy: absolute error normalized
    def accuracy(pred):
        return 1 - np.abs(np.array(pred) - true_dist) / np.maximum(true_dist, 1)

    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.title("Cells Explored vs Divergence")
    plt.plot(div, results["cells_edit"], label="Full DP")
    plt.plot(div, results["cells_banded"], label="Banded")
    plt.plot(div, results["cells_xdrop"], label="X-drop")
    plt.plot(div, results["cells_wfa"], label="WFA")
    plt.xlabel("Divergence (%)")
    plt.ylabel("Cells Explored (proxy)")
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.title("Accuracy vs Divergence")
    plt.plot(div, accuracy(results["banded"]), label="Banded")
    plt.plot(div, accuracy(results["xdrop"]), label="X-drop")
    plt.plot(div, accuracy(results["wfa"]), label="WFA")
    plt.xlabel("Divergence (%)")
    plt.ylabel("Accuracy (1 - error/true)")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    random.seed(42)
    results = compare_methods(seq_len=1000)
    plot_results(results)
