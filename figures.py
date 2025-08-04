import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def index_space_transformation():
    # Figure 2: Index‑Space Transformation with Smith‑Waterman annotation
    m, n = 10, 20
    orig = np.ones((m, n))
    trans = np.full((m + n, n), np.nan)
    for i in range(m):
        for j in range(n):
            k = i + j
            trans[k, j] = 1

    fig, axs = plt.subplots(1, 2, figsize=(8, 3))
    axs[0].imshow(orig, origin="lower")
    axs[0].set_title("Original (i,j)")
    axs[1].imshow(trans, origin="lower")
    axs[1].set_title("Transformed (k=i+j, m=j)")

    # Remove ticks
    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])

    # Annotate transformed subplot with SW recurrence mapping
    annotation = (
        "Smith-Waterman reindexed:\n"
        "S[k, m] where k = i + j, m = j\n\n"
        "Original: S[i,j] = max(0, S[i-1,j-1] + s, S[i-1,j] − g, S[i,j-1] − g)\n"
        "Reindexed: S[k,m] = max(\n"
        "  0,\n"
        "  S[k-2,m-1] + s,\n"
        "  S[k-m-1,m] - g,\n"
        "  S[k-m,m-1] - g\n"
        ")"
    )
    axs[1].text(
        0.5,
        -0.3,
        annotation,
        ha="center",
        va="top",
        transform=axs[1].transAxes,
        fontsize=8,
    )

    plt.tight_layout()
    plt.show()

    ######


def true_edit_distance(P, T):
    m, n = len(P), len(T)

    # Create a (m+1) x (n+1) matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if P[i - 1] == T[j - 1]:
                cost = 0
            else:
                cost = 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,  # Deletion
                dp[i][j - 1] + 1,  # Insertion
                dp[i - 1][j - 1] + cost,  # Substitution
            )

    return dp


def plot_edit_distance_3d(P, T):
    """
    Create a 3D plot of the edit distance matrix.

    Args:
        P: First string
        T: Second string
    """
    # Get the edit distance matrix
    dp = true_edit_distance(P, T)
    m, n = len(P), len(T)

    # Create coordinate grids
    i_coords = np.arange(m + 1)
    j_coords = np.arange(n + 1)
    I, J = np.meshgrid(i_coords, j_coords, indexing="ij")

    # Convert dp matrix to numpy array for easier handling
    scores = np.array(dp)

    # Create 3D plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Create surface plot
    surf = ax.plot_surface(I, J, scores, cmap="viridis", alpha=0.8)

    # Add wireframe for better visibility
    ax.plot_wireframe(I, J, scores, color="black", alpha=0.3, linewidth=0.5)

    # Customize the plot
    ax.set_xlabel("i (Position in P)")
    ax.set_ylabel("j (Position in T)")
    ax.set_zlabel("Edit Distance Score")
    ax.set_title(f'Edit Distance Matrix 3D Visualization\nP: "{P}", T: "{T}"')

    # Add colorbar
    fig.colorbar(surf, shrink=0.5, aspect=5)

    # Show the plot
    plt.tight_layout()
    plt.show()

    return fig, ax


def plot_2d_diag(P, T):
    """
    Create a 2D plot of edit distance scores vs diagonal offset (j - i).

    Args:
        P: First string
        T: Second string
    """
    # Get the edit distance matrix
    dp = true_edit_distance(P, T)
    m, n = len(P), len(T)

    # Collect all scores and their diagonal offsets
    diagonal_offsets = []
    scores = []
    positions = []

    for i in range(m + 1):
        for j in range(n + 1):
            diagonal_offset = j - i
            diagonal_offsets.append(diagonal_offset)
            scores.append(dp[i][j])
            positions.append((i, j))

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Scatter plot with color coding based on score
    scatter = ax.scatter(
        diagonal_offsets, scores, c=scores, cmap="viridis", alpha=0.7, s=30
    )

    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label("Edit Distance Score")

    # Customize the plot
    ax.set_xlabel("Diagonal Offset (j - i)")
    ax.set_ylabel("Edit Distance Score")
    ax.set_title(f'Edit Distance Scores vs Diagonal Offset\nP: "{P}", T: "{T}"')
    ax.grid(True, alpha=0.3)

    # Add some annotations for key points
    # Find min and max scores
    min_score = min(scores)
    max_score = max(scores)

    # Annotate a few interesting points
    for i, (offset, score) in enumerate(zip(diagonal_offsets, scores)):
        pos = positions[i]
        if (
            score == max_score
            or (pos[0] == 0 and pos[1] == 0)
            or (pos[0] == m and pos[1] == n)
        ):
            ax.annotate(
                f"({pos[0]},{pos[1]})",
                (offset, score),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                alpha=0.8,
            )

    plt.tight_layout()
    plt.show()

    return fig, ax


def plot_layers(P, T):
    """
    Create a layered plot showing i x j positions for each score value.
    Each subplot represents a "layer" of the 3D plot at a specific score.

    Args:
        P: First string
        T: Second string
    """
    # Get the edit distance matrix
    dp = true_edit_distance(P, T)
    m, n = len(P), len(T)

    # Find all unique scores and group positions by score
    score_positions = {}
    all_scores = []

    for i in range(m + 1):
        for j in range(n + 1):
            score = dp[i][j]
            all_scores.append(score)
            if score not in score_positions:
                score_positions[score] = []
            score_positions[score].append((i, j))

    # Sort scores for consistent ordering
    unique_scores = sorted(score_positions.keys())
    max_score = max(unique_scores)

    # Calculate subplot layout
    n_scores = len(unique_scores)
    cols = min(4, n_scores)  # Max 4 columns
    rows = (n_scores + cols - 1) // cols  # Ceiling division

    # Create subplots
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows))
    if n_scores == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.reshape(1, -1)

    # Flatten axes for easier indexing
    axes_flat = axes.flatten() if n_scores > 1 else axes

    for idx, score in enumerate(unique_scores):
        if idx >= len(axes_flat):
            break

        ax = axes_flat[idx]

        # Create a grid to show positions with this score
        grid = np.zeros((m + 1, n + 1))

        # Mark positions with this score
        for i, j in score_positions[score]:
            grid[i, j] = 1

        # Plot as heatmap
        im = ax.imshow(grid, cmap="Blues", origin="lower", alpha=0.8)

        # Add position labels
        for i, j in score_positions[score]:
            ax.text(
                j, i, f"({i},{j})", ha="center", va="center", fontsize=8, weight="bold"
            )

        # Customize subplot
        ax.set_title(f"Score = {score}\n({len(score_positions[score])} positions)")
        ax.set_xlabel("j (Position in T)")
        ax.set_ylabel("i (Position in P)")

        # Set ticks
        ax.set_xticks(range(n + 1))
        ax.set_yticks(range(m + 1))
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n_scores, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    # Overall title
    fig.suptitle(
        f'Edit Distance Matrix Layers\nP: "{P}", T: "{T}"\nEach layer shows positions with the same score',
        fontsize=14,
        y=0.98,
    )

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)  # Make room for suptitle
    plt.show()

    return fig, axes


def demo_3d_plot():
    """
    Demonstrate the 3D edit distance visualization with example strings.
    """
    # Example strings
    P = "KITTEN"
    T = "SITTING"

    print(f"Creating 3D plot for strings:")
    print(f"P: '{P}'")
    print(f"T: '{T}'")

    plot_edit_distance_3d(P, T)


def demo_2d_diag():
    """
    Demonstrate the 2D diagonal plot with example strings.
    """
    # Example strings
    P = "KITTEN"
    T = "SITTING"

    print(f"Creating 2D diagonal plot for strings:")
    print(f"P: '{P}'")
    print(f"T: '{T}'")

    plot_2d_diag(P, T)


def demo_layers():
    """
    Demonstrate the layers plot with example strings.
    """
    # Example strings
    P = "KITTEN"
    T = "SITTING"

    print(f"Creating layers plot for strings:")
    print(f"P: '{P}'")
    print(f"T: '{T}'")

    plot_layers(P, T)


def demo_all_plots():
    """
    Demonstrate all three types of plots: 3D, 2D diagonal, and layers.
    """
    # Example strings
    P = "KITTEN"
    T = "SITTING"

    print(f"Creating all plots for strings:")
    print(f"P: '{P}'")
    print(f"T: '{T}'")
    print()

    print("Generating 3D plot...")
    plot_edit_distance_3d(P, T)

    print("Generating 2D diagonal plot...")
    plot_2d_diag(P, T)

    print("Generating layers plot...")
    plot_layers(P, T)


if __name__ == "__main__":
    # Run the demo when script is executed directly
    demo_layers()  # Show the new layers plot


# create 3d representation, i, j, score
