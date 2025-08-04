# import matplotlib.pyplot as plt

# # Matrix dimensions
# rows, cols = 5, 6

# # Create figure and axes
# fig, ax = plt.subplots(figsize=(cols, rows))

# # Draw the grid and fill in with coordinates
# # for i in range(rows):
# #     for j in range(cols):
# #         ax.text(
# #             j + 0.5, rows - i - 0.5, f"({i},{j})", va="center", ha="center", fontsize=12
# #         )

# # Draw grid lines
# for x in range(cols + 1):
#     ax.plot([x, x], [0, rows], color="black")
# for y in range(rows + 1):
#     ax.plot([0, cols], [y, y], color="black")

# # Set row and column labels
# for i in range(rows):
#     ax.text(
#         -0.5,
#         rows - i - 0.5,
#         f" i : {i}",
#         va="center",
#         ha="right",
#         fontsize=12,
#     )
# for j in range(cols):
#     ax.text(
#         j + 0.5,
#         rows + 0.1,
#         f"j : {j}",
#         va="bottom",
#         ha="center",
#         fontsize=12,
#     )

# # Hide axes
# ax.axis("off")
# ax.set_xlim(-1, cols)
# ax.set_ylim(0, rows + 1)

# plt.tight_layout()
# plt.show()


##### (i,j)


# import matplotlib.pyplot as plt
# import matplotlib.patches as patches

# # New dimensions
# rows, cols = 11, 6  # k in [0,10], m in [0,5]

# # Create figure
# fig, ax = plt.subplots(figsize=(cols, rows))

# # Fill cells using loop bounds: for k in range(rows), for m in range(cols) if valid
# for k in range(rows):
#     m_start = max(0, k - cols + 1)
#     m_end = min(k, cols - 1)
#     for m in range(m_start, m_end + 1):
#         # Color the valid region light blue
#         rect = patches.Rectangle(
#             (m, rows - k - 1),
#             1,
#             1,
#             linewidth=1,
#             edgecolor="black",
#             facecolor="lightblue",
#         )
#         ax.add_patch(rect)
#         ax.text(
#             m + 0.5, rows - k - 0.5, f"({k},{m})", ha="center", va="center", fontsize=9
#         )

# # Draw the grid for all cells
# for i in range(rows):
#     for j in range(cols):
#         ax.plot([j, j + 1], [rows - i - 1, rows - i - 1], color="black")
#         ax.plot([j, j + 1], [rows - i, rows - i], color="black")
#         ax.plot([j, j], [rows - i - 1, rows - i], color="black")
#         ax.plot([j + 1, j + 1], [rows - i - 1, rows - i], color="black")

# # Label rows and columns
# for i in range(rows):
#     ax.text(
#         -0.5,
#         rows - i - 0.5,
#         f"k={i}",
#         va="center",
#         ha="right",
#         fontsize=10,
#     )
# for j in range(cols):
#     ax.text(
#         j + 0.5,
#         rows + 0.1,
#         f"m={j}",
#         va="bottom",
#         ha="center",
#         fontsize=10,
#     )

# # Final formatting
# ax.set_xlim(-1, cols)
# ax.set_ylim(0, rows + 1)
# ax.axis("off")
# plt.tight_layout()
# plt.show()


########## draw dependency graph
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# import matplotlib.patheffects as pe

# # Matrix size
# rows, cols = 2, 2

# # Cell size
# cell_size = 1

# # Coordinates for dependencies
# dependencies = {(1, 1): [(1, 0), (0, 1), (0, 0)]}  # (i,j): [left, top, top-left]

# # Create figure and axes
# fig, ax = plt.subplots(figsize=(4, 4))

# # Draw the grid and label the cells
# for i in range(rows):
#     for j in range(cols):
#         rect = patches.Rectangle(
#             (j, rows - i - 1),
#             cell_size,
#             cell_size,
#             linewidth=1,
#             edgecolor="black",
#             facecolor="white",
#         )
#         ax.add_patch(rect)
#         ax.text(
#             j + 0.5, rows - i - 0.5, f"({i},{j})", va="center", ha="center", fontsize=12
#         )

# # Draw dependency arrows into (1,1)
# target_i, target_j = 1, 1
# target_x, target_y = target_j + 0.5, rows - target_i - 0.5

# for src_i, src_j in dependencies[(1, 1)]:
#     src_x, src_y = src_j + 0.5, rows - src_i - 0.5
#     ax.annotate(
#         "",
#         xy=(target_x, target_y),
#         xytext=(src_x, src_y),
#         arrowprops=dict(arrowstyle="->", color="blue", lw=2),
#         path_effects=[pe.withStroke(linewidth=3, foreground="white")],
#     )

# # Final adjustments
# ax.set_xlim(-0.5, cols + 0.5)
# ax.set_ylim(-0.5, rows + 0.5)
# ax.set_aspect("equal")
# ax.axis("off")
# plt.tight_layout()
# plt.show()


####### dependency graph with k,m

# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# import matplotlib.patheffects as pe

# # Matrix size
# rows, cols = 3, 3

# # Target cell (k, m) = (2, 2)
# target_k, target_m = 2, 2
# dependencies = [
#     (target_k - 2, target_m - 1),  # (0,1)
#     (target_k - 1, target_m),  # (1,2)
#     (target_k - 1, target_m - 1),  # (1,1)
# ]

# # Create figure
# fig, ax = plt.subplots(figsize=(5, 5))

# # Draw cells and labels
# for i in range(rows):
#     for j in range(cols):
#         rect = patches.Rectangle(
#             (j, rows - i - 1), 1, 1, edgecolor="black", facecolor="white", linewidth=1
#         )
#         ax.add_patch(rect)
#         ax.text(
#             j + 0.5, rows - i - 0.5, f"({i},{j})", ha="center", va="center", fontsize=11
#         )

# # Highlight target cell
# ax.add_patch(
#     patches.Rectangle(
#         (target_m, rows - target_k - 1),
#         1,
#         1,
#         edgecolor="red",
#         facecolor="mistyrose",
#         linewidth=2,
#     )
# )

# # Draw arrows from dependencies
# for dep_k, dep_m in dependencies:
#     start_x, start_y = dep_m + 0.5, rows - dep_k - 0.5
#     end_x, end_y = target_m + 0.5, rows - target_k - 0.5
#     ax.annotate(
#         "",
#         xy=(end_x, end_y),
#         xytext=(start_x, start_y),
#         arrowprops=dict(arrowstyle="->", color="blue", lw=2),
#         path_effects=[pe.withStroke(linewidth=3, foreground="white")],
#     )

# # Final layout
# ax.set_xlim(-0.5, cols + 0.5)
# ax.set_ylim(-0.5, rows + 0.5)
# ax.set_aspect("equal")
# ax.axis("off")
# plt.tight_layout()
# plt.show()


####
