import random.islpy as isl
import networkx as nx
import matplotlib.pyplot as plt

# Dimensions
m, n = 3, 3  # small example for simplicity

# Define context
ctx = isl.DEFAULT_CONTEXT

# Define iteration domain: 1 <= i <= m, 1 <= j <= n
domain_str = f"[m,n] -> {{ S[i,j] : 1 <= i <= {m} and 1 <= j <= {n} }}"
domain = isl.Set.read_from_str(ctx, domain_str)

# Define dependency relations: each S[i,j] depends on 3 neighbors
deps = []

dep_strs = [
    f"[m,n] -> {{ S[i,j] -> S[i-1,j-1] : 1 <= i <= {m} and 1 <= j <= {n} }}",
    f"[m,n] -> {{ S[i,j] -> S[i-1,j]   : 1 <= i <= {m} and 1 <= j <= {n} }}",
    f"[m,n] -> {{ S[i,j] -> S[i,j-1]   : 1 <= i <= {m} and 1 <= j <= {n} }}",
]

for s in dep_strs:
    deps.append(isl.Map.read_from_str(ctx, s))

# Combine all dependencies into one map
dep_map = deps[0]
for d in deps[1:]:
    dep_map = dep_map.union(d)

# Create DAG as networkx DiGraph
G = nx.DiGraph()

# Enumerate all source and destination pairs from the dependency map
dep_map = dep_map.coalesce()
for b in dep_map.get_basic_maps():
    src_space = b.get_space().domain()
    dst_space = b.get_space().range()

    # Get bounding boxes from source and destination
    src_box = b.domain().sample_point()
    dst_box = b.range().sample_point()

    # Scan over the domain (for small m, n)
    points = domain.get_basic_set_list()
    for i in range(points.n_basic_set()):
        bset = points.get_basic_set(i)
        bset_points = bset.sample_point()
        if not bset_points:
            continue

    # Manually enumerate domain points for this small m x n case
for i in range(1, m + 1):
    for j in range(1, n + 1):
        src = (i, j)
        if i > 1 and j > 1:
            G.add_edge((i - 1, j - 1), src)
        if i > 1:
            G.add_edge((i - 1, j), src)
        if j > 1:
            G.add_edge((i, j - 1), src)

# Print all edges
print("Dependency edges:")
for u, v in G.edges:
    print(f"{u} -> {v}")

# Optional: visualize
pos = {(i, j): (j, -i) for i in range(1, m + 1) for j in range(1, n + 1)}
nx.draw(G, pos, with_labels=True, node_size=700, node_color="lightblue", arrows=True)
plt.title("Needleman-Wunsch Dependency DAG")
plt.show()
