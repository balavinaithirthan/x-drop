def edit_distance_wavefront(a: str, b: str) -> int:
    """
    Wavefront algorithm expressed with a two-level 'frontier':
      • frontier_on  – set  of diagonals that are active this round
      • frontier_off – dict {k: i_max} giving the furthest row i
                       reached on each active diagonal k (i-j = k)

    The inner loop literally follows:

        for diag in frontier_on:
            (i, j) = (off[k], off[k] - k)
            while a[i] == b[j]: …      # extend
            if reached end: return d
            turn on / update neighbours in next_frontier
    """
    m, n = len(a), len(b)
    if m == 0:
        return n
    if n == 0:
        return m

    # --- d = 0 ----------------------------------------------------------
    frontier_on = {0}  # which diagonals are “on” right now
    frontier_off = {0: 0}  # furthest i on each on-diagonal
    max_d = m + n

    for d in range(max_d + 1):  # ← ➊ “for diag …” outer layer
        next_on, next_off = set(), {}

        for k in list(frontier_on):  # ← ➋ “for on_bit,(i,j) in frontier”
            i = frontier_off[k]
            j = i - k

            # extend along exact matches (while a[i] == b[j] …)
            while i < m and j < n and a[i] == b[j]:
                i += 1
                j += 1

            # reached goal?  (“if end … terminate”)
            if i == m and j == n:
                return d

            # --- generate neighbours, store in next_frontier -------------
            # deletion  → diagonal k+1, position (i+1, j)
            if i + 1 <= m:
                next_off[k + 1] = max(next_off.get(k + 1, -1), i + 1)
                next_on.add(k + 1)  # “frontier[i+1] = … & turn on”

            # insertion → diagonal k-1, position (i,   j+1)
            if j + 1 <= n:
                next_off[k - 1] = max(next_off.get(k - 1, -1), i)
                next_on.add(k - 1)  # “frontier[i-1] = … & turn on”

            # substitution → diagonal k, position (i+1, j+1)
            if i + 1 <= m and j + 1 <= n:
                next_off[k] = max(next_off.get(k, -1), i + 1)
                next_on.add(k)  # same diag remains on

        # hand-off to the next edit-distance layer
        frontier_on, frontier_off = next_on, next_off

    raise RuntimeError("Unreachable: max_d = m+n is an upper bound")
