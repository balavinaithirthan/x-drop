def wavefront_edit_distance(seq1, seq2):
    n, m = len(seq1), len(seq2)
    max_score = n + m  # worst-case score (all inserts or deletes)

    # Each wavefront at score s maps diagonal k to offset (furthest i along diagonal)
    wavefronts = [{} for _ in range(max_score + 1)]

    # Base case: alignment starts at (0, 0) → diagonal k=0, offset i=0
    wavefronts[0][0] = 0

    for s in range(max_score + 1):
        wf = wavefronts[s]
        for k in wf:
            i = wf[k]
            j = i - k

            # extend matches as far as possible
            while i < n and j < m and seq1[i] == seq2[j]:
                i += 1
                j += 1
            wf[k] = i  # update furthest reaching point

            # check if reached bottom-right
            if i >= n and j >= m:
                return s

        # prepare next wavefront (s+1)
        if s + 1 <= max_score:
            next_wf = wavefronts[s + 1]
            for k in wf:
                i = wf[k]

                # Insertion → (i, j+1) → diagonal k-1
                next_wf[k - 1] = max(next_wf.get(k - 1, 0), i)

                # Deletion → (i+1, j) → diagonal k+1
                next_wf[k + 1] = max(next_wf.get(k + 1, 0), i + 1)

                # Mismatch → (i+1, j+1) → same diagonal
                next_wf[k] = max(next_wf.get(k, 0), i + 1)

    return -1  # should never reach here if sequences are finite


# Example usage
A = "GATTACA"
B = "GCATGCU"
print("Edit distance:", wavefront_edit_distance(A, B))
