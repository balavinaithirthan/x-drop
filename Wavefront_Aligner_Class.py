from collections import defaultdict
import sys


class WavefrontAligner:
    def __init__(self, pattern: str, text: str):
        self.pattern = pattern
        self.text = text
        self.plen = len(pattern)
        self.tlen = len(text)
        # wavefronts[score] = { k: offset }
        self.wavefronts = defaultdict(dict)

    def extend_matches(self, wf: dict) -> dict:
        """
        Given a wavefront (dict of k→offset), extend each diagonal as far as
        the two strings match along that diagonal.
        """
        extended = {}
        for k, off in wf.items():
            h = off + k  # text index
            v = off  # pattern index
            # slide along diagonal
            while v < self.plen and h < self.tlen and self.pattern[v] == self.text[h]:
                v += 1
                h += 1
                off += 1
            extended[k] = off
        return extended

    def compute_edit_wavefront(self, prev_wf: dict) -> dict:
        """
        Build the next wavefront (score +1) from prev_wf using
        ins/del/mismatch rules (pure edit distance).
        """
        new_wf = {}
        ks = list(prev_wf.keys())
        lo, hi = min(ks) - 1, max(ks) + 1

        for k in range(lo, hi + 1):
            ins = prev_wf.get(k - 1, -sys.maxsize)  # insertion
            dels = prev_wf.get(k + 1, -sys.maxsize)  # deletion
            mis = prev_wf.get(k, -sys.maxsize) + 1  # mismatch
            best = max(ins, dels, mis)
            # prune out-of-bounds
            h = best + k
            v = best
            if h < 0 or h > self.tlen or v < 0 or v > self.plen:
                continue
            new_wf[k] = best
        return new_wf

    def align(self):
        # Score = 0: initialize wavefront at k=0 offset=0, then extend matches
        self.wavefronts[0] = self.extend_matches({0: 0})
        if self._done(self.wavefronts[0]):
            return 0

        # Main scoring loop
        score = 1
        while True:
            # 1) compute next wavefront at this score
            wf = self.compute_edit_wavefront(self.wavefronts[score - 1])
            # 2) extend matches along each diagonal
            wf = self.extend_matches(wf)
            self.wavefronts[score] = wf

            if self._done(wf):
                return score
            score += 1

    def _done(self, wf: dict) -> bool:
        """
        Check if any diagonal has reached the end of both strings.
        Equivalently, k == tlen - plen and offset == plen
        """
        for k, off in wf.items():
            h = off + k
            v = off
            if h == self.tlen and v == self.plen:
                return True
        return False


# Example usage:
if __name__ == "__main__":
    pat = "ACGT"
    txt = "TTACT"
    aligner = WavefrontAligner(pat, txt)
    distance = aligner.align()
    print(f"Edit distance = {distance}")
