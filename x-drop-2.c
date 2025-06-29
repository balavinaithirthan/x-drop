
void main() {
    int max_score = 0;
    int m = 100;
    int n = 100;
    int var_lo = 10000;
    int var_hi = -10000;
    int x_drop = 1000;
        for (int d = 0; d < ((m + n) + 1); ++d) {
            for (int i = std::max((((-var_lo) / 2) + (d / 2)), std::max(0, (d - n))); i < std::max((((-var_hi) / 2) + (d / 2)), std::min(m, d)); ++i) {
                if ((S[d, i] == 0)) {
                    S[d, i] = 10;
                    if ((S[d, i] < (max_score - x_drop))) {
                        S[d, i] = -10000;
                    }
                }
            }
            int last_good = 0;
            for (int i = std::max((((-c2) / 2) + (d / 2)), std::min(m, d)); i < std::min((((-c1) / 2) + (d / 2)), std::max(0, (d - n))); ++i) {
                if ((S[d, i] != -10000)) {
                    last_good = i;
                    break;
                }
            }
            if ((last_good != max((((-c2) / 2) + (d / 2)), min(m, d)))) {
                c2 = max(((d - last_good) - last_good), c2);
            }
            int last_good = 0;
            for (int i = std::min((((-c1) / 2) + (d / 2)), std::max(0, (d - n))); i < std::max((((-c2) / 2) + (d / 2)), std::min(m, d)); ++i) {
                if ((S[d, i] != -inf)) {
                    last_good = i;
                    break;
                }
            }
            if ((last_good != min((((-c1) / 2) + (d / 2)), max(0, (d - n))))) {
                c2 = max(((d - last_good) - last_good), c2);
            }
        }
    }
}