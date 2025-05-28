#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <random>

using namespace std;
typedef unsigned int uint;

using Word = uint64_t;
constexpr int W = 64;

/* Сложение с переносом через __int128 (portable, ≥ C++17). */
inline Word add64(Word a, Word b, Word &carry) {
    __uint128_t sum = static_cast<__uint128_t>(a) + b + carry;
    carry = static_cast<Word>(sum >> W);
    return static_cast<Word>(sum);
}

/* CIPR для любой m (память ↦ O(k), k = ceil((m+1)/64)). */
size_t lcs_cipr(const string &x, const string &y) {
    const size_t m = x.size();
    if (m == 0) return 0;

    const size_t k = (m + 1 + W - 1) / W; // words for m bits + sentinel
    const size_t s_word = m >> 6; // word of sentinel
    const int s_off = m & 63; // bit offset
    const Word LOWER = (s_off == 63 ? ~0ULL : ((1ULL << s_off) - 1));

    /* Таблица масок M[c][i]. */
    array<vector<Word>, 256> M;
    for (auto &v: M) v.assign(k, 0);

    for (size_t i = 0; i < m; ++i) {
        size_t w = i >> 6;
        int b = i & 63;
        M[static_cast<unsigned char>(x[i])][w] |= (1ULL << b);
    }

    /* ΔL'₀ = …111 (m единиц), sentinel = 0. */
    vector<Word> S(k, ~0ULL);
    S[s_word] &= LOWER;

    size_t len = 0;

    for (unsigned char c: y) {
        Word carry = 0;
        vector<Word> v(k);

        /* v = S + (S & M) */
        for (size_t i = 0; i < k; ++i) {
            Word u = S[i] & M[c][i];
            v[i] = add64(S[i], u, carry);
        }

        /* sentinel set?  → LCS+1 */
        if (v[s_word] & (1ULL << s_off)) ++len;

        /* S ← v | (S & ~M) */
        for (size_t i = 0; i < k; ++i)
            S[i] = v[i] | (S[i] & ~M[c][i]);

        /* zero sentinel + bits above m */
        S[s_word] &= LOWER;
        for (size_t i = s_word + 1; i < k; ++i) S[i] = 0;
    }
    return len;
}

string generate_random_binary_string(int n, mt19937& gen) {
    string result;
    uniform_int_distribution<> dist(0, 1);
    for (int i = 0; i < n; ++i) {
        result += dist(gen) + '0'; // convert 0/1 to '0'/'1'
    }
    return result;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 100;
    mt19937 gen(static_cast<unsigned int>(time(nullptr)));
    int lcs;
    double s;

    for (int i = 0; i < 1e4; i++) {
        string str1 = generate_random_binary_string(n, gen);
        string str2 = generate_random_binary_string(n, gen);

        lcs = lcs_cipr(str1, str2);
        s += (double) lcs / n;
    }

    cout << "Approximation: " << s / 1e4 << endl;
}
