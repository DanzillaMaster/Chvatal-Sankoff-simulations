#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <omp.h>

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

string generate_random_binary_string(int n, mt19937 &gen) {
    string result;
    uniform_int_distribution<> dist(0, 1);
    for (int i = 0; i < n; ++i) {
        result += dist(gen) + '0'; // convert 0/1 to '0'/'1'
    }
    return result;
}

struct Stats {
    long long n = 0; // число наблюдений
    double m = 0.0; // среднее
    double M2 = 0.0; // сумма квадратов отклонений
    // обновление одним наблюдением
    void update(double x) {
        ++n;
        double delta = x - m;
        m += delta / n;
        M2 += delta * (x - m);
    }
};

/* ---------- объединение двух статистик (формулы Чана) ---------- */
Stats merge(const Stats &a, const Stats &b) {
    if (a.n == 0) return b;
    if (b.n == 0) return a;
    Stats res;
    res.n = a.n + b.n;
    double delta = b.m - a.m;
    res.m = (a.n * a.m + b.n * b.m) / res.n;
    res.M2 = a.M2 + b.M2 + delta * delta * a.n * b.n / res.n;
    return res;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int n = 100'000; // длина строк
    const int N = 100'000; // число симуляций
    const int P = 12; // число потоков (предзадано)

    vector<Stats> local(P); // локальные статистики потоков

#pragma omp parallel num_threads(P)
    {
        int tid = omp_get_thread_num();
        mt19937 gen((unsigned) time(nullptr) + tid * 997); // уникальный seed

        for (int iter = tid; iter < N; iter += P) {
            string s1 = generate_random_binary_string(n, gen);
            string s2 = generate_random_binary_string(n, gen);
            int lcs = lcs_cipr(s1, s2);
            double g = static_cast<double>(lcs) / n; // γ-наблюдение
            local[tid].update(g);
        }
    }

    /* --- агрегируем всё в одну статистику --- */
    Stats global = local[0];
    for (int i = 1; i < P; ++i) global = merge(global, local[i]);

    double variance = global.M2 / (global.n - 1);
    double stddev = sqrt(variance);
    double se = stddev / sqrt(global.n);
    double delta = sqrt(2) / sqrt(n);
    double general_error = 2.575 * se + delta;
    double ci_low = global.m - general_error;
    double ci_high = global.m + general_error;


    cout.setf(ios::fixed);
    cout.precision(6);
    cout << "Approximation (gamma_2): " << global.m << '\n'
            << "Standard deviation:      " << stddev << '\n'
            << "Standard error:          " << se << '\n'
            << "Delta:                   " << delta << '\n'
            << "General error:           " << general_error << '\n'
            << "99% confidence interval: [" << ci_low
            << ", " << ci_high << "]\n";
    return 0;
}
