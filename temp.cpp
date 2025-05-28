#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <thread>

using namespace std;
using namespace chrono;

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

// #define N 6
// #define M (2*N+1)
// #define bpow(n) (1 << n)
// #define mask(n) (bpow(n)-1)

constexpr int N = 5, M = 2 * N + 1;
constexpr int bpow(const int n) { return 1 << n; }
constexpr int mask(const int n) { return bpow(n) - 1; }

// probability parameters for model initialisation

// model B(1/2)
/* constexpr double U = sqrt(2)-1, V = 1 - U;
constexpr double P[4] = { 0.5, 0.5, 0.5, 0.5 };*/

// model B(opt)
const double U = -1 + sqrt(7. / 3) - sqrt(1. / 6 * (23 - 5 * sqrt(21))), V = 1 - U;
// initial particle probability for v-site = 0.585786, u-site complementary

const double P[4] = {
    -8. / 3 + 49. / 6 * U - U * U - 1. / 2 * U * U * U,
    29. / 2 - 51 * U + 75. / 2 * U * U + 9 * U * U * U,
    -2. / 3 + 34. / 3 * U - 19 * U * U - 4 * U * U * U,
    P[0]
};

const double R[4] = {P[0], 1 - (1 - P[1]) * U * U / (V * V), 0, P[3]};

// configuration probabilities: cells (N bits) . sites (2*N+1 bits) . cell (N bits) -> probability;
double prob[bpow(N)][bpow(M)][bpow(N)], new_prob[bpow(N)][bpow(M)][bpow(N)];

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

// bit reversal (widths N and M with complement)
uint32_t bitr_N[bpow(N)], bitrnot_M[bpow(M)];

void init_bitr() {
    for (uint16_t cell = 0; cell < bpow(N); cell++) {
        bitr_N[cell] = 0;
        for (int k = 0; k < N; k++) bitr_N[cell] ^= ((cell >> (N - 1 - k)) & 1) << k;
    }
    for (uint32_t site = 0; site < bpow(M); site++) {
        bitrnot_M[site] = mask(M);
        for (int k = 0; k < M; k++) bitrnot_M[site] ^= ((site >> (M - 1 - k)) & 1) << k;
    }
}

// sorting step
uint32_t step[bpow(M + 1)][bpow(N + 1)];

void init_step() {
    for (uint32_t site = 0; site < bpow(M + 1); site++)
        for (uint16_t cell = 0; cell < bpow(N + 1); cell++) {
            uint32_t site1 = site;

            for (int k = 0; k < N + 1; k++)
                if (((cell >> k) & mask(1)) && (((site >> (2 * k)) & mask(2)) == 2))
                    site1 ^= (mask(2) << (2 * k));

            step[site][cell] = site1;
        }
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

// initialise configuration probabilities from model B
void init_prob() {
    // iterate over site configurations
    for (uint32_t site1 = 0; site1 < bpow(M); site1++) {
        // site configuration probability
        double s = 1;
        for (int k = 0; k < M; k += 2) s *= (((site1 >> k) & mask(1)) ? V : U);
        for (int k = 1; k < M; k += 2) s *= (((site1 >> k) & mask(1)) ? U : V);

        // reverse and forward cell probabilities conditioned on site configuration
        double r[N], p[N];

        for (int k = 0; k < N; k++) r[k] = R[(site1 >> (2 * k)) & mask(2)];
        for (int k = 0; k < N; k++) p[k] = P[(site1 >> (2 * k + 1)) & mask(2)];

        // iterate over reverse and forward cell configurations assigning probabilities
        for (uint16_t cell0 = 0; cell0 < bpow(N); cell0++)
            for (uint16_t cell1 = 0; cell1 < bpow(N); cell1++) {
                prob[cell0][site1][cell1] = s;

                for (int k = 0; k < N; k++) prob[cell0][site1][cell1] *= ((cell0 >> k) & mask(1) ? r[k] : 1 - r[k]);
                for (int k = 0; k < N; k++) prob[cell0][site1][cell1] *= ((cell1 >> k) & mask(1) ? p[k] : 1 - p[k]);
            }
    }
}

// update probability distribution for full configuration
void update_prob_range(const uint16_t from, const uint16_t to) {
    for (uint16_t cell1 = from; cell1 < to; cell1++)
        for (uint16_t cell0 = 0; cell0 < bpow(N); cell0++) {
            const uint16_t cell2 = cell0 ^ cell1 ^ (cell1 >> 1);

            for (uint32_t site1 = 0; site1 < bpow(M); site1++) {
                // check admissibility for cell0/site1 (for performance only)
                const uint32_t site1_ = site1 >> 1;

                if (step[site1_][cell0] == site1_) {
                    // perform cell action at time=1
                    const uint32_t site2 = step[site1][cell1];

                    const uint16_t cell1_ = cell1 >> 1;

                    // current probability of Markov core
                    const double p_core = (prob[cell0][site1_][cell1_] + prob[cell0][site1_ | bpow(M - 1)][cell1_]) * 2;
                    if (p_core == 0) continue;

                    // current probability for right extension of Markov core
                    const double p_right = prob[bitr_N[cell0]][bitrnot_M[site1 & mask(M)]][bitr_N[cell1 & mask(N)]];

                    // update current probability for right extension of Markov core

                    // extension cell = 0, extension site = any: no swap
                    new_prob[cell1][site2][cell2] += p_right / 2;

                    // extension cell = 1
                    if (site1_ & bpow(M - 2))

                    // leftmost core site = 1: no swap
                        new_prob[cell1][site2][cell2 ^ bpow(N - 1)] += p_right / 2;

                    else {
                        const double update_p = p_right / p_core * prob[cell0][site1_][cell1_ | bpow(N - 1)];

                        // leftmost core site = 0, extension site = 0, no swap
                        new_prob[cell1][site2][cell2 ^ bpow(N - 1)] += update_p;

                        // leftmost core site = 0, extension site = 1. swap
                        new_prob[cell1][site2 | bpow(M - 1)][cell2 ^ bpow(N - 1)] += p_right / 2 - update_p;
                    }
                }
            }
        }
}

// run one step of model to update probability distribution
void update_prob() {
    for (uint16_t cell0 = 0; cell0 < bpow(N); cell0++)
        for (uint32_t site1 = 0; site1 < bpow(M); site1++)
            for (uint16_t cell1 = 0; cell1 < bpow(N); cell1++)
                new_prob[cell0][site1][cell1] = 0;

    // loop through cell/site/cell configurations at time=0, 1
    update_prob_range(0, bpow(N));

    swap(new_prob, prob);
}

const int n_th = 16, size_th = bpow(N) / n_th;
thread th[n_th];

void update_prob_par() {
    for (uint16_t cell0 = 0; cell0 < bpow(N); cell0++)
        for (uint32_t site1 = 0; site1 < bpow(M); site1++)
            for (uint16_t cell1 = 0; cell1 < bpow(N); cell1++)
                new_prob[cell0][site1][cell1] = 0;

    // loop through cell/site/cell configurations at time=0, 1
    for (uint16_t i_th = 0; i_th < n_th; i_th++)
        th[i_th] = thread(update_prob_range, i_th * size_th,
                          (i_th + 1) * size_th);
    for (uint16_t i_th = 0; i_th < n_th; i_th++) th[i_th].join();

    swap(new_prob, prob);
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

void print_stat() {
    double v = 0, p_total = 0;
    for (uint16_t cell0 = 0; cell0 < bpow(N); cell0++)
        for (uint32_t site1 = 0; site1 < bpow(M); site1++)
            for (uint16_t cell1 = 0; cell1 < bpow(N); cell1++) {
                p_total += prob[cell0][site1][cell1];
                if (site1 & mask(1)) v += prob[cell0][site1][cell1];
            }

    cout << "u: " << fixed << setw(10) << setprecision(10) << 1 - v << "; ";
    cout << "v: " << fixed << setw(10) << setprecision(10) << v << "; ";
    cout << "gamma: " << fixed << setw(10) << setprecision(10) << 2 * (1 - v) << "; ";
    cout << "total: " << fixed << setw(10) << setprecision(10) << p_total << "\n";
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

// main function
int main(int argc, char *argv[]) {
    int iter_count;
    iter_count = 1500;

    cout << "U: " << fixed << setw(10) << setprecision(10) << U << " ";
    cout << "V: " << fixed << setw(10) << setprecision(10) << V << "\n";

    init_bitr();
    init_step();

    init_prob();
    print_stat();

    for (int iter = 0; iter <= iter_count; iter++) {
        update_prob_par();
        if (iter % 20 == 0) {
            cout << setw(4) << iter << " ";
            print_stat();
        }
    }

    /*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
}


/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
