#include <cstring>
#include <string>
#include <iostream>
#include <random>

using namespace std;
typedef unsigned int uint;

class Solution {
public:
    int lcs(std::string X, std::string Y){
        if ( X.size() < Y.size() ) swap(X,Y) ;
        int m = X.size() , n = Y.size();
        if (m == 0 || n == 0) return 0;
        int w = (m + 31) >> 5;
        uint S[256][w];
        memset(S, 0, sizeof(S));
        uint set = 1;
        for (int i = 0, j = 0; i < m; ++i) {
            S[X[i]][j] |= set;
            set <<= 1;
            if (!set) set++,j++;
        }
        uint L[w];
        memset(L, 0, sizeof(L));
        for (int i = 0; i < n; ++i) {
            uint b1 = 1;
            uint b2 = 0;
            for (int j = 0; j < w; ++j) {
                uint U  = L[j] | S[Y[i]][j];
                uint V  = U - (L[j] << 1 | b1+b2);
                b1 = L[j] >> 31;
                b2 = (V > U);
                L[j] = U & (~V);
            }
        }
        int res = 0;
        for (int i = 0; i < w; ++i)
            res += __builtin_popcount(L[i]);
        return res;
    }

    int longestCommonSubsequence(string text1, string text2) {
        return lcs(text1, text2);
    }
};

string generate_random_binary_string(int n, std::mt19937& gen) {
    string result;
    uniform_int_distribution<> dist(0, 1);
    for (int i = 0; i < n; ++i) {
        result += dist(gen) + '0'; // convert 0/1 to '0'/'1'
    }
    return result;
}

int main(int argc, char *argv[]) {
    int n = 10000;
    std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr)));
    int lcs;
    double s;

    for (int i = 0; i < 10000; i++) {
        std::string str1 = generate_random_binary_string(n, gen);
        std::string str2 = generate_random_binary_string(n, gen);

        lcs = Solution().longestCommonSubsequence(str1, str2);
        s += (double) lcs / n;
    }

    std::cout << "Approximation: " << s / 10000 << std::endl;
    return 0;
}