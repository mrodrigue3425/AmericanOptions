#include "cr_scheme.h"
#include <cmath>

namespace cr_space
{
    double cr_scheme(double S_0, double K, double T, double sigma, double r, double D, int cp)
    {
        float S_max = 2 * K;
        const int N = 200;
        int n;
        float ds = S_max / N;
        const int M = 1000;
        int m;
        float dt = T / M;
    }

}