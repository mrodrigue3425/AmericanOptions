#include "cn_scheme.h"
#include <cmath>

namespace cn_space
{
    double S_max;
    const int N = 100;
    double ds;
    const int M = 2000;
    double dt;

    double cn_scheme(double S_0, double K, double T, double sigma, double r, double D, int cp)
    {
        S_max = 2 * K;
        ds = S_max / N;
        dt = T / M;

        Matrix_double FDM_grid = define_grid(cp, N, M, K, T, sigma, r, D, ds, dt, S_max, S_0);

        double V = interpolate(FDM_grid, S_0);

        return V;
    }

    Matrix_double define_grid(int cp, int N, int M, double K, double T, double sigma, double r, double D, double ds, double dt, double S_max, double S_0)
    {
        /* Filling out values at payoff using boundary condition
        V_n^M = max(nds-E,0) for calls and max(E-nds,0) for puts*/

        double Nds;
        double Nds_K;
        double V_N_M;

        Matrix_double FDM_grid(N + 1, std::vector<double>(M + 1));

        int m = M;

        switch (cp)
        {
        case 1:
            for (int n = 0; n <= N; n++)
            {
                Nds = n * ds;
                Nds_K = Nds - K;
                V_N_M = std::max(Nds_K, 0.0);
                FDM_grid[n][M] = V_N_M;
            }
            for (int m = 0; m <= M; m++)
            {
                FDM_grid[0][m] = 0;
                FDM_grid[N][m] = N * ds * exp(-D * (T - m * dt)) - K * exp(-r * (T - m * dt));
            }
            break;
        case 2:
            for (int n = 0; n <= N; n++)
            {
                Nds = n * ds;
                Nds_K = K - Nds;
                V_N_M = std::max(Nds_K, 0.0);
                FDM_grid[n][M] = V_N_M;
            }
            for (int m = 0; m <= M; m++)
            {
                FDM_grid[0][m] = K * exp(-r * (T - m * dt));
                FDM_grid[N][m] = 0;
            }
            break;
        }

        double A = 0; // a value for explicit M
        double B = 0; // b value for explicit M
        double C = 0; // c value for explicit M

        std::vector<double> a(N + 1); // a value for implicit M
        std::vector<double> b(N + 1); // b value for implicit M
        std::vector<double> c(N + 1); // c value for implicit M

        std::vector<double> V_exp(N + 1);
        std::vector<double> V_imp(N + 1); // v (m-1) vector for Mv=w

        for (int m = M - 1; m >= 0; m--)
        {
            // dettermine V at t=m-0.5 using explicit scheme
            for (int n = 1; n <= N - 1; n++)
            {
                a[n] = -0.25 * (sigma * sigma * n * n - n * (r - D)) * dt;
                b[n] = 1 + 0.5 * (r + sigma * sigma * n * n) * dt;
                c[n] = -0.25 * (sigma * sigma * n * n + n * (r - D)) * dt;

                A = 0.25 * (sigma * sigma * n * n - n * (r - D)) * dt;
                B = 1 - 0.5 * (r + sigma * sigma * n * n) * dt;
                C = 0.25 * (sigma * sigma * n * n + n * (r - D)) * dt;

                V_exp[n] = A * FDM_grid[n - 1][m + 1] + B * FDM_grid[n][m + 1] + C * FDM_grid[n + 1][m + 1];
            }
            // fill in V_imp from the FDM grid
            for (int n = 0; n <= N; n++)
            {
                V_imp[n] = FDM_grid[n][m + 1];
            }
            // use gauss seidel to find v at t = m - 1
            for (int i = 1; i < 1000; i++)
            {
                for (int n = 1; n <= N - 1; n++)
                {
                    V_imp[n] = (1 / b[n]) * (V_exp[n] - a[n] * V_imp[n - 1] - c[n] * V_imp[n + 1]);
                }
            }
            // fill in data into the FDM grid
            for (int n = 0; n <= N; n++)
            {
                FDM_grid[n][m] = V_imp[n];
            }
        }
        return FDM_grid;
    }

    double interpolate(Matrix_double FDM_grid, double S_0)
    {
        float V_inter;
        float V;
        float V1;
        float S;
        float S1;
        float dist;
        float dist1;
        for (int n = N; n >= 0; n--)
        {
            V = FDM_grid[N - n][0];
            // cout << "V = " << V << endl;
            V1 = FDM_grid[N - n + 1][0];
            S = S_max - (N - n) * ds;
            S1 = S_max - (N - n + 1) * ds;
            dist = abs(S - S_0);
            dist1 = abs(S1 - S_0);
            if (dist <= dist1)
            {
                if (dist == 0)
                {
                    V_inter = V;
                }
                else
                {
                    V_inter = (V * (S1 - S_0) + V1 * (S_0 - S)) / (S1 - S);
                }
                break;
            }
        }
        return V_inter;
    }
}