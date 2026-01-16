#include "implicit_scheme.h"
#include <iostream>
#include <algorithm>
#include <cmath>


namespace implicit_space
{

    const int N = 100;
    const int M = 2000;

    double implicit_scheme(double S_0, double K, double T, double sigma, double r, double D, int cp)
    {
        const double S_max = 2 * K;
        const double ds = S_max / N;
        const double dt = T / M;

        Matrix_double FDM_grid = define_grid(cp, N, M, K, sigma, r, D, ds, dt, S_max, S_0);

        double V = interpolate(FDM_grid, N, S_max, S_0, ds);

        return V;
    }

    Matrix_double define_grid(int cp, int N, int M, double K, double sigma, double r, double D, double ds, double dt, double S_max, double S_0)
    {

        Matrix_double FDM_grid(N + 1, std::vector<double>(M + 1));

        double Nds;
        double Nds_K;
        double V_N_M;

        int m = M;

        switch (cp)
        {
        case 1:
            for (int n = 0; n <= N; n++)
            {
                Nds = (N - n) * ds;
                Nds_K = Nds - K;
                V_N_M = std::max(Nds_K, 0.0);
                FDM_grid[n][m] = V_N_M;
            }
            break;
        case 2:
            for (int n = 0; n <= N; n++)
            {
                Nds = (N - n) * ds;
                Nds_K = K - Nds;
                V_N_M = std::max(Nds_K, 0.0);
                FDM_grid[n][m] = V_N_M;
            }
            break;
        }

        Matrix_double FDM_M(N + 1, std::vector<double>(N + 1)); // M Matrix for Mv=w
        std::vector<double> FDM_w(N + 1);                       // M Matrix for Mv=w
        std::vector<double> FDM_v(N + 1);                       // M Matrix for Mv=w
        double a;                                               // a values for backward marching FDM
        double b;                                               // b values for backward marching FDM
        double c;                                               // c values for backward marching FDM

        // fill out values of a and b for j=0 and j=N
        b = 1 + r * dt;
        c = 0;
        FDM_M[0][0] = b;
        FDM_M[0][1] = c;

        a = N * (r - D) * dt;           // a_hat
        b = 1 - (N * (r - D) - r) * dt; // b_hat
        FDM_M[N][N - 1] = a;
        FDM_M[N][N] = b;

        // fill out M matrix to establish MV(m-1)=V(m)
        for (int j = 1; j <= N - 1; j++)
        {
            a = -0.5 * (sigma * sigma * j * j - j * (r - D)) * dt;
            b = 1 + (r + sigma * sigma * j * j) * dt;
            c = -0.5 * (sigma * sigma * j * j + j * (r - D)) * dt;

            FDM_M[j][j - 1] = a;
            FDM_M[j][j] = b;
            FDM_M[j][j + 1] = c;
        }

        for (int m = M; m > 0; --m)
        {

            // fill in w from the FDM grid
            for (int j = 0; j <= N; ++j)
            {
                FDM_w[j] = FDM_grid[N - j][m];
            }

            // use gauss seidel to find v at t = m - 1
            for (int i = 1; i < 1000; i++)
            {
                for (int j = 0; j <= N; j++)
                {
                    if (j == 0)
                    {
                        FDM_v[j] = (1 / FDM_M[j][j]) * (FDM_w[j] - FDM_M[j][j + 1] * FDM_v[j + 1]);
                    }
                    else
                    {
                        if (j == N)
                        {
                            FDM_v[j] = (1 / FDM_M[j][j]) * (FDM_w[j] - FDM_M[j][j - 1] * FDM_v[j - 1]);
                        }
                        else
                        {
                            FDM_v[j] = (1 / FDM_M[j][j]) * (FDM_w[j] - FDM_M[j][j - 1] * FDM_v[j - 1] - FDM_M[j][j + 1] * FDM_v[j + 1]);
                        }
                    }
                }
            }
            // fill in data into the FDM grid
            for (int j = 0; j <= N; j++)
            {
                FDM_grid[j][m - 1] = FDM_v[N - j];
            }
        }
        return FDM_grid;
    }

    double interpolate(const Matrix_double &FDM_grid, int N, double S_max, double S_0, double ds)
    {
        double V_inter = 0.0;

        for (int i = 0; i < N; ++i) // i and i+1 valid
        {
            double V = FDM_grid[i][0];
            double V1 = FDM_grid[i + 1][0];

            double S = S_max - i * ds;
            double S1 = S_max - (i + 1) * ds;

            double dist = std::abs(S - S_0);
            double dist1 = std::abs(S1 - S_0);

            if (dist <= dist1)
            {
                if (dist == 0.0)
                    V_inter = V;
                else
                    V_inter = (V * (S1 - S_0) + V1 * (S_0 - S)) / (S1 - S);
                break;
            }
        }
        return V_inter;
    }

} // namespace implicit_space