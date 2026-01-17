#include <iostream>
#include <vector>
#include "explicit_scheme.h"

namespace explicit_space
{

    double S_max;
    const int N = 100;
    double ds;
    const int M = 2000;
    double dt;

    double explicit_scheme(double S_0, double K, double T, double sigma, double r, double D, int cp)
    {
        S_max = 2 * K;
        ds = S_max / N;
        dt = T / M;

        double V;

        if (is_stable(dt, sigma, N))
        {
            Matrix_double FDM_grid = define_grid(cp, N, M, K, sigma, r, D, ds, dt, S_max, S_0);

            V = interpolate(FDM_grid, S_0);
        }
        else
        {
            std::cout << "The selected scheme is unstable under the input parameters. "
                      << "Please try a different scheme or a different set of parameters." << std::endl
                      << std::endl;
            return -1;
        }

        return V;
    }

    bool is_stable(double dt, double sigma, int N)
    {
        return dt <= 1.0 / (sigma * sigma * N * N);
    }

    Matrix_double define_grid(int cp, int N, int M, double K, double sigma, double r, double D, double ds, double dt, double S_max, double S_0)
    {

        Matrix_double FDM_grid(N + 1, std::vector<double>(M + 1));

        int Nds;
        int Nds_K;
        float V_N_M;

        /* Filling out values at payoff using boundary condition
            V_n^M = max(nds-E,0) for calls and max(E-nds,0) for puts*/
        int m = M;

        switch (cp)
        {
        case 1:
            for (int n = 0; n <= N; n++)
            {
                Nds = (N - n) * ds;
                Nds_K = Nds - K;
                V_N_M = std::max(Nds_K, 0);
                FDM_grid[n][m] = V_N_M;
            }
            break;
        case 2:
            for (int n = 0; n <= N; n++)
            {
                Nds = (N - n) * ds;
                Nds_K = K - Nds;
                V_N_M = std::max(Nds_K, 0);
                FDM_grid[n][m] = V_N_M;
            }
            break;
        }

        float A_N_hat;
        float B_N_hat;

        float A_n;
        float B_n;
        float C_n;

        // Filling out the rest of the grid
        float V_n_m_1;
        for (int m = M; m > 0; m--)
        {
            /* solving at N(S = S_max)
            using V_N^(m-1) = A_hat_N*V_(N-1)^m + B_hat_N*V_N^m
            */
            A_N_hat = -N * (r - D) * dt;
            B_N_hat = 1 + (N * (r - D) - r) * dt;
            FDM_grid[0][m - 1] = A_N_hat * FDM_grid[1][m] + B_N_hat * FDM_grid[0][m];

            /*Fill out the rest of the column using
            V_n^(m-1) = A_n*V_(n-1)^m +B_n etc etc*/
            for (int n = N - 1; n > 0; n--)
            {
                A_n = 0.5 * (sigma * sigma * n * n - n * (r - D)) * dt;
                B_n = 1 - (r + sigma * sigma * n * n) * dt;
                C_n = 0.5 * (sigma * sigma * n * n + n * (r - D)) * dt;

                V_n_m_1 = A_n * FDM_grid[N - n + 1][m] + B_n * FDM_grid[N - n][m] + C_n * FDM_grid[N - n - 1][m];

                FDM_grid[N - n][m - 1] = V_n_m_1;
            }
        }
        return FDM_grid;
    }

    double interpolate(Matrix_double FDM_grid, double S_0)
    {
        /*finding the S nearest to S0 and corresponding V via
        linear interpolation
        */
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
                // std::cout << "V = " << V_inter << std::endl;
                break;
            }
        }
        return V_inter;
        // delete FDM_grid;
    }

} // namespace explicit_scheme
