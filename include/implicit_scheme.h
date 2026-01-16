#pragma once

#include <vector>

namespace implicit_space
{
    using Matrix_double = std::vector<std::vector<double>>;

    double implicit_scheme(double S_0, double K, double T,
                           double sigma, double r, double D, int cp);

    Matrix_double define_grid(int cp, int N, int M, double K, double sigma,
                              double r, double D, double ds, double dt,
                              double S_max, double S_0);

    double interpolate(const Matrix_double &FDM_grid,
                       int N, double S_max, double S_0, double ds);
} // namespace implicit_space