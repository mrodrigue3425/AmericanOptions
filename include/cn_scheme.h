#pragma once

#include <iostream>
#include <vector>

namespace cr_space
{
    using Matrix_double = std::vector<std::vector<double>>;
    double cr_scheme(double S_0, double K, double T, double sigma, double r, double D, int cp);
}