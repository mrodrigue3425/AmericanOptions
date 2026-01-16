#include <vector>

namespace explicit_space
{

    using Matrix_double = std::vector<std::vector<double>>;

    double explicit_scheme(double S_0, double K, double T, double sigma, double r, double D, int cp);
    bool is_stable(double dt, double sigma, int N);
    Matrix_double define_grid(int cp, int N, int M, double K, double sigma, double r, double D, double ds, double dt, double S_max, double S_0);
    double interpolate(Matrix_double FDM_grid, double S_0);
}
