#include "gridcalc.hpp"
#include <fstream>
#include <fmt/ostream.h>
#include <cmath>

extern "C" double ewa_wp_p_(double& x, double& Q2, int& pol);
extern "C" double ewa_wm_p_(double& x, double& Q2, int& pol);
extern "C" double ewa_z_p_(double& x, double& Q2, int& pol);

void test_x(double xlow, double xupp, double Q, int pol)
{
    int n = 100;
    double ylow = log(xlow);
    double yupp = log(xupp);
    std::ofstream out("comp.dat");
    double Q2 = Q * Q;
    int type = 23;
    for (int i = 0; i < n; ++i)
    {
        double y = ylow + ((yupp - ylow) * i) / (n - 1);
        double x = exp(y);
        double res = fx(x, Q, pol, type);
        double gridres;
        switch (type)
        {
            case 24:
                gridres = ewa_wp_p_(x, Q2, pol);
                break;
            case -24:
                gridres = ewa_wm_p_(x, Q2, pol);
                break;
            case 23:
                gridres = ewa_z_p_(x, Q2, pol);
                break;
        }
        fmt::print(out, "{:.9} {} {}\n", x, res, gridres);
    }
}

int main()
{
    test_x(1e-7, 1 - 1e-7, 1300, -1);
    return 0;
}
