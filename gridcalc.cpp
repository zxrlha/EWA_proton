#include "gridcalc.hpp"
#include <LHAPDF/LHAPDF.h>

using namespace std;

typedef double (*func)(double&);

extern "C" double dgauss_(func f, double& a, double& b, double& eps);
extern "C" double ewa_wx_bypid_(double& x, double& q2max, int& pol, int& ppid);
extern "C" double ewa_zx_bypid_(double& x, double& q2max, int& pol, int& ppid);

namespace
{
    double x_internal;
    double Q_internal;
    int pol_internal;
    int type_internal;
}

double integrand(double& t)
{
	//PDF name in LHAPDF
    static auto pdf = LHAPDF::mkPDF("PDF4LHC15_nlo_100_pdfas", 0);
    double Q = Q_internal;
    double Q2 = Q * Q;
    double y = pow(x_internal, t);
    double z = pow(x_internal, 1 - t);
    //For W+: positive charged (anti-)quarks: dbar,u,sbar,c,bbar
    std::vector<double> Qp{-1, 2, -3, 4, -5};
    //For W-: negative charged (anti-)quarks: d,ubar,s,cbar,b
    std::vector<double> Qm{1, -2, 3, -4, 5};
    //For Z: all quarks
    std::vector<double> Qpm{1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
    std::vector<double> Qv;
    switch (type_internal)
    {
        //W+
        case 24:
            Qv = Qp;
            break;
        //W-
        case -24:
            Qv = Qm;
            break;
        //Z
        case 23:
            Qv = Qpm;
            break;
        default:
            exit(1);
    };
    double res = 0;
    for (auto& i : Qp)
    {
        double v1 = pdf->xfxQ(i, y, Q) / y;
        int j = i;
        double v2;
        switch (type_internal)
        {
            case 24:
            case -24:
                v2 = ewa_wx_bypid_(z, Q2, pol_internal, j);
                break;
            case 23:
                v2 = ewa_zx_bypid_(z, Q2, pol_internal, j);
                break;
            default:
                exit(1);
        }
        res += v1 * v2;
    }
    return res * (-log(x_internal));
}

double fx(double x, double Q, int pol, int type)
{
    ::x_internal = x;
    ::Q_internal = Q;
    ::pol_internal = pol;
    ::type_internal = type;
    double a = 0.0;
    double b = 1.0;
    double eps = 1e-6;
    double res = dgauss_(integrand, a, b, eps);
    return res;
}
