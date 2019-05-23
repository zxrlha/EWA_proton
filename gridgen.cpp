#include <iostream>
#include <cstdio>
#include <LHAPDF/LHAPDF.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include "gridcalc.hpp"

using namespace std;

void write_a_function(std::ostream& out, const std::vector<double>& vy, const std::vector<double>& vz, const std::string& name, int pol, int type)
{
    assert(vy.size() % 4 == 0);
    assert(vz.size() % 4 == 0);
    std::vector<std::vector<double>> vg(vy.size(), std::vector<double>(vz.size()));
    for (int i = 0; i < vy.size(); ++i)
    {
        for (int j = 0; j < vz.size(); ++j)
        {
            vg[i][j] = exp(vy[i]) * fx(exp(vy[i]), exp(vz[j]), pol, type) / vy[i] / vy[i];
        }
    }
    fmt::print(out, "      function {}(y,z)\n", name);
    fmt::print(out, "      implicit none\n");
    fmt::print(out, "      real*8 {},y,z\n", name);
    fmt::print(out, "      integer narg,nny,nnz\n");
    fmt::print(out, "      parameter (narg=2)\n");
    fmt::print(out, "      parameter (nny={})\n", vy.size());
    fmt::print(out, "      parameter (nnz={})\n", vz.size());
    fmt::print(out, "      integer iny,inz,nent(narg)\n");
    fmt::print(out, "      real*8 tmp,dfint,ymap,zmap\n");
    fmt::print(out, "      real*8 arg(narg),ent(nny+nnz)\n");
    fmt::print(out, "      real*8 yv(nny),zv(nnz),gridv(nny,nnz)\n");
    fmt::print(out, "      logical firsttime\n");
    fmt::print(out, "      external dfint,ymap,zmap\n");
    fmt::print(out, "      data yv/\n");
    for (int i = 0; i < vy.size() / 4; ++i)
    {
        fmt::print(out, "     #{:15.7e},{:15.7e},{:15.7e},{:15.7e}{}\n", vy[4 * i], vy[4 * i + 1], vy[4 * i + 2], vy[4 * i + 3], 4 * (i + 1) == vy.size() ? '/' : ',');
    }
    fmt::print(out, "      data zv/\n");
    for (int i = 0; i < vz.size() / 4; ++i)
    {
        fmt::print(out, "     #{:15.7e},{:15.7e},{:15.7e},{:15.7e}{}\n", vz[4 * i], vz[4 * i + 1], vz[4 * i + 2], vz[4 * i + 3], 4 * (i + 1) == vz.size() ? '/' : ',');
    }
    for (int j = 0; j < vz.size(); ++j)
    {
        fmt::print(out, "      data (gridv(iny, {}),iny=1,{})/\n", j + 1, vy.size());
        for (int i = 0; i < vy.size() / 4; ++i)
        {
            fmt::print(out, "     #{:15.7e},{:15.7e},{:15.7e},{:15.7e}{}\n", vg[4 * i][j], vg[4 * i + 1][j], vg[4 * i + 2][j], vg[4 * i + 3][j], 4 * (i + 1) == vg.size() ? '/' : ',');
        }
    }
    fmt::print(out, "      data firsttime/.true./\n");
    fmt::print(out, "      save\n");
    fmt::print(out, "      if(firsttime)then\n");
    fmt::print(out, "        firsttime=.false.\n");
    fmt::print(out, "        nent(1)=nny\n");
    fmt::print(out, "        nent(2)=nnz\n");
    fmt::print(out, "        do iny=1,nny\n");
    fmt::print(out, "          ent(iny)=(yv(iny))\n");
    fmt::print(out, "        enddo\n");
    fmt::print(out, "        do inz=1,nnz\n");
    fmt::print(out, "          ent(nny+inz)=(zv(inz))\n");
    fmt::print(out, "        enddo\n");
    fmt::print(out, "      endif\n");
    fmt::print(out, "      arg(1)=(y)\n");
    fmt::print(out, "      arg(2)=(z)\n");
    fmt::print(out, "      tmp=dfint(narg,arg,nent,ent,gridv)\n");
    fmt::print(out, "      {}=tmp\n", name);
    fmt::print(out, "      return\n");
    fmt::print(out, "      end\n");
}

void write_output(std::ostream& out, const std::string& name, int type)
{
    //ranges
    double Qmin = 80.385;
    double Qmax = 1000000;
    double xmin = 1e-7;
    double xmax = 1 - 1e-7;
    double ylow = log(xmin);
    double yupp = log(xmax);
    double zlow = log(Qmin);
    double zupp = log(Qmax);
    const std::string spaces = "      ";
    fmt::print(out, "      function {}(x,Q2,pol)\n", name);
    fmt::print(out, "      implicit none\n");
    fmt::print(out, "      real*8 {}\n", name);
    fmt::print(out, "      real*8 x,Q2\n");
    fmt::print(out, "      integer pol\n");
    fmt::print(out, "      real*8 tmp\n");
    fmt::print(out, "      real*8 y,z\n");
    fmt::print(out, "      real*8 ylow,yupp,zlow,zupp\n");
    fmt::print(out, "      parameter (ylow={},yupp={})\n", ylow, yupp);
    fmt::print(out, "      parameter (zlow={},zupp={})\n", zlow, zupp);
    fmt::print(out, "      real*8 {}_1p\n", name);
    fmt::print(out, "      real*8 {}_1m\n", name);
    fmt::print(out, "      real*8 {}_0\n", name);
    fmt::print(out, "      y=log(x)\n");
    fmt::print(out, "      z=0.5*log(Q2)\n");
    fmt::print(out, "      if(pol.eq.1)then\n");
    fmt::print(out, "        tmp={}_1p(y,z)\n", name);
    fmt::print(out, "      else if(pol.eq.-1)then\n");
    fmt::print(out, "        tmp={}_1m(y,z)\n", name);
    fmt::print(out, "      else\n");
    fmt::print(out, "        tmp={}_0(y,z)\n", name);
    fmt::print(out, "      endif\n");
    fmt::print(out, "      {}=tmp/x*y*y\n", name);
    fmt::print(out, "      end\n");
    //preparing points
    std::vector<double> vy;
    std::vector<double> vz;
    //equal space
    int nny = 32;
    for (int i = 0; i < nny; ++i)
    {
        double tmp = yupp - (yupp - ylow) * (nny - i - 1) * (nny - i - 1) / (nny - 1) / (nny - 1);
        vy.push_back(tmp);
    }
    int nnz = 8;
    for (int j = 0; j < nnz; ++j)
    {
        double tmp = zlow + (zupp - zlow) * j / (nnz - 1);
        vz.push_back(tmp);
    }
    write_a_function(out, vy, vz, name + "_1p", 1, type);
    write_a_function(out, vy, vz, name + "_1m", -1, type);
    write_a_function(out, vy, vz, name + "_0", 0, type);
}


int main()
{
    std::ofstream out("ewa_pp.f");
    write_output(out, "ewa_wp_p", 24);
    write_output(out, "ewa_wm_p", -24);
    write_output(out, "ewa_z_p", 23);
    return 0;
}
