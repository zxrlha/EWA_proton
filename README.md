# Effective W approximation for proton

This program calculates the W/Z PDF inside a proton according to the Effective W approximation(EWA) (see e.g. [S. Dawson, Nucl. Phys. B249 (1985) 42](https://doi.org/10.1016/0550-3213(85)90038-0)).
It co-evolves the EWA with (anti-)quark inside a proton, and creates a grid
as well as Fortran code for obtaining the PDF values from such grid.

## Installation and Usage
### External dependency:
1. [LHAPDF](https://lhapdf.hepforge.org):
  For accessing the (anti-)quark PDF of a proton.
  Please follow the manual of LHAPDF for installation.
  After installing LHAPDF, please modify the variable "LHAPDFCONFIG" in **Makefile** into the path to "lhapdf-config" if your "lhapdf-config" is not present in $PATH.
2. [fmt](https://fmt.dev)
  For formatting output code.
  For GNU/Linux distributions, usually it can be installed through the package manager.
  If you install it by yourself, please also add corresponding path to headers and libraries in **Makefile*.

With those packages installed, you can simple type
```Shell
make
```
to compile and run the code.
A file called **ewa_pp.f** will be generated and it includes the grid as well as interface to the grid.
Two extra file **dfint.f** and **kerset.f** are necessary to use the grid.
An executable **gridtest** will be generated and it will compare the results from grid interpolation and direct calculation, to validate the grid can reproduce the results to good accuracy.

### Usage:
To use the generated W/Z PDF in your code, just includes the file **dfint.f**, **kerset.f** and the generated **ewa_pp.f** in your code.
It provides the following Fortran routine:
```Fortran
real*8 function ewa_wp_p(x,Q2,pol)
real*8 function ewa_wm_p(x,Q2,pol)
real*8 function ewa_z_p(x,Q2,pol)
```
Which will return the W<sup>+</sup>,W<sup>-</sup>,Z boson PDF correspondingly.
Parameters:

- real*8 x: the energy fraction of the W/Z boson.
- real*8 Q2: the square of factorisation scale Q.
- integer pol: the helicity of the W/Z boson.

## Customize:
Note that you need to type **make** again to re-generate the grid after modifications.

- To modify the range of **x** and **Q** of the generated grid, please open the file **gridgen.cpp* and modify the variable **xmin**,**xmax**,**Qmin**,**Qmax** in the first several lines of function **write_output**. Note that outside such range the grid will also return a value through extrapolation, which is not as reliable as interpolation inside the range but usually acceptable.
- To modify the (anti-)quark PDF adopted, please open the file **gridcalc.cpp** and modify the call to **LHAPDF::mkPDF**. The first argument is the name of PDF in LHAPDF and the second argument is the set number of such PDF. You should have the corresponding PDF installed through LHAPDF.

## Brief explanation for the code:

- **ElectroweakFlux.f** and **ElectroweakFlux.inc** : the EWA approximation from a fermion.
- **gridcalc.cpp** and **gridcalc.hpp** : co-evolve the EWA function with parton PDF inside a proton
- **gridgen.cpp** : generate the grid and Fortran code
- **ewa_pp.f** : the generated Fortran code which includes the grid.
- **dgauss.f** : the one-dimensional integration routine, which utilizes 8-points and 16-points Gauss-Legendre quardrature rules and compare.
- **gridtest.cpp** : comparing the results from grid interpolation and direct calculation.
- **dfint.f** and **kerset.f** : routines for grid interpolation