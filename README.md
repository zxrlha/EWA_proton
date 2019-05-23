# Effective W approximation for proton

This program calculates the W/Z PDF inside a proton according to the Effective W approximation(EWA).
It co-evolves the EWA with (anti-)quark inside a proton, and creates a grid
as well as FORTRAN code for obtaining the PDF values from such grid.

## INSTALL
### External dependency:
1. [LHAPDF](https://lhapdf.hepforge.org):
  For accessing the (anti-)quark PDF of a proton.
  Please follow the manual of LHAPDF for installation.
2. [fmt](https://fmt.dev)
  For formatting output code.
  For GNU/Linux distributions, usually it can be installed through the package manager.

With those packages installed, you can simple type
```Shell
make
```
to compile and run the code.
