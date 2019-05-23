	double precision pi, pi2, sqrt2, eps
	double precision ewa_mw2, ewa_mz2, ewa_gf
	double precision ewa_qu,ewa_qd,ewa_qv,ewa_ql,ewa_t3u,ewa_t3d
	double precision ewa_cw2,ewa_sw2,ewa_gW2
	double precision ewa_gV,ewa_gVu,ewa_gVd,ewa_gVv,ewa_gVl
	double precision ewa_gA,ewa_gAu,ewa_gAd,ewa_gAv,ewa_gAl
c Math
	parameter (pi=3.141592653589793d0)
	parameter (pi2=9.869604401089358d0)
	parameter (sqrt2=1.4142135623730951d0)
	parameter (eps=1d-16)
c Input EW parameters
	parameter (ewa_mw2 = 80.379d0**2)  ! [GeV]^2 PDG 2018
	parameter (ewa_mz2 = 91.1876d0**2) ! [GeV]^2 PDG 2018
	parameter (ewa_gf  = 1.1663787d-5) ! [GeV]^-2 PDG 2018
	parameter (ewa_qu  =  2d0/3d0)
	parameter (ewa_qd  = -1d0/3d0)
	parameter (ewa_qv  =  0d0)
	parameter (ewa_ql  = -1d0)
	parameter (ewa_t3u =  0.5d0)
	parameter (ewa_t3d = -0.5d0)
c Derived EW parameters
  	parameter (ewa_cw2 = ewa_mw2/ewa_mz2)
	parameter (ewa_sw2 = 1d0 - ewa_mw2/ewa_mz2)
	parameter (ewa_gW2 = sqrt2*4d0*ewa_gf*ewa_mw2)
c Weak Vector Couplings
	parameter (ewa_gV  = 1d0)
	parameter (ewa_gVu = 0.5d0*ewa_t3u - ewa_qu*ewa_sw2)
	parameter (ewa_gVd = 0.5d0*ewa_t3d - ewa_qd*ewa_sw2)
	parameter (ewa_gVv = 0.5d0*ewa_t3u)
	parameter (ewa_gVl = 0.5d0*ewa_t3d - ewa_ql*ewa_sw2)
c Weak Axialvector Couplings
	parameter (ewa_gA  = -1d0)
	parameter (ewa_gAu = -0.5d0*ewa_t3u)
	parameter (ewa_gAd = -0.5d0*ewa_t3d)
	parameter (ewa_gAv = -0.5d0*ewa_t3u)
	parameter (ewa_gAl = -0.5d0*ewa_t3d)
c