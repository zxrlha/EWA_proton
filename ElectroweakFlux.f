c/* ********************************************************* */
c/*      Effective W Approximation Splitting Functions        */
c/*      Pulled from Barger and Phillips, pg 276              */
c/* ********************************************************* */
c
c List of functions:
c ewa_wX_ByPID,ewa_wX,ewa_wT,ewa_w0
c ewa_zX_ByPID,ewa_zT,ewa_z0
c
c Input for ewa_vX_ByPID:
c x: longitudinal momentum fraction carried by boson v
c q2max: scale^2 of v boson pdf
c pol: polarization of v boson (pol = +1,-1,0)
c ppid: pid of v boson's parent fermion
c/* ********************************************************* */

c/* ********************************************************* */
c return ewa splitting function for w boson by parent PID
      double precision function ewa_wX_ByPID(x,q2max,pol,ppid)
      implicit none
      integer pol,ppid
      double precision x,q2max
      double precision ewa_wX
      external ewa_wX

      include 'ElectroweakFlux.inc'

      if(q2max.lt.ewa_mw2) then
c         write (*,*) 'ERROR: q2max below MW2',q2max,ewa_mw2 
         ewa_wX_ByPID = 0d0
         return
      endif
      if(x.lt.eps .or. x.gt.(1.d0-eps)) then
c         write (*,*) 'ERROR: x out of range',x
         ewa_wX_ByPID =0d0
         return
      endif

      ewa_wX_ByPID = ewa_wX(x,q2max,pol)

      return
      end
c/* ********************************************************* */
c return ewa splitting function for w boson
      double precision function ewa_wX(x,q2max,pol)
      implicit none
      integer pol
      double precision x,q2max
      double precision ewa_wT,ewa_w0
      external ewa_wT,ewa_w0

      if(abs(pol).gt.0) then
         ewa_wX = ewa_wT(x,q2max,pol)
      else
         ewa_wX = ewa_w0(x)
      endif
      return
      end
c/* ********************************************************* */
c return ewa splitting function for z boson by parent PID
      double precision function ewa_zX_ByPID(x,q2max,pol,ppid)
      implicit none
      integer pol,ppid
      double precision x,q2max
      double precision ewa_zT,ewa_z0
      external ewa_zT,ewa_z0

      include 'ElectroweakFlux.inc'

      if(q2max.lt.ewa_mz2) then
c         write (*,*) 'ERROR: q2max below MZ2',q2max,ewa_mz2
         ewa_zX_ByPID = 0d0
         return
      endif
      if(x.lt.eps .or. x.gt.(1.d0-eps))then
c         write (*,*) 'ERROR: x out of range',x
         ewa_zX_ByPID = 0d0
         return
      endif

c is parent a charged lepton?
      if(     abs(ppid).eq.11
     & .or.   abs(ppid).eq.13
     & .or.   abs(ppid).eq.15) then
         if(abs(pol).gt.0) then
            ewa_zX_ByPID = ewa_zT(x,q2max,pol,ewa_gVl,ewa_gAl)
         else
            ewa_zX_ByPID = ewa_z0(x,ewa_gVl,ewa_gAl)
         endif
         return
c is parent a neutrino?
      elseif( abs(ppid).eq.12
     & .or.   abs(ppid).eq.14
     & .or.   abs(ppid).eq.16) then
         if(abs(pol).gt.0) then
            ewa_zX_ByPID = ewa_zT(x,q2max,pol,ewa_gVv,ewa_gAv)
         else
            ewa_zX_ByPID = ewa_z0(x,ewa_gVv,ewa_gAv)
            endif
         return
c is parent an up-type quark?
      elseif( abs(ppid).eq.2
     & .or.   abs(ppid).eq.4
     & .or.   abs(ppid).eq.6) then
         if(abs(pol).gt.0) then
            ewa_zX_ByPID = ewa_zT(x,q2max,pol,ewa_gVu,ewa_gAu)
         else
            ewa_zX_ByPID = ewa_z0(x,ewa_gVu,ewa_gAu)
         endif
         return
c is parent a down-type quark?
      elseif( abs(ppid).eq.1
     & .or.   abs(ppid).eq.3
     & .or.   abs(ppid).eq.5) then
         if(abs(pol).gt.0) then
            ewa_zX_ByPID = ewa_zT(x,q2max,pol,ewa_gVd,ewa_gAd)
         else
            ewa_zX_ByPID = ewa_z0(x,ewa_gVd,ewa_gAd)
         endif
         return
c or something else?
      else
         ewa_zX_ByPID = 0d0
         return
      endif
      end
c/* ********************************************************* */
c w boson splitting function: longitudinal polarization
      double precision function ewa_w0(x)
      implicit none
      double precision x
      double precision coup
      
      include 'ElectroweakFlux.inc'
      
c     P_W(x,lambda=0) = (gW/4pi)**2 (1-x)/x
      coup = ewa_gW2/(16d0*pi2)
      ewa_w0 = coup * (1d0-x)/x
      return
      end
c/* ********************************************************* */
c z boson splitting function: longitudinal polarization
      double precision function ewa_z0(x,gV,gA)
      implicit none
      double precision x,gV,gA
      double precision coup

      include 'ElectroweakFlux.inc'

c     P_Z(x,lambda=0) = (gW2/cw2*4pi2) * (gV2 + gA2) (1-x)/x
      coup = ewa_gW2/(ewa_cw2*4d0*pi2) * (gV**2 + gA**2)
      ewa_z0 = coup * (1d0-x)/x
      return
      end
c/* ********************************************************* */
c w boson splitting function: transverse polarization
      double precision function ewa_wT(x,q2max,pol)
      implicit none
      integer pol
      double precision x,q2max
      double precision coup,numer

      include 'ElectroweakFlux.inc'
      
c     P_W(x,lambda=\pm) = coup * [(gV \mp gA)^2 + (gV \pm gA)^2 (1-x)^2] * log(Q2/MV2)
      coup  = ewa_gW2/(8d0*16d0*pi2)
      numer = (ewa_gV-pol*ewa_gA)**2 
     &      + (ewa_gV+pol*ewa_gA)**2 * (1d0-x)**2
      ewa_wT = coup * numer * log(q2max/ewa_mw2) / x
      return
      end
c/* ********************************************************* */
c z boson splitting function: transverse polarization
      double precision function ewa_zT(x,q2max,pol,gV,gA)
      implicit none
      integer pol
      double precision x,q2max,gV,gA
      double precision coup,numer

      include 'ElectroweakFlux.inc'
      
c     P_Z(x,lambda=\pm) = coup * [(gV \mp gA)^2 + (gV \pm gA)^2 (1-x)^2] * log(Q2/MV2)
      coup = ewa_gW2/(ewa_cw2*16d0*pi2)
      numer = (gV-pol*gA)**2
     &      + (gV+pol*gA)**2 * (1d0-x)**2
      ewa_zT = coup * numer * log(q2max/ewa_mz2) / x
      return
      end
