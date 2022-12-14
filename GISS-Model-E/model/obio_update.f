#include "rundeck_opts.h"

      subroutine obio_update(vrbos,kmax,i,j,nstep)
 
c  Performs updating of biological particles, uses mid-point
c  leap frog method.

      USE obio_dim
      USE obio_com, only: P_tend,obio_deltat,D_tend,C_tend
     .                   ,obio_P,det,car,dp1d,p1d
#ifdef TRACERS_Alkalinity
     .                   ,A_tend,alk1d
#ifdef TOPAZ_params
     .                   ,Ca_tend,ca_det_calc1d
#endif
#endif
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     .                    ,O_tend,o21d,errchk1,errchk2 !@PLdbg added test for NaNs in O2
#endif 
#ifdef TRACERS_abio_O2
     .                   ,Abo_tend,abo21d
#endif
#endif
#ifdef TRACERS_degC
     .                   ,ndegC1d,Ndeg_tend
#endif

      implicit none

      integer :: i,j,k
 
      integer :: nt,kmax,nstep !@PLdb nstep added for test for NaNs
      real    :: Pnew,Dnew,Cnew,Anew,Canew
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     .           ,O2new
#endif
#ifdef TRACERS_abio_O2
     .           ,Abo2new
#endif
#endif
#ifdef TRACERS_degC
     .           ,Ndegnew !@PL
#endif
      logical :: vrbos
 
c  Loop to update
c   indexes are mixed because new H has already been computed
c   in update.F, but P has not been updated yet

      do 1000 k = 1,kmax

        do nt = 1,ntyp
         Pnew = (obio_P(k ,nt) +  P_tend(k,nt)*obio_deltat)
         obio_P(k,nt) = max(0.d0,pnew)
        enddo
 
        do nt = 1,ndet
         Dnew = (det(k,nt)     +  D_tend(k,nt)*obio_deltat)
          det(k,nt) = max(0.d0,Dnew)
        enddo
#ifdef TRACERS_degC
         Ndegnew = (ndegC1d(k) + Ndeg_tend(k)*obio_deltat) !@PL
          ndegC1d(k) = max(0.d0,Ndegnew)
#endif
        do nt = 1,ncar
         Cnew = (car(k,nt) +  C_tend(k,nt)*obio_deltat)
         car(k,nt) = max(0.d0,Cnew) 
        enddo

#ifdef TRACERS_Alkalinity
         nt=1
         Anew = (alk1d(k) +  A_tend(k)*obio_deltat)
         alk1d(k) = max(0.d0,Anew)
#ifdef TOPAZ_params
         nt = 2
         Canew = (ca_det_calc1d(k) +  Ca_tend(k)*obio_deltat)
         ca_det_calc1d(k) = Canew
#endif
#endif

#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
         O2new = (o21d(k) +  O_tend(k)*obio_deltat)
#ifdef prescribe_o2sf
         if (k.eq.1) then
         O2new = o21d(k)
         endif
#endif
         o21d(k) = O2new
!@PLdebg
!       if (vrbos) then
!          write(6,'(a,3i7,3e12.4)')'obio_o2(postbio):',
!     .      nstep,i,j,p1d(k),o21d(k),O_tend(k)
!        endif

!      if (ISNAN(o21d(k))) then
!          errchk2=1
!          write(6,'(a,4i7,5e12.4)')'obio_o2(postbionan):',
!     .      nstep,i,j,errchk2,p1d(k),o21d(k),O_tend(k)
!     .      ,car(k,nt),C_tend(k,nt)

!        endif
!@PLdbg
#endif
#ifdef TRACERS_abio_O2
         Abo2new = (abo21d(k) + Abo_tend(k)*obio_deltat)
         abo21d(k) = Abo2new
#endif
#endif


 1000 continue


      return
      end
