#include "rundeck_opts.h"

      subroutine obio_alkalinity(kmax,i,j,nstep)

!@sum  online computation of alkalinity
!@auth Natassa Romanou

! compute source/sink term for alkalinity
! zc     compensation depth: photosynthesis = respiration
! or more precisely p1d(kzc)
! J_PO4  source/sinks of phosphate
! J_Ca   source/sinks of calcium carbonate
! Jprod  production of organic phosphorus which is linked to particulate
!        organic carbon POC production by 
!        J_POC=(1-sigma_Ca)*J_DOC=(1-sigma_Ca)*rC:P*Jprod
!        to estimate J_DOC use total productivity pp
!        Really, production of particulate organic carbon (or phosphorus)
!        above the compensation depth
! F_Ca   downward flux of CaCO3
! Fc     instantaneous downward flux of particulate organic phosphorus at 
!        compensation depth
! sigma_Ca fraction of P converted to DOP, or fraction of C converted to DOC
! 1-sigma_Ca fraction of P NOT converted to DOP, or fraction of C NOT converted to DOC
! alk tendency is also computed under ice.  this is because 
! alk changes due to sal and temp changes.
! alkalinity units should be umol/kg (uE/kg to go into obio_carbon)
! alk tendency computed here is in umol/m3/s, convert to umol/kg
! *** a note on units: I assume here that units are umol,N/m3
!-------------------------------------------------------------------------
! J_ALK = -rN:P* J_PO4 + 2*J_Ca
! J_PO4 = J_NO3/rN:P, tendency of phosphates
! J_NO3 = tendency of nitrates, P_tend(:,1)
! J_Ca = - R* rC:P* (1-sigma_Ca)* Jprod, for  z<zc  ***NOTE sign was wrong in OCMIP documentation
! J_Ca = - dF_Ca/dz,   for z>zc
! Jprod = pp, for z<zc 
! Jprod = 0,  for z>zc
! Distinguish two cases: pp based on all species
!                        pp based on coccolithophores only
! F_Ca = R * rC:P * Fc * exp(-(z-zc)/d)
!   Fc = (1-sigma_Ca) * integral_0_zc (Jprod*dz)
!   ! Fc = (1-sigma_Ca) * integral_0_zc (J_DOC*dz)
!
! note sign errors in equations 28a and 31 of notes.
!-------------------------------------------------------------------------

      USE obio_dim
      USE obio_incom, only: rain_ratio,cpratio,sigma_Ca,d_Ca,
     .      npratio,uMtomgm3,cnratio,bn,zc,mgchltouMC
      USE obio_com, only: P_tend,p1d,pp2_1d,dp1d,A_tend,
     .      rhs,alk1d,caexp,kzc,sday,rho_water

      implicit none

      integer, intent(in) :: nstep

      integer nt,k,kmax,nchl1,nchl2,i,j
      real*8 J_PO4(kmax),pp,Jprod(kmax),Jprod_sum,Fc,zz,F_Ca(kmax+1),
     .       J_Ca(kmax),term,term1,term2,DOP,offterm



!--------------------------------------------------------------------------
!only compute tendency terms if total depth greater than conpensation depth
      if (p1d(kmax+1) .lt. p1d(kzc)) then
         A_tend(:) = 0.
         go to 100
      endif

!pp and therefore Jprod are defined at dp points (mid level)
!F_Ca is defined at pressure interfaces
!J_Ca is therefore defined mid-layer (dp points)

!compute sources/sinks of phosphate
!J_PO4 units uM/hr
! now uM/s  July 2016
      do k=1,kmax
!     J_PO4(k) =  P_tend(k,1)           !approximate by nitrate conc tendency
!                                       !NO3/PO4 ratio from Conkright et al, 1994
!                                       !expressed as nitrates
!                                       ! in notes npratio that multiplies tendency term cancells 
!                                       ! out with npratio that divides to get units of nitrate
!     J_PO4(k) =  rhs(k,1,5)+rhs(k,1,6)+rhs(k,1,7)+rhs(k,1,8)    !uptake of nitrate
!     J_PO4(k) =((rhs(k,1,5)+rhs(k,1,6)+rhs(k,1,7)+rhs(k,1,8))/bn
!    .         +  rhs(k,14,6)/mgchltouMC
!    .         +(rhs(k,2,5)+rhs(k,2,6)
!    .         + rhs(k,2,7)+rhs(k,2,8))/bn)/cnratio

      J_PO4(k) =  rhs(k,1,5)+rhs(k,1,6)+rhs(k,1,7)+rhs(k,1,8)+rhs(k,1,9)
     .          + rhs(k,1,10)+rhs(k,1,11)+rhs(k,1,14)+rhs(k,1,15)

      term = -1.d0* J_PO4(k)            !uM,N/s= mili-mol,N/m3/s  July 2016
#ifdef alk_adj
      term = 0.2d0*term
#endif
#ifdef alk_adj2
      term = -0.2d0*P_tend(k,1)
#endif
#ifdef alk_adj3
      term = -1.d0*(P_tend(k,1)+P_tend(k,2))
#endif
#ifdef alk_adj4
      term = -0.15d0*(P_tend(k,1)+P_tend(k,2))
#endif
      rhs(k,15,1) = term
      A_tend(k)= term 
      enddo

      call cadet_topaz(kmax)

!distinguish two cases: OCMIP uses total pp for Jprod, whereas here we also
!consider the case where Jprod only includes coccolithophores
#ifdef Jprod_based_on_pp
      nchl1=1
      nchl2=nchl
#endif
#ifdef Jprod_based_on_cocc
      nchl1=4
      nchl2=4   
#endif

      if (nchl1.le.0 .or. nchl2.le.0) then
          print*, nchl1, nchl2, 
     .   'MUST SET Jprod_based_on_pp or Jprod_based_on_cocc in rundeck'
          stop  
      endif

!compute total primary production
!and the Jprod term (production of organic phosphorus)
      do k=1,kmax
        if (p1d(k) .le. p1d(kzc))  then
          pp=0.
          do nt=nchl1,nchl2
             pp=pp+pp2_1d(k,nt)/sday/dp1d(k)   ![pp2_1d]=mgC/m2/day -> [pp]=mgC/m3/s  July 2016
          enddo
          Jprod(k) = pp/cpratio  !mgP/m3/s    July 2016
        else
          Jprod(k)=0.
        endif

!     write(*,'(a,5i5,5e12.4)')'obio_alkalinity1:',
!    .            nstep,i,j,k,kzc,
!    .            p1d(k),p1d(k+1),p1d(kzc),pp,Jprod(k)
      enddo

!integrate net primary production down to zc
      Jprod_sum = 0.
      do k=1,kmax
      if (p1d(k) .le. p1d(kzc))  then
       Jprod_sum = Jprod_sum + Jprod(k)*dp1d(k)   ! mgP/m2/s   July 2016
      endif
      enddo

      Fc = (1.d0-sigma_Ca)* Jprod_sum     !mgP/m2/s  July 2016

!compute downward flux of CaCO3
!only below the euphotic zone (the compensation layer)
      F_Ca = 0.d0
      do k=kzc,kmax+1
           F_Ca(k) = rain_ratio*cpratio*Fc
     .             * exp(-1.d0*(p1d(k)-p1d(kzc))/d_Ca)    !mgC/m2/s   July 2016
      enddo

      !p1d(kzc) is really the compensation depth
      !F_Ca(kzc) is the CaCO3 export
!     write(*,'(a,i8,3i5,7e12.4)')'CaCO3 downward flux:',
!    .        nstep,i,j,kzc,p1d(kzc),zc,rain_ratio,cpratio,Fc,
!    .        exp(-1.d0*(p1d(kzc)-p1d(kzc))/d_Ca),F_Ca(kzc)


      caexp = F_Ca(kzc)        !mili-gC/m2/s   July 2016

!compute sources/sinks of CaCO3
      offterm= 0.d0
      do k=1,kmax
         if (p1d(k) .le. p1d(kzc))  then
             !formation of calcium carbonate above compensation depth
             J_Ca(k) = -1.d0*rain_ratio*cpratio*(1.-sigma_Ca)*Jprod(k)   !mgC/m3/s  July 2016
         else
             !dissolution of calcium carbonate below compensation depth
             J_Ca(k) = -1.d0* (F_Ca(k+1)-F_Ca(k)) / dp1d(k)   !mgC/m3/s   July 2016
         endif
       term = 2.d0* J_Ca(k)/cnratio  ! mgC/m3/s -> mili-mol,N/m3/s  
                                     ! no need to multiply here by mol.weight
                                     ! because already in cnratio (see obio_init)
       rhs(k,15,5) = term
       A_tend(k) = A_tend(k) + term      

       offterm = offterm + term * dp1d(k)

!    .write(*,'(a,4i5,3e12.4)')'obio_alkalinity2:',
!    .   nstep,i,j,k,term,offterm,rhs(k,15,5)


#ifdef no_offtermalk
#else
      !bottom boundary condition adjust bottom layer
      if (k.eq.kmax) then
         rhs(kmax,15,5) = rhs(kmax,15,5) - offterm /dp1d(kmax)
         A_tend(kmax) = A_tend(kmax) - offterm / dp1d(kmax)
      endif
#endif

!     if(vrbos)
!    .write(*,'(a,4i5,3e12.4)')'obio_alkalinity3:',
!    .   nstep,i,j,k,term,offterm,rhs(k,15,5)

      enddo

      !for consistency, keep term that goes into rhs table in uM/hr = mili-mol,N/m3/hr
      !convert A_tend terms into uE/kg/hr, the actual units of alkalinity 
      A_tend = A_tend /rho_water *1.d3     ! mili-mol,N/m3/s -> umol/m3/s -> umol/kg/s

!!!!!!!!!! NEED TO ADD BOTTOM BOUNDARY CONDITIONS 

!     do k=1,kmax
!     k=1
!     write(*,'(a,4i5,12e12.4)')'obio_alkalinity; ',
!    .    nstep,i,j,k
!    .   ,p1d(k),p1d(kzc),J_PO4(k),pp,Jprod_sum,Fc
!    .   ,F_Ca(k),J_Ca(k),alk1d(k),A_tend(k),term1,term2
!     enddo

 100  continue
      end subroutine obio_alkalinity

      subroutine cadet_topaz(kmax,vrbos)

      USE MODEL_COM, only: dtsrc
      USE obio_dim
      USE obio_incom, only: cnratio,mgchltouMC
      USE obio_com, only: co3_conc,temp1d,saln1d,obio_P,p1d
     .                   ,wsdet

      implicit none

      integer k,kmax
      real*8 T,Salt,TK,PKSPC,wsink, gamma_cadet_calc
      real*8 co3_calc_sol,Omega_calc,J_cadet_calc
      real*8 lambda0,KEppley,P_star,P_min,CaN_calc_ratio,Omega_satmax
      real*8 N_cyanob,jgraz_cyanob_N,P_insitu,J_prod_cadet_calc
      real*8 Ca_det_calc
      logical vrbos

      do k=1,kmax
      !rate of constant dissolution of cadet_calc
      !wsink=1.d0
      wsink =  wsdet(k,1)    !wsdet for carbon????
      gamma_cadet_calc = wsink/1343.d0    ! in s-1      ! July 2016
!     gamma_cadet_calc = gamma_cadet_calc * 3600.d0     ! in hr-1

      T = temp1d(k)
      Salt = saln1d(k)
      TK = T + 273.15
      P_insitu = 0.1016*p1d(k)+1.013
      PKSPC = 171.9065
     .      + 0.077993 * TK
     .      - 2903.293 / TK
     .      - 71.595* dlog10(TK)
     .      - (-0.77712 + 2.8426d-3 *TK + 178.34/TK)*Salt**(1./2.)
     .      + 0.07711*salt
     .      - 4.1249d-3 * salt**(3./2.)
     .      - 0.02
     .      - (48.76 - 0.5304*T) * (P_insitu - 1.013) / (191.46*TK)
     .      + (1.d-3 *(11.76 -0.3692 *T))
     .      * (P_insitu -1.013)*(P_insitu-1.013)/(382.92*TK)
      co3_calc_sol = 10**(-PKSPC) / (2.937d-4 *max(5.d0, Salt))
      !??? co3_conc, how do I apply it over all depths?
      Omega_calc = co3_conc / co3_calc_sol

      !how do i compute ca_det_calc ????
!     Ca_det_calc =
      ! dissolution rate in units of mol kg-1 s-1 calculated in every layer,
      ! and has to be multiplied by the thickness to get a layer flux.
!     J_cadet_calc =
!    .     gamma_cadet_calc*max(0.d0,1.d0-Omega_calc)* Ca_det_calc

      !grazing rate constant at 0C
!     lambda0 = 0.19/86400. *3600.   !hr^-1
      lambda0 = 0.19/86400.   !s^-1
      !Temp. coeff for growth
      KEppley = 0.063   ! C^-1 (Eppley 1972)
      !nitrogen in small phytoplankton (=here cyanobacteria)
      N_cyanob = obio_P(k,3) *mgchltouMC/cnratio
      !pivot phyto. conc. for grazing allometry
      P_star = 1.9d-6 * 16/106   !mol,N/kg
      !min. phyto. conc. threshold for grazing
      P_min = 1.d-10   !mol,N/kg (added for stability)
      jgraz_cyanob_N = min(1.d0/dtsrc,
     .                     lambda0 * exp(KEppley*T)
     .      *(N_cyanob * N_cyanob /
     .       (P_star *(N_cyanob + P_min)))) * N_cyanob
      !Calcite CaCO3 to nitrogen uptake ratio
      CaN_calc_ratio =  0.005 * 106/16    !mol,Ca/mol,N *** tunable to give right caexp
      !maximum saturation state
      Omega_satmax = 10.  !dimensionless; to limit possible extreme values
      J_prod_cadet_calc = jgraz_cyanob_N * CaN_calc_ratio
     .                  * exp(-0.0539*T)
     .                  * min(Omega_satmax,max(0.d0,Omega_calc - 1.d0))

!     A_tend() = 2.d0 * (J_cadet_calc - J_prod_cadet_calc)

      enddo

      end subroutine cadet_topaz
