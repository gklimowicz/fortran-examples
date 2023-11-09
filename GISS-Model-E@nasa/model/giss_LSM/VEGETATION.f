#include "rundeck_opts.h"

      module vegetation

      use constant, only : zero,tfrz=>tf,gasc
      use ghy_h, only : ngm, imt

      implicit none
      save

      private

c**** public functions:
      public veg_conductance, update_veg_locals
      public t_vegcell

ccc t_vegcell is a hack to pass vegetation-related data to 
ccc ground hydrology. Basically it is an equivalent of entcell.
      type t_vegcell
        real*8 fr(ngm)
        real*8 snowm
        real*8 ws_can
        real*8 shc_can
        real*8 alaie,rs,nm,nf,alai,vh
     &     ,alait(11),vfraction(11)
     &     ,fdir,parinc,vegalbedo,sbeta,Ci,Qf ! added by adf, nyk
     &     ,Ca     ! nyk
     &     ,cond_scheme,vegCO2X_off  !nyk
     &     ,cnc
      end type t_vegcell

ccc   rundeck parameters  5/1/03 nyk
!@dbparam cond_scheme selects vegetation conductance scheme:
!@+       = 1     default, original conductance scheme.
!@+       = 2     conductance scheme of Andrew Friend
      integer, public :: cond_scheme = 2
!@dbparam vegCO2X_off turns off CO2X doubling for vegetation.
!@+       = 0     default, vegetation sees both ghg_yr and CO2X.
!@+       = 1     turn off doubling of CO2 for veg by parameter CO2X.
      integer, public :: vegCO2X_off = 0

!input from driver:
!from veg_set_cell:
      real*8, public :: alaie,rs,alai,nm,nf,vh,CO2atm
!nyk alait is lai by functional type.  Must multiply by vfraction
!     to get fraction cover per grid cell
      real*8, dimension(11), public :: alait, vfraction
!from earth:
      real*8, public :: srht,pres,ch,vsm

!@var parinc Incident photosynthetically active (visible solar, dir+dif)
!@+   radiation on the surface (W/m2) (nyk)
      real*8, public :: parinc
!@var fdir Fraction of surface visible radiation that is direct (adf)
      real*8, public :: fdir
!@var parinc Incident photosynthetically active (visible solar, dir+dif)
!     radiation on the surface (W/m2) (nyk)
!veg      real*8, public :: parinc
!@var vegalbedo  Vegetation canopy albedo averaged over the grid cell.
!     This is a temporary hack to keep the prescribed albedos until
!     a canopy radiation scheme is in place. (nyk)
      real*8, public :: vegalbedo
!@var sbeta Sine of solar zenith angle (rad).
      real*8, public :: sbeta
!@var Ci Original internal foliage CO2 (mol/m3) (adf)
      real*8, public :: Ci
!@var Qf Original foliage surface mixing ratio (kg/kg) (adf)
      real*8, public :: Qf
!****************************************************
!NEED TO REPLACE Ca WITH CO2 FROM THE CLIMATE MODEL TRACERS
!@var Ca Atmospheric CO2 concentration at surface height (mol/m3).
      real*8, public :: Ca !!=.0127D0 !Now is a variable, nyk, dummy initial
!****************************************************

!out:
      real*8, public :: CNC

!input from ghy:
      real*8 betad,tcan,qv
!@var qv Canopy saturated mixing ratio (kg/kg)

!out:
!      real*8 GPP 

      contains

      subroutine veg_conductance(
     &     vegcell
     &     ,cnc_out 
     &     ,gpp_out 
     &     ,trans_sw_out
     &     ,betad_in ! evaporation efficiency
     &     ,tcan_in ! canopy temperature C
     &     ,qv_in
     $     ,dt_in  ! GHY time step
     &     )
      implicit none
      type(t_vegcell) :: vegcell
      real*8, intent(out) :: cnc_out
      real*8, intent(out) :: gpp_out
      real*8, intent(out) :: trans_sw_out
      real*8, intent(in) :: betad_in,tcan_in,qv_in,dt_in


      alaie         = vegcell%alaie       
      rs            = vegcell%rs          
      nm            = vegcell%nm          
      nf            = vegcell%nf          
      alai          = vegcell%alai        
      vh            = vegcell%vh          
      alait         = vegcell%alait       
      vfraction(:)  = vegcell%vfraction(:)
      fdir          = vegcell%fdir        
      parinc        = vegcell%parinc      
      vegalbedo     = vegcell%vegalbedo   
      sbeta         = vegcell%sbeta       
      Ci            = vegcell%Ci          
      Qf            = vegcell%Qf          
      Ca            = vegcell%Ca          
cddd      cond_scheme   = vegcell%cond_scheme 
cddd      vegCO2X_off   = vegcell%vegCO2X_off 
      cnc           = vegcell%cnc         

      ! call stop_model("veg_conductance called",255)

      !print *, "VVV ", parinc, fdir*parinc, sbeta

      betad = betad_in
      tcan = tcan_in
      qv = qv_in

      if ( cond_scheme.eq.1 ) then !if switch added by nyk 5/1/03
        gpp_out=0               ! dummy value for GPP if old cond_scheme
        trans_sw_out=0          ! dummy value for old cond_scheme
        call cond
      else         ! cond_scheme=2 or anything else
        call veg(dt_in, gpp_out, trans_sw_out)
      endif
      cnc_out = CNC

      vegcell%alaie         = alaie       
      vegcell%rs            = rs          
      vegcell%nm            = nm          
      vegcell%nf            = nf          
      vegcell%alai          = alai        
      vegcell%vh            = vh          
      vegcell%alait         = alait       
      vegcell%vfraction(:)  = vfraction(:)
      vegcell%fdir          = fdir        
      vegcell%parinc        = parinc      
      vegcell%vegalbedo     = vegalbedo   
      vegcell%sbeta         = sbeta       
      vegcell%Ci            = Ci          
      vegcell%Qf            = Qf          
      vegcell%Ca            = Ca          
      vegcell%cond_scheme   = cond_scheme 
      vegcell%vegCO2X_off   = vegCO2X_off 
      vegcell%cnc           = cnc         

      end subroutine veg_conductance


      subroutine cond
c**** calculates the canopy conductance
c**** input:
c**** betad - beta due to roots
c**** alaie - effective leaf area index
c**** rs - minimum stomatal resistance, m s-1
c**** xinc - incoming solar radiation, w m-2
c**** tp - temperature of canopy, c
c**** tfrz - freezing point of water, 0 c in k
c**** output:
c**** cnc - canopy conductance, m s-1
c**** adjust canopy conductance for soil water potential
!@var c1 canopy conductance related parameter
      real*8, parameter :: c1 = 90.d0
      real*8 srht0
      CNC=betad*alaie/rs
c**** adjust canopy conductance for incoming solar radiation
      srht0=max(srht,zero)
      CNC=CNC*(srht0/c1)/(1.d0+srht0/c1)
      CNC=CNC/(1.d0+((tcan+tfrz-296.d0)/15.d0)**4)
      return
      end subroutine cond

!----------------------------------------------------------------------!
      subroutine veg( dt, GPP, TRANS_SW )
!----------------------------------------------------------------------!
! Vegetation dynamics routine (adf). Currently (April, 2003)
! calculates canopy stomatal conductance (CNC, m/s) using a
! biologically-based formulation, along with canopy CO2 fluxes (GPP;
! kg[C]/m2/s). Vegetation growth will be included next.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
      real*8, intent(in) :: dt
      real*8, intent(out) :: GPP
      real*8, intent(out) :: TRANS_SW

!@var sigma Leaf scattering coefficient (?unitless).
      real*8, parameter :: sigma=0.2D0
      real*8  temp
!@var kdf Canopy extinction coeff. for diffuse radiation (unitless).
      real*8, parameter :: kdf=0.71D0
!@var Oi Internal foliage O2 partial pressure (kPa).
      real*8, parameter :: Oi=20.9D0
!@var k Canopy nitrogen extinction coefficient (?unitless).
      real*8, parameter :: k=0.11D0
!@var PAR Incident photosynthetically active radiation (umol/m2/s).
      real*8 PAR
!@var tk Canopy temperature (K).
      real*8 tk
!@var prpa Atmospheric pressure (Pa).
      real*8 prpa
!@var I0dr Direct beam PAR incident on canopy (umol/m2/s).
      real*8 I0dr
!@var I0df Diffuse PAR incident on canopy (umol/m2/s).
      real*8 I0df
!------Full canopy radiation transmittance, nyk------------------------
!@var absdf Absorbance of diffuse SW radiation (fraction)
      real*8 absdf
!@var absdr Absorbance of direct SW radiation (fraction)
      real*8 absdr
!@var absdrdr Absorbance of direct SW that remains direct (fraction)
      real*8 absdrdr
!@var abssh Absorbance of SW through shaded foliage (fraction)
      real*8 abssh
!@var abssl Absorbance of SW through sunlit foliage (fraction)
      real*8 abssl
!@var fracsl Fraction of leaves in canopy that are sunlit (fraction)
      real*8 fracsl
!Output variable trans_sw Total canopy transmittance of shortwave (fraction)
!      real*8 TRANS_SW  !Subroutine var parameter
!----------------------------------------------------------------------
!@var rho Canopy reflectivity (?unitless).
      real*8 rhor
!@var kbl Canopy extinction coeff. for direct radiation (unitless).
      real*8 kbl
!@var CiPa Internal foliage CO2 partial pressure (Pa).
      real*8 CiPa
!@var pcp Photorespiratory compensation point (Pa).
      real*8 pcp
!@var Kc Rubisco Michaelis-Menten constant for CO2 (Pa).
      real*8 Kc
!@var Ko Rubisco Michaelis-Menten constant for O2 (Pa).
      real*8 Ko
!@var n1 Proportionality between Jmax and N (umol[CO2]/mmol[N]/s).
      real*8 n1
!@var n2 Proportionality between Vcmax and N (umol[CO2]/mmol[N]/s).
      real*8 n2
!@var m1 CO2 dependence of light-limited photosynthesis (?units).
      real*8 m1
!@var m2 CO2 dependence of light-saturated photosynthesis (?units).
      real*8 m2
!@var msat Nitrogen dependence of photosynthesis (?units).
      real*8 msat
!@var Ntot Total canopy nitrogen (g/m[ground]2).
      real*8 Ntot
!@var N0 Foliage nitrogen at top of canopy (mmol/m2).
      real*8 N0
!@var Acan Canopy A (umol[CO2]/m[ground]2/s).
      real*8 Acan
!@var Acan Canopy net A (umol[CO2]/m[ground]2/s).
      real*8 Anet
!@var Acan Canopy A at saturating CiPa (umol[CO2]/m[ground]2/s).
      real*8 Amax
!@var Acan Canopy net A at saturating CiPa (umol[CO2]/m[ground]2/s).
      real*8 Anet_max
!@var Rcan Canopy mitochondrial respiration (umol/m2/s).
      real*8 Rcan
!@var dqs Humidity deficit across canopy surface (kg/kg).
      real*8 dQs
!@var CNCN New equilibrium canopy conductance to water (m/s).
      real*8 CNCN
!nu@var dCNC Change in canopy conductance to reach equilibrium (m/s).
      real*8 dCNC
!nu@var dCNC_max Maximum possible change in CNC over timestep (m/s).
      real*8 dCNC_max
!@var gt Conductance from inside foliage to surface height at 30m (m/s).
      real*8 gt
!@var EPS to stop dbzs
      real*8, parameter :: EPS = 1.d-8
!@var GPP Gross primary productivity (kg[C]/m2/s).
!      real*8 GPP    !moved to public, nyk 04/25/03
!----------------------------------------------------------------------!
! Make sure is some leaf area (m2/m2).
      if(alai.le.EPS)alai=EPS
! Convert radiation from W/m2 (total shortwave) to umol/m2/s (PAR).
! Really want incident here, but I think srht is absorbed. Based on
! Hansen et al. (1983), assume vegetation albedo is about 0.08, and
! so allow for this here.  2.3 umol/J for conversion from shortwave to
! fraction that is PAR (Monteith & Unsworth).
!nu    PAR=2.3d0*max(srht,zero)/(1.0D0-0.08D0)
! Replaced back-calculation with actual incident PAR at surface.
! nyk  For 400-700 nm, energy is 3.3-6 umol photons/J.
!      Top-of-atmosphere flux-weighted average of energy is 4.54 umol/J.
!      Univ. Maryland, Dept. of Met., PAR Project, suggests nominal
!      485 nm for conversion, which gives 4.05 umol/J.
      PAR=4.05d0*parinc          !W/m2 to umol/m2/s
!      write (99,*) 'PAR umol/m2/s',PAR,'PARsrht',
!     $     2.3d0*max(srht,zero)/(1.0D0-0.08D0)
      !DEBUG
! Convert canopy temperature from oC to K.
      tk=tcan+tfrz
! Convert atmospheric pressure from mbar to Pa.
      prpa=100.0D0*pres
! Internal foliage CO2 partial pressure from concentration (Pa).
      CiPa=Ci*(gasc*tk)
! Direct beam PAR incident on canopy (umol/m2/s).
      I0dr=fdir*PAR
! Diffuse PAR incident on canopy (umol/m2/s).
      I0df=(1.0D0-fdir)*PAR
! Canopy reflectivity depends on sun angle. Or take prescribed albedos.
      temp=sqrt(1.0D0-sigma)
      !Hack canopy reflectivity
      if (vegalbedo.eq.0) then !because radiation gets initialized late
        call stop_model("vegalbedo == 0", 255)
        rhor=((1.0D0-temp)/(1.0D0+temp))*(2.0D0/(1.0D0+1.6D0*sbeta))
        !write (99,*) 'rhor', rhor
      else
       !Temporary: keep prescribed veg albedos until have canopy scheme.
        rhor = vegalbedo
        !write (99,*) 'vegalbedo', vegalbedo
      end if
! Canopy extinction coefficient for direct beam radiation depends on
! sun angle (unitless).
      kbl=0.5D0*kdf/(0.8D0*temp*sbeta+EPS)
! Photorespiratory compensation point (Pa).
      pcp=exp(19.02D0-37.83D3/(gasc*tk))*prpa/1.0D6
! Rubisco Michaelis-Menten constant for CO2 (Pa).
      Kc=exp(38.05D0-79.43D3/(gasc*tk))*prpa/1.0D6
! Rubisco Michaelis-Menten constant for O2 (Pa).
      Ko=exp(20.30D0-36.38D3/(gasc*tk))*prpa/1.0D6
! Proportionality coefficient between Jmax and N (umol[CO2]/mmol[N]/s).
      n1=nf*0.12D0*0.95D+14*exp(-79.5D3/(gasc*tk))/
     1    (1.0D0+exp((650.0D0*tk-199.0D3)/(gasc*tk)))
! Proportionality coefficient between Vcmax and N (umol[CO2]/mmol[N]/s).
      n2=nf*0.23D0*exp(26.35D0-65.33D3/(gasc*tk))
! CO2 dependence of light-limited photosynthesis (?units).
      m1=CiPa/(CiPa+2.0D0*pcp)
! CO2 dependence of light-saturated photosynthesis (?units).
      m2=CiPa/(CiPa+kc*(1.0D0+Oi/ko))
! Nitrogen dependence of photosynthesis (?units).
      msat=min(m1*n1,m2*n2)
! Total canopy nitrogen (g/m[ground]2).
      Ntot=nm*alai
! Top of canopy nitrogen (mmol/m[foliage]2).
      N0=(1000.0D0/14.0D0)*k*Ntot/(1.0D0-exp(-k*alai))
! Canopy photosynthesis (Acan: umol/m2/s).
      call qsimp(alai,Acan,I0dr,I0df,CiPa,tk,prpa,nf,pcp,Kc,Oi,Ko,N0,
     1k,n1,n2,m1,msat,sigma,temp,rhor,kdf,kbl)
!      write(99,*)'sbeta,Acan',sbeta,Acan
!----------------------------------------------------------------------!
! Transmission of shortwave radiation through canopy.
! Diffuse shortwave in canopy (umol/m[ground]2/s).
      if(sbeta.gt.zero)then
        absdf=(1-fdir)*(1.0D0-rhor)*kdf*exp(-kdf*alai)
! Direct shortwave in canopy (umol/m[ground]2/s).
        absdr=fdir*(1.0D0-rhor)*temp*kbl*exp(-temp*kbl*alai)
! Direct shortwave that remains direct (umol/m[ground]2/s).
        absdrdr=absdr*(1.0D0-sigma)*kbl*exp(-temp*kbl*alai)
! Shortwave penetrating shaded foliage (umol/m[foliage]2/s).
        abssh=absdf + (absdr-absdrdr) 
! Shortwave penetrating sunlit foliage (umol/m[foliage]2/s).
        abssl=abssh+(1.0D0-sigma)*kbL*fdir
        if(abssh.lt.0.0D0) abssh=0.0D0
        if(abssl.lt.0.0D0) abssl=0.0D0
        if(abssh.gt.1.0D0) abssh=1.0D0
        if(abssl.gt.0.0D0) abssl=1.0D0
        fracsl=exp(-kbl*alai)
      else
        abssh=0.001D0
        abssl=0.0D0
        fracsl=0
      endif
!      TRANS_SW = 1-((1-fracsl)*abssh + fracsl*abssl)
      TRANS_SW = 1-((1-fracsl)*abssh + fracsl*abssl)
!      write(110,*) TRANS_SW,I0dr,I0df,abssh, abssl, sbeta, sigma, kbl
!----------------------------------------------------------------------!
! Saturating Ci to calculate canopy conductance (mol/m3).
      CiPa=1.0D6
! CO2 dependence of light-limited photosynthesis (?units).
      m1=CiPa/(CiPa+2.0D0*pcp)
! CO2 dependence of light-saturated photosynthesis (?units).
      m2=CiPa/(CiPa+kc*(1.0D0+Oi/ko))
! Nitrogen dependence of photosynthesis (?units).
      msat=min(m1*n1,m2*n2)
! Canopy photosynthesis at saturating CiPa (Amax: umol/m[ground]2/s).
      call qsimp(alai,Amax,I0dr,I0df,CiPa,tk,prpa,nf,pcp,Kc,Oi,Ko,N0,
     1k,n1,n2,m1,msat,sigma,temp,rhor,kdf,kbl)
!----------------------------------------------------------------------!
! Canopy mitochondrial respiration (umol/m[ground]2/s).
      Rcan=0.20D0*Ntot*exp(18.72D0-46.39D3/(gasc*tk))
! Net canopy photosynthesis (umol/m[ground]2/s).
      Anet=Acan-Rcan
! Net canopy photosynthesis at saturating CiPa (umol/m[ground]2/s).
      Anet_max=Amax-Rcan
!----------------------------------------------------------------------!
! Humidity deficit across canopy surface (kg/kg).
      dQs=Qv-Qf
      if(dQs.lt.zero)dQs=zero
!----------------------------------------------------------------------!
! New equilibrium canopy conductance to moisture (m/s). 
      CNCN=betad*(1.0D0-0.0075D0*vh)*650.0D-6*Anet_max*
     &   ((Ci+0.004D0)/(Ci+EPS))*2.8D0**(-80.0D0*dQs)
! Required change in canopy conductance to reach equilibrium (m/s).
      dCNC=CNCN-CNC
!nu Limit CNC change over timestep because of guard cell mechanics (m/s)
      dCNC_max=dt*alai*(0.006D0-0.00006D0)/1800.0D0
      if( dCNC.gt.dCNC_max)CNCN=CNC+dCNC_max
      IF(-dCNC.gt.dCNC_max)CNCN=CNC-dCNC_max
! Biological limits of absolute CNC (m/s).
      if(CNCN.gt.0.006*alai)CNCN=0.006*alai
!NOTE:  Water balance issue due to not modeling canopy water content
!       explicitly.  This needs to be considered at some point. -nyk
      if(CNCN.lt.0.00006*alai)CNCN=0.00006*alai
!----------------------------------------------------------------------!
! Update Ci.
! Total conductance from inside foliage to surface height at 30m (m/s),
! where CO2 is assumed at Ca.
      gt=1.0D0/(1.42D0/CNCN+1.65D0/(ch*vsm+EPS))
! Foliage internal CO2 concentration for next timestep (mol/m3).
      Ci=Ca-1.0D-6*Anet/gt
! Limit Cin to physical realism (mol/m3). It is possible that
! oscilliations could occur due the Ci<>CNC feedback. Also, something
! to watch out for is that setting Ci like this does not conserve CO2.
      if(Ci.lt.EPS)Ci=EPS
! NOTE:  Cannot update Qf here, because requires evap_tot calculated by
!        ground hydrology.  Therefore updated through update_veg_locals.
!----------------------------------------------------------------------!
! OUTPUTS:
! Transmission of shortwave through canopy.
!!!  trying something simple while the main algorithm is debugged
!      TRANS_SW = (1.d0-2.d0/3.14d0)**alai
      if ( TRANS_SW < 0.d0 .or. TRANS_SW > 1.d0 )
     &     call stop_model("veg: unphysical TRANS_SW", 255)
! Gross primary productivity (kg[C]/m2/s).
      GPP=0.012D-6*Anet   ! should be dependent on conductance ??
! Canopy conductance for next timestep (m/s).
      CNC=CNCN
      return
      end subroutine veg
!----------------------------------------------------------------------!
      subroutine qsimp(B,S,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,
     &                 m1,msat,sigma,temp,rhor,kdf,kbl)
!----------------------------------------------------------------------!
! qsimp calculates canopy photosynthesis by increasing the number of
! layers in the integration until the result (S) changes by less than
! 0.1 umol[CO2]/m2/s.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Passed parameters
      real*8 B,S,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,msat
      real*8 sigma,temp,rhor,kdf,kbl
!----------------------------------------------------------------------!
!Local variables
!nu   real*8, parameter :: EPS=1.D-3
      integer, parameter :: MAXIT=6
      integer IT, canopylayers
      real*8 A, OST,OS,ST,ERROR
!----------------------------------------------------------------------!
      A=0.0D0
      OST=-1.D30
      OS= -1.D30
      canopylayers=1
      do 11 IT=1,MAXIT
         CALL TRAPZD(A,B,ST,IT,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,
     &               n2,m1,msat,sigma,temp,rhor,kdf,kbl,canopylayers)
         S=(4.D0*ST-OST)/3.D0
         ERROR=ABS(S-OS)
!nu      IF (ERROR.lt.EPS*ABS(OS)) RETURN
         IF (ERROR.lt.0.1D0) RETURN
         OS=S
         OST=ST
         if(IT.gt.1) canopylayers=canopylayers*2
   11 enddo
!      write(99,*) 'Too many steps.'
!      write(99,*) S,ERROR,100.0D0*ERROR/S
      return
      end subroutine qsimp
!----------------------------------------------------------------------!
      subroutine trapzd(A,B,S,N,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,
     &n1,n2,m1,msat,sigma,temp,rhor,kdf,kbl,canopylayers)
!----------------------------------------------------------------------!
! Integrates canopy photosynthesis over canopy layers using Simpson's
! Rule (Press et al., 19??).
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
      integer N,canopylayers, L
      real*8 A,B,S,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,msat
      real*8 sigma,temp,rhor,kdf,kbl
      real*8 DEL,X,SUM,RCL
!@var func1 Mean net photosynthesis at Lc (umol[CO2]/m2/s).
      real*8 func1
!@var func2 Mean net photosynthesis at Lc (umol[CO2]/m2/s).
      real*8 func2
      if(N.eq.1)then
         call phot(A,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,
     &             msat,sigma,temp,rhor,kdf,kbl,B,func1)
         call phot(B,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,
     &             msat,sigma,temp,rhor,kdf,kbl,B,func2)
         S=0.5D0*(B-A)*(func1+func2)
      else
         RCL=canopylayers  ! convert to real*8
         DEL=(B-A)/RCL
         X=A+0.5D0*DEL
         SUM=0.D0
         do 11 L=1,canopylayers
           call phot(X,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,
     &               msat,sigma,temp,rhor,kdf,kbl,B,func1)
           SUM=SUM+func1
           X=X+DEL
   11    continue
         S=0.5D0*(S+(B-A)*SUM/RCL)
      endif
      return
      end subroutine trapzd
!----------------------------------------------------------------------!
      subroutine phot(Lc,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,
     &                m1,msat,sigma,temp,rhor,kdf,kbl,alai,func)
!----------------------------------------------------------------------!
! Calculate mean leaf photosynthesis at cumulative leaf area index Lc
! (from top of canopy), returning result as func (umol[CO2]/m2/s).
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
      real*8 func,I0dr,I0df,Ci,T,P,nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,msat
      real*8 sigma,temp,rhor,kdf,kbl,alai
!----------------------------------------------------------------------!
!@var alpha Intrinsic quantum efficiency (?units).
      real*8, parameter :: alpha=0.08D0
!@var ka Chlorophyll extinction coefficient (?units).
      real*8, parameter :: ka=0.005D0
!@var n3 Ratio of foliage chlorophyll to N (?units).
      real*8 n3
!@var Lc Cumulative LAI from top of canopy (m2/m2).
      real*8 Lc
!@var Np Foliage nitrogen content (mmol/m2).
      real*8 Np
!@var Idfa Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idfa
!@var Idra Direct PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idra
!@var Idrdra Direct PAR at Lc that remains direct (umol/m[ground]2/s).
      real*8 Idrdra
!@var Isha PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 Isha
!@var Isla PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 Isla
!@var fsl Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 fsl
!@var Nlim Theoretical cumulative N for light-saturated A (mmol/m2).
      real*8 Nlim
!@var Nsat Actual cumulative N for light-saturated A (mmol/m2).
      real*8 Nsat
!@var fabs Absorbed radiation in light-limited chloroplasts (fraction).
      real*8 fabs
!@var FUNCsl Photosynthesis in sunlit foliage (umol/m2/s).
      real*8 FUNCsl
!@var FUNCsh Photosynthesis in shaded foliage (umol/m2/s).
      real*8 FUNCsh
!@var To stop dbzs.
      real*8, parameter :: EPS = 1.d-8
!----------------------------------------------------------------------!
! Check following OK. n3 is ratio of chlorophyll to foliage N (units?).
      n3=6.0D0-3.6D0*exp(-0.7D0*Lc)
! Foliage N at canopy depth Lc (mmol/m2).
      Np=N0*exp(-k*Lc)
! Following conditional required as canopy radiation profile only works
! when sun is above the horizon.
      if(sbeta.gt.zero)then
! Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
        Idfa=(1.0D0-rhor)*I0df*kdf*exp(-kdf*Lc)
! Direct PAR in canopy at Lc (umol/m[ground]2/s).
        Idra=(1.0D0-rhor)*I0dr*temp*kbl*exp(-temp*kbl*Lc)
! Direct PAR at Lc that remains direct (umol/m[ground]2/s).
        Idrdra=(1.0D0-sigma)*I0dr*kbl*exp(-temp*kbl*Lc)
! PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
        Isha=Idfa+(Idra-Idrdra)
! PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s).
        Isla=Isha+(1.0D0-sigma)*kbL*I0dr
! Fraction of sunlit foliage in layer at Lc (unitless).
        fsl=exp(-kbL*Lc)
        if(Isha.lt.0.1D0)Isha=0.1D0
        if(Isla.lt.0.1D0)Isla=0.1D0
        if(fsl.lt.zero)fsl=zero
      else
        Isha=0.1D0
        Isla=0.1D0
        fsl=zero
      endif
!      write(99,*) I0dr, I0df,Isha, Isla, kbl
!----------------------------------------------------------------------!
! Cumulative nitrogen concentration at which photosynthesis becomes
! light-limited in sunlit foliage (mmol/m2).
      Nlim=-dlog(msat/(alpha*Isla*ka*n3*m1+EPS))/(ka*n3)
! Impose limits on Nsat based on actual foliage nitrogen (mmol/m2).
      if(Nlim.lt.zero)then
        Nsat=zero
      elseif(Nlim.gt.Np)then
        Nsat=Np
      else
        Nsat=Nlim
      endif
! Absorbed radiation in light-limited chloroplasts (fraction).
      fabs=exp(-ka*n3*Nsat)-exp(-ka*n3*Np)
! Photosynthesis in sunlit foliage (umol/m2/s).
      FUNCsl=(1.0D0-PCP/(Ci+EPS))*(msat*Nsat+alpha*m1*Isla*fabs)
!----------------------------------------------------------------------!
! Cumulative nitrogen concentration at which photosynthesis becomes
! light-limited in shaded foliage (mmol/m2).
      Nlim=-log(msat/(alpha*Isha*ka*n3*m1+EPS))/(ka*n3)
! Impose limits on Nsat based on actual foliage nitrogen (mmol/m2).
      if(Nlim.lt.0.0D0)then
        Nsat=0.0D0
      elseif(Nlim.gt.Np)then
        Nsat=Np
      else
        Nsat=Nlim
      endif
! Absorbed radiation in light-limited chloroplasts (fraction).
      fabs=exp(-ka*n3*Nsat)-exp(-ka*n3*Np)
! Photosynthesis in shaded foliage (umol/m2/s).
      FUNCsh=(1.0D0-PCP/(Ci+EPS))*(msat*Nsat+alpha*m1*Isha*fabs)
!----------------------------------------------------------------------!
! Mean photosynthesis in layer (umol/m2/s).
      func=fsl*FUNCsl+(1.0D0-fsl)*FUNCsh
!----------------------------------------------------------------------!
      return
      end subroutine phot
!----------------------------------------------------------------------!

      !----------------------------------------------------------------
      subroutine update_veg_locals(evap_tot2, rho, rhow, ch, vsm,qs
     &     ,qf_out)
      !Update vegetation input variables that require values external
      !to vegetation module to update.  For values that change within
      !the smaller time step of the GHY modules.

      real*8,intent(in):: evap_tot2
      real*8,intent(in):: rho, rhow, ch, vsm, qs
      real*8, intent(out) :: qf_out
      real*8 rho3,cna
      
      !adf Get new Qf
      rho3=rho/rhow
      cna=ch*vsm
      Qf=evap_tot2/(rho3*cna)+qs  ! New mixing ratio for next timestep (adf)
      qf_out = Qf
      end subroutine update_veg_locals
      !----------------------------------------------------------------


      subroutine veg_accm

      ! agpp = agpp + gpp*(1.d0-fr_snow(2)*fm)*fv*dts

      entry veg_accm0
      
      !agpp=0.d0   ! new accumulator, nyk 4/25/03
      

      end subroutine veg_accm


      end module vegetation








