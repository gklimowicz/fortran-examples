#include "rundeck_opts.h"

      subroutine calc_flammability(t,p,r,v,flam,burnt_area,
     &                             fearth_axyp)
!@sum calculated the flammability of vegetation based on model
!@+ variables for temperature, precipitation, relative humidity,
!@+ and an index of vegetation density.
!
!@auth Keren Mezuman and Greg Faluvegi in part based on work by Olga Pechony
!@var T the surface air temperature passed in Kelvin for current I,J
!@var P the precipitation rate passed in mm/day for current I,J
!@var R the relative humidity passed as fraction for current I,J
!@var V vegetative density (unitless 0 to ~1) for current I,J
!@var burned_area area burned in a grid box (m^2)
!@var Z component of the Goff-Gratch saturation vapor pressure equation
!@var tsbyt = reciprocal of the temperature times ts
!@var flam the flammability returned for current I,J
!@param a,b,s,d,f,h,ts,cr coefficients for the parameterization
!@+ e.g. from Goff-Gratch saturation vapor pressure equation
!
      use constant, only: tf!freezing point of water at 1 atm [K]

      implicit none

      real*8, parameter :: a=-7.90298d0,d=11.344d0,c=-1.3816d-7,
     & b=5.02808,f=8.1328d-3,h=-3.49149d0,ts=tf+100.d0,cr=-2.d0
      real*8, intent(in) :: t,p,r,v,fearth_axyp,burnt_area
      real*8, intent(out) :: flam
      real*8 :: z,tsbyt,ba_frac

      tsbyt=ts/t

      z= a*(tsbyt-1.d0) + b*log10(tsbyt) + 
     &   c*(10.d0**(d*(1.d0-tsbyt))-1.d0) +
     &   f*(10.d0**(h*(tsbyt-1.d0))-1.d0)

      if (fearth_axyp <= 0.d0) then
        ba_frac = 0.d0
      else
        ba_frac = burnt_area/fearth_axyp
      endif
      if (ba_frac > 1) then
        ba_frac=1.d0
      endif
      flam=min( (10.d0**z*(1.d0-r)) * exp(cr*p) * v *
     &           (1.d0-ba_frac) , 1.d0)

      return
      end subroutine calc_flammability


      subroutine prec_running_average(p,avg,iH,iD,i0,first,HRA,DRA,PRS)
!@sum prec_running_average keeps a running average of the model
!@+ precipitation variable for use in the flammability model.
!@+ In practice, this does hourly and daily running averages and
!@+ uses those to get the period-long running average, to avoid
!@+ saving a huge array.
!@auth Greg Faluvegi
      use flammability_com, only: nday=>nday_prec,nmax=>maxHR_prec
 
      implicit none
      
!@var p model variable which will be used in the average (precip)
!@var temp just for holding current day average for use in avg_prec
!@var nmax number of accumulations in one day (for checking)
!@var bynmax reciprocal of nmax
!@var bynday reciprocal of nday
!@var i0 see i0fl(i,j) (day in period marker)
!@var iH see iHfl(i,j) (hour in day index)
!@var iD see iDfl(i,j) (day in period index)
!@var first see first_prec(i,j) whether in first period
!@var HRA see HRAfl(time,i,j) (hourly average)
!@var DRA see DRAfl(time,i,j) (daily average)
!@var PRS see PRSfl(i,j) (period running sum)
!@var avg see ravg_prec(i,j) the running average prec returned
      real*8, intent(IN) :: p
      real*8, dimension(nday) :: DRA
      real*8, dimension(nmax) :: HRA
      real*8 :: temp, bynmax, PRS, avg, iH, iD, i0, first, bynday
      integer :: n

      if(nint(iH) < 0 .or. nint(iH) > nmax) then
        write(6,*) 'iH maxHR_prec=',iH,nint(iH),nmax
        call stop_model('iHfl or maxHR_prec problem',255)
      endif
      bynmax=1.d0/real(nmax)
      bynday=1.d0/real(nday)
      iH = iH + 1.d0
      HRA(nint(iH)) = p
      ! do no more, unless it is the end of the day:

      if(nint(iH) == nmax) then ! end of "day":
        iH = 0.d0
        if(nint(first) == 1) then ! first averaging period only
          iD = iD + 1.d0
          do n=1,nmax
            DRA(nint(iD)) = DRA(nint(iD)) + HRA(n)
          end do
          DRA(nint(iD)) = DRA(nint(iD))*bynmax
          if(nint(iD) == nday) then ! end first period
            PRS = 0.d0
            do n=1,nday
              PRS = PRS + DRA(n)
            end do
            avg = PRS * bynday
            first=0.d0
            iD=0.d0
            i0=0.d0
          end if
        else ! not first averaging period: update the running average
          i0 = i0 + 1.d0 ! move marker
          if(nint(i0) == nday+1) i0=1.d0 ! reset marker
          temp=0.d0
          do n=1,nmax
            temp = temp + HRA(n)
          end do
          temp = temp * bynmax ! i.e. today's average
          PRS = PRS - DRA(nint(i0))
          DRA(nint(i0)) = temp
          PRS = PRS + DRA(nint(i0))
          avg = PRS * bynday
        end if
      end if

      end subroutine prec_running_average

      subroutine step_ba(burnt_area,RH1,wsurf,saveFireCount,
     &                    pvt,fearth_axyp,i,j)
!@sum step_ba calculates the burnt area of the model at every time step
!@+ burnt area variable for use in the flammability model.
!@auth Keren Mezuman 
      use ent_const, only : N_COVERTYPES
      use ent_pfts, only: ent_cover_names
      use MathematicalConstants_mod, only: PI
      use TimeConstants_mod, only: SECONDS_PER_DAY,SECONDS_PER_YEAR
      use model_com, only : DTsrc
      use diag_com, only: ij_a_tree,ij_a_shrub,ij_a_grass,aij=>aij_loc
 
      implicit none
      
      real*8 :: inst_FC
      real*8 :: RHlow,RHup,CRH,Cm,Cb,Lb,Hb,gW,g0,nv,up,a,tau,N_steps
      real*8 :: new_burnt_area,previous_BA,pvt_area
      real*8, Dimension(N_COVERTYPES), intent(INOUT) :: burnt_area 
      real*8, intent(IN) :: RH1,wsurf,saveFireCount,fearth_axyp
      real*8, Dimension(N_COVERTYPES), intent(IN) :: pvt 
      real*8 :: umax
      integer, intent(in) :: i,j
      !if there are no fires to create new BA and no old BA to recover
      if ((saveFireCount <= 0.d0) .AND. (sum(burnt_area) <= 0.d0)) 
     &  return
!@var RHlow lower bound of RH for fire spread (fraction)
!@var RHup upper bound of RH for fire spread (fraction)
!@var CRH response of fuel combustibility to real-time 
!@+   climate conditions (unitless)
!@var Cb root zone soil wetness (unitless)
!@var g0 dependance of fire spread perpendicular to wind 
!@var tau average fire duration (seconds/fire)
      !+direction (unitless)
      RHlow=0.3d0
      RHup=0.7d0
      Cb=0.5d0
      g0=0.05d0
      tau=SECONDS_PER_DAY
      N_steps=SECONDS_PER_DAY/DTsrc

!@var saveFireCount fire count rate (fire/m2/s)
!@var DTSRC source time step (s) 
!@var inst_FC number of fires in a time step in 
      !+ a square meter (#fires/m^2)
      inst_FC = saveFireCount*DTsrc
      if(RH1 <= RHlow) then
        CRH=1.d0
      else if ((RH1 > RHlow) .AND. (RH1 < RHup)) then
        CRH = (RHup - RH1) / (RHup - RHlow)
      else
        CRH = 0.d0
      end if
!@var Cm dependence of downwind on fuel wetness (Li et al. 2018)
      Cm=Cb*CRH
!@var Lb length-to-breadth ratio (unitless)
!@var wsurf surface wind velocity (m/s)
      Lb=1. + 10. * (1. - (EXP(-0.06 * wsurf)))
!@var Hb head-to- back ratio (unitless)
      Hb= (Lb + (Lb**2 - 1.)**0.5) / (Lb - (Lb**2 - 1.)**0.5)
!@var the dependence of fire spread on wind speed
      gW=2*Lb/(1.+1./Hb)*g0
      do nv=1,N_COVERTYPES
        !if there are no fires to create new BA and no old BA to recover
        if ((saveFireCount <= 0.d0) .AND. (burnt_area(nv) <= 0.d0)) 
     &    cycle
        !umax: average maximum fire spread rate [m/s]
      select case(ent_cover_names(nv))
        case('arid_shrub')
          umax=0.17d0!(m/s)
        case('c3_grass_ann')
          umax=0.2d0!(m/s)
        case('c3_grass_arct')
          umax=0.2d0!(m/s)
        case('c3_grass_per')
          umax=0.2d0!(m/s)
        case('c4_grass')
          umax=0.2d0!(m/s)
        case('cold_br_late')
          umax=0.11d0!(m/s)
        case('cold_shrub')
          umax=0.17d0!(m/s)
        case('decid_nd')
          umax=0.15d0!(m/s)
        case('drought_br')
          umax=0.11d0!(m/s)
        case('ever_br_late')
          umax=0.11d0!(m/s)
        case('ever_nd_late')
          umax=0.15d0!(m/s)
        case default
          umax=0.d0!(m/s)
      end select
!@var up fire spread rate in the downwind direction (m/s)
        up=umax*Cm*gW
        if (Lb == 0.d0) then
          a=0.d0!(m)
        else
!@var a average fire spread area (m^2)
          a = PI * up**2 * tau**2 / (4*Lb) * (1+1/Hb)**2 / N_steps
        endif
        pvt_area = pvt(nv) * fearth_axyp
        new_burnt_area = inst_FC * a * pvt_area * 250.d0
        burnt_area(nv) =  min(new_burnt_area , pvt_area) 
      select case(ent_cover_names(nv))
        case('arid_shrub')
          aij(i,j,ij_a_shrub) = aij(i,j,ij_a_shrub) + a 
        case('c3_grass_ann')
          aij(i,j,ij_a_grass) = aij(i,j,ij_a_grass) + a 
        case('c3_grass_arct')
          aij(i,j,ij_a_grass) = aij(i,j,ij_a_grass) + a 
        case('c3_grass_per')
          aij(i,j,ij_a_grass) = aij(i,j,ij_a_grass) + a 
        case('c4_grass')
          aij(i,j,ij_a_grass) = aij(i,j,ij_a_grass) + a 
        case('cold_br_late')
          aij(i,j,ij_a_tree) = aij(i,j,ij_a_tree) + a 
        case('cold_shrub')
          aij(i,j,ij_a_shrub) = aij(i,j,ij_a_shrub) + a 
        case('decid_nd')
          aij(i,j,ij_a_tree) = aij(i,j,ij_a_tree) + a 
        case('drought_br')
          aij(i,j,ij_a_tree) = aij(i,j,ij_a_tree) + a 
        case('ever_br_late')
          aij(i,j,ij_a_tree) = aij(i,j,ij_a_tree) + a 
        case('ever_nd_late')
          aij(i,j,ij_a_tree) = aij(i,j,ij_a_tree) + a 
        case default
      end select
      end do

      end subroutine step_ba 

      subroutine lai_running_average(lai,avg,iH,iD,i0,first,HRA,DRA,PRS)
!@sum lai_running_average keeps a running average of the Ent leaf area 
!@+ index variable for use in the flammability model.
!@+ In practice, this does hourly and daily running averages and
!@+ uses those to get the period-long running average, to avoid
!@+ saving a huge array.
!@auth Greg Faluvegi
      use flammability_com, only: nday=>nday_lai,nmax=>maxHR_lai
      use constant, only: undef

      implicit none

!@var lai model variable which will be used in the average
!@var temp just for holding current day average for use in avg_prec
!@var nmax number of accumulations in one day (for checking)
!@var bynmax reciprocal of nmax
!@var bynday reciprocal of nday
!@var i0 see i0lai(i,j) (day in period marker)
!@var iH see iHlai(i,j) (hour in day index)
!@var iD see iDlai(i,j) (day in period index)
!@var first see first_lai(i,j) whether in first period
!@var HRA see HRAlai(time,i,j) (hourly average)
!@var DRA see DRAlai(time,i,j) (daily average)
!@var PRS see PRSlai(i,j) (period running sum)
!@var avg see ravg_lai(i,j) the running average prec returned
      real*8, intent(IN) :: lai
      real*8, dimension(nday) :: DRA
      real*8, dimension(nmax) :: HRA
      real*8 :: temp, bynmax, PRS, avg, iH, iD, i0, first, bynday
      integer :: n
      
      if(lai == undef) then
        write(6,*) 'undefined LAI found in lai_running_average.'
        call stop_model('undef LAI in lai_running_average',255)
      endif
      if(nint(iH) < 0 .or. nint(iH) > nmax) then
        write(6,*) 'iH maxHR_lai=',iH,nint(iH),nmax
        call stop_model('iHlai or maxHR_lai problem',255)
      endif
      bynmax=1.d0/real(nmax)
      bynday=1.d0/real(nday)
      iH = iH + 1.d0
      HRA(nint(iH)) = lai
      ! do no more, unless it is the end of the day:

      if(nint(iH) == nmax) then ! end of "day":
        iH = 0.d0
        if(nint(first) == 1) then ! first averaging period only
          iD = iD + 1.d0
          do n=1,nmax
            DRA(nint(iD)) = DRA(nint(iD)) + HRA(n)
          end do
          DRA(nint(iD)) = DRA(nint(iD))*bynmax
          if(nint(iD) == nday) then ! end first period
            PRS = 0.d0
            do n=1,nday
              PRS = PRS + DRA(n)
            end do
            avg = PRS * bynday
            first=0.d0
            iD=0.d0
            i0=0.d0
          end if
        else ! not first averaging period: update the running average
          i0 = i0 + 1.d0 ! move marker
          if(nint(i0) == nday+1) i0=1.d0 ! reset marker
          temp=0.d0
          do n=1,nmax
            temp = temp + HRA(n)
          end do
          temp = temp * bynmax ! i.e. today's average
          PRS = PRS - DRA(nint(i0))
          DRA(nint(i0)) = temp
          PRS = PRS + DRA(nint(i0))
          avg = PRS * bynday
        end if
      end if

      end subroutine lai_running_average


#if defined DYNAMIC_BIOMASS_BURNING && defined CALCULATE_FLAMMABILITY

      subroutine calculate_fire_count
!@sum calculate_fire_count calculated the #fires rate for the
!@+ dynamic biomass burning sources.
!@auth Greg Faluvegi based on direction from Olga Pechony
!@+ later modified by Keren Mezuman
      use model_com, only : DTsrc
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR,DAYS_PER_YEAR,
     & SECONDS_PER_DAY
      use domain_decomp_atm,only: grid, getDomainBounds
      use constant, only: undef
      use flammability_com, only: mfcc,flammability,first_prec,
     & saveFireCount
      use diag_com, only: ij_fireC,aij=>aij_loc
#ifdef ANTHROPOGENIC_FIRE_MODEL
      use lightning, only : CG_DENS 
      use flammability_com, only: populationDensity, nonSuppressFrac
      use diag_com, only: ij_cgign,ij_humanign
#endif
#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_ngroups,subdd_type
     &     ,inc_subdd,find_groups
#endif
      implicit none

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups),k
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
#endif
      integer :: J_0S, J_1S, I_0H, I_1H, i, j
      real*8 :: CtoG, humanIgn, tuneToMODIS, conv,
     & monthPerSecond,yearsPerSecond
!@var CtoG local copy of cloud-to-ground lightning strikes
!@var humanIgn the human-induced fire ignition rate (before 
!@+ supression in units of #/m2/sec)
!@var nonSuppressFrac the fraction of fire ignitions not supressed 
!@var tuneToMODIS a tuning factor of the fire count to MODIS obs.

! Care should be taken if one were to move the call to this routine.
! For example: I think the flammability is calculated after
! the tracer 2D sources are set in the main loop. This is probably 
! OK, since that shouldn't change very quickly, but the fire count is a
! function of the lightning, so in the main loop, it should be called
! after CONDSE but before the tracer 2D sources are set. Similarly, if
! the fire count were filled after it is used in the main loop, then
! one would have to save it to the rsf file, like flammability( ) is.

      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &                          I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      yearsPerSecond=1.d0/(SECONDS_PER_DAY*DAYS_PER_YEAR)
      monthPerSecond=real(INT_MONTHS_PER_YEAR,kind=8)*yearsPerSecond

      do j=J_0S,J_1S
        do i=I_0H,I_1H

          ! only do calculation after enough precip averaging done,
          ! and where flammability is defined:
          if(first_prec(i,j)==0 .and. flammability(i,j)/=undef) then
#ifdef ANTHROPOGENIC_FIRE_MODEL
            ! Anthropogenic/lightning fire model ignition/supression, based on Olga's  
            ! document: "Anthropogenic ignitions and supression.docx" Nov 2012.
            ! First, the lightning-induced portion, where CtoG is the total 
            ! cloud-to-ground lightning flashes per m2 per second, so that the
            ! fire count will also be in units of #fire. Starting with 
            ! CG_DENS in flashes/m2/s:
            ! a factor of 0.1 is necessary to reduce overestimation compared to WWLLN data
            CtoG=CG_DENS(i,j)*0.1d0 ! #/m2/s

            ! Human ingition portion: The population density units are humans/km2, 
            ! and the formula then puts the human ingition rate in #/km2/month. Thus
            ! we need a conversion factor (conv) to go to #/m2/s:
            ! 1.d-6 = 1km/1000m * 1km/1000m 
            conv=1.d-6*monthPerSecond !(km^2 * M)/(m^2 * s)
            humanIgn=conv*0.2d0*populationDensity(i,j)**(0.4) ! #/m2/s


            ! Olga says: "While theoretically, this is supposed to give the absolute number
            ! of fire counts, this is not actually so, since (A) We don't know the real 
            ! lightning flash rate and we don't know the real number of fire human-caused
            ! ignitions and so we can not really calibrate the number of fire ignitions 
            ! to reflect the true absolute values.  (B) We don't know the real fire counts
            ! and we don't know if MODIS fire counts are any closer to them in absolute values. 
            ! And since we want to use the EPFCs derived relying on the MODIS fire counts, we'll
            ! still need to calibrate the modeled fire counts to be comparable with the MODIS
            ! absolute values." The following is that tuning parameter:
            tuneToMODIS=30.d0

            ! Putting that all together to get the fire count rate (fire/m2/s):
            ! change variable name of saveFireCount to something better
            saveFireCount(i,j)=tuneToMODIS*
     &       flammability(i,j)*(CtoG+humanIgn)*nonSuppressFrac(i,j)

            ! Save a daignostic for the portion that is human-caused. (1.0-this) is the
            ! portion that is lightning-caused, so no reason to save that. Also save the
            ! fire count:
            aij(i,j,ij_humanign)=aij(i,j,ij_humanign)+humanIgn
            aij(i,j,ij_cgign)=aij(i,j,ij_cgign)+CtoG

#else /* ubiquitous only */

            ! Ubiquitous fire model. Note on units:
            ! flammability*mfcc*yearsPerSecond = [#fire/m2]
            ! mfcc*yearsPerSecond includes DTsrc in it
            ! yearsPerSecond = [yr/s]
            ! saveFireCount = [#fire/m2]
            saveFireCount(i,j)=
     &      flammability(i,j)*mfcc*yearsPerSecond
#endif /* anthro vs ubiquitous */
          else ! flammability not ready yet or undefined here:
            saveFireCount(i,j)=0.d0
          end if
          aij(i,j,ij_fireC)=aij(i,j,ij_fireC)+saveFireCount(i,j)
        end do ! i
      end do   ! j
#ifdef CACHED_SUBDD
****
**** Collect some high-frequency outputs
      call find_groups('aijh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('fireCount')
            sddarr2d =saveFireCount(:,:)
            call inc_subdd(subdd,k,sddarr2d)
          case ('f_ignCG')
            sddarr2d =CG_DENS(:,:)*0.1d0
            call inc_subdd(subdd,k,sddarr2d)
          end select
        enddo
      enddo
#endif

      end subroutine calculate_fire_count


      subroutine dynamic_biomass_burning(n,ns)
!@sum dynamic_biomass_burning fills in the surface source ns for
!@+ tracer n with biomass burning based on flammability, offline
!@+ correlations of emissions with observed fire counts by vegetation
!@+ type, and the GCM's online vegetation. For now, this is mapped
!@+ onto the traditional VDATA( ) types, even if Ent is on.
!@auth Greg Faluvegi based on direction from Olga Pechony, Igor A.

      use ent_const, only : N_COVERTYPES
      use domain_decomp_atm,only: grid, getDomainBounds
      use flammability_com, only: flammability,first_prec,
     & saveFireCount,EPFCByVegType
      use constant, only: undef
      use tracer_com, only: sfc_src
      use ghy_com, only: fearth
      use geom, only : axyp
      use ent_com, only: entcells
      use ent_mod, only: ent_get_exports
     &                   ,n_covertypes !YKIM-temp hack

      implicit none
   
      integer :: J_0S, J_1S, I_0H, I_1H, i, j, nv
      integer, intent(in) :: n,ns
!@var emisPerFire emission per fire count, generally kg/m2/fire
      real*8 :: emisPerFire
!@var pvt fraction vegetation type for N_COVERTYPES (fraction)
!@var N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER = 16 + 2 + 0
      real*8, dimension(N_COVERTYPES):: pvt, EPFBVT

      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &              I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

!@var EPFCByVegType emission factor (kg/#fire)
      EPFBVT=EPFCByVegType(:,n)

      do j=J_0S,J_1S
        do i=I_0H,I_1H

          sfc_src(i,j,n,ns)=0.d0
          if(fearth(i,j)==0.d0) cycle
          ! only do calculation after enough precip averaging done,
          ! and where flammability is defined:
          if(first_prec(i,j)==0 .and. flammability(i,j)/=undef) then
            ! Obtain the vegetation types in the box:
            ! For now, the same way RAD_DRV does it, as per Greg F.'s 
            ! e-mails with Igor A. Mar-Apr,2010:
              call ent_get_exports(entcells(i,j),
     &                             vegetation_fractions=pvt)
            ! Notes on units:
            ! sfc_src = [kg/m2/s] gridbox
            ! emisPerFire = [kg/#fire]
            ! EPFBVT = [kg/#fire]
            ! saveFireCount = [#fire/m2/sec] gridbox
    
            ! construct emisPerFire from EPFCByVegType:
            emisPerFire = 0.d0
            do nv=1,N_COVERTYPES
              emisPerFire = emisPerFire + pvt(nv)*fearth(i,j)*EPFBVT(nv)
            end do
            sfc_src(i,j,n,ns) = emisPerFire*saveFireCount(i,j)
          end if
        end do ! i
      end do   ! j
    
      end subroutine dynamic_biomass_burning

#ifdef DETAILED_FIRE_OUTPUT
      subroutine accumulateVegTypesDiag
!@sum accumulateVegTypesDiag for updating the optional vegetation cover
!@+ by vegetation types diagnostic for the fire model.
!@auth Greg Faluvegi

      use ent_const, only : N_COVERTYPES
      use domain_decomp_atm,only: grid, getDomainBounds
      use flammability_com, only: ij_flamV
      use ghy_com, only: fearth

      use diag_com, only: ij_flam,aij=>aij_loc

      use ent_com, only: entcells
      use ent_mod, only: ent_get_exports
     &                   ,n_covertypes !YKIM-temp hack
      use ent_pfts, only: ent_cover_names

      implicit none

      integer :: J_0S, J_1S, I_0H, I_1H, i, j, nv
!@var pvt percent vegetation type for 12 VDATA types (per ice-free land)
      real*8, dimension(N_COVERTYPES):: pvt
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &                           I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      do j=J_0S,J_1S
        do i=I_0H,I_1H
          ! Obtain the vegetation types in the box:
          ! For now, the same way RAD_DRV does it, as per Greg F.'s 
          ! e-mails with Igor A. Mar-Apr,2010:
          if(fearth(i,j)>0.d0) then
            call ent_get_exports(entcells(i,j), 
     &           vegetation_fractions=pvt)
          else
            pvt(:) = 0.d0
          end if
          do nv=1,N_COVERTYPES
          select case(ent_cover_names(nv))
            case('arid_shrub')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('cold_shrub')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('c3_grass_ann')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('c3_grass_arct')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('c3_grass_per')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('c4_grass')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('cold_br_late')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('decid_nd')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('drought_br')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('ever_br_late')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
            case('ever_nd_late')
              aij(i,j,ij_flamV)=aij(i,j,ij_flamV)+
     &        fearth(i,j)*pvt(nv)
          end select
          end do
        end do ! i
      end do   ! j

      end subroutine accumulateVegTypesDiag

#endif /* DETAILED_FIRE_OUTPUT */
#endif /* defined DYNAMIC_BIOMASS_BURNING && defined CALCULATE_FLAMMABILITY */

#if defined DYNAMIC_BIOMASS_BURNING && defined ANTHROPOGENIC_FIRE_MODEL

      subroutine readFlamPopDens(xyear,xday)
!@sum reads 2D human population density for flammability purposes
!@auth Greg Faluvegi
!@+ Later modified by Keren Mezuman 
      use geom, only : lon_to_i,lat_to_j
      use domain_decomp_atm, only: GRID,getDomainBounds,readt_parallel, 
     & write_parallel,rewind_parallel
      use filemanager, only: openunit,closeunit,nameunit !,is_fbsa
      use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
      use timestream_mod, only : init_stream,read_stream
      use flammability_com, only: populationDensity,flamPopB,flamPopA,
     & firstFlamPop,flamPopYearStart,flamPopYearEnd,flamPopDelYear,
     & nonSuppressFrac
      use diag_com, only: ij_nsuppress,aij=>aij_loc
 
      implicit none

      integer :: iuPopDen,k,ipos,kx,kstep=10,half
      character*80 :: fname,title
      character(len=300) :: out_line
      real*8 :: alpha
      integer, intent(IN) :: xyear, xday

      SAVE iuPopDen ! flamPopB,flamPopA, etc. saved in module

      integer :: J_1, J_0, J_0H, J_1H, I_0, I_1, j, i

      fname='FLAMPOPDEN'

      half=NINT((EARTH_DAYS_PER_YEAR+1.)/2.)

!     if(.not.is_fbsa(fname)) then
!       call stop_model('population netCDF input not implemented.',255) 
!     else

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)


      ! If first time through after restart, open file and read its header:
      if(firstFlamPop) then
        call openunit(fname,iuPopDen,.true.)
        read(iuPopDen)
     &  title,flamPopYearStart,flamPopYearEnd,flamPopDelYear
      end if

      ! now read the data: (should be in humans/km2, please)

! -------------- non-transient case ----------------------------!
      if(flamPopYearStart==flamPopYearEnd)then
        if(firstFlamPop) then 
          call readt_parallel
     &    (grid,iuPopDen,fname,populationDensity(:,:),1)
          write(out_line,*)'Population Density for Fire Model Read.'
          call write_parallel(trim(out_line))
          firstFlamPop = .false.
          call closeunit(iuPopDen)
        end if

! --------------- transient case -------------------------------!
      else ! only annual files allowed; only read first time + new steps
        kstep=flamPopDelYear
        ipos=1
        alpha=0.d0 ! before start year, use start year value
        kx=flamPopYearStart 
        if(xyear>flamPopYearEnd .or.
     &    (xyear==flamPopYearEnd.and.xday>=half))then
          alpha=1.d0 ! after end year, use end year value     
          ipos=(flamPopYearEnd-flamPopYearStart)/kstep
          kx=flamPopYearEnd-kstep
        endif
        do k=flamPopYearStart,flamPopYearEnd-kstep,kstep
          if(xyear==k .and. xday==half)firstFlamPop=.true. ! tell model need to read
          if(xyear>k .or. (xyear==k.and.xday>=half)) then
            if(xyear<k+kstep.or.(xyear==k+kstep.and.xday<half))then
              ipos=1+(k-flamPopYearStart)/kstep ! (integer artithmatic)
              alpha=(EARTH_DAYS_PER_YEAR*(0.5+real(xyear-1-k))+xday) /
     &              (EARTH_DAYS_PER_YEAR*real(kstep))
              kx=k
              exit
            endif
          endif
        enddo
        ! if time to read, update interpolation arrays:
        if(firstFlamPop) then
          call rewind_parallel(iuPopDen)
          call readt_parallel(grid,iuPopDen,fname,flamPopA(:,:),ipos+1) ! +1 for header
          call readt_parallel(grid,iuPopDen,fname,flamPopB(:,:),1)
          firstFlamPop = .false.
        end if
        populationDensity(I_0:I_1,J_0:J_1)=flamPopA(I_0:I_1,J_0:J_1)*
     &   (1.d0-alpha)+flamPopB(I_0:I_1,J_0:J_1)*alpha

        write(out_line,'("FireModel PopDens at ",F9.4,a16,I4,a8,I4)')
     &  100.d0*alpha,'% of period mid ',kx,' to mid ',kx+kstep
        call write_parallel(trim(out_line))
      end if ! tras/non-trans
!calculate default fns 
!loop over i,j 
!in the loop if in the regions overwrite 
            ! Fraction not supressed by humans (unitless):
      do j=J_0,J_1
        do i=I_0,I_1
          !Temperate North America (USA) 
          if ((j>lat_to_j(29.d0)) .AND. (j<lat_to_j(49.d0)) .AND.
     &      (i>lon_to_i(-128.75d0)) .AND. (i<lon_to_i(-66.25d0))) then
            nonSuppressFrac(i,j)=0.2d0*exp(-0.05*populationDensity(i,j))
          !Middle East 
          else if ((j>lat_to_j(21.d0)) .AND. (j<lat_to_j(37.d0)) .AND.
     &      (i>lon_to_i(-18.75d0)) .AND. (i<lon_to_i(28.75d0))) then
            nonSuppressFrac(i,j)=0.2d0*exp(-0.05*populationDensity(i,j))
          else if ((j>lat_to_j(15.d0)) .AND. (j<lat_to_j(43.d0)) .AND.
     &      (i>lon_to_i(26.25d0)) .AND. (i<lon_to_i(71.25d0))) then
            nonSuppressFrac(i,j)=0.2d0*exp(-0.05*populationDensity(i,j))
          !Northern Hemisphere Africa 
          else if ((j>lat_to_j(-1.d0)) .AND. (j<lat_to_j(23.d0)) .AND.
     &      (i>lon_to_i(-21.25d0)) .AND. (i<lon_to_i(51.25d0))) then
            nonSuppressFrac(i,j)=1.d0
          !Southern Hemisphere Africa 
          else if ((j>lat_to_j(-37.d0)) .AND. (j<lat_to_j(1.d0)) .AND.
     &      (i>lon_to_i(0.d0)) .AND. (i<lon_to_i(53.75d0))) then
            nonSuppressFrac(i,j)=1.d0
          else
            nonSuppressFrac(i,j)=
     &       0.05d0+0.9d0*exp(-0.05*populationDensity(i,j))
          end if
          aij(i,j,ij_nsuppress)=nonSuppressFrac(i,j)
        enddo
      enddo

!     endif ! fbsa format or not
      return
      ! this keeps transient files open. maybe study TracerSurfaceSource.F90/
      ! subroutine readSurfaceSource on how to not do that. 
      end subroutine readFlamPopDens

#endif /* DYNAMIC_BIOMASS_BURNING && defined ANTHROPOGENIC_FIRE_MODEL */

