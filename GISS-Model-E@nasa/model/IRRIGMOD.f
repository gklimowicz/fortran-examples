#include "rundeck_opts.h" 
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
      module irrigmod

!@sum  Module irrigmod contains the arrays/subroutines needed to prescribe
!@+    irrigation rates from input files.
!@auth M.J. Puma
      use timestream_mod, only : timestream
      implicit none
      save

!@dbparam irrig_yr year for irrigation data (if 0: time var)
      INTEGER :: irrig_yr = 2000 !placeholder value, should be read in

!@dbparam year to start adding "fossil" groundwater to system
      INTEGER :: yr_fossilGW = 1900 

!@var potential monthly (gross) irrigation water requirements (m/s)
      real*8, dimension(:,:), allocatable :: irrig_water_pot

!@var IRRIGstream interface for reading and time-interpolating irrigation files
!@+   See general usage notes in timestream_mod.
      type(timestream) :: IRRIGstream

      logical :: irrig_exists=.true.

#ifdef TRACERS_WATER
      real*8, dimension(:), allocatable :: tr_conc_in_aquifer
#endif

      contains

      subroutine init_irrigmod
!@sum init_irrigmod initializes the IRRIGstream object
      use Dictionary_mod, only : sync_param,get_param
      use domain_decomp_atm, only : grid,getDomainBounds
      USE model_com,ONLY : modelEclock,master_yr
      use timestream_mod, only : init_stream
      use dictionary_mod, only : get_param,is_set_param
      use filemanager, only : file_exists
      implicit none

      integer :: i_0h,i_1h,j_0h,j_1h,ier
      integer :: jyear,jday,year_start, year_end
      logical :: cyclic

      irrig_exists = file_exists('IRRIG')
      if(.not.irrig_exists) return


      call getDomainBounds(grid,j_strt_halo=j_0h,j_stop_halo=j_1h)
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      allocate(irrig_water_pot(i_0h:i_1h,j_0h:j_1h))


      ! If parameter irrig_yr exists, IRRIG data from that year is
      ! selected  (if 0: time var, else cyclic)
      ! Otherwise, irrig_yr is set to master_yr.
      call get_param( "master_yr", master_yr )
      call get_param( "irrig_yr", irrig_yr, default=master_yr )

      cyclic = irrig_yr /= 0 ! irrig_yr==0 implies transient mode.
      irrig_yr = abs(irrig_yr)


      call modelEclock%get(year=jyear, dayOfYear=jday)

      if(cyclic)jyear = irrig_yr
         
      if (jyear > 2100 .or. jyear < 1848) then 
         call stop_model("No irrigation for that yr;turn off",255)
      endif
!------
      ! Check if irrigation data is available all years in the run
!      if(cyclic)then 
!         jyear = irrig_yr
!         if (jyear > 2100 .or. jyear < 1848) then 
!            call stop_model("No irrigation for that yr;turn off",255)
!         endif
!
!      else
!         call modelEclock%get(YEARI=year_start, YEARE=year_end)
!
!         if (year_end > 2100 .or. year_start < 1848) then 
!            call stop_model("No irrigation for yr range;turn off",255)
!         endif
!      endif
!-----------
      call init_stream(grid,IRRIGstream,'IRRIG','irrigation_per_m2',
     &              0d0,1000d0,'linm2m',jyear,jday,cyclic=cyclic)      

      call read_irrig(.false.)

#ifdef TRACERS_WATER
      call init_irrigmod_tracers
#endif

      end subroutine init_irrigmod


#ifdef TRACERS_WATER
      subroutine init_irrigmod_tracers
      USE TRACER_COM, only : ntm, n_H2O18, n_H2O17, n_HDO, n_HTO, 
     &                       n_Water
#ifdef TRACERS_SPECIAL_O18
      real*8, external :: water_iso_conc_in_aquifer
#endif
      integer :: n

      allocate(tr_conc_in_aquifer(ntm))

      ! by defaut nothing is present in aquifer water except regular water
      tr_conc_in_aquifer(:) = 0.0
      if ( n_Water > 0 ) then
        tr_conc_in_aquifer(n_Water) = 1.d0
      endif

#ifdef TRACERS_SPECIAL_O18
      do n=1,ntm
        if ( n==n_H2O18 .or. n==n_H2O17 .or. 
     &       n==n_HDO .or. n==n_HTO ) then
          tr_conc_in_aquifer(n) = water_iso_conc_in_aquifer(n)
        endif
      enddo
#endif

      end subroutine init_irrigmod_tracers
#endif


      subroutine read_irrig(end_of_day)
!@sum read_irrig invokes procedures to read irrigation rates from
!@+   input files and perform time interpolation
      use Dictionary_mod, only : sync_param
      use dictionary_mod, only : get_param,is_set_param
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use domain_decomp_atm, only : getDomainBounds,grid
      use model_com,only : itime,itimei,master_yr
      use model_com, only :  modelEclock
      use resolution, only : im,jm
      use geom, only : imaxj
      use timestream_mod, only : read_stream
      use DIAG_COM, only : aij=>aij_loc, ij_irrW_tot
      implicit none

!@var end_of_day flag to determine whether this is an initial or daily call
      logical, INTENT (IN) :: end_of_day ! true if called from daily
      real*8 :: tfo
      integer i,j
      integer :: jyear,jday

      integer :: j_0,j_1, i_0,i_1
      logical :: have_north_pole, have_south_pole
      logical :: cyclic

      if(.not.irrig_exists) return

      call getDomainBounds(grid,i_strt=i_0,i_stop=i_1,
     &         j_strt=j_0,j_stop=j_1,
     &         have_south_pole=have_south_pole,
     &         have_north_pole=have_north_pole)


      ! If parameter irrig_yr exists, IRRIG data from that year is
      ! selected  (if 0: time var, else cyclic)
      ! Otherwise, irrig_yr is set to master_yr.
      call get_param( "master_yr", master_yr )
      call get_param( "irrig_yr", irrig_yr, default=master_yr )
 
      cyclic = irrig_yr /= 0 ! irrig_yr==0 implies transient mode.
      irrig_yr = abs(irrig_yr)


      call modelEclock%get(year=jyear, dayOfYear=jday)
      if(cyclic) jyear = irrig_yr

      if (jyear > 2100 .or. jyear < 1848) then 
         call stop_model("No irrigation for that yr;turn irrig off",255)
      endif

      call read_stream(grid,IRRIGstream,jyear,jday,irrig_water_pot)


!**** Replicate values at pole (not relevant for present-day Earth)
      if(have_north_pole) then
        irrig_water_pot(2:,jm)=irrig_water_pot(1,jm)
      end if
      if(have_south_pole) then
        irrig_water_pot(2:,1)=irrig_water_pot(1,1)
      end if

!**** Make sure no negative irrigation and update diagnostic
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        if (irrig_water_pot(i,j).lt.0) irrig_water_pot(i,j)=0

!**** Diagnostic
        if (end_of_day) aij(i,j,ij_irrW_tot)=aij(i,j,ij_irrW_tot)+
     &          irrig_water_pot(i,j)*SECONDS_PER_DAY*1000.d0

      end do
      end do

      return
      end subroutine read_irrig


      subroutine irrigate_extract(i,j,MWL,GML,MLDLK,tlake,flake
     *     ,hlake_min,MWL_to_irrig,GML_to_irrig,irrig_gw,irrig_gw_energy
     *     ,irrig_water_act,irrig_energy_act
#ifdef TRACERS_WATER
     *     ,TRML,TRML_to_irrig,irrig_tracer_act,irrig_gw_tracer
#endif
     *       )
!@sum irrigate_extract calculate water used for irrigation
!@auth Michael Puma
      USE CONSTANT, only : rhow,teeny,shw
      USE MODEL_COM, only : dtsrc
      use model_com, only :  modelEclock
!      USE sle001, only : tp ! tp is not saved for each gridpoint
      USE GHY_COM, only : tearth
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
      USE FLUXES, only : atmlnd
#endif
! fixed i,j arrays - feed in from call?
      USE GEOM, only : axyp
      USE GHY_COM, only : fearth
      use dictionary_mod, only : get_param
      implicit none
      integer,intent(in):: i,j
      integer :: jyear,master_yr
      logical :: cyclic

C**** Inputs
!@var MWL,GML,MLDLK,tlake,flake local versions of lake variables
      REAL*8, INTENT(IN):: MWL,GML,MLDLK,tlake,flake,hlake_min
#ifdef TRACERS_WATER
     *     ,TRML(NTM,2)
#endif
C**** Outputs
!@var MWL_to_irrig (kg), GML_to_irrig (J) mass/energy changes
      REAL*8, INTENT(OUT) :: MWL_to_irrig,GML_to_irrig
#ifdef TRACERS_WATER
!@var TRML_to_irrig (kg) tracer changes
     *     ,TRML_to_irrig(NTM,2)
#endif
!@var irrig_gw (m/s), irrig_gw_energy (W/m2) ground water diagnostics 
!@var irrig_water_act (m/s), irrig_energy_act (W/m2) actual irrigation diagnostics 
      REAL*8, INTENT(OUT) :: irrig_gw,irrig_gw_energy,irrig_water_act
     *     ,irrig_energy_act
#ifdef TRACERS_WATER
!@var irrig_tracer_act (kg/s) actual irrigation tracer diagnostics 
!@var irrig_gw_tracer_act (kg/s) implied ground water tracer diagnostics 
     *     ,irrig_tracer_act(NTM),irrig_gw_tracer(NTM)
#endif

!@var Mass of irrigation water withdrawn from rivers/groundwater [kg]
      real*8 :: m_irr_pot
!@var Mass of water available for irrigation [kg]
      real*8 :: m_avail
!@var Temperature of abstracted irrigation water [deg C]
      real*8 :: T_irr, T_irr_g, T_irr2
!@var Vegetated and bare fractions of the grid cell
      real*8 :: fv, fb
!@param Conversion factor from m^3 water/(m^2 grid cell) to kg water
      real*8 :: m_to_kg
!@param Flag for external irrigation source (from deep aquifers - not modeled) 
      integer :: flag_irrig_grndwat  ! =1 use ground water

C**** set default output
      irrig_water_act = 0.d0 ; irrig_energy_act = 0.d0
      irrig_gw = 0.d0 ; irrig_gw_energy = 0.d0
      MWL_to_irrig = 0 ; GML_to_irrig = 0
#ifdef TRACERS_WATER
      irrig_tracer_act = 0 
      irrig_gw_tracer = 0 
      TRML_to_irrig = 0
#endif

! ->  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ->  begin block copied here by MK for correctness/consistency
! ->  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If parameter irrig_yr exists, IRRIG data from that year is
      ! selected  (if 0: time var, else cyclic)
      ! Otherwise, irrig_yr is set to master_yr.
      call get_param( "master_yr", master_yr )
      call get_param( "irrig_yr", irrig_yr, default=master_yr )
      cyclic = irrig_yr /= 0 ! irrig_yr==0 implies transient mode.
      irrig_yr = abs(irrig_yr)
      call modelEclock%get(year=jyear)
      if(cyclic)jyear = irrig_yr
! ->  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ->  end block copied here.  Todo: set effective year info in 1 place rather than 3!
! ->  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Determine if we add water to the system (not conserving water) 
      if (jyear <yr_fossilGW) then
         flag_irrig_grndwat = 0 ! Don't add 'fossil' groundwater
      else
         flag_irrig_grndwat = 1 ! Add 'fossil' groundwater
      endif

!    Unit conversion factor
      m_to_kg = rhow*axyp(i,j)

      call get_fb_fv( fb, fv, i, j )

!***  Irrigation based on potential values in m/s per area grid cell 
      if ( (irrig_water_pot(i,j) > teeny)  .and. 
     &     (fearth(i,j)          > teeny)  .and.
     &     (fv                   > teeny)) then

         m_irr_pot = irrig_water_pot(i,j)* m_to_kg *dtsrc

         if (flake > 0.d0) then
            m_avail = mwl - hlake_min*flake*m_to_kg
            m_avail = max(m_avail, 0.d0)
            T_irr = tlake
            T_irr2 =0.
            if (mwl.gt.flake*mldlk*m_to_kg+teeny) T_irr2 = 
     *           (gml-mldlk*m_to_kg*flake*tlake*shw)/
     *           (mwl-mldlk*m_to_kg*flake+teeny)/shw 
         else
            m_avail = mwl
            T_irr = gml/(mwl*shw+teeny)
            T_irr2 = T_irr
         endif
!        Check these limits
         T_irr = max(T_irr, 0.d0)
         T_irr2 = max(T_irr2, 0.d0)
! need to reconstuct local tp(1,2) using ground hydrology code
!         T_irr_g = max(tp(1,2),0.d0)
         T_irr_g = max(tearth(i,j),0.d0)

!***     Set actual irrigation rates and update mwl and gml (if necessary)
         if (m_avail <= teeny) then

            if(flag_irrig_grndwat > 0) then
               irrig_water_act = irrig_water_pot(i,j)
               irrig_energy_act= irrig_water_act*shw*
     &                                 T_irr_g*rhow
               irrig_gw        = irrig_water_act
               irrig_gw_energy = irrig_energy_act
#ifdef TRACERS_WATER
               irrig_tracer_act(:) =
     &              irrig_water_act*tr_conc_in_aquifer(:)
               irrig_gw_tracer(:) = irrig_gw*tr_conc_in_aquifer(:)
#endif
            else
               irrig_water_act = 0.d0
               irrig_energy_act= 0.d0
               irrig_gw        = 0.d0
               irrig_gw_energy = 0.d0
#ifdef TRACERS_WATER
               irrig_tracer_act = 0.d0
               irrig_gw_tracer  = 0.d0
#endif

            endif

         elseif (m_avail >= m_irr_pot) then

            mwl_to_irrig = irrig_water_pot(i,j)*m_to_kg*dtsrc
            if (flake.gt.0 .and. mwl_to_irrig .gt. mldlk*m_to_kg*flake) 
     *           then           ! need layer 2 water
              gml_to_irrig = mldlk*m_to_kg*flake*shw*T_irr + 
     *             (mwl_to_irrig-mldlk*m_to_kg*flake)*shw*T_irr2
#ifdef TRACERS_WATER
              trml_to_irrig(:,1)=trml(:,1)
              trml_to_irrig(:,2)=trml(:,2)/(mwl-mldlk*m_to_kg*flake)
     *             *(mwl_to_irrig-mldlk*m_to_kg*flake)
#endif
            else
              gml_to_irrig = mwl_to_irrig*shw*T_irr
#ifdef TRACERS_WATER
              if (flake.gt.0) then
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/
     *                (mldlk*m_to_kg*flake)
              else
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/mwl
              endif                 
              trml_to_irrig(:,2)=0.
#endif
            end if
            irrig_water_act = irrig_water_pot(i,j)
            irrig_energy_act= gml_to_irrig*rhow / (m_to_kg*dtsrc) 
#ifdef TRACERS_WATER
            irrig_tracer_act = (trml_to_irrig(:,1)+trml_to_irrig(:,2))/
     *           (m_to_kg*dtsrc)
#endif
            
         else !!! (m_avail < m_irr_pot)

            mwl_to_irrig = m_avail
            if (flake.gt.0 .and. mwl_to_irrig .gt. mldlk*m_to_kg*flake) 
     *           then           ! need layer 2 water 
              gml_to_irrig = mldlk*m_to_kg*shw*T_irr*flake + 
     *             (mwl_to_irrig-mldlk*m_to_kg*flake)*shw*T_irr2
#ifdef TRACERS_WATER
              trml_to_irrig(:,1)=trml(:,1)
              trml_to_irrig(:,2)=trml(:,2)/(mwl-mldlk*m_to_kg*flake)
     *             *(mwl_to_irrig-mldlk*m_to_kg*flake)
#endif
            else
              gml_to_irrig = m_avail*shw*T_irr
#ifdef TRACERS_WATER
              if (flake.gt.0) then
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/
     *                (mldlk*m_to_kg*flake)
              else
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/mwl
              end if
              trml_to_irrig(:,2)=0.
#endif
            end if


            if (flag_irrig_grndwat > 0) then
               irrig_water_act  = irrig_water_pot(i,j)
               irrig_energy_act = irrig_water_act*shw*T_irr*rhow
               irrig_gw = irrig_water_act - m_avail / (m_to_kg*dtsrc)
               irrig_gw_energy  = irrig_energy_act - gml_to_irrig*rhow /
     *              (m_to_kg*dtsrc)  
#ifdef TRACERS_WATER
               irrig_gw_tracer(:) = irrig_gw*tr_conc_in_aquifer(:)
               irrig_tracer_act(:) =  (trml_to_irrig(:,1)
     *              + trml_to_irrig(:,2)) / (m_to_kg*dtsrc)
     &              + irrig_gw_tracer(:)
#endif
             else
               irrig_water_act  = m_avail / (m_to_kg*dtsrc) 
               irrig_energy_act = gml_to_irrig*rhow / (m_to_kg*dtsrc) 
#ifdef TRACERS_WATER
               irrig_tracer_act=(trml_to_irrig(:,1)+trml_to_irrig(:,2))/
     *              (m_to_kg*dtsrc)
#endif
             endif

         endif

      endif

!!!!! Test (HACK... revisit)
      if( (irrig_water_act-irrig_water_pot(i,j)) > 1.d-30 ) then
         write(6,*) 'error in irrigation: act> pot'
         write(6,*) 'i,j,irrig_water_act, irrig_water_pot(i,j):',
     &              i,j,irrig_water_act, irrig_water_pot(i,j)
         write(6,*) 'm_avail, m_irr_pot:', m_avail, m_irr_pot
         irrig_water_act=irrig_water_pot(i,j)
      endif
!!!!!
      end subroutine irrigate_extract

      end module irrigmod

