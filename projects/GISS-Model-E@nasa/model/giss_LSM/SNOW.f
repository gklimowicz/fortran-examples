!@sum contains code related to 3 layers snow model
!@auth I.Aleinov

#define DEBUG_SNOW

      MODULE SNOW_MODEL
!@sum SNOW_MODEL does column physics for 3 layers snow model
!@auth I.Aleinov
      USE FILEMANAGER, only: openunit
      USE CONSTANT, only :
     $     rho_water => rhow,   ! (kg m-3)
     $     rho_ice => rhoi,     ! (kg m-3)
     $     lhm_kg => lhm,       ! (J kg-1)
     $     lhe_kg => lhe,       ! (J kg-1)
     $     shi_kg => shi,       ! (J kg-1 C-1)
     $     tfrz => tf,          ! (K)
     $     sigma => stbo        ! (W m-2 K-4)
      USE TRIDIAG_MOD, only :  TRIDIAG

      IMPLICIT NONE
      SAVE
      PRIVATE

      PUBLIC snow_adv, snow_redistr, snow_fraction, i_earth, j_earth
      public MIN_SNOW_THICKNESS, MIN_FRACT_COVER, rho_water,
     &     rho_fresh_snow, minSnowTemperature

ccc physical parameters
!@var rho_fresh_snow density of fresh snow (kg m-3)
      real*8, parameter :: rho_fresh_snow =  150.d0 ! (kg m-3)
!@var lat_fusion volumetric latent heat of fusion for water  (J m-3)
      real*8, parameter :: lat_fusion = lhm_kg * rho_water ! (J m-3)
!@var lat_evap volumetric latent heat of evaporation for water  (J m-3)
      real*8, parameter :: lat_evap = lhe_kg * rho_water ! (J m-3)
!@var max_fract_water max fraction of free water in a wet snow (1)
      real*8, parameter :: max_fract_water = .055d0

ccc model parameters
!@var EPS small number (used as cutoff parameter) (1)
      real*8, parameter :: EPS = 1.d-8  ! was 1.d-12
!@var MAX_NL maximal number of snow layers (only 3 are used, really) (1)
      integer, parameter :: MAX_NL = 16
!@var TOTAL_NL maximal number of snow layers (actual number used) (1)
ccc actually MAX_NL can be set = TOTAL_NL but should test first ...
      integer, parameter :: TOTAL_NL=3
!@var MIN_SNOW_THICKNESS minimal thickness of snow (m)
ccc trying to increase MIN_SNOW_THICKNESS for stability
ccc      real*8, parameter :: MIN_SNOW_THICKNESS =  0.01d0  ! was 0.09d0
      real*8, parameter :: MIN_SNOW_THICKNESS =  0.1d0
!@var MIN_FRACT_COVER minimal snow cover (cutoff param) (1)
      real*8, parameter :: MIN_FRACT_COVER = 0.0001d0

!@var DEB_CH channel for debug output
      integer :: DEB_CH = 0
!@var i_earth, j_earth coordinate of current point (for debugging)
      integer i_earth, j_earth

!@var minSnowTemperature minimum allowed snow temperature (C)
!@+   (in GCM runs is reset in GHY_DRV to minGroundTemperature, 
!@+    which is a rundeck parameter)
      real*8 :: minSnowTemperature = -120.d0

      CONTAINS

      subroutine pass_water( wsn, hsn, dz, nl,
     &                       water_down, heat_down,
     &                       lat_fusion, max_fract_water,
     &                       rho_water, rho_fresh_snow)
!@sum  pass extra water to the lower layers
!@auth I.Aleinov
      implicit none
      integer nl
      real*8  wsn(nl), hsn(nl), dz(nl)
      real*8 water_down, heat_down, lat_fusion, max_fract_water
      real*8 rho_water, rho_fresh_snow
      integer n
      real*8 ice, free_water, ice_old
      integer pass_water_flag

      pass_water_flag = 0
      do n=1,nl
        ice_old = min(wsn(n),-hsn(n)/lat_fusion)
        ice = ice_old
        wsn(n) = wsn(n) + water_down
        hsn(n) = hsn(n) + heat_down
        water_down = 0.d0
        heat_down = 0.d0
        if( hsn(n) .ge. 0.d0 .or. wsn(n) .le. 0.d0 ) then ! all melted
          water_down = wsn(n)
          heat_down = hsn(n)
          wsn(n) = 0.d0
          hsn(n) = 0.d0
          dz(n) = 0.d0
          pass_water_flag = 1
        else if (hsn(n).gt.-wsn(n)*lat_fusion) then !/* some melted */
          ice = -hsn(n)/lat_fusion
          free_water = wsn(n) - ice
c!!!                may be ice*max_fract_water ??
          water_down = max(0d0, free_water - ice*max_fract_water)
          wsn(n) = wsn(n) - water_down
          !/* should I reduce dz here ??? */
          !/* dz(n) = dz(n) * min( ice/ice_old, 1.d0 )  */
          dz(n) = min( dz(n), ice*rho_water/rho_fresh_snow )
          pass_water_flag = 2
        else !/* all frozen */
          if (wsn(n)+EPS .lt. ice_old) dz(n) = dz(n)*wsn(n)/ice_old
          dz(n) = min( dz(n), wsn(n)*rho_water/rho_fresh_snow )
        endif
        dz(n) = max( dz(n), wsn(n)*rho_water/rho_ice )
      enddo
      return
      end subroutine pass_water


      subroutine snow_fraction(
     &     dz, nl, prsnow, dt, fract_cover, fract_cover_new)
!@sum computes new snow fraction and returns it in fract_cover_new
!@auth I.Aleinov
      implicit none
      integer nl
      !real*8 dz(TOTAL_NL+1), prsnow, dt
      real*8 dz(:), prsnow, dt
      real*8 fract_cover, fract_cover_new
      real*8 dz_aver, fresh_snow

      fresh_snow = rho_water/rho_fresh_snow * prsnow*dt

      dz_aver = sum(dz(1:nl))*fract_cover + fresh_snow

      !!fract_cover_new = min( 1.d0, dz_aver/MIN_SNOW_THICKNESS )
      !! changing max snow fraction to 95%
      fract_cover_new = min( .95d0, dz_aver/MIN_SNOW_THICKNESS )
      if ( fract_cover_new < MIN_FRACT_COVER ) fract_cover_new = 0.d0

      return
      end subroutine snow_fraction


      subroutine snow_redistr(
     &     dz, wsn, hsn, nl, fract_cover_ratio, tr_flux, dt )
!@sum  redistributes snow between the layers
!@auth I.Aleinov
      implicit none
      integer nl
      !real*8 dz(TOTAL_NL+1), wsn(TOTAL_NL), hsn(TOTAL_NL)
      real*8 dz(:), wsn(:), hsn(:)
      real*8 fract_cover_ratio
      integer nlo
      real*8 dzo(TOTAL_NL+1), wsno(TOTAL_NL), hsno(TOTAL_NL)
      integer n, no
      real*8 total_dz, ddz, delta, fract
      real*8, optional :: tr_flux(0:TOTAL_NL), dt
      integer i

cc      if( fract_cover_new .lt. MIN_FRACT_COVER ) then
cc        print *, 'snow_redistr: fract_cover_new < MIN_FRACT_COVER'
cc        print *, 'fract_cover_new = ', fract_cover_new
cc        call stop_model(
cc     &       'snow_redistr: fract_cover_new < MIN_FRACT_COVER',255)
cc      endif

      if ( dz(1) == 0.d0 ) return  ! no snow - nothing to redistribute

      dzo(1:nl)  = dz(1:nl)
      wsno(1:nl) = wsn(1:nl)
      hsno(1:nl) = hsn(1:nl)
      nlo = nl

      total_dz = sum( dzo(1:nl) )

      !fract_cover_ratio = fract_cover/fract_cover_new
      total_dz = total_dz*fract_cover_ratio

      if ( total_dz .gt. MIN_SNOW_THICKNESS*1.5d0 ) then
        nl = TOTAL_NL
        dz(1) = MIN_SNOW_THICKNESS
        dz(2:nl) = (total_dz-dz(1))/(nl-1)
      else
        nl = 1
        dz(1) = total_dz
      endif

      wsn(1:TOTAL_NL) = 0.d0
      hsn(1:TOTAL_NL) = 0.d0

      no = 0
      delta = 0.d0
      do n=1,nl
        do while ( delta.lt.dz(n) .and. no.lt.nlo )
          no = no + 1
          delta = delta + dzo(no)*fract_cover_ratio
          wsn(n) = wsn(n) + wsno(no)*fract_cover_ratio
          hsn(n) = hsn(n) + hsno(no)*fract_cover_ratio
          enddo
        ddz = delta - dz(n)
        fract = ddz/(dzo(no)*fract_cover_ratio)
ccc the following is just for check
          if ( fract.lt.-EPS .or. fract.gt.1.d0+EPS ) then
            print *, 'fract= ', fract
            call stop_model('snow_redistr: internal error 3',255)
            endif
        wsn(n) = wsn(n) - fract*wsno(no)*fract_cover_ratio
        hsn(n) = hsn(n) - fract*hsno(no)*fract_cover_ratio
        if ( n.lt.nl ) then
          wsn(n+1) = wsn(n+1) + fract*wsno(no)*fract_cover_ratio
          hsn(n+1) = hsn(n+1) + fract*hsno(no)*fract_cover_ratio
          delta = ddz
        else
          if ( abs(fract).gt.EPS ) then
            print *, 'fract= ', fract
            call stop_model('snow_redistr: internal error 1',255)
            endif
          endif
        enddo

ccc for tracers
      if ( present(tr_flux) ) then
        if ( nlo < TOTAL_NL ) wsno(nlo+1:TOTAL_NL) = 0.d0
        tr_flux(0) = 0.d0
        do i=1,TOTAL_NL
          tr_flux(i) =
     &         - (wsn(i) - wsno(i)*fract_cover_ratio)/dt
     &         + tr_flux(i-1)
        enddo

        if ( abs(tr_flux(TOTAL_NL)) > 1.d-12 )
     &       call stop_model("SNOW: snow redistr flux error",255)
      endif

      return
      end subroutine snow_redistr


      subroutine snow_adv(dz, wsn, hsn, nl,
     &    srht, trht, snht, htpr, evaporation, pr, dt,
     &    t_ground, dz_ground,
     &    water_to_ground, heat_to_ground,
     &    radiation_out, snsh_dt, evap_dt, evap_min,! fb_or_fv,
     &     tr_flux )
      implicit none
!@sum  a wrapper that calles real snow_adv (introduced for debugging)
!@auth I.Aleinov
ccc input:
      integer nl
      real*8 srht, trht, snht, htpr, evaporation
      real*8 pr, dt, t_ground, dz_ground
      real*8 snsh_dt, evap_dt, evap_min !, fb_or_fv

ccc output:
      real*8 water_to_ground, heat_to_ground
      real*8 radiation_out
ccc data arrays
      real*8 dz(TOTAL_NL+1), wsn(TOTAL_NL), hsn(TOTAL_NL)

ccc tracer variables
!@var tr_flux flux of water between snow layers (>0 is down) (m/s)
      real*8 tr_flux(0:TOTAL_NL)
      real*8 wsn_o(MAX_NL), evap_o
      integer nl_o

ccc for debug
      real*8 total_energy, total_water  !, lat_evap
      integer i, retcode
      integer, save :: counter=0

      counter = counter + 1
      !rint *, counter

ccc for tracers
      wsn_o(:) = 0.d0
      nl_o = nl
      wsn_o(1:nl) = wsn(1:nl)
      evap_o = evaporation

ccc checking if the model conserves energy (part 1) (for debugging)
      total_energy = 0.d0
      total_water = 0.d0
      do i=1,nl
        total_energy = total_energy - hsn(i)
        total_water = total_water - wsn(i)
      enddo

      call snow_adv_1(dz, wsn, hsn, nl,
     &    srht, trht, snht, htpr, evaporation, pr, dt,
     &    t_ground, dz_ground,
     &    water_to_ground, heat_to_ground,
     &    radiation_out, snsh_dt, evap_dt, evap_min, retcode )

!      if (fb_or_fv .le. 0.) return
ccc checking if the model conserves energy (part 2) (for debugging)
      do i=1,nl
        total_energy = total_energy + hsn(i)
        total_water = total_water + wsn(i)
      enddo
      total_energy = total_energy -
     &    (srht+trht-snht-lat_evap*evaporation+htpr)*dt
      total_energy = total_energy + heat_to_ground*dt + radiation_out*dt
      if ( abs(total_energy) .gt. 1.d0 ) then
        print*, "total energy error",i_earth, j_earth,total_energy,
     *       heat_to_ground*dt,radiation_out*dt
        call stop_model('snow_adv: total energy error',255)
      end if
      total_water = total_water
     &     - (pr - evaporation - water_to_ground)*dt
      if ( abs(total_water)/dt .gt. 1.d-15 ) then
c        print*, "water cons error",i_earth, j_earth, total_water/dt
        if ( abs(total_water)/dt .gt. 1.d-10 )
     &    call stop_model('snow_adv: water conservation error',255)
      end if

cc    if( fr_type .lt. 1.d-6 .and. abs(total_energy) .gt. 1.d-6 ) then
cc      print*, "total energy error",i_earth, j_earth,total_energy
cc      call abort
cc    endif

ccc for tracers
      !!!tr_flux(0) = tr_flux(0) + pr - evaporation
      tr_flux(0) = pr - evaporation
      do i=1,TOTAL_NL
        tr_flux(i) = -(wsn(i) - wsn_o(i))/dt
     &       + tr_flux(i-1)
      enddo
ccc checking if preserve water
      if ( abs( tr_flux(TOTAL_NL) - water_to_ground ) > 1.d-15 ) then
        if ( DEB_CH == 0 )
     $       call openunit("snow_debug", DEB_CH, .false., .false.)
        write(DEB_CH,*) "snow_adv: H2O error "
     &       , abs( tr_flux(TOTAL_NL) - water_to_ground )
     &       , tr_flux(TOTAL_NL) , water_to_ground
     &       , nl, nl_o, counter, evaporation, evap_o
        !call abort
        !stop 'snow_adv: H2O error'
      endif

      return
      end subroutine snow_adv


      subroutine snow_adv_1(dz, wsn, hsn, nl,
     &    srht, trht, snht, htpr, evaporation, pr, dt,
     &    t_ground, dz_ground,
     &    water_to_ground, heat_to_ground,
     &    radiation_out, snsh_dt, evap_dt, evap_min,
     &     retcode )
      implicit none
!@sum main program that does column snow physics
!@auth I.Aleinov
ccc input:
      integer nl
      real*8 srht, trht, snht, htpr, evaporation
      real*8 pr, dt, t_ground, dz_ground
      real*8 snsh_dt, evap_dt, evap_min

ccc constants: (now defined as global params)
      real*8 k_ground, c_ground
ccc output:
      real*8 water_to_ground, heat_to_ground
      real*8 radiation_out
      integer retcode

ccc main parameters: layer thickness, water equivalent, heat content
      real*8 dz(TOTAL_NL+1), wsn(TOTAL_NL), hsn(TOTAL_NL)
ccc!!! I wonder if the following arrays should have dim (nl) instead
ccc    of (MAX_NL) to force the allocation on a stack (for OpenMP) ?
      real*8 tsn(MAX_NL+1), csn(MAX_NL), ksn(MAX_NL+1), isn(MAX_NL)
      real*8 wsn_o(MAX_NL), hsn_o(MAX_NL), dz_o(MAX_NL+1)

      real*8 water_down, heat_down, fresh_snow, rho_snow, flux_in
      real*8 mass_layer, mass_above, scale_rho
      real*8 flux_in_deriv, flux_corr
      real*8 delta_tsn_impl ! shft of temperature due to implicit method
      real*8 delta_evap, evap_corr
      integer n, nl_o
c      real*8 dz_he(MAX_NL+1)

ccc!!!  check if lat_evap shoud be replaced by heat of sublimation

ccc !!! ground properties should be passed as formal parameters !!!
      k_ground =        3.4d0    !/* W K-1 m */    /* --??? */
      c_ground =        1.d5     !/* J m-3 K-1 */  /* -- ??? */

ccc will need initial values if all snow melted and for tracers
      wsn_o(1:nl) = wsn(1:nl)
      hsn_o(1:nl) = hsn(1:nl)
      dz_o(1:nl) = dz(1:nl)
      nl_o = nl

      !!! just to deceive the optimizer
      if ( dz_o(1) > 1.d6 ) call stop_model("snow_adv_1: what?!!",255)

ccc just in case, may need reasonable tsn(1) if all melted
      tsn(1) = 0.d0

ccc fluxes in/out
      heat_to_ground = 0.d0
      water_to_ground = 0.d0

#ifdef DEBUG_SNOW
      call check_rho_snow( dz, wsn, nl )
#endif


ccc compute amount of fresh snow
ccc !!! insert evaporation into computation of the amount of fresh snow
ccc use fresh_snow only to change the depth of the first layer
      fresh_snow = rho_water/rho_fresh_snow
     &                * min(pr*dt-evaporation*dt, -htpr*dt/lat_fusion)
      if(fresh_snow > 0.d0) then
        dz(1) = dz(1) + fresh_snow
        nl = max(nl,1) ! make sure that we create a layer if first snowfall
      else
        if ( wsn(1) < EPS ) goto 1000
        !goto 1000
      endif

ccc water and heat to pass through the snow pack
      water_down = (pr - evaporation)*dt
      heat_down = htpr*dt
ccc will include heat of evaporation later
c!!      heat_down = heat_down - lat_evap*evaporation*dt

!#ifdef DEBUG_SNOW
!      call check_rho_snow( dz, wsn, nl )
!#endif

      call pass_water( wsn, hsn, dz, nl, water_down, heat_down,
     &                       lat_fusion, max_fract_water,
     &                       rho_water, rho_fresh_snow)
      heat_to_ground = heat_to_ground + heat_down/dt
      water_to_ground = water_to_ground + water_down/dt

      if ( sum( wsn(1:nl) ) < EPS ) goto 1000

#ifdef DEBUG_SNOW
      call check_rho_snow( dz, wsn, nl )
#endif

      call snow_redistr(dz, wsn, hsn, nl, 1.d0)
      dz(nl+1) = dz_ground

ccc compute spec. heat and thermal conductivity
      do n=1,nl
        rho_snow = wsn(n)*rho_water/dz(n)
        !csn(n) = (1.9d6) * rho_snow / rho_ice
        csn(n) = shi_kg * rho_snow
        ksn(n) = (3.22d-6) * rho_snow**2
      enddo
      ksn(nl+1) = k_ground

ccc compute temperature of the layers (and amount of ice)
      do n=1,nl
        if ( hsn(n) .gt. 0.d0 ) then
          print *, 'OOPS, No snow in layer ', n
          call stop_model('snow_adv_1: empty snow layer found (1)',255)
        else if ( hsn(n) .gt. -wsn(n)*lat_fusion ) then
          tsn(n) = 0.d0
          isn(n) = -hsn(n)/lat_fusion
        else
          tsn(n) = (hsn(n)+wsn(n)*lat_fusion)/(csn(n)*dz(n))
          isn(n) = wsn(n)
        endif
      enddo
      tsn(nl+1) = t_ground

c!!! this is for debugging
      if(tsn(1).lt.minSnowTemperature)
     &     call stop_model("SNOW:tsn<minSnowTemperature",255)

ccc compute incomming heat flux (from atm.)
ccc include all fluxes except htpr (which is already included)
      flux_in =
     &      srht
     &      +trht - sigma*(tsn(1)+tfrz)**4
     &      -lat_evap*evaporation
     &      -snht

      flux_in_deriv =
     &     -4.d0*sigma*(tsn(1)+tfrz)**3 - lat_evap*evap_dt - snsh_dt

      radiation_out = sigma*(tsn(1)+tfrz)**4

ccc this is a hack to keep heat_eq stable: force dz >= 0.1
      !dz_he(1:nl+1) = dz(1:nl+1)
      !dz_he(1) = max( dz_he(1), MIN_SNOW_THICKNESS )

ccc solve heat equation
      call heat_eq(
     &     dz, tsn, hsn, csn, ksn, nl,
     &     flux_in, flux_in_deriv, flux_corr, dt )

c!!! this is for debugging
      if(tsn(1).lt.minSnowTemperature)
     &     call stop_model("SNOW:he:tsn<minSnowTemperature",255)

      heat_to_ground = heat_to_ground + flux_in

ccc now redistribute extra flux proportionallly to input derivatives
ccc snht_out_cor, delta_tsn_impl
      delta_tsn_impl = flux_corr / flux_in_deriv
      radiation_out = radiation_out -
     &     (-4.d0*sigma*(tsn(1)+tfrz)**3)*delta_tsn_impl
      snht = snht + snsh_dt * delta_tsn_impl
      !evaporation = evaporation + evap_dt * delta_tsn_impl
      delta_evap = evap_dt * delta_tsn_impl
      if ( evaporation + delta_evap < evap_min ) then !hack to restrict the dew
        evap_corr = evap_min - (evaporation + delta_evap)
        delta_evap = delta_evap + evap_corr
        snht = snht - evap_corr*lat_evap
        !write(0,*) "SHOW: evap_corr = ", evap_corr
      endif
      evaporation = evaporation + delta_evap

ccc and now remove (add) water due to extra evaporation.
c!!! this may make wsn(1) negative, the only way I see now to prevent it
c!!! is to keep minimal thickness of snow big enough

      water_down = - delta_evap * dt ! ??

#ifdef DEBUG_SNOW
      call check_rho_snow( dz, wsn, nl )
#endif

ccc pass extra water down
ccc      water_down = 0.d0
      heat_down = 0.d0
      call pass_water( wsn, hsn, dz, nl, water_down, heat_down,
     &                       lat_fusion, max_fract_water,
     &                       rho_water, rho_fresh_snow)
      heat_to_ground = heat_to_ground + heat_down/dt
      water_to_ground = water_to_ground + water_down/dt

      if ( sum( wsn(1:nl) ) < EPS ) goto 1000

      call snow_redistr(dz, wsn, hsn, nl, 1.d0)
      dz(nl+1) = dz_ground

ccc compute temperature of the layers
      do n=1,nl
        if ( hsn(n) .gt. 0.d0 ) then
          print *, 'OOPS, No snow in layer ', n
          call stop_model('snow_adv_1: empty snow layer found (2)',255)
        else if ( hsn(n) .gt. -wsn(n)*lat_fusion ) then
          tsn(n) = 0.d0
        else
          tsn(n) = (hsn(n)+wsn(n)*lat_fusion)/(csn(n)*dz(n))
        endif
      enddo
      tsn(nl+1) = t_ground

ccc repack the layers
      mass_above = 0.d0
      do n=1,nl
        if( dz(n) .gt. EPS ) then
          mass_layer = wsn(n)*rho_water
          mass_above = mass_above + .5d0 * mass_layer
          scale_rho = .5d-7 * 9.8d0 * mass_above
     &     * exp(14.643d0 - 4000.d0/(tsn(n)+tfrz)
     &                  - .02d0 * mass_layer/dz(n))
     &     * dt
          scale_rho = 1.d0 + scale_rho
          dz(n) = dz(n) / scale_rho
          if( dz(n).lt.mass_layer/rho_ice )
     &                   dz(n) = mass_layer/rho_ice
          mass_above = mass_above + .5d0 * mass_layer
        endif
      enddo

#ifdef DEBUG_SNOW
      call check_rho_snow( dz, wsn, nl )
#endif

c!!! this is for debugging
      if(tsn(1).lt.minSnowTemperature) then
        print*,"tsn error",i_earth, j_earth,1,tsn(1:nl)
        call stop_model('snow_adv_1: tsn error',255)
      end if

      if ( radiation_out > 500.d0 ) call stop_model("SNOW:T>0",255)

      retcode = 0
      return

 1000 continue ! all melted
ccc!!! may create water_to_ground<0 if evaporation is not properly
ccc    limited. Check later.
      water_to_ground = sum( wsn_o(1:nl_o) )/dt + pr - evaporation
      heat_to_ground  = sum( hsn_o(1:nl_o) )/dt + htpr
     &     - lat_evap*evaporation
     &     - snht + srht + trht - sigma*(tsn(1)+tfrz)**4
      radiation_out = sigma*(tsn(1)+tfrz)**4
      if ( radiation_out > 500.d0 ) call stop_model("SNOW:T>0",255)
      wsn(1:nl) = 0.d0
      hsn(1:nl) = 0.d0
      dz(1:nl) = 0.d0
      nl = 1
      return
      end subroutine snow_adv_1


      subroutine heat_eq(
     &     dz, tsn, hsn, csn, ksn, nl,
     &     flux_in, flux_in_deriv, flux_corr, dt)
      implicit none
!@sum solves heat transport equation for snow layers
!@auth I.Aleinov
!@alg This solver is currently using a half-implicit scheme.
!@+   Implicitness is controlled by the parameters:
!@+       gamma - for upper boundary
!@+       eta(MAX_NL+1) - for the rest of the boundaries
!@+   In general we are trying to use half-implicit method (gamma = .5)
!@+   at the surface. But this may introduce a systematic error when
!@+   temperature of the snow is close to 0 C.
!@+   When option DO_EXPLIC_0 is enabled the program checks for possible
!@+   error and if necessary recomputes the solution with reduced gamma,
!@+   i.e. in more explicit way.
!@+   DO_EXPLIC_0 is enabled by default. One may try to turn it off if
!@+   there are serious stability problems.
!@calls :TRIDIAG
      integer nl
      real*8 dz(nl+1), tsn(nl+1), hsn(nl)
      real*8 csn(nl), ksn(nl+1)
ccc n+1 corresponds to the first layer of the ground
      real*8 flux_in, flux_in_deriv, flux_corr, dt, syst_flux_err

      real*8 eta(MAX_NL+1), tnew(MAX_NL), gamma
ccc arrays for coefficients:
      real*8 a(MAX_NL), b(MAX_NL), c(MAX_NL), f(MAX_NL)
      real*8 left, right, dt_to_cdz
      integer n, iter
      integer :: itermax=2

ccc setting implicit/explicit choice for each layer
ccc eta=1 - implicit, eta=0 - explicit
ccc for now using explicit method for layers with T=0.d0

      do n=1,nl
         eta(n) = .5d0
ccc        eta(n) = 1.d0
c!!        if( tsn(n) .ge. 0.d0 ) eta(n) = 0.d0
        enddo
      eta(nl+1) = 0.d0
ccc      gamma = .5d0
c!! I am trying to make it nearly expicit for a test !!!
      gamma = .5d0

!!! testing: trying to improve stability
ccc switching to explicit method for very thin snow
      if ( dz(1) < MIN_SNOW_THICKNESS * .5d0 ) then
        eta(1) = 1.d0
        gamma = 1.d0
      endif

ccc In general we are trying to use half-implicit method (gamma = .5d0)
ccc on the interface. But this introduces a systematic error when
ccc tnew(1) > 0. So in the case DO_EXPLIC_0 we move to explicit method
ccc (i.e. reduce gamma)

ccc DO_EXPLIC_0 enables correction of the error which is introduced
ccc by implicit scheme when tsn(1) == 0 C
ccc (of course one has to use preprocessing to do it ...)
ccc it is enabled by default

#define DO_EXPLIC_0
#ifdef DO_EXPLIC_0
      do iter=1,itermax
#endif
ccc equation for upper snow layer
        n = 1
        dt_to_cdz = dt/(csn(n)*dz(n))
        right = 2.d0*dt_to_cdz/( dz(n)/ksn(n) + dz(n+1)/ksn(n+1) )
        a(n) = 1.d0 + right*eta(n) - dt_to_cdz*flux_in_deriv*gamma
        b(n) = 0.d0
        c(n) = -right*eta(n+1)
        f(n) = tsn(n)*( 1.d0 - right*(1.d0-eta(n))
     &                  - dt_to_cdz*flux_in_deriv*gamma )
     &         + tsn(n+1)*right*(1.d0-eta(n+1))
     &         + dt_to_cdz * flux_in

ccc all other snow layers including the lower one
        do n=2,nl

          dt_to_cdz = dt/(csn(n)*dz(n))
          right = 2.d0*dt_to_cdz/( dz(n)/ksn(n) + dz(n+1)/ksn(n+1) )
          left  = 2.d0*dt_to_cdz/( dz(n)/ksn(n) + dz(n-1)/ksn(n-1) )

          a(n) = 1.d0 + ( left + right )*eta(n)
          b(n) = -left*eta(n-1)
          c(n) = -right*eta(n+1)

          f(n) = tsn(n)*( 1.d0 - ( left + right )*(1.d0-eta(n)) )
     &         + tsn(n-1)*left*(1.d0-eta(n-1))
     &         + tsn(n+1)*right*(1.d0-eta(n+1))
          enddo

c        call sweep3diag(b, c, a, f, tnew, nl)
          call TRIDIAG(b,a,c,f,tnew,nl)

ccc flux_corr is the energy wich should be returned to the atmosphere
        flux_corr = flux_in_deriv*( tnew(1) - tsn(1) )*gamma

#ifdef DO_EXPLIC_0
        syst_flux_err = flux_in_deriv*( tnew(1) - 0.d0 )*gamma
        ! if( iter/=2 .and. tnew(1)>0.d0 .and. flux_in_deriv<0.d0 ) then
        ! back to 77 :-L
        if ( iter.ne.itermax .and.
     &         tnew(1).gt.0.d0 .and. flux_in_deriv.lt.0.d0 ) then
          syst_flux_err = flux_in_deriv*( tnew(1) - 0.d0 )*gamma
          gamma = (1.d0 - syst_flux_err/flux_corr) *gamma
        else
          ! exit
          ! back to 77 :-L
          goto 77
        endif
      enddo
 77   continue
#endif

      do n=1,nl
        hsn(n) = hsn(n) + (tnew(n)-tsn(n))*csn(n)*dz(n)
        enddo

ccc flux to the ground :
      n = nl
      flux_in = - ( tsn(n+1)-tnew(n)*eta(n)-tsn(n)*(1.d0-eta(n)) ) *
     &           2.d0/( dz(n)/ksn(n) + dz(n+1)/ksn(n+1) )

      return
      end subroutine heat_eq

      subroutine check_rho_snow( dz, wsn, nl )
      integer nl
      real*8 dz(nl+1), wsn(nl)
      integer n
      do n=1,nl
        if( wsn(n) < EPS ) cycle
        if ( wsn(n)/dz(n)*rho_water+EPS < rho_fresh_snow) then
          print *,"rho_snow < rho_fresh_snow at i,j=",i_earth, j_earth
          print *,"n, wsn, rho_sn, rho_fresh= ", n,wsn(n),
     *         wsn(n)/dz(n)*rho_water ,rho_fresh_snow
          call stop_model('check_rho_snow: rho < rho_fresh_snow',255)
        end if
      enddo

      do n=1,nl
        if( wsn(n) < EPS ) cycle
        if ( wsn(n)/dz(n)*rho_water-EPS > rho_ice) then
          print *,"rho_snow > rho_ice at i,j=",i_earth, j_earth
          print *,"n, wsn, rho_sn, rho_ice= ", n,wsn(n),
     *         wsn(n)/dz(n)*rho_water ,rho_ice
          call stop_model('check_rho_snow: rho > rho_ice',255)
        end if
      enddo


      end subroutine check_rho_snow

      END MODULE SNOW_MODEL

