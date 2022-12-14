#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  2010/10/13
!@cont init_LI,PRECIP_LI,GROUND_LI,daily_LI

      SUBROUTINE init_LI(istart)
!@sum  init_ice initialises landice arrays
!@auth Original Development Team
      USE CONSTANT, only : lhm,tf
      use TimeConstants_mod, only: SECONDS_PER_DAY,EARTH_DAYS_PER_YEAR
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only : dtsrc
      USE FLUXES, only : flice,focean
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      USE GEOM, only : axyp,imaxj,lat2d
      USE LANDICE, only: ace1li,ace2li,glmelt_on,glmelt_fac_nh
     *     ,glmelt_fac_sh,fwarea_sh,fwarea_nh,accpda,accpdg,eaccpda
     *     ,eaccpdg,snmin,micbimp,eicbimp
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,traccpda,traccpdg
#endif
#endif               /* TNL: inserted */
      Use LANDICE_COM, Only: NHC,FHC,TLANDI,SNOWLI, MDWNIMP,EDWNIMP,
     *                       FSHGLM,FNHGLM
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trdwnimp
#endif
      USE FLUXES, only : atmgla,atmglas,atmocn
#ifdef GLINT2
      use fluxes, only : atmglas_hp, flice
      use hp2hc
      use landice_com, only : fhp_approx
#endif
#ifdef TRACERS_OCEAN
#ifdef TRACERS_OCEAN_INDEP
#else
      use OldTracer_mod, only: trglac, trw0
#endif
#endif
      USE DIAG_COM, only : npts,icon_MLI,icon_HLI,title_con,conpt0
     *     ,icon_MICB,icon_HICB
      USE Dictionary_mod
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds, 
     &     GLOBALSUM,AM_I_ROOT
      USE EXCHANGE_TYPES, only : avg_patches_srfstate_exports
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER, INTENT(IN) :: istart
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     *                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     r8mask ! dummy array to read in 0/1 mask
      Logical*4,Dimension(GRID%I_STRT_HALO:grid%I_STOP_HALO,
     *                    GRID%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     LOC_GLM  !  mask for GLMELT around Antarctica and Greenland
      LOGICAL :: do_glmelt = .false.
      INTEGER I,J,N, IHC, fid
      INTEGER :: I_0,I_1, J_0,J_1
      Real*8  :: DTSRCzY,byoarea

! Change the variable we'll initialize depending on GLINT2 existence
#ifdef GLINT2
#define ATMGLAX atmglas_hp
#else
#define ATMGLAX atmglas
#endif

      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      ! Make sure fhc is set, even when running with EC's
      atmglas(1)%fhc(:,:) = fhc(:,:,1)

      ! temporary
      do ihc=1+lbound(atmglas,1),ubound(atmglas,1)
        do j=j_0,j_1
        do i=i_0,i_1
#ifdef GLINT2
          ! Set height point "area" equal to amount that total value
          ! will increase when we increase this height point value by
          ! 1.  This was obtained in glint2_modele_init_landice_com()
          ! by summing over the hp_to_hc matrix, combined with fhc
          ! The "area" here is used by certain diagnostics.
          atmglas_hp(ihc)%fhc(i,j) = fhp_approx(i,j,ihc)
          atmglas_hp(ihc)%ftype(i,j) =
     &         flice(i,j)*atmglas_hp(ihc)%fhc(i,j)
#endif
          atmglas(ihc)%fhc(i,j) = fhc(i,j,ihc)
          atmglas(ihc)%ftype(i,j) =
     &         flice(i,j)*atmglas(ihc)%fhc(i,j)
        enddo
        enddo
      enddo

C**** set GTEMP array for landice
      DO IHC=1,NHC
      DO J=J_0,J_1
        DO I=I_0,I_1
          IF (FLICE(I,J).gt.0) THEN
            ATMGLAX(ihc)%SNOW(I,J)=SNOWLI(I,J,IHC)
            ATMGLAX(ihc)%GTEMP(I,J)=TLANDI(1,I,J,IHC)
            ATMGLAX(ihc)%GTEMP2(I,J)=TLANDI(2,I,J,IHC)
            ATMGLAX(ihc)%GTEMPR(I,J)=TLANDI(1,I,J,IHC)+TF

#ifdef SCM
            if( SCMopt%Tskin )then
              ATMGLAX(ihc)%GTEMP(I,J) = SCMin%Tskin - TF
              ATMGLAX(ihc)%GTEMP2(I,J) = SCMin%Tskin - TF
              ATMGLAX(ihc)%GTEMPR(I,J) = SCMin%Tskin
            endif
#endif
#ifdef TRACERS_WATER
!TODO This should depend on whether itime>itime_tr0(n), not on istart
!TODO and may be different for each individual tracer (reto) 
            if (istart.ge.9) then ! ok if all tracers start at beg.of run
            IF (SNOWLI(I,J,IHC).gt.SNMIN) THEN
              ATMGLAX(ihc)%GTRACER(:,I,J)=
     &             TRSNOWLI(:,I,J,IHC)/SNOWLI(I,J,IHC)
            ELSE
              ATMGLAX(ihc)%GTRACER(:,I,J)=
     &             TRLNDI(:,I,J,IHC)/(ACE1LI+ACE2LI)
            END IF
            end if
#endif
          END IF
        END DO
      END DO
      END DO ! ihc
#undef ATMGLAX
#ifdef GLINT2
      call bundle_hp_to_hc(bundle_init_li,
     &     atmglas_hp(1:), atmglas(1:))
#endif
      call avg_patches_srfstate_exports(grid,
c     &     atmglas,atmgla, ! gfortran prob. if passed as class() args
     &     atmglas(:)%atmsrf_xchng_vars,atmgla%atmsrf_xchng_vars,
     &     rel=.true.)

C**** Calculate (fixed) iceberg melt terms from Antarctica and Greenland
      call sync_param("glmelt_on",glmelt_on)
      call sync_param("glmelt_fac_nh",glmelt_fac_nh)
      call sync_param("glmelt_fac_sh",glmelt_fac_sh)

C**** Read in GLMELT file to distribute glacial melt
      IF (JM.gt.24) THEN ! for finer than old 8x10

      do_glmelt = file_exists('GLMELT')

      if(do_glmelt) then
        fid = par_open(grid,'GLMELT','read')
        call read_dist_data(grid,fid,'mask',r8mask)
        call par_close(grid,fid)
        loc_glm(i_0:i_1,j_0:j_1) = r8mask(i_0:i_1,j_0:j_1).eq.1d0
      else
        loc_glm(:,:) = .false.
      endif
C**** Calculate hemispheric fractions (weighted by area and landmask)
C**** This could be extended by a different number in GLMELT to
C**** give ice sheet (rather than hemispheric) dependent areas
      FSHGLM(:,:) = 0
      FNHGLM(:,:) = 0
      DO J=J_0,J_1
        DO I=I_0,I_1
          IF (LOC_GLM(I,J) .and. FOCEAN(I,J).gt.0 ) THEN
            IF(LAT2D(I,J).LT.0.) THEN
              FSHGLM(I,J) = AXYP(I,J)*FOCEAN(I,J)
            ELSE
              FNHGLM(I,J) = AXYP(I,J)*FOCEAN(I,J)
            ENDIF
          END IF
        END DO
      END DO
      Call GLOBALSUM (GRID, FSHGLM, FWAREA_SH, ALL=.True.)
      Call GLOBALSUM (GRID, FNHGLM, FWAREA_NH, ALL=.True.)
      If (FWAREA_SH > 0)
     *  FSHGLM(I_0:I_1,J_0:J_1) = FSHGLM(I_0:I_1,J_0:J_1) / FWAREA_SH
      If (FWAREA_NH > 0)
     *  FNHGLM(I_0:I_1,J_0:J_1) = FNHGLM(I_0:I_1,J_0:J_1) / FWAREA_NH

      END IF

C**** Intialise gmelt fluxes
      atmocn%GMELT = 0. ; atmocn%EGMELT = 0.
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
      atmocn%TRGMELT = 0.
#endif
#endif /* TNL: inserted */

      if (do_glmelt .and. glmelt_on > 0) then
C****  Note that water goes in with as if it is ice at 0 deg.
C****  Possibly this should be a function of in-situ freezing temp?
C****
C**** gmelt_fac_nh/sh should be used to match accumulation in control
C**** run so that global mean sea level is nearly constant. We can't
C**** get perfect values, but the drift should hopefully be less
C**** than 1 mm/year. This will minimise shock once interactive values
C**** are calculated.
C****
C**** Initial average mass fluxes for Antarctica/Greenland

      if (istart.le.8 .and. ACCPDA.eq.0) then
C**** Initiallise total mass/energy fluxes (only at start of run)
C**** The net accumulation from IPCC2 report is 2016x10**12 kg/year
C**** for Antarctica and for Greenland it is 316x10**12 kg/year
        ACCPDA = glmelt_fac_sh*2016d12 ! kg/year
        ACCPDG = glmelt_fac_nh*316d12 ! kg/year
        EACCPDA = -LHM*ACCPDA ; EACCPDG = -LHM*ACCPDG ! J/year
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
        TRACCPDA(:) = trglac()*ACCPDA ! kg/year
        TRACCPDG(:) = trglac()*ACCPDG ! kg/year
#endif
#endif /* TNL: inserted */

C**** initiallise implicit accumulators (note that these arrays are
C**** not used until at least one full year has passed)
        MDWNIMP=0.  ;  EDWNIMP=0.
        MICBIMP=0.  ;  EICBIMP=0.
#ifdef TRACERS_WATER
        TRDWNIMP=0.
c        TRICBIMP=0.
#endif
      end if

C**** Set GMELT (kg per source time step)
C****   from total hemispere ACCPDA (kg per year)
      DTSRCzY = DTSRC / (SECONDS_PER_DAY*EARTH_DAYS_PER_YEAR)
#ifndef SCM
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          IF (LOC_GLM(I,J) .and. lat2D(I,J).lt.0 ) THEN ! SH
            atmocn% GMELT(I,J) =  ACCPDA*FSHGLM(I,J)*DTSRCzY  !  kg
            atmocn%EGMELT(I,J) = EACCPDA*FSHGLM(I,J)*DTSRCzY  !  J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            atmocn%TRGMELT(:,I,J) = TRACCPDA(:)*FSHGLM(I,J)*DTSRCzY  !  kg
#endif
#endif /* TNL: inserted */
          END IF
          IF (LOC_GLM(I,J) .and. lat2D(I,J).gt.0 ) THEN ! NH
            atmocn% GMELT(I,J) =  ACCPDG*FNHGLM(I,J)*DTSRCzY  !  kg
            atmocn%EGMELT(I,J) = EACCPDG*FNHGLM(I,J)*DTSRCzY  !  J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            atmocn%TRGMELT(:,I,J) = TRACCPDG(:)*FNHGLM(I,J)*DTSRCzY  !  kg
#endif
#endif /* TNL: inserted */
          END IF
          IF (LOC_GLM(I,J) .and. FOCEAN(I,J).GT.0.) THEN
            byoarea = 1.d0/(AXYP(I,J)*FOCEAN(I,J))
            atmocn% GMELT(I,J) = atmocn% GMELT(I,J)*byoarea
            atmocn%EGMELT(I,J) = atmocn%EGMELT(I,J)*byoarea
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            atmocn%TRGMELT(:,I,J) = atmocn%TRGMELT(:,I,J)*byoarea
#endif
#endif /* TNL: inserted */
          ENDIF
        END DO
      END DO
#endif

      end if

C**** Set conservation diagnostics for land ice mass, energy
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LAND ICE","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_MLI)
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LNDI ENR","(10**6 J/M^2)   ",
     *     "(10**-3 W/M^2)  ",1d-6,1d3,icon_HLI)

C**** Set conservation diagnostics for implicit iceberg mass, energy
      QCON=(/ F, F, F, F, F, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"ICEBERG ","(KG/M^2)        ",
     *     "(10**-7 KG/SM^2)",1d0,1d7,icon_MICB)
      QCON=(/ F, F, F, F, F, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"ICBRG EN","(10**6 J/M^2)   ",
     *     "(10**-1 W/M^2)  ",1d-6,1d1,icon_HICB)
C****
      END SUBROUTINE init_LI

      subroutine reset_glaacc
      use landice_com, only : ijhc
      implicit none
      ijhc = 0.
      return
      end subroutine reset_glaacc

      subroutine ijhc_defs
      use CONSTANT, only : bygrav,sha,rgas
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : dtsrc
      use LANDICE_COM, only : kijhc,ia_ijhc,sname_ijhc,lname_ijhc
     &     ,units_ijhc,denom_ijhc,scale_ijhc,cdl_ijhc
      use LANDICE_COM, only : nhc
      use LANDICE_COM, only :
     &     ijhc_frac,ijhc_fhc,ijhc_one,
     %     IJHC_SRFP,
     &     IJHC_PRECLI,  ! done
     &     IJHC_RUNLI,   ! done
     &     IJHC_EVAPLI,  ! done
     &     IJHC_F0LI,    ! done
     &     IJHC_TSLI,    ! done GROUND_LI
     &     IJHC_SHDTLI,   ! done
     &     IJHC_EVHDT,    ! done
     &     IJHC_TRHDT,    ! done
     &     IJHC_IMPMLI,		! done GROUND_LI
     &     IJHC_IMPHLI		! done GROUND_LI

      use DIAG_COM, only : ia_src,ia_srf,cdl_ij_template
      use MDIAG_COM, only : make_timeaxis
#ifdef CUBED_SPHERE
      use LANDICE_COM, only : cdl_ijhc_latlon
      use DIAG_COM, only : cdl_ij_latlon_template
#endif
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use cdl_mod
      USE FLUXES, only : nisurf

      implicit none
      integer :: k,kk
      character(len=32) :: dimstr,lldimstr
      logical :: set_miss
c
      do k=1,kijhc
         write(sname_ijhc(k),'(a4,i3.3)') 'IJHC',k
         lname_ijhc(k) = 'no output'
         units_ijhc(k) = 'unused'
         scale_ijhc(k) = 1.
         denom_ijhc(k) = 0
         ia_ijhc(k) = ia_src
      enddo

c
      k=0
c
      k=k+1				! ijhc
      ijhc_frac = k
      sname_ijhc(k) = 'frac'
      lname_ijhc(k) = 'area fraction'
      units_ijhc(k) = '1'
c
      k=k+1				! ijhc
      ijhc_fhc = k
      sname_ijhc(k) = 'fhc'
      lname_ijhc(k) = 'area fraction (ice only)'
      units_ijhc(k) = '1'
      denom_ijhc(k) = 0  ! not ijhc_frac
c
      k=k+1				! ijhc
      ijhc_one = k
      sname_ijhc(k) = 'one'
      lname_ijhc(k) = 'test diagnostic, should == 1'
      units_ijhc(k) = 'K'
      scale_ijhc(k) = 1d0/DTsrc ! to cancel acc factor of dtsurf
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1				! ijhc
      ijhc_srfp = k
      sname_ijhc(k) = 'srfp'
      lname_ijhc(k) = 'surface air pressure'
      units_ijhc(k) = 'hPa'
      scale_ijhc(k) = 1d0/DTsrc ! [s-1] to cancel acc factor of dtsurf
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_PRECLI = k ! PREC OVER LAND ICE [mm day-1]       1 CN
      lname_ijhc(k) = 'PRECIPITATION OVER LAND ICE (HC)'
      units_ijhc(k) = 'mm day-1'
      sname_ijhc(k) = 'pr_lndice'
      scale_ijhc(k) = SECONDS_PER_DAY/DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_RUNLI = k ! RUN1 OVER LAND ICE  [kg m-2] (NO PRT)    1 PG
      lname_ijhc(k) = 'SURFACE RUNOFF OVER LAND ICE (HC)'
      units_ijhc(k) = 'mm day-1'
      sname_ijhc(k) = 'runoff_lndice'
      scale_ijhc(k) = SECONDS_PER_DAY/DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_EVAPLI = k ! EVAP OVER LAND ICE  [kg m-2]          1 GD
      lname_ijhc(k) = 'LAND ICE EVAPORATION (HC)'
      units_ijhc(k) = 'mm day-1'
      sname_ijhc(k) = 'evap_lndice'
      scale_ijhc(k) = SECONDS_PER_DAY/DTsrc
c     iw built-in
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_F0LI = k ! F0DT, NET HEAT AT Z0 OVER LAND ICE [J m-2] 1 GD
      lname_ijhc(k) = 'NET HEAT INTO LAND ICE (HC)'
      units_ijhc(k) = 'W m-2'
      sname_ijhc(k) = 'netht_lndice'
      scale_ijhc(k) = 1./DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_TSLI = k ! SURF AIR TEMP OVER LAND ICE  [K]  NISURF*1 SF
      lname_ijhc(k) = 'SURF AIR TEMP OVER LAND ICE (HC)'
      units_ijhc(k) = 'K'
      sname_ijhc(k) = 'tsurf_lndice'
      scale_ijhc(k) = 1.d0/DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_SHDTLI = k ! SHDT OVER LAND ICE  [J m-2]           1 SF
      lname_ijhc(k) = 'SENS HEAT FLUX OVER LAND ICE (HC)'
      units_ijhc(k) = 'W m-2'
      sname_ijhc(k) = 'sensht_lndice'
      scale_ijhc(k) = 1./DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_EVHDT = k ! EVHDT OVER LAND ICE  [J m-2]           1 SF
      lname_ijhc(k) = 'LATENT HEAT FLUX OVER LAND ICE (HC)'
      units_ijhc(k) = 'W m-2'
      sname_ijhc(k) = 'latht_lndice'
      scale_ijhc(k) = 1./DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      IJHC_TRHDT = k ! TRHDT OVER LAND ICE  [J m-2]           1 SF
      lname_ijhc(k) = 'NET THERMAL RADIATION INTO LAND ICE (HC)'
      units_ijhc(k) = 'W m-2'
      sname_ijhc(k) = 'trht_lndice'
      scale_ijhc(k) = 1./DTsrc
      denom_ijhc(k) = ijhc_fhc
c
      k=k+1 ! ijhc
      IJHC_IMPMLI = k ! IMPLICIT MASS FLUX over LAND ICE [kg m-2]
      lname_ijhc(k) = 'IMPLICIT MASS FLUX over LAND ICE (HC)'
      units_ijhc(k) = 'kg m-2 s-1'
      sname_ijhc(k) = 'impm_lndice'
      scale_ijhc(k) = 1./DTsrc
      denom_ijhc(k) = ijhc_fhc

c
      k=k+1 ! ijhc
      IJHC_IMPHLI = k ! IMPLICIT HEAT FLUX over LAND ICE [J m-2]
      lname_ijhc(k) = 'IMPLICIT HEAT FLUX over LAND ICE (HC)'
      units_ijhc(k) = 'W m-2'
      sname_ijhc(k) = 'imph_lndice'
      scale_ijhc(k) = 1./DTsrc
      denom_ijhc(k) = ijhc_fhc

      if (k .gt. kijhc) then
        if(am_i_root())
     &       write (6,*) 'ijhc_defs: Increase kijhc=',kijhc,' to ',k
        call stop_model( 'kijhc too small', 255 )
      end if
      if(AM_I_ROOT())
     &     write (6,*) 'Number of IJHC diagnostics defined: kijhcmax=',k

c
c Declare the dimensions and metadata of IJHC output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      cdl_ijhc = cdl_ij_template
      call add_dim(cdl_ijhc,'nhc',nhc)

      lldimstr='(nhc,lat,lon) ;'
#ifdef CUBED_SPHERE
      cdl_ijhc_latlon = cdl_ij_latlon_template
      call add_dim(cdl_ijhc,'nhc',nhc)
      dimstr='(tile,nhc,y,x) ;'
#else
      dimstr=lldimstr
#endif
      do k=1,kijhc
        if(trim(units_ijhc(k)).eq.'unused') cycle
        set_miss = denom_ijhc(k).ne.0
        call add_var(cdl_ijhc,
     &       'float '//trim(sname_ijhc(k))//trim(dimstr),
     &       units=trim(units_ijhc(k)),
     &       long_name=trim(lname_ijhc(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#ifdef CUBED_SPHERE
        call add_var(cdl_ijhc_latlon,
     &       'float '//trim(sname_ijhc(k))//trim(lldimstr),
     &       units=trim(units_ijhc(k)),
     &       long_name=trim(lname_ijhc(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#endif
      enddo

      return
      end subroutine ijhc_defs

      SUBROUTINE PRECIP_LI(atmgla,ihc)
!@sum  PRECIP_LI driver for applying precipitation to land ice fraction
!@auth Original Development team
!@calls LANDICE:PRECLI
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      USE CONSTANT, only : tf
      USE GEOM, only : imaxj
      USE LANDICE, only: ace1li,ace2li,precli,snmin
      USE LANDICE_COM, only : snowli,tlandi
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi
      USE TRACER_COM, only : NTM
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      USE EXCHANGE_TYPES
      USE LANDICE_COM, only : ijhc,ijhc_frac,ijhc_fhc,
     &       IJHC_PRECLI,IJHC_RUNLI
#ifdef GLINT2
      use fluxes, only : flice
      use landice_com, only : usedhp
#endif
      IMPLICIT NONE
      type(atmgla_xchng_vars) :: atmgla
      integer :: ihc

      real*8, dimension(:,:), pointer :: prec,eprec

      REAL*8 SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
      real*8, dimension(:,:,:), pointer :: trprec
#endif
C**** Get useful grid parameters
      INTEGER :: I, J
      INTEGER :: J_0, J_1, J_0H, J_1H ,I_0,I_1

      call getDomainBounds(GRID,J_STRT=J_0      , J_STOP=J_1      ,
     &              J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
      PRCP=atmgla%prec(i,j)   ! [kg m-2]
      atmgla%RUNO(I,J)=0
#ifdef TRACERS_WATER
      atmgla%TRUNO(:,I,J)=0.
#endif
      atmgla%IMPLM(I,J)=0.
      atmgla%IMPLH(I,J)=0.
#ifdef TRACERS_WATER
      atmgla%IMPLT(:,I,J)=0.
#endif
      ! demo diagnostic
      ijhc(i,j,ihc,ijhc_frac) = ijhc(i,j,ihc,ijhc_frac) +
     &       atmgla%ftype(i,j)
      ijhc(i,j,ihc,ijhc_fhc) = ijhc(i,j,ihc,ijhc_fhc) +
     &       atmgla%fhc(i,j)

#ifdef GLINT2
      IF (usedhp(i,j,ihc) /= 0 .and. PRCP.gt.0) THEN
#else
      IF (atmgla%ftype(i,j).gt.0 .and. PRCP.gt.0) THEN
#endif
        ENRGP=atmgla%eprec(I,J)      ! energy of precipitation
!      if (i==53.and.j==83.and.ihc==1)
!     &       print *,'AZZ PRECLI in',PRCP,ENRGP/PRCP

        SNOW=SNOWLI(I,J,IHC)
        TG1=TLANDI(1,I,J,IHC)
        TG2=TLANDI(2,I,J,IHC)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J,IHC)
        TRSNOW(:)=TRSNOWLI(:,I,J,IHC)
        TRPRCP(:)=atmgla%trprec(:,I,J)
#endif
        ! PRECLI() uses nothing, no globals involved.
        CALL PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *       TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *       EDIFS,DIFS,ERUN2,RUN0)

!      if (i==53.and.j==83.and.ihc==1)
!     &       print *,'AZZ PRECLI out',RUN0,ERUN2/RUN0

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J,IHC)=SNOW
        TLANDI(1,I,J,IHC)=TG1
        TLANDI(2,I,J,IHC)=TG2
        atmgla%RUNO(I,J)  =RUN0
        atmgla%GTEMP(I,J)=TLANDI(1,I,J,IHC)
        atmgla%GTEMP2(I,J)=TLANDI(2,I,J,IHC)
        atmgla%GTEMPR(I,J)   =TLANDI(1,I,J,IHC)+TF
#ifdef SCM
        if( SCMopt%Tskin )then
          atmgla%GTEMP(I,J) = SCMin%Tskin - TF
          atmgla%GTEMP2(I,J) = SCMin%Tskin - TF
          atmgla%GTEMPR(I,J) = SCMin%Tskin
        endif
#endif
#ifdef TRACERS_WATER
        TRLNDI(:,I,J,IHC)=TRLI(:)
        TRSNOWLI(:,I,J,IHC)=TRSNOW(:)
        atmgla%TRUNO(:,I,J)=TRUN0(:)
        IF (SNOW.gt.SNMIN) THEN
          atmgla%GTRACER(:,I,J)=TRSNOW(:)/SNOW
        ELSE
          atmgla%GTRACER(:,I,J)=TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
        atmgla%IMPLM(I,J)=DIFS
        atmgla%IMPLH(I,J)=ERUN2
#ifdef TRACERS_WATER
        atmgla%IMPLT(:,I,J)=atmgla%IMPLT(:,I,J)+TRDIFS(:)
#endif
        atmgla%E1(I,J)=EDIFS

        ! demo diagnostic
        ijhc(i,j,ihc,IJHC_PRECLI) = ijhc(i,j,ihc,IJHC_PRECLI)+prcp  ! [kg m-2]
        ijhc(i,j,ihc,IJHC_RUNLI)=ijhc(i,j,ihc,IJHC_RUNLI) +
     &        atmgla%RUNO(i,j)  ! [kg m-2]
      ELSE
        atmgla%IMPLM(I,J)=0.
        atmgla%IMPLH(I,J)=0.
#ifdef TRACERS_WATER
        atmgla%IMPLT(:,I,J)=0.
#endif
        atmgla%E1(I,J)=0.
      END IF
      END DO
      END DO

      END SUBROUTINE PRECIP_LI

      SUBROUTINE GROUND_LI(atmgla,ihc)
!@sum  GROUND_LI driver for applying surface fluxes to land ice fraction
!@auth Original Development team
!@ver  2010/10/06
!@calls LANDICE:LNDICE
      USE CONSTANT, only : tf
      USE MODEL_COM, only : dtsrc
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      USE GEOM, only : imaxj
      USE LANDICE, only : lndice,ace1li,ace2li,snmin
      USE SEAICE, only : rhos
      USE LANDICE_COM, only : snowli,tlandi
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi
      USE TRACER_COM, only : NTM
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      USE EXCHANGE_TYPES
      USE LANDICE_COM, only : ijhc
      USE LANDICE_COM, only : IJHC_SRFP
      USE LANDICE_COM, only : IJHC_SHDTLI,IJHC_EVHDT,IJHC_TRHDT
      USE LANDICE_COM, only : IJHC_F0LI,IJHC_EVAPLI,IJHC_TSLI
      USE LANDICE_COM, only : IJHC_IMPMLI,IJHC_IMPHLI, IJHC_RUNLI
#ifdef GLINT2
      use fluxes, only : flice
      use landice_com, only : usedhp
#endif

      IMPLICIT NONE
      type(atmgla_xchng_vars) :: atmgla
      integer :: ihc

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0
     *     ,SCOVLI
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
      INTEGER I,J,N
      INTEGER :: J_0,J_1, J_0H, J_1H ,I_0,I_1

      call startTimer('GROUND_LI()')
      call getDomainBounds(GRID,J_STRT=J_0      ,J_STOP=J_1
     &             ,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
      atmgla%RUNO(I,J)=0.
#ifdef GLINT2
      if (usedhp(i,j,ihc) /= 0) then
#else
      if (atmgla%ftype(i,j) > 0) then
#endif
        SNOW=SNOWLI(I,J,IHC)
        TG1=TLANDI(1,I,J,IHC)
        TG2=TLANDI(2,I,J,IHC)
        F0DT=atmgla%E0(I,J)
        F1DT=atmgla%E1(I,J)
        EVAP=atmgla%EVAPOR(I,J)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J,IHC)
        TRSNOW(:)=TRSNOWLI(:,I,J,IHC)
        TREVAP(:)=atmgla%TREVAPOR(:,I,J)
#ifdef TRACERS_DRYDEP
     *       -atmgla%trdrydep(:,i,j)
#endif
#endif
        CALL LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)
C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J,IHC)=SNOW
        TLANDI(1,I,J,IHC)=TG1
        TLANDI(2,I,J,IHC)=TG2
        atmgla%RUNO(I,J) = RUN0
        atmgla%GTEMP(I,J)=TLANDI(1,I,J,IHC)
        atmgla%GTEMP2(I,J)=TLANDI(2,I,J,IHC)
        atmgla%GTEMPR(I,J)   =TLANDI(1,I,J,IHC)+TF
#ifdef SCM
        if( SCMopt%Tskin )then
          atmgla%GTEMP(I,J) = SCMin%Tskin - TF
          atmgla%GTEMP2(I,J) = SCMin%Tskin - TF
          atmgla%GTEMPR(I,J) = SCMin%Tskin
        endif
#endif
#ifdef TRACERS_WATER
        TRLNDI(:,I,J,IHC)=TRLI(:)
        TRSNOWLI(:,I,J,IHC)=TRSNOW(:)
        atmgla%TRUNO(:,I,J)=TRUN0(:)
        IF (SNOW.gt.SNMIN) THEN
          atmgla%GTRACER(:,I,J)=TRSNOW(:)/SNOW
        ELSE
          atmgla%GTRACER(:,I,J)=TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
        atmgla%IMPLM(I,J)=atmgla%IMPLM(I,J)+DIFS
        atmgla%IMPLH(I,J)=atmgla%IMPLH(I,J)+EDIFS
#ifdef TRACERS_WATER
        atmgla%IMPLT(:,I,J)=atmgla%IMPLT(:,I,J)+TRDIFS(:)
#endif
        atmgla%E1(I,J)=EDIFS+F1DT
        SCOVLI=0
        IF (SNOWLI(I,J,IHC).GT.0.) SCOVLI=1.
        atmgla%snow(i,j) = snow
        atmgla%snowfr(i,j) = scovli
        atmgla%snowdp(i,j) = snow/rhos

        ! Put height-classified diagnostics here.
        ! Put PBL-related diagnostics in SURFACE_LANDICE.f
        ! This accumulation is run once per Model Timestep (dtsrc=3600s)
        ijhc(i,j,ihc,IJHC_SHDTLI)=ijhc(i,j,ihc,IJHC_SHDTLI) +
     &       atmgla%SENSHT(I,J)
        ijhc(i,j,ihc,IJHC_EVHDT)=ijhc(i,j,ihc,IJHC_EVHDT) +
     &       atmgla%LATHT(I,J)
        ijhc(i,j,ihc,IJHC_TRHDT)=ijhc(i,j,ihc,IJHC_TRHDT) +
     &       atmgla%TRHEAT(I,J)
        ijhc(i,j,ihc,IJHC_F0LI)=ijhc(i,j,ihc,IJHC_F0LI) +
     &       atmgla%E0(I,J) + atmgla%EPREC(I,J)
        ijhc(i,j,ihc,IJHC_EVAPLI)=ijhc(i,j,ihc,IJHC_EVAPLI) +
     &       atmgla%EVAPOR(I,J)

!        ijhc(i,j,ihc,IJHC_TSLI)=ijhc(i,j,ihc,IJHC_TSLI) +
!     &       atmgla%TSAVG(I,J)
        ijhc(i,j,ihc,IJHC_IMPMLI)=ijhc(i,j,ihc,IJHC_IMPMLI) +
     &       atmgla%IMPLM(I,J)
        ijhc(i,j,ihc,IJHC_IMPHLI)=ijhc(i,j,ihc,IJHC_IMPHLI) +
     &       atmgla%IMPLH(I,J)

        ijhc(i,j,ihc,IJHC_RUNLI)=ijhc(i,j,ihc,IJHC_RUNLI) +
     &        atmgla%RUNO(i,j)    ! [kg m-2]

        ijhc(i,j,ihc,IJHC_SRFP)=ijhc(i,j,ihc,IJHC_SRFP) +
     &        atmgla%SRFP(i,j) * dtsrc   ! [hPa s]

      END IF

C****
      END DO
      END DO

      call stopTimer('GROUND_LI()')
      END SUBROUTINE GROUND_LI

      SUBROUTINE conserv_MLI(ICE)
!@sum  conserv_MLI calculates total amount of snow and ice in land ice
!@auth Gavin Schmidt
      USE RESOLUTION, only : im,jm
      USE FLUXES, only : flice,atmgla
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : fhc,snowli
      USE LANDICE, only : lndice,ace1li,ace2li
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      IMPLICIT NONE
!@var ICE total land ice snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: ICE
      INTEGER I,J,IHC
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ICE(I,J)=FLICE(I,J)*
     &    (ACE1LI+ACE2LI+SUM(FHC(I,J,:)*SNOWLI(I,J,:)))
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) ICE(2:im,1) =ICE(1,1)
      IF(HAVE_NORTH_POLE) ICE(2:im,JM)=ICE(1,JM)

      RETURN
C****
      END SUBROUTINE conserv_MLI

      SUBROUTINE conserv_HLI(EICE)
!@sum  conserv_HLI calculates total land ice energy
!@auth Gavin Schmidt
      USE CONSTANT, only : shi,lhm
      USE RESOLUTION, only : im,jm
      USE FLUXES, only : flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : nhc,fhc,snowli,tlandi
      USE LANDICE, only : ace1li,ace2li
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      IMPLICIT NONE
!@var EICE total land ice energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: EICE
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      REAL*8 :: EICEIJ(NHC)

      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EICEIJ = ((TLANDI(1,I,J,:)*SHI-LHM)*(ACE1LI
     *         +SNOWLI(I,J,:))+(TLANDI(2,I,J,:)*SHI-LHM)*ACE2LI)
        EICE(I,J)=FLICE(I,J)*SUM(FHC(I,J,:)*EICEIJ)
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) EICE(2:im,1) =EICE(1,1)
      IF(HAVE_NORTH_POLE) EICE(2:im,JM)=EICE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_HLI

      Subroutine CONSERV_MICB (MICB)
!@sum  CONSERV_MICB calculates ice berg mass.  MICBIMP(1) + MDWNIMP_SH
!@+      is distributed spatially into MICB (kg/m^2) according to FSHGLM
!@auth Gavin Schmidt, Gary Russell
!@ver  2010/10/14
      USE RESOLUTION, only : im,jm
      Use GEOM,              Only: LAT2D,AXYP
      Use LANDICE_COM,       Only: FSHGLM,FNHGLM, MDWNIMP
      Use LANDICE,           Only: MICBIMP
      Use DOMAIN_DECOMP_ATM, Only: GRID,getDomainBounds, GLOBALSUM
      Implicit None
!@var MICB implicit mass of iceberg (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:grid%J_STOP_HALO) ::
     &        MICB, ARR_S,ARR_N
      Real*8  MDWNIMP_SH,MDWNIMP_NH
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      Logical :: QSP,QNP

      call getDomainBounds(GRID, J_STRT=J_0, J_STOP=J_1,
     &          HAVE_SOUTH_POLE=QSP, HAVE_NORTH_POLE=QNP)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Sum MDWNIMP (kg) into hemispheric values
      Do J=J_0,J_1  ;  Do I=I_0,I_1
         If (LAT2D(I,J) < 0)
     *      Then  ;  ARR_S(I,J) = MDWNIMP(I,J)
                     ARR_N(I,J) = 0
            Else  ;  ARR_S(I,J) = 0
                     ARR_N(I,J) = MDWNIMP(I,J)  ;  EndIf
      EndDo  ;  EndDo
      If (QSP)  ARR_S(2:IM,J_0) = ARR_S(1,J_0)
      If (QNP)  ARR_N(2:IM,J_1) = ARR_N(1,J_1)
      Call GLOBALSUM (GRID, ARR_S, MDWNIMP_SH, ALL=.True.)
      Call GLOBALSUM (GRID, ARR_N, MDWNIMP_NH, ALL=.True.)

C**** Distribute summed hemispheric values MICBIMP(:)+MDWNIMP_?H into
C**** array MICB(I,J) acording to FSHGLM and FNHGLM
      Do J=J_0,J_1  ;  Do I=I_0,I_1
         MICB(I,J) = (FSHGLM(I,J)*(MICBIMP(1)+MDWNIMP_SH) +
     +                FNHGLM(I,J)*(MICBIMP(2)+MDWNIMP_NH)) / AXYP(I,J)
      EndDo  ;  EndDo
      RETURN
      END SUBROUTINE conserv_MICB

      Subroutine CONSERV_HICB (HICB)
!@sum  CONSERV_HICB calculates ice berg mass.  EICBIMP(1) + EDWNIMP_SH
!@+      is distributed spatially into HICB (J/m^2) according to FSHGLM
!@auth Gavin Schmidt, Gary Russell
!@ver  2010/10/14
      USE RESOLUTION, only : im,jm
      Use GEOM,              Only: LAT2D,AXYP
      Use LANDICE_COM,       Only: FSHGLM,FNHGLM, EDWNIMP
      Use LANDICE,           Only: EICBIMP
      Use DOMAIN_DECOMP_ATM, Only: GRID,getDomainBounds, GLOBALSUM
      Implicit None
!@var EICB impicit iceberg energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &        HICB, ARR_S,ARR_N
      Real*8  EDWNIMP_SH,EDWNIMP_NH
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      Logical :: QSP,QNP

      call getDomainBounds(GRID, J_STRT=J_0, J_STOP=J_1,
     &          HAVE_SOUTH_POLE=QSP, HAVE_NORTH_POLE=QNP)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Sum EDWNIMP (J) into hemispheric values
      Do J=J_0,J_1  ;  Do I=I_0,I_1
         If (LAT2D(I,J) < 0)
     *      Then  ;  ARR_S(I,J) = EDWNIMP(I,J)
                     ARR_N(I,J) = 0
            Else  ;  ARR_S(I,J) = 0
                     ARR_N(I,J) = EDWNIMP(I,J)  ;  EndIf
      EndDo  ;  EndDo
      If (QSP)  ARR_S(2:IM,J_0) = ARR_S(1,J_0)
      If (QNP)  ARR_N(2:IM,J_1) = ARR_N(1,J_1)
      Call GLOBALSUM (GRID, ARR_S, EDWNIMP_SH, ALL=.True.)
      Call GLOBALSUM (GRID, ARR_N, EDWNIMP_NH, ALL=.True.)

C**** Distribute summed hemispheric values EICBIMP(:)+EDWNIMP_?H into
C**** array HICB(I,J) acording to FSHGLM and FNHGLM
      Do J=J_0,J_1  ;  Do I=I_0,I_1
         HICB(I,J) = (FSHGLM(I,J)*(EICBIMP(1)+EDWNIMP_SH) +
     +                FNHGLM(I,J)*(EICBIMP(2)+EDWNIMP_NH)) / AXYP(I,J)
      EndDo  ;  EndDo
      RETURN
      END SUBROUTINE conserv_HICB

      SUBROUTINE daily_LI
!@sum  daily_ice does daily landice things
!@auth Gavin Schmidt
!@ver  2010/10/13
      USE CONSTANT, only : lhm,shi
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only : dtsrc, modelEclock, modelEclockI
     &     , itime, itimei, nday, calendar
      use TimeConstants_mod, only: SECONDS_PER_DAY, EARTH_DAYS_PER_YEAR,
     &                             INT_DAYS_PER_YEAR
      USE GEOM, only : axyp,imaxj,lat2d
      USE LANDICE, only: ace1li,ace2li,glmelt_on,glmelt_fac_nh
     *     ,glmelt_fac_sh,fwarea_sh,fwarea_nh,accpda,accpdg,eaccpda
     *     ,eaccpdg,micbimp,eicbimp
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,traccpda,traccpdg   ! tricbimp
#endif
      USE TRDIAG_COM, only : to_per_mil
      use OldTracer_mod, only :  trw0, nWater, trname
#endif    /* TNL: inserted */
      Use LANDICE_COM, Only: TLANDI,SNOWLI, MDWNIMP,EDWNIMP,
     *                       FSHGLM,FNHGLM
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trdwnimp
      USE TRACER_COM, only : NTM
#endif
      USE FLUXES, only : atmocn,flice,focean
      USE Dictionary_mod
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, 
     &     GLOBALSUM, AM_I_ROOT
      use Time_mod, only: Time
      use Rational_mod

      IMPLICIT NONE
!@var gm_relax Glacial Melt relaxation parameter (1/year)
      REAL*8, PARAMETER :: gm_relax = 0.1d0  ! 10 year relaxation

      Real*8 MDWNIMP_SH,MDWNIMP_NH, EDWNIMP_SH,EDWNIMP_NH, DTSRCzY
     &     ,byoarea
#ifdef TRACERS_WATER
      REAL*8 trdwnimp_SH(NTM),trdwnimp_NH(NTM)
#endif
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO) :: mask_s,arr_s,arr_n
      INTEGER :: J_0,J_1,I_0,I_1,I,J,ITM
      LOGICAL :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE
      type (Time) :: now, startTime
      type (Rational) :: oneYear
      logical :: atLeastOneYearHasPassed

      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Keep track of accumlated implicit iceberg amount
C**** hemispheric masking (why isn't this in GEOM?)
      do j=j_0,j_1; do i=i_0,i_1
        if(lat2d(i,j).lt.0.) then
          mask_s(i,j) = 1.
        else
          mask_s(i,j) = 0.
        endif
      enddo; enddo

C**** Subtract daily ice berg mass that is added to ocean in DAILY_OC
      MICBIMP(:) = MICBIMP(:)-(/  ACCPDA,  ACCPDG /)/EARTH_DAYS_PER_YEAR
      EICBIMP(:) = EICBIMP(:)-(/ EACCPDA, EACCPDG /)/EARTH_DAYS_PER_YEAR
c#ifdef TRACERS_WATER
C     DO N=1,NTM
C       TRICBIMP(N,:) = TRCIBIMP(N,:) -
C    -                  (/ TRACCPDA(:), TRACCPDG(:) /) / EARTH_DAYS_PER_YEAR
C     END DO
c#endif

C**** Every year update gmelt factors in order to balance downward
C**** implicit fluxes. If this is used, then ice sheets/snow are FORCED to
C**** be in balance. This may not be appropriate for transient runs but
C**** we aren't getting that right anyway.

      If (modelEclock%getDayOfYear()==1)  Then  !  Jan 1 only, EndIf at 500

C**** Calculate mass/energy/tracer accumulation for the past year
        do j=j_0,j_1; do i=i_0,i_1
          arr_s(i,j) = mdwnimp(i,j)*mask_s(i,j)
          arr_n(i,j) = mdwnimp(i,j)*(1.-mask_s(i,j))
        enddo; enddo
        If (have_south_pole) ARR_S(2:IM,J_0) = ARR_S(1,J_0)
        If (have_north_pole) ARR_N(2:IM,J_1) = ARR_N(1,J_1)
        CALL GLOBALSUM(grid, arr_s, mdwnimp_SH ,ALL=.TRUE.)
        CALL GLOBALSUM(grid, arr_n, mdwnimp_NH ,ALL=.TRUE.)
        do j=j_0,j_1; do i=i_0,i_1
          arr_s(i,j) = edwnimp(i,j)*mask_s(i,j)
          arr_n(i,j) = edwnimp(i,j)*(1.-mask_s(i,j))
        enddo; enddo
        If (have_south_pole) ARR_S(2:IM,J_0) = ARR_S(1,J_0)
        If (have_north_pole) ARR_N(2:IM,J_1) = ARR_N(1,J_1)
        CALL GLOBALSUM(grid, arr_s, edwnimp_SH ,ALL=.TRUE.)
        CALL GLOBALSUM(grid, arr_n, edwnimp_NH ,ALL=.TRUE.)
        MICBIMP(:) = MICBIMP(:) + (/ MDWNIMP_SH, MDWNIMP_NH /)
        EICBIMP(:) = EICBIMP(:) + (/ EDWNIMP_SH, EDWNIMP_NH /)
C****   MICBIMP and EICBIMP are now only correct for root processor
#ifdef TRACERS_WATER
        DO ITM=1,NTM
          do j=j_0,j_1; do i=i_0,i_1
            arr_s(i,j) = trdwnimp(itm,i,j)*mask_s(i,j)
            arr_n(i,j) = trdwnimp(itm,i,j)*(1.-mask_s(i,j))
          enddo; enddo
          If (have_south_pole) ARR_S(2:IM,J_0) = ARR_S(1,J_0)
          If (have_north_pole) ARR_N(2:IM,J_1) = ARR_N(1,J_1)
          CALL GLOBALSUM(grid, arr_s, trdwnimp_SH(itm) ,ALL=.TRUE.)
          CALL GLOBALSUM(grid, arr_n, trdwnimp_NH(itm) ,ALL=.TRUE.)
C         TRICBIMP(ITM,:) = TRCIBIMP(ITM,:) +
C    +                      (/ TRDWNIMP_SH(:), TRDWNIMP_NH(:) /)
        END DO
#endif

! only adjust after at least one full year
        now = modelEclock%getCurrentTime()
        startTime = modelEclockI%getCurrentTime()
        oneYear = 
     &       calendar%getSecondsPerDay() * calendar%getMaxDaysInYear()
        atLeastOneYearHasPassed = (now >= startTime + oneYear)

        if (atLeastOneYearHasPassed .and. GLMELT_ON==1) then
                                                        !  EndIf at 400
C*** prevent iceberg sucking
          if(mdwnimp_NH.lt.0) then
            write(99,*) "Limiting NH icesheet replacement mass/energy",
     *           mdwnimp_NH,edwnimp_NH
            mdwnimp_NH=0.
            edwnimp_NH=0.
#ifdef TRACERS_WATER
            trdwnimp_NH=0.
#endif
          endif

          if(mdwnimp_SH.lt.0) then
            write(99,*) "Limiting SH icesheet replacement mass/energy",
     *           mdwnimp_SH,edwnimp_SH
            mdwnimp_SH=0.
            edwnimp_SH=0.
#ifdef TRACERS_WATER
            trdwnimp_SH=0.
#endif
          endif

#ifdef TRACERS_WATER
          DO ITM=1,NTM
            if (to_per_mil(itm).gt.0) then
              if (mdwnimp_NH .gt. 0) print*,'Tracers(NH) (permil):',itm
     *             ,(trdwnimp_NH(ITM)/mdwnimp_NH/trw0(itm)-1.)*1000
     *             ,trdwnimp_NH(ITM),mdwnimp_NH
              if (mdwnimp_SH .gt. 0)  print*,'Tracers(SH) (permil):',itm
     *             ,(trdwnimp_SH(ITM)/mdwnimp_SH/trw0(itm)-1.)*1000
     *             ,trdwnimp_SH(ITM),mdwnimp_SH
            else
              print*,'Tracers(NH) (kg):',itm,trdwnimp_NH(ITM),mdwnimp_NH
              print*,'Tracers(SH) (kg):',itm,trdwnimp_SH(ITM),mdwnimp_SH
            end if
          END DO
#endif

C**** adjust hemispheric mean glacial melt amounts (only on root processor)
      if (AM_I_ROOT()) THEN
          write(6,*) "Adjusting glacial melt: ", 
     &         modelEclock%getDayOfYear()
     *        ,modelEclock%getYear()
          write(6,*) "Mass (before): ",accpda,accpdg,mdwnimp_SH
     *         ,mdwnimp_NH
          write(6,*) "Temp (before): ",(eaccpda/accpda+lhm)/shi
     *         ,(eaccpdg/accpdg+lhm)/shi
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
          do itm=1,ntm
            if (to_per_mil(itm).gt.0) then
              write(6,*) trim(trname(itm))," (before) (permil)",1000
     *             *(traccpda(itm)/accpda/trw0(itm)-1.)
            else
              write(6,*) trim(trname(itm))," (before) (kg)",traccpda(itm
     *             ),accpda
            end if
          end do
#endif
#endif   /* TNL: inserted */
        EndIf

         accpda =  accpda + gm_relax*(mdwnimp_SH -  accpda)
         accpdg =  accpdg + gm_relax*(mdwnimp_NH -  accpdg)
        eaccpda = eaccpda + gm_relax*(edwnimp_SH - eaccpda)
        eaccpdg = eaccpdg + gm_relax*(edwnimp_NH - eaccpdg)
        If (AM_I_ROOT())  Then
          write(6,*) "Mass (after): ",accpda,accpdg
          write(6,*) "Temp (after): ",(eaccpda/accpda+lhm)/shi,
     *         (eaccpdg/accpdg+lhm)/shi
        EndIf

#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
        traccpda(:)=traccpda(:)+gm_relax*(trdwnimp_SH(:)-traccpda(:))
        traccpdg(:)=traccpdg(:)+gm_relax*(trdwnimp_NH(:)-traccpdg(:))
        If (AM_I_ROOT())  Then
          do itm=1,ntm
            if (to_per_mil(itm).gt.0) then
              write(6,*) trim(trname(itm))," (after) (permil)",1000
     *             *(traccpda(itm)/accpda/trw0(itm)-1.)
            else
              write(6,*) trim(trname(itm))," (after) (kg)",traccpda(itm)
     *             ,accpda
            end if
          end do
        ENDIF
#endif
#endif   /* TNL: inserted */

C**** Set GMELT (kg of ice per source time step
C****   from total hemisphere ACCPDA or ACCPDG (kg per year)
      DTSRCzY = DTSRC / (SECONDS_PER_DAY*EARTH_DAYS_PER_YEAR)
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          If (LAT2D(I,J) < 0)  Then  !  SH
             atmocn% GMELT(I,J) =  ACCPDA*FSHGLM(I,J)*DTSRCzY  !  kg
             atmocn%EGMELT(I,J) = EACCPDA*FSHGLM(I,J)*DTSRCzY  !  J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
              atmocn%TRGMELT(:,I,J)=TRACCPDA(:)*FSHGLM(I,J)*DTSRCzY  !  kg
#endif
#endif    /* TNL: inserted */
          Else  !  NH
             atmocn% GMELT(I,J) =  ACCPDG*FNHGLM(I,J)*DTSRCzY  !  kg
             atmocn%EGMELT(I,J) = EACCPDG*FNHGLM(I,J)*DTSRCzY  !  J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
             atmocn%TRGMELT(:,I,J) = TRACCPDG(:)*FNHGLM(I,J)*DTSRCzY  !  kg
#endif
#endif    /* TNL: inserted */
          END IF
          IF (FOCEAN(I,J).GT.0.) THEN
            byoarea = 1.d0/(AXYP(I,J)*FOCEAN(I,J))
            atmocn% GMELT(I,J) = atmocn% GMELT(I,J)*byoarea
            atmocn%EGMELT(I,J) = atmocn%EGMELT(I,J)*byoarea
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            atmocn%TRGMELT(:,I,J) = atmocn%TRGMELT(:,I,J)*byoarea
#endif
#endif /* TNL: inserted */
          ENDIF
        END DO
      END DO

  400 EndIf  !  run not in its first year

C**** Add MDWNIMP to MICBIMP and reset implicit accumulators
      MDWNIMP=0.
      EDWNIMP=0.
#ifdef TRACERS_WATER
      TRDWNIMP=0.
#endif

  500 EndIf  !  start of new year

      END SUBROUTINE daily_LI

      SUBROUTINE CHECKLI (SUBR)
!@sum  CHECKLI checks whether the land ice variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : teeny
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only : qcheck
      USE FLUXES, only : flice
      USE DOMAIN_DECOMP_ATM, only : HALO_UPDATE, getDomainBounds, GRID
      USE GEOM, only : imaxj
#ifdef TRACERS_WATER
      use OldTracer_mod, only: trname
      USE TRACER_COM, only : NTM
#endif
      USE LANDICE, only : ace1li, ace2li
      USE LANDICE_COM, only : nhc,tlandi,snowli,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trdwnimp
#endif
      IMPLICIT NONE
      INTEGER :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0H,I_1H,I_0,I_1,njpol
      INTEGER I,J,N,IHC !@var I,J loop variables
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKL
#ifdef TRACERS_WATER
      integer :: imax,jmax
      real*8 relerr,errmax
#endif
      call getDomainBounds(grid, J_STRT=J_0,      J_STOP=J_1,
     *               J_STRT_HALO=J_0H,J_STOP_HALO=J_1H,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in land ice data
      CALL CHECK4B(tlandi(:,I_0:I_1,J_0:J_1,:),2,I_0,I_1,J_0,J_1,NJPOL,
     &     NHC,SUBR,'tli')
      CALL CHECK3B(snowli(I_0:I_1,J_0:J_1,:) ,I_0,I_1,J_0,J_1,NJPOL,NHC,
     &     SUBR,'sli')
      CALL CHECK3B(mdwnimp(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'mdw')
      CALL CHECK3B(edwnimp(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'edw')

      QCHECKL = .FALSE.
#ifdef TRACERS_WATER
      do ihc=1,nhc
      do n=1,ntm
C**** Check conservation of water tracers in land ice
        if (trname(n).eq.'Water') then
          errmax = 0. ; imax=I_0 ; jmax=J_0
          do j=J_0, J_1
          do i=I_0,imaxj(j)
            if (flice(i,j).gt.0) then
              relerr=max(
     *             abs(trlndi(n,i,j,ihc)/(ace1li+ace2li)-1.)
     *             ,abs((trsnowli(n,i,j,ihc)-snowli(i,j,ihc))/
     *             (snowli(i,j,ihc)+teeny
     *             )),abs((trdwnimp(n,i,j)-mdwnimp(i,j))/(mdwnimp(i,j)
     *             +teeny)))
            else
              relerr=abs(trdwnimp(n,i,j)-mdwnimp(i,j))/(mdwnimp(i,j)
     *             +teeny)
            end if
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          print*,"Relative error in land ice mass after ",trim(subr),":"
     *         ,imax,jmax,errmax,trlndi(n,imax,jmax,ihc),ace1li+ace2li
     *         ,trsnowli(n,imax,jmax,ihc),snowli(imax,jmax,ihc)
     *         ,trdwnimp(n,imax,jmax),mdwnimp(imax,jmax)
        end if
      end do
      end do
#endif

      IF (QCHECKL)
     &call stop_model('CHECKLI: Land Ice variables out of bounds',255)
      RETURN
C****
      END SUBROUTINE CHECKLI
