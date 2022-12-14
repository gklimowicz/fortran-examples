!@sum  Routines to calculates solar zenith angles, weighted by time/sunlight
!@auth Original Development Team
!@auth Rewritten for non-latlon grids by M. Kelley
      module RAD_COSZ0
      USE CONSTANT, only : twopi,pi,teeny
      USE GEOM, only : lon2d,sinlat2d,coslat2d
      USE DOMAIN_DECOMP_ATM, ONLY: grid
      implicit none
      save
      private
      public :: coszs,coszt,cosz_init,daily_cosz

!@var SIND,COSD sin,cos of solar declination angle
!@+   (these are local copies that are set by daily_cosz)
      real*8 :: sind,cosd

c
c local variables common to coszs and coszt routines
c
      REAL*8 :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2,LT1_SV
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8, PARAMETER :: ZERO1=1.D-2
      INTEGER I,J
      REAL*8 ZERO2,DUSK,DAWN,CJCD,SJSD,DROT
      INTEGER :: I_0, I_1, J_0, J_1
      logical :: lwrap

      real*8, dimension(:,:), allocatable :: duskij,sinlatij,coslatij

      logical :: use_const_cosz
      real*8 :: const_cosz

      contains

      subroutine cosz_init(cosz_const)
      implicit none
      real*8, optional :: cosz_const
      if(present(cosz_const)) then
        const_cosz = cosz_const
        use_const_cosz = .true.
      else
        use_const_cosz = .false.
      endif
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      allocate(duskij(i_0:i_1,j_0:j_1))
      allocate(sinlatij(i_0:i_1,j_0:j_1))
      allocate(coslatij(i_0:i_1,j_0:j_1))
c use center lons/lats for now 
      do j=j_0,j_1
      do i=i_0,i_1
        sinlatij(i,j) = sinlat2d(i,j)
        coslatij(i,j) = coslat2d(i,j)
      enddo
      enddo
      return
      end subroutine cosz_init

      subroutine COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****

      IMPLICIT NONE
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ

      if(use_const_cosz) then
        cosz = const_cosz
        return
      endif

      DROT=ROT2-ROT1

C****
C**** LOOP OVER GRID BOXES
C****
      DO J=J_0,J_1
      DO I=I_0,I_1

      cosz(i,j) = 0d0

      DUSK=DUSKIJ(I,J)

      if(dusk.eq.0d0) cycle ! constant night at this latitude

      DAWN=-DUSK

      LT1=ROT1+LON2D(I,J)
      LT2=ROT2+LON2D(I,J)

      call adjust_angles

      if(dawn.gt.lt2 .or. dusk.lt.lt1) then ! constant darkness in this cell
c        cosz(i,j) = 0d0
      else
        SJSD=SINLATIJ(I,J)*SIND
        CJCD=COSLATIJ(I,J)*COSD
        lt1 = max(lt1,dawn)
        lt2 = min(lt2,dusk)
        slt1 = sin(lt1)
        slt2 = sin(lt2)
        COSZ(I,J)=SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
        if(lwrap) then ! second contribution: lt1 to dusk
          lt1 = lt1_sv
          lt2 = dusk
          slt1 = sin(lt1)
          slt2 = sin(lt2)
          COSZ(I,J)=COSZ(I,J)+SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
        endif
        COSZ(I,J)=COSZ(I,J)/DROT
      endif

      ENDDO ! I
      ENDDO ! J

      RETURN
      END subroutine COSZT

      subroutine COSZS (ROT1,ROT2,COSZ,COSZA)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.
C****
      IMPLICIT NONE
      REAL*8 ROT1,ROT2
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ,COSZA
      REAL*8 ECOSZ,ECOSQZ,ECOSZ1

      if(use_const_cosz) then
        cosz = const_cosz
        cosza = const_cosz
        return
      endif

      DROT=ROT2-ROT1

C****
C**** LOOP OVER GRID BOXES
C****
      DO J=J_0,J_1
      DO I=I_0,I_1

      cosz(i,j) = 0d0
      cosza(i,j) = 0d0

      DUSK=DUSKIJ(I,J)

      if(dusk.eq.0d0) cycle ! constant night at this latitude

      DAWN=-DUSK

      LT1=ROT1+LON2D(I,J)
      LT2=ROT2+LON2D(I,J)

      call adjust_angles

      if(dawn.gt.lt2.or.dusk.lt.lt1) then ! constant darkness in this cell
c        cosz(i,j) = 0d0
c        cosza(i,j) = 0d0
      else
        SJSD=SINLATIJ(I,J)*SIND
        CJCD=COSLATIJ(I,J)*COSD
        lt1 = max(lt1,dawn)
        lt2 = min(lt2,dusk)
        slt1 = sin(lt1)
        slt2 = sin(lt2)
        S2LT1 = 2.*SLT1*cos(lt1)
        S2LT2 = 2.*SLT2*cos(lt2)
        ECOSZ=SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2-SLT1)+
     *       .5*CJCD*(LT2-LT1+.5*(S2LT2-S2LT1)))
        if(lwrap) then ! second contribution: lt1 to dusk
          lt1 = lt1_sv
          lt2 = dusk
          slt1 = sin(lt1)
          slt2 = sin(lt2)
          S2LT1=2.*SLT1*cos(lt1)
          S2LT2=2.*SLT2*cos(lt2)
          ECOSZ1=SJSD*(LT2-LT1)+CJCD*(SLT2-SLT1)
          ECOSQZ = ECOSQZ + SJSD*ECOSZ1+CJCD*(SJSD*(SLT2-SLT1)+
     *         .5*CJCD*(LT2-LT1+.5*(S2LT2-S2LT1)))
          ECOSZ = ECOSZ + ECOSZ1
        endif
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      endif

      ENDDO ! I
      ENDDO ! J

      return
      end subroutine COSZS

      subroutine daily_cosz(SIND_in,COSD_in, COSZ_day,DUSK)
c Resets parameters needed for zenith angle calculations, and calculates
c daily-average cosz and the hour of dusk
      implicit none
      REAL*8 :: SIND_in,COSD_in
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     COSZ_day, ! average of cos(zenith angle) for the current day
     &     DUSK      ! hour of dusk, radians from local noon
      REAL*8 :: CJCD_SDUSK
      SIND = SIND_in
      COSD = COSD_in
C**** COMPUTE DUSK FOR EACH GRID CELL AND DAILY-AVERAGE COSZ
      DO J=J_0,J_1
      DO I=I_0,I_1
        SJSD=SINLATIJ(I,J)*SIND
        CJCD=COSLATIJ(I,J)*COSD
        IF(CJCD-SJSD.LT.0.) THEN     ! constant daylight
          DUSKIJ(I,J) = PI
          COSZ_day(I,J) = SJSD
        ELSEIF(CJCD+SJSD.LT.0.) THEN ! constant darkness
          DUSKIJ(I,J) = 0D0
          COSZ_day(I,J) = 0D0
        ELSE
          DUSKIJ(I,J) = ACOS(-SJSD/CJCD)
          CJCD_SDUSK = SQRT((CJCD-SJSD)*(CJCD+SJSD))
          COSZ_day(I,J) = max(0d0,(SJSD*DUSKIJ(I,J)+CJCD_SDUSK)/PI)
        ENDIF
        DUSK(I,J) = DUSKIJ(I,J)
      ENDDO
      ENDDO
      return
      end subroutine daily_cosz

      subroutine adjust_angles
c
c domain is lt=-pi:+pi, noon is at lt=0
c
      implicit none
      do while(lt1.gt.pi)
        lt1 = lt1 - twopi
      enddo
      do while(lt2.gt.pi)
        lt2 = lt2 - twopi
      enddo
      lwrap = .false.
      if(lt1.gt.lt2) then
        if(lt1.le.dusk .and. lt2.gt.dawn) then
c
c daylight at lt1 and lt2, brief darkness in between.
c computation will be divided into two parts:
c 1. dawn to lt2
c 2. lt1 to dusk
c
          lt1_sv = lt1
          lt1 = dawn
          lwrap = .true.
        elseif(lt1.lt.dusk .and. lt2.lt.dawn) then
          lt2 = dusk
        elseif(lt2.gt.dawn .and. lt1.gt.dusk) then
          lt1 = dawn
        endif
      endif
      end subroutine adjust_angles

      end module RAD_COSZ0
