#include "rundeck_opts.h"
      SUBROUTINE dust_emission_constraints(itype,wsgcm,pbl_args)
!@sum  local constrainsts for dust tracer emission valid for all dust bins
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      USE model_com,ONLY : dtsrc
      use fluxes, only : nisurf
      USE socpbl,ONLY : t_pbl_args
      USE trdust_mod,ONLY : imDust,lim,ljm,lkm,table,x1,x2,x3
      use trdust_drv, only: wsgInterp
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: itype
      REAL*8,INTENT(IN) :: wsgcm

      type(t_pbl_args),INTENT(INOUT) :: pbl_args

      REAL*8 :: snowe,vtrsh
      REAL*8 :: dsteve1,dsteve2
      REAL*8 :: soilvtrsh
      LOGICAL :: qdust
      REAL*8 :: dryhr,hbaij,pprec,pevap,ricntd
      REAL*8 :: hbaijd,hbaijold
      LOGICAL :: pmei
      REAL*8 :: wearth,aiearth,wfcs,pdfint,wsubtke,wsubwd,wsubwm
      REAL(KIND=8) :: soilwet,sigma,ans,dy,workij1,workij2,wsgcm1,mcfrac
      CHARACTER(17) :: fname='WARNING_in_TRDUST'
      CHARACTER(25) :: subr='dust_emission_constraints'
      CHARACTER(5) :: vname1='wsgcm',vname2='sigma'
      CHARACTER(9) :: vname3='soilvtrsh'

c**** input
      snowe=pbl_args%snow
      vtrsh=pbl_args%vtrsh
      dryhr=pbl_args%dryhr
      pprec=pbl_args%pprec
      pevap=pbl_args%pevap
      hbaij=pbl_args%hbaij
      ricntd=pbl_args%ricntd
      wearth=pbl_args%wearth
      aiearth=pbl_args%aiearth
      wfcs=pbl_args%wfcs
      wsubtke=pbl_args%wsubtke
      wsubwd=pbl_args%wsubwd
      wsubwm=pbl_args%wsubwm
      mcfrac=pbl_args%mcfrac

      qdust = .false.
      dsteve1 = 0.D0
      dsteve2 = 0.D0
      soilvtrsh = 0.D0
      pdfint = 0.D0

      IF (imDUST == 2) THEN

c**** legacy dust emission
c     Checking if accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission

      hbaijold=hbaij
      hbaij=hbaijold+pprec/nisurf-pevap
      hbaijd=hbaij-hbaijold
      IF (itype == 4 .AND. hbaijd <= 0.D0) THEN
        ricntd=ricntd+Dtsrc/3600./nisurf
        IF (ricntd >= dryhr .AND. dryhr /= 0.D0) THEN
          pmei=.TRUE.
        ELSE
          pmei=.FALSE.
        END IF
      ELSE
        ricntd=0.D0
        pmei=.FALSE.
      END IF

      IF (vtrsh > 0.D0 .AND. wsgcm > vtrsh) THEN
        dsteve2=1.D0
      ELSE
        dsteve2=0.D0
      END IF
      IF (pmei .AND. snowe <= 1 .AND. vtrsh > 0.D0 .AND. wsgcm > vtrsh)
     &     THEN
        dsteve1=1.D0
        qdust=.TRUE.
      ELSE
        dsteve1=0.D0
        qdust=.FALSE.
      END IF

      ELSE IF ( imDUST == 0 .or. imDust >= 3 ) THEN

c**** dust emission using probability density function of wind speed
      IF (itype == 4 .AND. snowe <= 1.D0) THEN
        qdust=.TRUE.
      ELSE
        qdust=.FALSE.
        dsteve1=0.D0
        dsteve2=0.D0
        soilvtrsh=0.D0
      END IF

      IF (qdust) THEN

        soilwet=(wearth+aiearth)/(wfcs+1.D-20)
        if (soilwet.gt.1.D0) soilwet=1.d0
        soilvtrsh=8.d0*(exp(0.7d0*soilwet))

        pdfint=0.d0
        workij1=0.d0
        workij2=0.d0
        wsgcm1=wsgcm

c     There is no moist convection, sigma is composed of TKE and DRY
c     convective velocity scale
        IF (wsubwm == 0.D0) THEN
          sigma=wsubtke+wsubwd
c     No need to calculate the emission below these values since
c     the emission is zero
          IF (sigma > 0.1D0 .OR. wsgcm1 > 1.D0) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005D0 .AND. wsgcm1 > 1.D0) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                pdfint=(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                pdfint=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1D0 .AND. wsgcm1 < 0.0005D0) THEN
              wsgcm1=0.0005d0
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              ans = wsgInterp%interpolate3dlin([wsgcm1, sigma, 
     &          soilvtrsh])
              pdfint=exp(ans)
            ELSE
              CALL check_upper_limit(wsgcm1,x1(Lim),fname,subr,vname1)
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              ans = wsgInterp%interpolate3dlin([wsgcm1, sigma, 
     &          soilvtrsh])
              pdfint=exp(ans)
            END IF
          END IF

        ELSE

c     When there is moist convection, the sigma is the combination of
c     all three subgrid scale parameters (i.e. independent or dependent)
c     Takes into account that the moist convective velocity scale acts
c     only over the area with downdrafts (mcfrac).

          sigma=wsubtke+wsubwd+wsubwm
c     No need to calculate the emission below these values since
c     the emission is Zero
          IF (sigma > 0.1D0 .OR. wsgcm1 > 1.D0) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005D0 .AND. wsgcm1 > 1.D0) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                workij1=mcfrac*(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                workij1=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1D0 .AND. wsgcm1 < 0.0005D0) THEN
              wsgcm1=0.0005d0
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              ans = wsgInterp%interpolate3dlin([wsgcm1, sigma, 
     &          soilvtrsh])
              workij1=mcfrac*exp(ans)
            ELSE
              CALL check_upper_limit(wsgcm1,x1(Lim),fname,subr,vname1)
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              ans = wsgInterp%interpolate3dlin([wsgcm1, sigma, 
     &          soilvtrsh])
              workij1=mcfrac*exp(ans)
            END IF
          END IF

          sigma=wsubtke+wsubwd
c     No need to calculate the emission below these values since
c     the emission is Zero
          IF (sigma > 0.1D0 .OR. wsgcm1 > 1.D0) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005D0 .AND. wsgcm1 > 1.D0) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                workij2=(1.d0-mcfrac)*(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                workij2=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1D0 .AND. wsgcm1 < 0.0005D0) THEN
              wsgcm1=0.0005d0
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              ans = wsgInterp%interpolate3dlin([wsgcm1, sigma, 
     &          soilvtrsh])
              workij2=(1.d0-mcfrac)*exp(ans)
            ELSE
              CALL check_upper_limit(wsgcm1,x1(Lim),fname,subr,vname1)
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              ans = wsgInterp%interpolate3dlin([wsgcm1, sigma, 
     &          soilvtrsh])
              workij2=(1.d0-mcfrac)*exp(ans)
            END IF
          END IF
          pdfint=workij1+workij2
        END IF

        IF (sigma == 0.D0) THEN
          IF (wsgcm1 > soilvtrsh) THEN
            pdfint=(wsgcm1-soilvtrsh)*wsgcm1**2.D0
          ELSE
            pdfint=0.D0
          END IF
        END IF

        IF (pdfint > 0.D0) THEN
          dsteve1=1.D0
        ELSE
          dsteve1=0.D0
        END IF

        IF (vtrsh > 0.D0 .AND. wsgcm1 > vtrsh) THEN
          dsteve2=1.D0
        ELSE
          dsteve2=0.D0
        END IF

      END IF

      ELSE IF (imDUST == 1) THEN
        IF (itype == 4) THEN
          qdust=.TRUE.
        ELSE
          qdust=.FALSE.
        END IF
        soilvtrsh=0.D0
      END IF

c**** output
      pbl_args%dust_event1=dsteve1
      pbl_args%dust_event2=dsteve2
      pbl_args%qdust=qdust
      pbl_args%hbaij=hbaij
      pbl_args%ricntd=ricntd
      pbl_args%wtrsh=soilvtrsh
      pbl_args%pdfint=pdfint

#endif

      RETURN
      END SUBROUTINE dust_emission_constraints

      SUBROUTINE local_dust_emission(n,wsgcm,pbl_args,dsrcflx,
     &     dsrcflx2)
!@sum  selects routine for calculating local dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use RunTimeControls_mod, only : tracers_minerals
      USE socpbl,ONLY : t_pbl_args
      use tracer_com, only: Ntm_dust
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
     &     , n_soildust
#endif
      use OldTracer_mod, only: trname
      use trdust_mod,only : nDustBins, CWiCub, FClWiCub, FSiWiCub,
     &     CWiPdf, scaleDustEmission, fracClayPDFscheme,
     &     fracSiltPDFscheme, imDust

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: n
      REAL*8,INTENT(IN) :: wsgcm
      TYPE(t_pbl_args),INTENT(IN) :: pbl_args

      REAL*8,INTENT(OUT) :: dsrcflx,dsrcflx2

      integer :: n1, n_bin
      REAL*8 :: vtrsh
      real( kind=8 ) :: d_dust( nDustBins )
      REAL*8 :: frtrac
      LOGICAL :: qdust
      REAL*8 :: frclay,frsilt
      real(kind=8) :: ers_data,dustSourceFunction,soilvtrsh,pdfint
      real(kind=8) :: mineralFractions( max( nDustBins, ntm_dust ) )
      real( kind=8 ) :: zsum

c**** input
      qdust=pbl_args%qdust
      vtrsh=pbl_args%vtrsh
      IF ( imDust == 1 .or. imDust == 3 .or. imDust == 5 ) d_dust( : ) =
     &     pbl_args%d_dust( : )
      frclay=pbl_args%frclay
      frsilt=pbl_args%frsilt
      ers_data=pbl_args%ers_data
      dustSourceFunction = pbl_args%dustSourceFunction
      soilvtrsh=pbl_args%wtrsh
      pdfint=pbl_args%pdfint
      mineralFractions( : ) = pbl_args%mineralFractions( : )

c**** initialize
      dsrcflx=0.D0
      dsrcflx2=0.D0
      frtrac = scaleDustEmission
      IF (qdust) THEN

      IF (imDUST /= 1) THEN
c**** Interactive dust emission

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)

        SELECT CASE(trname(n))
        case('Clay','ClayIlli','ClayKaol','ClaySmec','ClayCalc'
     &         ,'ClayQuar','ClayFeld','ClayHema','ClayGyps','ClayIlHe'
     &         ,'ClayKaHe','ClaySmHe','ClayCaHe','ClayQuHe','ClayFeHe'
     &         ,'ClayGyHe')
        IF ( imDust == 0 .or. imDust >= 3 ) THEN
            frtrac = FracClayPDFscheme*frtrac
          ELSE IF (imDust == 2) THEN
            frtrac=FClWiCub*frclay*frtrac
          END IF
        case('Silt1','Silt2','Silt3','Silt4','Silt5','Sil1Quar'
     &         ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &         ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &         ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &         ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &         ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &         ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &         ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &         ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &         ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &         ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &         ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &         ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &         ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &         ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &         ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')
        IF ( imDust == 0 .or. imDust >= 3 ) THEN
            frtrac = FracSiltPDFscheme*frtrac
          ELSE IF (imDust == 2) THEN
            frtrac=FSiWiCub*frsilt*frtrac
          END IF
        case default
          return
        END SELECT

        if ( imDust == 3 .or. imDust == 5 ) then
c*****    shape size distribution of dust emission flux according to
c         AeroCom size distribution (imDust=3) or use AeroCom
c         distribution as additional mask for dust emission (imDust=5)
          SELECT CASE(trname(n))

          case('Clay','ClayIlli' ,'ClayKaol','ClaySmec','ClayCalc'
     &         ,'ClayQuar','ClayFeld','ClayHema','ClayGyps','ClayIlHe'
     &         ,'ClayKaHe','ClaySmHe','ClayCaHe','ClayQuHe','ClayFeHe'
     &         ,'ClayGyHe')
          n_bin = 1

          case('Silt1','Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema'
     &         ,'Sil1Gyps','Sil1Illi','Sil1Kaol','Sil1Smec','Sil1QuHe'
     &         ,'Sil1FeHe','Sil1CaHe','Sil1GyHe','Sil1IlHe','Sil1KaHe'
     &         ,'Sil1SmHe')
          n_bin = 2

          case('Silt2','Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema'
     &         ,'Sil2Gyps','Sil2Illi','Sil2Kaol','Sil2Smec','Sil2QuHe'
     &         ,'Sil2FeHe','Sil2CaHe','Sil2GyHe','Sil2IlHe','Sil2KaHe'
     &         ,'Sil2SmHe')
          n_bin = 3

          case('Silt3','Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema'
     &         ,'Sil3Gyps','Sil3Illi','Sil3Kaol','Sil3Smec','Sil3QuHe'
     &         ,'Sil3FeHe','Sil3CaHe','Sil3GyHe','Sil3IlHe','Sil3KaHe'
     &         ,'Sil3SmHe')
          n_bin = 4

          case('Silt4','Sil4Quar','Sil4Feld','Sil4Calc','Sil4Hema'
     &         ,'Sil4Gyps','Sil4Illi','Sil4Kaol','Sil4Smec','Sil4QuHe'
     &         ,'Sil4FeHe','Sil4CaHe','Sil4GyHe','Sil4IlHe','Sil4KaHe'
     &         ,'Sil4SmHe')
          n_bin = 5

          case('Silt5','Sil5Quar','Sil5Feld','Sil5Calc','Sil5Hema'
     &         ,'Sil5Gyps','Sil5Illi','Sil5Kaol','Sil5Smec','Sil5QuHe'
     &         ,'Sil5FeHe','Sil5CaHe','Sil5GyHe','Sil5IlHe','Sil5KaHe'
     &         ,'Sil5SmHe')
          n_bin = 6

          case default

          n_bin = 0

          END SELECT

          if ( n_bin > 0 ) then

            select case ( imDust )
            case ( 3 )
              zsum = sum( d_dust( : ) )
              if ( zsum > 0.d0 ) then
                frtrac = d_dust( n_bin ) / zsum * frtrac
              else
                frtrac = 0.d0
              end if
            case ( 5 )
              if ( d_dust( n_bin ) <= 0.d0 ) frtrac = 0.d0
            end select

          end if

        end if

c**** mineral fractions of dust aerosols
        if ( tracers_minerals .or. imDust == 4 .or. imDust == 5 ) frtrac
     &       = frtrac * mineralFractions( n - n_soildust + 1 )

#else /* TRACERS_DUST || TRACERS_MINERALS */
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)

        SELECT CASE (n)
        CASE (1)
          IF ( imDust == 0 .or. imDust >= 3 ) THEN
            frtrac = FracClayPDFscheme*frtrac
          ELSE IF (imDust == 2) THEN
            frtrac=FClWiCub*frclay*frtrac
          END IF
        CASE (2,3,4)
          IF ( imDust == 0 .or. imDust >= 3 ) THEN
            frtrac = FracSiltPDFscheme*frtrac
          ELSE IF (imDust == 2) THEN
            frtrac=FSiWiCub*frsilt*frtrac
          END IF
        END SELECT

        select case ( imDust )
        case ( 3 )
          zsum = sum( d_dust( : ) )
          if ( zsum > 0.d0 ) then
            frtrac = d_dust( n ) / zsum * frtrac
          else
            frtrac = 0.d0
          end if
        case ( 5 )
          if ( d_dust( n ) <= 0.d0 ) frtrac = 0.d0
        end select

        if ( imDust == 4 .or. imDust == 5 ) frtrac = mineralFractions(
     &       n ) * frtrac

#endif /* TRACERS_AMP || TRACERS_TOMAS */

#endif

c**** legacy dust emission scheme
        IF (imDust == 2) THEN
          dsrcflx=CWiCub*frtrac*(wsgcm-vtrsh)*wsgcm**2

c**** default case
c ..........
c dust emission above threshold from sub grid scale wind fluctuations
c ..........
        ELSE IF ( imDust == 0 .or. imDust >= 3 ) THEN
          dsrcflx = CWiPdf*frtrac*ers_data*dustSourceFunction*pdfint
c ..........
c emission according to cubic scheme, but with pdf scheme parameters
c (only used as diagnostic variable)
c ..........
          IF (soilvtrsh > 0. .AND. wsgcm > soilvtrsh) THEN
            dsrcflx2 = CWiPdf*frtrac*ers_data*dustSourceFunction
     &           *(wsgcm-soilvtrsh)*wsgcm**2
          END IF
        END IF

      ELSE IF (imDUST == 1) THEN
c**** prescribed (AeroCom) dust emission

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)

        SELECT CASE(trname(n))

        case('Clay','ClayIlli' ,'ClayKaol','ClaySmec','ClayCalc'
     &         ,'ClayQuar','ClayFeld','ClayHema','ClayGyps','ClayIlHe'
     &         ,'ClayKaHe','ClaySmHe','ClayCaHe','ClayQuHe','ClayFeHe'
     &         ,'ClayGyHe')
          n_bin = 1

        case('Silt1','Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema'
     &         ,'Sil1Gyps','Sil1Illi','Sil1Kaol','Sil1Smec','Sil1QuHe'
     &         ,'Sil1FeHe','Sil1CaHe','Sil1GyHe','Sil1IlHe','Sil1KaHe'
     &         ,'Sil1SmHe')
          n_bin = 2

        case('Silt2','Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema'
     &         ,'Sil2Gyps','Sil2Illi','Sil2Kaol','Sil2Smec','Sil2QuHe'
     &         ,'Sil2FeHe','Sil2CaHe','Sil2GyHe','Sil2IlHe','Sil2KaHe'
     &         ,'Sil2SmHe')
          n_bin = 3

        case('Silt3','Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema'
     &         ,'Sil3Gyps','Sil3Illi','Sil3Kaol','Sil3Smec','Sil3QuHe'
     &         ,'Sil3FeHe','Sil3CaHe','Sil3GyHe','Sil3IlHe','Sil3KaHe'
     &         ,'Sil3SmHe')
          n_bin = 4

        case default

          n_bin = 0

        END SELECT

        if ( n_bin > 0 ) dsrcflx = frtrac * d_dust( n_bin )

        if ( tracers_minerals ) dsrcflx = dsrcflx * mineralFractions( n
     &       - n_soildust + 1 )

#else /* TRACERS_DUST || TRACERS_MINERALS */

#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
        dsrcflx = frtrac * d_dust( n )
#endif

#endif

      END IF

      END IF

#endif
            
      return

      END SUBROUTINE local_dust_emission

c****
c**** This is legacy code for the old wet deposition scheme.
c**** To use the code the comiler directive TRACERS_WATER must
c**** not be defined and WET_DEPO_Ina must be defined in rundeck.
c****
      SUBROUTINE dust_wet(i,j)
!@sum  Simple scheme for wet deposition of dust/mineral tracers
!@auth Ina Tegen, Reha Cakmur, Jan Perlwitz

#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE constant,ONLY : Grav
      USE resolution,ONLY : Jm,Lm
      USE atm_com,ONLY : zatmo,gz
      USE fluxes,ONLY : prec
      USE clouds,ONLY : tm_dust,tmom_dust,trprc_dust
      USE tracer_com,ONLY : Ntm_dust
      USE trdust_mod,ONLY : prelay

      IMPLICIT NONE

      REAL*8,PARAMETER :: Z=700.

      INTEGER,INTENT(IN) :: i,j

      INTEGER :: l,layer,n
      INTEGER,DIMENSION(Jm) :: lwdep
      REAL*8,DIMENSION(Jm) :: h
      REAL*8 :: y
      REAL*8 :: height
      REAL*8,DIMENSION(Lm,Ntm_dust) :: tmold

#ifdef WET_DEPO_Ina
      SELECT CASE(j)
      CASE (1:6)
        lwdep(j)=3
        h(j)=2800
        lwdep(jm+1-j)=3
        h(jm+1-j)=2800
      CASE (7:12)
        lwdep(j)=4
        h(j)=4900
        lwdep(jm+1-j)=4
        h(jm+1-j)=4900
      CASE (13:16)
        LWDEP(J)=5
        h(j)=7400
        lwdep(jm+1-j)=5
        h(jm+1-j)=7400
      CASE (17:23)
        lwdep(j)=6
        h(j)=10300
        lwdep(jm+1-j)=6
        h(jm+1-j)=10300
      END SELECT
      
      y = Z*(prec(i,j)/h(j))
      IF (y > 1.) y=1.

#else

      layer=0
      DO l=Lm,1,-1
        IF (prelay(i,j,l) /= 0.) THEN 
          layer=l
          EXIT
        END IF
      END DO

      IF (layer == 1) THEN
        height=(gz(i,j,layer)-zatmo(i,j))/Grav
      ELSE IF (layer /= 0) THEN
        height=(gz(i,j,layer)-gz(i,j,1))/Grav
      END IF
      IF (layer /= 0) THEN
        y=Z*(prec(i,j)/height)
        IF (y > 1.) y=1.
      ELSE
        y=0.D0
      END IF

#endif

c**** Wet Deposition

      tmold=tm_dust
      DO n=1,Ntm_dust
        DO l=1,layer
          tm_dust(l,n)=tm_dust(l,n)*(1-y)
          tmom_dust(:,l,n)=tmom_dust(:,l,n)*(1-y)
        END DO
      END DO
      trprc_dust=tmold-tm_dust

#endif
#endif

      RETURN
      END SUBROUTINE dust_wet
