#include "rundeck_opts.h"

!@sum TRACER_PRT Tracer diagnostics.  These routines are generic for
!@+    all tracers
!@+   TRACEA and DIAGTCA are called throughout the day
!@+   The other routines are for printing
!@vers 2013/03/28
!@auth Jean Lerner (with ideas stolen from G.Schmidt, R. Ruedy, etc.)

#ifdef TRACERS_ON
      SUBROUTINE TRACEA
!@sum TRACEA accumulates tracer concentration diagnostics (IJL, JL)
!@auth J.Lerner
      USE CONSTANT, only : rgas
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only: itime
      USE ATM_COM, only: qcl,qci,t
      USE ATM_COM, only: pmid,pk
      USE DIAG_COM, only: jl_dpasrc,jl_dwasrc
      USE GEOM, only: imaxj,axyp,byaxyp
      USE SOMTQ_COM, only: mz
      use OldTracer_mod, only: itime_tr0, dowetdep
      USE TRACER_COM, only: ntm, trm
#ifdef TRACERS_WATER
      USE TRACER_COM, only: trwm
#endif
      USE TRDIAG_COM, only : taijln => taijln_loc, taijn  => taijn_loc,
     *     tij_mass, tij_conc, jlnt_conc, jlnt_mass, tajln => tajln_loc,
     &     taijls => taijls_loc,
     $     to_conc, ijlt_airmass
#ifdef SAVE_AEROSOL_3DMASS_FOR_NINT
     *     , taijls => taijls_loc, ijlt_3Dmass
#endif
#ifdef TRACERS_WATER
     *     ,jlnt_cldh2o
#endif
      USE ATM_COM, only: MA,byMA
      use oldtracer_mod, only: src_dist_index
      implicit none

      integer i,j,l,n
      real*8 tsum,asum

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C****
C**** Accumulate concentration for all tracers
C****

C**** save some basic model diags for weighting
      do l=1,lm
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            n = ijlt_airmass
            taijls(i,j,l,n) = taijls(i,j,l,n) + ma(l,i,j)
            call inc_ajl2(i,j,l,jl_dpasrc,axyp(i,j)*MA(l,i,j))
            call inc_ajl2(i,j,l,jl_dwasrc,axyp(i,j)*MA(l,i,j)*
     &        (qcl(i,j,l)+qci(i,j,l)))
          end do
        end do
      end do

      do 600 n=1,NTM
      IF (itime.lt.itime_tr0(n)) cycle
C**** Latitude-longitude by layer concentration
      if (to_conc(n).eq.1) then ! kg/m3
        do l=1,lm
          taijln(:,J_0:J_1,l,n) = taijln(:,J_0:J_1,l,n) + trm(:,J_0:J_1
     $          ,l,n)*byMA(l,:,J_0:J_1)*1d2*pmid(l,:,J_0:J_1)/(rgas*t(:
     $          ,J_0:J_1,L)*pk(L,:,J_0:J_1))
        end do
      else ! mixing ratio
        do l=1,lm
          taijln(:,J_0:J_1,l,n) = taijln(:,J_0:J_1,l,n) + trm(:,J_0:J_1
     $          ,l,n)*byMA(l,:,J_0:J_1)
        end do
!$OMP END PARALLEL DO
#ifdef SAVE_AEROSOL_3DMASS_FOR_NINT
!$OMP PARALLEL DO PRIVATE (L)
      if (ijlt_3Dmass(n).gt.0) then ! Ron: 3D mass distribution
         do l=1,lm
            taijls(:,J_0:J_1,l,ijlt_3Dmass(n)) =
     *      taijls(:,J_0:J_1,l,ijlt_3Dmass(n)) +
     *      trm   (:,J_0:J_1,l,            n )*byaxyp(:,J_0:J_1)
           end do
        endif
!$OMP END PARALLEL DO
# endif /* accumulate aerosol 3Dmass (Ron) */
      end if
C**** Average concentration; surface concentration; total mass
      if (src_dist_index(n)==0) then
        do j=J_0,J_1
        do i=I_0,I_1
          tsum = sum(trm(i,j,:,n))*byaxyp(i,j)  !sum over l
          asum = sum(MA(:,i,j))     !sum over l
          taijn(i,j,tij_mass,n) = taijn(i,j,tij_mass,n)+tsum  !MASS
          taijn(i,j,tij_conc,n) = taijn(i,j,tij_conc,n)+tsum/asum
        enddo; enddo
C**** Zonal mean concentration and mass
        do l=1,lm
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            call inc_tajln(i,j,l,jlnt_mass,n,trm(i,j,l,n))
          end do
        enddo; enddo

#ifdef TRACERS_WATER
C**** Zonal mean cloud water concentration
        if (dowetdep(n)) then
        do l=1,lm
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            call inc_tajln(i,j,l,jlnt_cldh2o,n,trwm(i,j,l,n))
          end do
        enddo; enddo
        end if
#endif
      endif ! src_dist_index(n)==0

  600 continue
      return
      end SUBROUTINE TRACEA
#endif

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      SUBROUTINE DIAGTCA (M,NT)
!@sum  DIAGTCA Keeps track of the conservation properties of tracers
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
      USE DIAG_COM, only : jm_budg,j_budg, j_0b, j_1b
      USE TRDIAG_COM, only: tconsrv=>tconsrv_loc,nofmt,title_tcon
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use oldtracer_mod, only: src_dist_index
      IMPLICIT NONE
!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT index denoting tracer number
      INTEGER, INTENT(IN) :: nt
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     *                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: total
      REAL*8, DIMENSION(JM_BUDG) :: TOTALJ
      INTEGER :: nm,ni
      INTEGER :: I, J, I_0, I_1, J_0, J_1

      if (src_dist_index(nt)/=0) return
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.
C**** NOFMT(1,NT) is the index for the instantaneous value.
      if (nofmt(m,nt).gt.0) then
C**** Calculate current value TOTAL
        call consrv_tr(nt,total)

        nm=nofmt(m,nt)
        ni=nofmt(1,nt)

C**** Calculate zonal sums
        totalj(j_0b:j_1b)=0.
        do j=j_0,j_1
        do i=i_0,i_1
          totalj(j_budg(i,j)) = totalj(j_budg(i,j)) + total(i,j)
        end do
        end do

c**** Accumulate difference from last time in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(J_0b:J_1b,nm,nt) =
     &         tconsrv(J_0b:J_1b,nm,nt)+(totalj(J_0b:J_1b)
     *         -tconsrv(J_0b:J_1b,ni,nt))
        end if
C**** Save current value in TCONSRV(NI)
        tconsrv(J_0b:J_1b,ni,nt)=totalj(J_0b:J_1b)
      end if
      return
      end subroutine diagtca

      subroutine consrv_tr(nt,total)
!@sum consrv_tr calculate total zonal tracer amount (kg)
!@auth Gavin Schmidt
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use resolution, only : ls1=>ls1_nominal
      use resolution, only : lm,jm,im
      use geom, only : imaxj
      use OldTracer_mod, only: trname
      use tracer_com, only : trm
#ifdef TRACERS_WATER
     *     ,trwm
#endif
      implicit none
      integer, intent(in) :: nt
!@var total = zonal total of tracer (kg)
      real*8, intent(out), dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &     GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: total
      integer :: i,j,l,ltop

      INTEGER :: I_0, I_1, J_1, J_0
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

      ltop=lm
#ifdef TRACERS_SPECIAL_Shindell
      if(trname(nt).eq.'Ox'.or.trname(nt).eq.'NOx') ltop=LS1-1
#endif

      total(:,:) = 0.

      do l=1,ltop
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            total(i,j) = total(i,j) + trm(i,j,l,nt)
#ifdef TRACERS_WATER
     *           +trwm(i,j,l,nt)
#endif
          end do
        end do
      end do

      IF (HAVE_SOUTH_POLE) total(2:im,1) = total(1,1)
      IF (HAVE_NORTH_POLE) total(2:im,jm)= total(1,jm)
      return
      end subroutine consrv_tr

      SUBROUTINE INC_DIAGTCB(I,J,DTRACER,M,NT)
!@sum  INC_DIAGTCB Keeps track of the conservation properties of tracers
!@+    This routine takes an already calculated difference
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
      USE DIAG_COM, only : j_budg
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE RESOLUTION, only: im,jm
      USE TRDIAG_COM, only: tconsrv=>tconsrv_loc,nofmt
      IMPLICIT NONE
      real*8, parameter :: fim=im
!@var I, J indices denoting gird box
      INTEGER, INTENT(IN) :: I,J
!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT index denoting tracer number
      INTEGER, INTENT(IN) :: nt
!@var DTRACER change of conserved quantity at this time
      REAL*8  :: DTRACER
      INTEGER :: nm

      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.

      nm=nofmt(m,nt)
      if (nm .gt.0) then
C**** Calculate latitudinal mean of change DTRACER
        IF (HAVE_SOUTH_POLE .AND. J.EQ.1) dtracer = fim*dtracer
        IF (HAVE_NORTH_POLE .AND. J.EQ.JM) dtracer = fim*dtracer

C**** Accumulate difference in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(J_BUDG(I,J),nm,nt)=tconsrv(J_BUDG(I,J),nm,nt)+dtracer
        end if
C**** No need to save current value
      end if
      return
      end subroutine inc_diagtcb

      SUBROUTINE DIAGTCB (DTRACER,M,NT)
!@sum  DIAGTCB Keeps track of the conservation properties of tracers
!@+    This routine takes an already calculated difference
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
      USE DIAG_COM, only : jm_budg,j_budg, j_0b, j_1b
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE RESOLUTION, only: jm,im
      USE TRDIAG_COM, only: tconsrv=>tconsrv_loc,nofmt,title_tcon
      IMPLICIT NONE

!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT index denoting tracer number
      INTEGER, INTENT(IN) :: nt
!@var DTRACER change of conserved quantity at this time
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: DTRACER
      REAL*8, DIMENSION(JM_BUDG) :: DTJ
      INTEGER :: j,nm

      INTEGER :: I, I_0, I_1, J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE,
     &               J_STRT = J_0, J_STOP = J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.

      if (nofmt(m,nt).gt.0) then
C**** Calculate latitudinal mean of change DTRACER
        IF (HAVE_SOUTH_POLE) dtracer(2:im,1) = dtracer(1,1)
        IF (HAVE_NORTH_POLE) dtracer(2:im,jm)= dtracer(1,jm)

C**** Calculate zonal sums
        DTJ(J_0B:J_1B)=0.
        DO J=J_0,J_1
          DO I=I_0,I_1
            DTJ(J_BUDG(I,J)) = DTJ(J_BUDG(I,J)) + DTRACER(I,J)
          END DO
        END DO

        nm=nofmt(m,nt)
C**** Accumulate difference in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(J_0b:J_1b,nm,nt)=tconsrv(J_0b:J_1b,nm,nt)
     *         +dtj(J_0b:J_1b)
        end if
C**** No need to save current value
      end if
      return
      end subroutine diagtcb
#endif /* TRACERS_ON or TRACERS_OCEAN */

#ifndef CUBED_SPHERE /* skip prt routines */

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      SUBROUTINE DIAGTCP
!@sum  DIAGCP produces tables of the conservation diagnostics
!@auth Gary Russell/Gavin Schmidt
!@ESMF This subroutine should only be called from a serial region.
!      It is NOT parallelized.

      USE CONSTANT, only: teeny, twopi
      use model_com, only: modelEclock
      USE MODEL_COM, only:
     &     idacc,jhour0,jdate0,amon,amon0,
     &     jyear0,nday,itime,itime0,xlabel,lrunid
      USE CONSTANT, only: areag
      use OldTracer_mod, only: itime_tr0
      USE TRACER_COM, only: NTM
      USE TRDIAG_COM, only:
     &     TCONSRV,ktcon,scale_tcon,title_tcon,nsum_tcon,ia_tcon,nofmt,
     &     lname_tconsrv,name_tconsrv,units_tconsrv,
     &     natmtrcons,nocntrcons
      USE DIAG_COM, only: inc=>incj,kdiag,qdiag
     &     ,jm=>jm_budg,dxyp_budg,lat_budg
      USE GC_COM, only : jeq
      USE DIAG_ZONAL, only : xwon
      USE MDIAG_COM, only: acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MAREA
      REAL*8, DIMENSION(KTCON) :: FGLOB
      REAL*8, DIMENSION(2,KTCON) :: FHEM
c GISS-ESMF EXCEPTIONAL CASE
c   JM+3 used for Lat, N-Hemi, S-Hemi, and Global Sums Storage
      REAL*8, DIMENSION(JM+3,KTCON) :: CNSLAT
      CHARACTER*4, PARAMETER :: HEMIS(2) = (/' SH ',' NH '/),
     *     DASH = ('----')
      INTEGER :: j,jhemi,jnh,jp1,jpm,jsh,jx,n,k,KTCON_max
      REAL*8 :: aglob,ahem,fnh,fsh,days
C**** Arrays needed for full output and pdE
      CHARACTER*38, DIMENSION(KTCON) :: TITLEO
      CHARACTER(len=lname_strlen), DIMENSION(KTCON) :: LNAMEO
      CHARACTER(len=sname_strlen), DIMENSION(KTCON) :: SNAMEO
      CHARACTER(len=units_strlen), DIMENSION(KTCON) :: UNITSO
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      if (kdiag(8).ge.2) return
C**** CALCULATE SCALING FACTORS
      IF (IDACC(12).LT.1) IDACC(12)=1
C**** Calculate areas
      DO J=1,JM
        MAREA(J)=1.D-10*XWON*DXYP_BUDG(J)+.5
      END DO
      AGLOB=1.D-10*AREAG*XWON
      AHEM=1.D-10*(.5*AREAG)*XWON
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF (QDIAG)
     *     call open_jc(trim(acc_period)//'.jct'//XLABEL(1:LRUNID),
     *     jm,lat_budg)
C****
C**** Outer loop over tracers
C****
      ktcon_max=0
      DO 900 N=1,natmtrcons+nocntrcons
C**** CALCULATE SUM OF CHANGES
C**** LOOP BACKWARDS SO THAT INITIALIZATION IS DONE BEFORE SUMMATION!
      DO J=1,JM
        DO K=KTCON,1,-1
          IF (NSUM_TCON(K,N).eq.0) THEN
            TCONSRV(J,K,N)=0.
          ELSEIF (NSUM_TCON(K,N).gt.0) THEN
            TCONSRV(J,NSUM_TCON(K,N),N)=
     *      TCONSRV(J,NSUM_TCON(K,N),N)+TCONSRV(J,K,N)
     *           *SCALE_TCON(K,N)*IDACC(12)/(IDACC(IA_TCON(K,N))+teeny)
            KTCON_max = NSUM_TCON(K,N)
          END IF
        END DO
      END DO

C**** CALCULATE ALL OTHER CONSERVED QUANTITIES ON TRACER GRID
      DO K=1,KTCON_max
        FGLOB(K)=0.
        FHEM(1,K)=0.
        FHEM(2,K)=0.

c GISS-ESMF EXCEPTIONAL CASE
c     GLobal/Hemi Sums
c     J-loop usage N-Hemi vs S-Hemi

        DO JSH=1,JEQ-1
          JNH=1+JM-JSH
          FSH=TCONSRV(JSH,K,N)*SCALE_TCON(K,N)/(IDACC(IA_TCON(K,N))
     *         +teeny)
          FNH=TCONSRV(JNH,K,N)*SCALE_TCON(K,N)/(IDACC(IA_TCON(K,N))
     *         +teeny)
          FGLOB (K)=FGLOB (K)+(FSH+FNH)
          FHEM(1,K)=FHEM(1,K)+FSH
          FHEM(2,K)=FHEM(2,K)+FNH
          CNSLAT(JSH,K)=FSH/(DXYP_BUDG(JSH))
          CNSLAT(JNH,K)=FNH/(DXYP_BUDG(JNH))
        END DO
        FGLOB (K)=FGLOB (K)/AREAG
        FHEM(1,K)=FHEM(1,K)/(.5*AREAG)
        FHEM(2,K)=FHEM(2,K)/(.5*AREAG)
c GISS-ESMF EXCEPTIONAL CASE
c     Storage of Hemi/Global Sums
        CNSLAT(JM+1,K)=FHEM(1,K)
        CNSLAT(JM+2,K)=FHEM(2,K)
        CNSLAT(JM+3,K)=FGLOB (K)
          titleo(k)=title_Tcon(k,n)
          lnameo(k)=lname_Tconsrv(k,n)
          snameo(k)=name_Tconsrv(k,n)
          unitso(k)=units_Tconsrv(k,n)
      END DO
C**** LOOP OVER HEMISPHERES
      DAYS=(Itime-Itime0)/FLOAT(nday)
      WRITE (6,'(''1'')')
      DO JHEMI=2,1,-1
        WRITE (6,'(''0'',A)') XLABEL
        if (n.le.natmtrcons) then
           WRITE (6,901) JYEAR0,AMON0,JDATE0,JHOUR0,
     *          YEAR,AMON,DATE,HOUR,ITIME,DAYS
        else
           WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *          YEAR,AMON,DATE,HOUR,ITIME,DAYS
        end if
        JP1=1+(JHEMI-1)*(JEQ-1)
        JPM=JHEMI*(JEQ-1)
C**** PRODUCE TABLES
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_BUDG(JX)),JX=JPM,JP1,-INC)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        DO K=1,KTCON_max
        WRITE (6,905) TITLE_TCON(K,N),FGLOB(K),FHEM(JHEMI,K),
     *         (NINT(MAX(-1d5,MIN(CNSLAT(JX,K),1d5))),JX=JPM,JP1,-INC)
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
      END DO
      IF (QDIAG) CALL POUT_JC(TITLEO,SNAMEO,LNAMEO,UNITSO,cnslat,
     *   ktcon_max)
  900 CONTINUE
      if(qdiag) call close_JC
      RETURN
C****
  901 FORMAT ('0Conservation Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  902 FORMAT ('0Ocean Consrv Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  903 FORMAT (1X,28('--'),13(A4,'--'))
  904 FORMAT (41X,'GLOBAL',A8,2X,13I6)
  905 FORMAT (A38,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10^10 M^2)',F30.1,F9.1,1X,13I6)
      END SUBROUTINE DIAGTCP
#endif /* TRACERS_ON or TRACERS_OCEAN */

#ifdef TRACERS_ON
      MODULE BDjlt
!@sum  stores info for outputting lat-sig/pressure diags for tracers
!@auth J. Lerner
      IMPLICIT NONE
      SAVE
!@param names of derived JLt/JLs tracer output fields
      INTEGER jlnt_nt_eddy,jlnt_vt_eddy

      END MODULE BDjlt


      SUBROUTINE JLt_TITLEX
!@sum JLt_TITLEX sets up titles, etc. for composite JL output (tracers)
!@auth J. Lerner
      use OldTracer_mod, only: trname, ntm_power
      USE TRACER_COM, only: ntm
      USE TRDIAG_COM
      USE BDjlt
      IMPLICIT NONE
      character*50 :: unit_string
      integer k,n,kpmax

!@var jlnt_xx Names for TAJL diagnostics
      do n=1,NTM
      k = ktajl+1
        jlnt_nt_eddy = k
        sname_jln(k,n) = 'tr_nt_eddy_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY EDDIES'
        jlq_power(k) = 10.
        jgrid_jlq(k) = 2
      k = k+1
        jlnt_vt_eddy = k
        sname_jln(k,n) = 'tr_vt_eddy_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY EDDIES'
        jlq_power(k) = 10.

      kpmax = k

      if (k .gt. ktajlx) then
        write (6,*) 'JLt_TITLEX: Increase ktajlx=',ktajlx,
     &         ' to at least ',k
        call stop_model('ktajlx too small',255)
      end if

C**** Construct UNITS string for output
      do k=ktajl+1,kpmax
      units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),' kg/s')
      end do

      end do

      END SUBROUTINE JLt_TITLEX


      SUBROUTINE DIAGJLT
!@sum DIAGJLT controls calls to the pressure/lat print routine for
!@+      tracers
!@auth J. Lerner
!@calls open_jl, JLt_TITLEX, JLMAP_t
C****
C**** THIS ROUTINE PRODUCES LATITUDE BY LAYER TABLES OF TRACERS
C****
!@ESMF This routine should only be called from a serial region.
!@     It is NOT parallelized.
      USE CONSTANT, only : undef,teeny
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE RESOLUTION, only : ls1=>ls1_nominal
      USE RESOLUTION, only : jm,lm
      USE MODEL_COM, only: itime,idacc,xlabel,lrunid
      USE DYNAMICS, only : dsig
      USE GEOM, only: bydxyp,dxyp,lat_dg
      use OldTracer_mod, only: ntm_power, dowetdep, trw0
#ifdef TRACERS_SPECIAL_Lerner
      use OldTracer_mod, only: trname
#endif
      USE TRACER_COM, only: ntm, n_water
#ifdef TRACERS_SPECIAL_O18
      use tracer_com, only: n_h2o18, n_hdo, n_h2o17
#endif
      USE DIAG_COM, only: linect,plm,qdiag,lm_req,ia_dga,ajl
     *     ,jl_dpa,jl_dpasrc,jl_dwasrc,fim
      USE MDIAG_COM, only: acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE ATM_COM, only : pednl00,pdsigl00
      USE TRDIAG_COM, only : PDSIGJL, tajln, tajls, lname_jln, sname_jln
     *     , units_jln,  scale_jln, lname_jls, sname_jls, units_jls,
     *     scale_jls, jls_power, jls_ltop, ia_jls, jwt_jls, jgrid_jls,
     *     jls_3Dsource, jlnt_conc, jlnt_mass, jlnt_nt_tot, jlnt_nt_mm,
     *     jlnt_lscond,  jlnt_turb,  jlnt_vt_tot, jlnt_vt_mm, jlnt_mc,
     *     jgrid_jlq, ia_jlq, scale_jlq, jlq_power, ktajls, jls_source
#ifdef TRACERS_WATER
     *     ,jlnt_cldh2o
#endif
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      USE TRDIAG_COM, only : to_per_mil
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRDIAG_COM, only : jls_H2Omr, jls_day
#endif
#ifdef TRACERS_SPECIAL_Lerner
      use TRACER_COM, only: n_CH4, n_O3
#endif
      USE BDJLT
      use oldtracer_mod, only: src_dist_index
      IMPLICIT NONE

      REAL*8, DIMENSION(1:JM)   :: ONESPO, BYAPO
      REAL*8, DIMENSION(JM+LM) :: ONES
      REAL*8, DIMENSION(JM,LM) :: A
      REAL*8, DIMENSION(LM) :: PM
      REAL*8 :: scalet
      INTEGER :: J,L,N,K,jtpow,n1,n2
      character(len=lname_strlen) :: lname
      character(len=sname_strlen) :: sname
      character(len=units_strlen) :: units

      REAL*8 :: dD, d18O, d17O

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_jl(trim(acc_period)//'.jlt'//XLABEL(1:LRUNID)
     *  ,jm,lm,0,lat_dg)
!?   *  ,jm,lm,lm_req,lat_dg)
C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      call JLt_TITLEX

      do l=1,lm  ! one level offset? why?
        pm(l)=pednl00(l+1)  !psfmpt*sige(l+1)+ptop
      end do
      onespo(1:JM)  = 1.d0
      onespo(1)  = fim
      onespo(jm) = fim
      ones(:) = 1.d0
      byapo(1:JM)=bydxyp(1:JM)*onespo(1:JM)/fim
!     do j=1,jm
!       ap=onespo(j)*APJ(j,1)/(fim*IDACC(ia_dga)+teeny)
!       call CALC_VERT_AMP (ap,lm,PL,MA,PDSIG,PEDN,PMID)
!       pdsigjl(j,l)=PDSIG(l)
!     end do
      do L=1,LS1-1
c        pdsigjl(1:JM,L)=dsig(l)*onespo(1:JM)*APJ(1:JM,1)/
c     *     (fim*IDACC(ia_dga)+teeny)
        do j=1,jm
          pdsigjl(J,L)=dsig(l)*onespo(J)*SUM(AJL(J,1:LS1-1,JL_DPA))/
     *         (IDACC(ia_dga)+teeny)
        end do
      end do
      do L=LS1,LM
        pdsigjl(1:JM,L)=pdsigl00(l)   ! psfmpt*dsig(l)
      end do

      linect = 65
C****
C**** LOOP OVER TRACERS
C****

C**** Note: why are jtpow defined here AND in units definition above?
C**** there is needless scope for inconsistency....

      DO 400 N=1,NTM

      if (src_dist_index(n)/=0) cycle

c     IF (itime.LT.itime_tr0(N)) cycle
C****
C**** TRACER CONCENTRATION
C****
      k = jlnt_conc

#ifdef TRACERS_WATER
      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
      scalet = 1.
      jtpow = 0.
      do l=1,lm
      do j=1,JM
        if (tajln(j,l,jlnt_mass,n_water).gt.0) then
          a(j,l)=1d3*(tajln(j,l,jlnt_mass,n)/(trw0(n)*tajln(j,l
     *         ,jlnt_mass,n_water))-1.)
        else
          a(j,l)=undef
        end if
      end do
      end do
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *     plm,a,scalet,ones,ones,lm,2,jgrid_jlq(k))
      else
#endif

      scalet = scale_jln(n)*scale_jlq(k)  !/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      do l=1,lm
      do j=1,JM
        if (ajl(j,l,jl_dpasrc).gt.0) then
          a(j,l)=tajln(j,l,jlnt_mass,n)/ajl(j,l,jl_dpasrc)
        else
          a(j,l)=undef
        end if
      end do
      end do
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *     plm,a,scalet,ones,ones,lm,2,jgrid_jlq(k))

#ifdef TRACERS_WATER
      end if
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AEROSOLS_SEASALT)
C****
C**** Mass diagnostic (this is saved for everyone, but only output
C**** for Dorothy and Drew for the time being)
C****
      k=jlnt_mass
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AEROSOLS_SEASALT)
      jtpow = ntm_power(n)+jlq_power(k)+13
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *     plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,jgrid_jlq(k))
#else
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *     plm,tajln(1,1,k,n),scalet,byapo,ones,lm,2,jgrid_jlq(k))
#endif
#endif

#ifdef TRACERS_WATER
C****
C**** TRACER CLOUD WATER CONCENTRATION
C****
      if (dowetdep(n)) then
      k = jlnt_cldh2o

      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
        scalet = 1.
        jtpow = 0.
        do l=1,lm
        do j=1,JM
        if (tajln(j,l,k,n_water).gt.0) then
          a(j,l)=1d3*(tajln(j,l,k,n)/(trw0(n)*tajln(j,l,k,n_water))-1.)
        else
          a(j,l)=undef
        end if
        end do
        end do
        CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *       plm,a,scalet,ones,ones,lm,2,jgrid_jlq(k))
      else
        scalet = scale_jlq(k)  !/idacc(ia_jlq(k))
        jtpow = ntm_power(n)+jlq_power(k)
        scalet = scalet*10.**(-jtpow)
        do l=1,lm
          do j=1,JM
            if (ajl(j,l,jl_dwasrc).gt.0) then
              a(j,l)=tajln(j,l,k,n)/ajl(j,l,jl_dwasrc)
            else
              a(j,l)=undef
            end if
          end do
        end do
        CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *       plm,a,scalet,ones,ones,lm,2,jgrid_jlq(k))
c        CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
c     *       plm,tajln(1,1,k,n),scalet,bydxyp,ones,lm,2,jgrid_jlq(k))
      end if
      end if
#endif
C****
C**** NORTHWARD TRANSPORTS: Total and eddies
C****
      a(:,:) = 0.
      k = jlnt_nt_tot
      do 205 l=1,lm
      do 205 j=1,JM-1
  205 a(j+1,l) = tajln(j,l,k,n)
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,a,scalet,ones,ones,lm,1,jgrid_jlq(k))
      do 210 l=1,lm
      do 210 j=1,JM-1
  210 a(j+1,l) = tajln(j,l,k,n) -tajln(j,l,jlnt_nt_mm,n)
      k = jlnt_nt_eddy
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,a,scalet,ones,ones,lm,1,jgrid_jlq(k))
C****
C**** VERTICAL TRANSPORTS: Total and eddies
C****
      k = jlnt_vt_tot
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  pm,tajln(1,1,k,n),scalet,ones,ones,lm-1,1,jgrid_jlq(k))
      a(:,:) = 0.
      do 260 l=1,lm-1
      do 260 j=1,JM
  260 a(j,l) = tajln(j,l,k,n)-tajln(j,l,jlnt_vt_mm,n)
      k = jlnt_vt_eddy
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  pm,a,scalet,ones,ones,lm-1,1,jgrid_jlq(k))
C****
C**** Convective Processes
C****
C**** MOIST CONVECTION
      k = jlnt_mc
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,onespo,ones,lm,1,Jgrid_jlq(k))
C**** LARGE-SCALE CONDENSATION
      k = jlnt_lscond
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,onespo,ones,lm,1,jgrid_jlq(k))
C**** TURBULENCE (or Dry Convection)
      k = jlnt_turb
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,onespo,ones,lm,1,jgrid_jlq(k))
  400 CONTINUE
C****
C**** JL Specials (incl. Sources and sinks)
C**** Partial move towards correct units (kg/(mb m^2 s)).
C**** Plot depends on jwt_jls.
C**** Note that only jwt_jls=3 is resolution independent.
C****
      do k=1,ktajls
        if (sname_jls(k).eq."daylight" .or. sname_jls(k).eq."H2O_mr"
     *       .or. lname_jls(k).eq."unused") cycle
        scalet = scale_jls(k)*10.**(-jls_power(k))/idacc(ia_jls(k))
        select case (jwt_jls(k))
        case (1)   !  simple sum (like kg/s),
          CALL JLMAP_t (lname_jls(k),sname_jls(k),units_jls(k),plm,
     *         tajls(1,1,k),scalet,onespo,ones,jls_ltop(k),jwt_jls(k)
     *         ,jgrid_jls(k))

        case (2)   !  area weighting (like kg/m^2 s)
          CALL JLMAP_t (lname_jls(k),sname_jls(k),units_jls(k),plm
     *         ,tajls(1,1,k),scalet,bydxyp,ones,jls_ltop(k),jwt_jls(k)
     *         ,jgrid_jls(k))

        case (3)   !  area + pressure weighting (like kg/mb m^2 s)
          CALL JLMAP_t (lname_jls(k),sname_jls(k),units_jls(k),plm
     *         ,tajls(1,1,k),scalet,byapo,ones,jls_ltop(k),jwt_jls(k)
     *         ,jgrid_jls(k))
        end select
        end do

#ifdef TRACERS_SPECIAL_Lerner
C**** some special combination diagnostics

C**** total chemical change for CH4
      if (n_CH4.gt.0) then
        k=jls_3Dsource(1,n_CH4)
        a(:,:) = tajls(:,:,jls_3Dsource(1,n_CH4))
     *         + tajls(:,:,jls_3Dsource(2,n_CH4))
        sname = 'Total_Chem_change'//trname(n_CH4)
        lname = 'TOTAL CHANGE OF CH4 BY CHEMISTRY'
        scalet = scale_jls(k)*10.**(-jls_power(k))/idacc(ia_jls(k))
        CALL JLMAP_t (lname,sname,units_jls(k),plm,a
     *     ,scalet,onespo,ones,jls_ltop(k),jwt_jls(k),jgrid_jls(k))
      end if
C**** total chemical change for O3
      if (n_O3.gt.0) then
        k=jls_3Dsource(1,n_O3)
        a(:,:) = tajls(:,:,jls_3Dsource(1,n_O3))
     *         + tajls(:,:,jls_3Dsource(2,n_O3))
     *         + tajls(:,:,jls_3Dsource(3,n_O3))
     *         + tajls(:,:,jls_source(1,n_O3))
        sname = 'Total_change_chem+depo'//trname(n_O3)
        lname = 'Total Change of O3 by Chemistry and deposition'
        scalet = scale_jls(k)*10.**(-jls_power(k))/idacc(ia_jls(k))
        CALL JLMAP_t (lname,sname,units_jls(k),plm,a
     *     ,scalet,onespo,ones,jls_ltop(k),jwt_jls(k),jgrid_jls(k))
      end if
#endif

#ifdef TRACERS_COSMO
C**** ratios : Be7/Pb210 and Be10/Be7
      if (n_Be7.gt.0 .and. n_Be10.gt.0 .and. n_Pb210.gt.0) then
        scalet = 1.
        jtpow = 0.
C*** ratio Be10/Be7
        k=jlnt_conc
        do l=1,lm
          do j=1,JM
            if (tajln(j,l,k,n_Be7).gt.0) then
              a(j,l)=tajln(j,l,k,n_Be10)/tajln(j,l,k,n_Be7)
            else
              a(j,l)=undef
            end if
          end do
        end do
        lname="Be10 to Be7 ratio"
        sname="be10be7"
        units=" "
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k))

C*** ratio Be7/Pb210
        k=jlnt_conc
        do l=1,lm
          do j=1,JM
            if (tajln(j,l,k,n_Pb210).gt.0) then
              a(j,l)=tajln(j,l,k,n_Be7)/tajln(j,l,k,n_Pb210)
C*** scale by (Be7decay/mm_Be7)/(Pb210decay/mm_Pb210) to convert to mBq
              a(j,l) = a(j,l)*trdecay(n_Be7)*tr_mm(n_Pb210)
     *             /trdecay(n_Pb210)/tr_mm(n_Be7)
            else
              a(j,l)=undef
            end if
          end do
        end do

        lname="Be7 to Pb210 ratio"
        sname="be7pb210"
        units="mBq/mBq"        !be sure this is 1/scalet
        scalet = 1.d0;
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k))
      end if

#endif

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-8*d18O)
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
        scalet = 1.
        jtpow = 0.

C**** Concentration in water vapour
        k=jlnt_mass
        do l=1,lm
          do j=1,jm
            if (tajln(j,l,k,n_water).gt.0) then
              d18O=1d3*(tajln(j,l,k,n1)/(trw0(n1)*tajln(j,l,k,n_water))
     *             -1.)
              dD=1d3*(tajln(j,l,k,n2)/(trw0(n2)*tajln(j,l,k,n_water))
     *             -1.)
              a(j,l)=dD-8.*d18O
            else
              a(j,l)=undef
            end if
          end do
        end do
        lname="Deuterium excess"
        sname="dexcess"
        units="per mil"
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k))

C**** Concentration in cloud water
        k=jlnt_cldh2o
        do l=1,lm
          do j=1,jm
            if (tajln(j,l,k,n_water).gt.0) then
              d18O=1d3*(tajln(j,l,k,n1)/(trw0(n1)*tajln(j,l,k,n_water))
     *             -1.)
              dD=1d3*(tajln(j,l,k,n2)/(trw0(n2)*tajln(j,l,k,n_water))
     *             -1.)
              a(j,l)=dD-8.*d18O
            else
              a(j,l)=undef
            end if
          end do
        end do
        lname="Deuterium excess in cloud water"
        sname="dexcess_cldh2o"
        units="per mil"
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k))

      end if

C****
C**** Calculations of 17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
        scalet = 1.
        jtpow = 0.

C**** Concentration in water vapour
        k=jlnt_conc
        do l=1,lm
          do j=1,jm
            if (tajln(j,l,k,n_water).gt.0) then
              d18O=tajln(j,l,k,n1)/(trw0(n1)*tajln(j,l,k,n_water))
              d17O=tajln(j,l,k,n2)/(trw0(n2)*tajln(j,l,k,n_water))
              a(j,l)=1d6*(log(d17O)-0.529d0*log(d18O))
            else
              a(j,l)=undef
            end if
          end do
        end do
        lname="D17O excess"
        sname="D17O_excess"
        units="per meg"
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k))

C**** Concentration in cloud water
        k=jlnt_cldh2o
        do l=1,lm
          do j=1,jm
            if (tajln(j,l,k,n_water).gt.0) then
              d18O=tajln(j,l,k,n1)/(trw0(n1)*tajln(j,l,k,n_water))
              d17O=tajln(j,l,k,n2)/(trw0(n2)*tajln(j,l,k,n_water))
              a(j,l)=1d6*(log(d17O)-0.529d0*log(d18O))
            else
              a(j,l)=undef
            end if
          end do
        end do
        lname="D17O excess in cloud water"
        sname="D17O_excess_cldh2o"
        units="per meg"
        CALL JLMAP_t (lname,sname,units,plm,a,scalet,ones,ones,lm,2
     *       ,jgrid_jlq(k))

      end if
#endif
      if (qdiag) call close_jl

      RETURN
      END SUBROUTINE DIAGJLT

      SUBROUTINE JLMAP_t (LNAME,SNAME,UNITS,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,JG)
C****
C**** THIS ROUTINE PRODUCES LAYER BY LATITUDE TABLES ON THE LINE
C**** PRINTER.  THE INTERIOR NUMBERS OF THE TABLE ARE CALCULATED AS
C****               AX * SCALET * SCALEJ * SCALEL.
C**** WHEN JWT=1, THE INSIDE NUMBERS ARE NOT AREA WEIGHTED AND THE
C****    HEMISPHERIC AND GLOBAL NUMBERS ARE SUMMATIONS.
C**** WHEN JWT=2, ALL NUMBERS ARE PER UNIT AREA.
C**** WHEN JWT=3, ALL NUMBERS ARE PER UNIT AREA AND PRESSURE, THE
C****    VERTICAL INTEGRAL GIVES TOTAL
C**** JG INDICATES PRIMARY OR SECONDARY GRID.
C**** THE BOTTOM LINE IS CALCULATED AS THE SUMMATION OF DSIG TIMES THE
C**** NUMBERS ABOVE
C****
      USE CONSTANT, only : undef
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, GLOBALSUM
      USE RESOLUTION, only : jm,lm
      use model_com, only: modelEclock
      USE MODEL_COM, only: jdate0,amon,amon0,jyear0,xlabel
      USE DYNAMICS, only : dsig,sige
      USE GEOM, only: wtj,jrange_hemi,lat_dg
      USE DIAG_COM, only: qdiag,inc=>incj,linect,jmby2,lm_req
      USE MDIAG_COM, only: acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE TRDIAG_COM, only : pdsigjl
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=units_strlen) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=lname_strlen) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=sname_strlen) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      INTEGER :: JG,JWT,LMAX
      REAL*8 :: SCALET
      REAL*8, DIMENSION(JM,LM) :: AX
      REAL*8, DIMENSION(JM) :: SCALEJ
      REAL*8, DIMENSION(LM) :: SCALEL,PL

      CHARACTER*4 :: DASH = '----',
     *     WORD(4) = (/ 'SUM ','MEAN','MEAN','.1* '/),
     *     BWORD(4)*6= (/ ' SUM  ',' MEAN ','INT/1K',' .1*  '/)

      INTEGER, DIMENSION(JM) :: MLAT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8, DIMENSION(2) :: FHEM,HSUM,PJSUM
      INTEGER :: J,JHEMI,L
      REAL*8 :: FGLOB,FLATJ,GSUM,SDSIG

!?    REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      REAL*8, DIMENSION(JM+3,LM+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
      integer :: year, date

      call modelEclock%get(year=year,date=date)

C form title string
      title = trim(lname)//' ('//trim(units)//')'

C****
C**** WRITE XLABEL ON THE TOP OF EACH OUTPUT PAGE
C****
      LINECT = LINECT+LMAX+7
      IF(LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,DATE,AMON,YEAR
      LINECT = LMAX+8
   20 CONTINUE
      WRITE (6,901) TITLE,(DASH,J=JG,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,JG)),J=JM,JG,-INC)
      WRITE (6,905) (DASH,J=JG,JM,INC)
C****
C**** CALCULATE TABLE NUMBERS AND WRITE THEM TO THE LINE PRINTER
C****
      XJL(:,:) = undef
      MLAT(:) = -1d5
      SDSIG = 1.-SIGE(LMAX+1)
      ASUM(:) = 0.
      HSUM(:) = 0.
      GSUM = 0.
      DO L=LMAX,1,-1
      FGLOB = 0.

c GISS-ESMF EXCEPTIONAL CASE
c   Hemisphere specific Loops and I/O issues, plus
c   N-Hemi, S-Hemi, and Global Sums

      DO JHEMI=1,2
        FHEM(JHEMI) = 0.
        PJSUM(JHEMI)= 0.
        DO J=JRANGE_HEMI(1,JHEMI,JG),JRANGE_HEMI(2,JHEMI,JG)
          if (AX(J,L).ne.undef) then
            FLATJ = AX(J,L)*SCALET*SCALEJ(J)*SCALEL(L)
            IF (JWT.eq.3) FLATJ=FLATJ/PDSIGJL(J,L)
            XJL(J,L) = FLATJ
            MLAT(J) = NINT(MAX(-1d5,MIN(FLATJ,1d5))) ! prevent integer overflow
            IF (JWT.EQ.1) THEN
              ASUM(J) = ASUM(J)+FLATJ !!!!most
              FHEM(JHEMI) = FHEM(JHEMI)+FLATJ*WTJ(J,JWT,JG)
            ELSEIF (JWT.EQ.2) THEN
              ASUM(J) = ASUM(J)+FLATJ*DSIG(L)/SDSIG !!!! concentration
              FHEM(JHEMI) = FHEM(JHEMI)+FLATJ*WTJ(J,JWT,JG)
            ELSEIF (JWT.EQ.3) THEN !! mass and area weighting
              ASUM(J) = ASUM(J)+FLATJ*PDSIGJL(J,L)
              FHEM(JHEMI)=FHEM(JHEMI)+FLATJ*WTJ(J,2,JG)*PDSIGJL(J,L)
              PJSUM(JHEMI)=PJSUM(JHEMI)+WTJ(J,2,JG)*PDSIGJL(J,L)
            ENDIF
          end if
        END DO                  ! loop over J
        FGLOB = FGLOB+FHEM(JHEMI)/REAL(MIN(JWT,2),KIND=8)
        IF (JWT.eq.3) FHEM(JHEMI) = FHEM(JHEMI) / PJSUM(JHEMI)
      END DO                    ! loop over hemisphere

      IF (JWT.EQ.1) THEN
        HSUM(1) = HSUM(1)+FHEM(1)
        HSUM(2) = HSUM(2)+FHEM(2)
        GSUM    = GSUM   +FGLOB
      ELSEIF (JWT.EQ.2) THEN
        HSUM(1) = HSUM(1)+FHEM(1)*DSIG(L)/SDSIG
        HSUM(2) = HSUM(2)+FHEM(2)*DSIG(L)/SDSIG
        GSUM    = GSUM   +FGLOB  *DSIG(L)/SDSIG
      ELSE
        HSUM(1) = HSUM(1)+FHEM(1)*PJSUM(1)
        HSUM(2) = HSUM(2)+FHEM(2)*PJSUM(2)
        GSUM    = GSUM   +FGLOB
        FGLOB=2.*FGLOB/(PJSUM(1)+PJSUM(2))
      ENDIF
C**** Output for each layer
         XJL(JM+3,L)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L)=FGLOB     ! GLOBAL
      WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),(MLAT(J),J=JM,JG,-INC)
      END DO                    ! loop over Layer
      WRITE (6,905) (DASH,J=JG,JM,INC)
      IF (JWT.eq.3) THEN ! scale integrated sums to look neater
        ASUM(:)= ASUM(:)*1d-3
        HSUM(:)= HSUM(:)*1d-3
        GSUM   = GSUM   *1d-3
      END IF
      ASUM(jmby2+1) = ASUM(jmby2+1)/JG
         DO J=JG,JM
           XJL(J   ,LM+1)=ASUM(J)
         END DO
         XJL(JM+3,LM+1)=HSUM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+1)=HSUM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+1)=GSUM      ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QDIAG) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS,
     *        JG,LMAX,XJL,PL,CLAT,CPRES)
      IF (LMAX.EQ.1) THEN
         LINECT = LINECT-1
         RETURN
      ENDIF
      IF (LMAX.EQ.1) RETURN
      WRITE (6,903) BWORD(JWT),GSUM,HSUM(2),HSUM(1),
     *    (NINT(MAX(-1d5,MIN(ASUM(J),1d5))),J=JM,JG,-INC)
      RETURN
C****
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (2X,A6,1X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)    ',A4,' G     NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END

      SUBROUTINE DIAGIJt
!@sum  DIAGIJt produces lat-lon fields as maplets (6/page) or full-page
!@+    digital maps, and binary (netcdf etc) files (if qdiag=true)
!@auth Jean Lerner (adapted from work of G. Russell,R. Ruedy)
!@ESMF This routine should only be called from a serial region.
!@     It is NOT parallelized.
      USE RESOLUTION, only : im,jm,lm
      use model_com, only: modelEclock
      USE MODEL_COM, only: jhour0,jdate0,amon,amon0
     *     ,jyear0,nday,itime,itime0,xlabel,lrunid,idacc
      use OldTracer_mod, only: dodrydep, dowetdep, trname, trw0,
     &   src_dist_index
      USE TRACER_COM, only: ntm, n_water
#ifdef TRACERS_SPECIAL_O18
      use tracer_com, only: n_h2o18, n_hdo, n_h2o17
#endif
      USE DIAG_COM
      USE DIAG_COM_RAD, only : ij_cldcv
      USE TRDIAG_COM, only : taijn, taijs, sname_tij, lname_tij,
     *     units_tij, scale_tij, tij_mass, lname_ijts,  sname_ijts,
     *     units_ijts,  scale_ijts,  ia_ijts, ktaij, ktaijs,
     *     tij_drydep, tij_gsdep, tij_surf, tij_grnd, tij_prec,
     *     tij_uflx, tij_vflx, ijts_HasArea, denom_ijts, ijts_clrsky,
     *     ijts_pocean, denom_tij, dname_tij
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
     &     ,to_per_mil
#endif
      USE DIAG_SERIAL, only : MAPTXT
      USE CONSTANT, only : teeny
      USE MDIAG_COM, only: acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

!TODO fix kludge
      integer, parameter :: MAXNTM = 1000 ! kludge - NTM is now dynamic

      integer, parameter :: ktmax = ktaij*MAXNTM+ktaijs
#ifdef TRACERS_SPECIAL_O18
     *     + 4                  ! include dexcess + D17O diags
#endif
#ifdef TRACERS_DRYDEP
     *     + MAXNTM                ! include dry dep % diags
#endif

!@var Qk: if Qk(k)=.true. field k still has to be processed
      logical, dimension (ktmax) :: Qk
!@var Iord: Index array, fields are processed in order Iord(k), k=1,2,..
!@+     only important for fields 1->nmaplets which appear in printout
!@+     Iord(k)=0 indicates that a blank space replaces a maplet
      INTEGER nmaplets
      INTEGER, DIMENSION(ktmax) :: ijtype,Iord,irange,iacc
      CHARACTER(len=lname_strlen), DIMENSION(ktmax) :: lname
      CHARACTER(len=sname_strlen), DIMENSION(ktmax) :: name
      CHARACTER(len=units_strlen), DIMENSION(ktmax) :: units
      REAL*8, DIMENSION(ktmax) :: scale
      REAL*8, DIMENSION(:,:,:), allocatable :: aij1,aij2
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      CHARACTER xlb*32,title*48
!@var LINE virtual half page (with room for overstrikes)
      CHARACTER*133 LINE(53)
      INTEGER ::  I,J,K,kx,kd,L,N,nd,kcolmn,nlines,jgrid,n1,n2
      REAL*8 :: DAYS,gm
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      if (kdiag(8).ge.1) return

      allocate(aij1(IM,JM,ktmax),aij2(IM,JM,ktmax))
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG)call open_ij(trim(acc_period)//'.ijt'//XLABEL(1:LRUNID)
     *     ,im,jm)

c fill in some denominators if needed
      if(ijts_clrsky.gt.0) then
        taijs(:,:,ijts_clrsky) = real(idacc(ia_rad))-aij(:,:,ij_cldcv)
      endif
      if(ijts_pocean.gt.0) then
        taijs(:,:,ijts_pocean) = aij(:,:,ij_pocean)
      endif

c**** always skip unused fields
      Qk = .true.
      k = 0

C**** Fill in the undefined pole box duplicates
      do i=2,im
        taijn(i,1,:,:) = taijn(1,1,:,:)
        taijn(i,jm,:,:) = taijn(1,jm,:,:)
        taijs (i,1,:) = taijs(1,1,:)
        taijs (i,jm,:) = taijs(1,jm,:)
      end do

C**** Fill in maplet indices for tracer sums/means and ground conc
      do n=1,NTM
      if (src_dist_index(n)/=0) cycle
      do kx=1,ktaij
        if (index(lname_tij(kx,n),'unused').gt.0) cycle
        k = k+1
        iord(k) = kx
        name(k) = sname_tij(kx,n)
        lname(k) = lname_tij(kx,n)
        units(k) = units_tij(kx,n)
        irange(k) = ir_log2
        scale(k) = scale_tij(kx,n)
        iacc(k) = ia_src
        ijtype(k) = 2
        aij1(:,:,k) = taijn(:,:,kx,n)
        aij2(:,:,k) = 1.
!@auth Kelley postprocessing decisions use predeclared metadata
        nd = denom_tij(kx,n)
        if(nd.gt.0) then
          ijtype(k)=3  ! ratio
          aij2(:,:,k) = taijn(:,:,kx,nd)
        endif
        if(dname_tij(kx,n).eq.'oicefr') then
          ijtype(k)=3  ! ratio
          scale(k) = scale(k)/idacc(iacc(k)) ! ijt_mapk workaround
          aij2(:,:,k)=aij(:,:,ij_rsoi)/(idacc(ia_ij(ij_rsoi))+teeny)
        endif
#ifdef TRACERS_WATER
        if (index(units(k),'er mil').gt.0) then
          aij1(:,:,k)=1d3*(aij1(:,:,k)/trw0(n)-taijn(:,:,kx,n_water))
          scale(k) = 1
        end if
#endif
      end do

#if (defined TRACERS_WATER) && (defined TRACERS_DRYDEP)
C****
C**** Calculation of dry deposition as percent
C****
      if (dodrydep(n).and.dowetdep(n)) then
        k=k+1
        ijtype(k) = 3
        name(k) = trim(trname(n))//"_pc_dry_dep"
        lname(k) = trim(trname(n))//" Percent Dry Deposition"
        units(k) = "%"
        irange(k) = ir_pct
        iacc(k) = ia_src
        iord(k) = 1
        aij1(:,:,k) = 100.*(taijn(:,:,tij_drydep,n)+
     *       taijn(:,:,tij_gsdep,n))
        aij2(:,:,k) = taijn(:,:,tij_prec,n)+taijn(:,:,tij_drydep,n)
     *       +taijn(:,:,tij_gsdep,n)
        scale(k) = 1.
      end if
#endif

      end do

C**** Fill in maplet indices for sources and sinks
      do kx=1,ktaijs
        if (index(lname_ijts(kx),'unused').gt.0) cycle
        k = k+1
        iord(k) = kx
        ijtype(k) = 1
        name(k) = sname_ijts(kx)
        lname(k) = lname_ijts(kx)
        units(k) = units_ijts(kx)
        irange(k) = ir_log2    ! should be the correct default
        iacc(k) = ia_ijts(kx)
        aij1(:,:,k) = taijs(:,:,kx)
        aij2(:,:,k) = 1.
        scale(k) = scale_ijts(kx)

        if (name(k)=='NO2_1030c' .or. name(k)=='NO2_1330c')then
          scale(k)=real(idacc(iacc(k)))+teeny
        endif

!@auth Kelley postprocessing decisions use predeclared metadata
        if(.not.ijts_HasArea(kx)) ijtype(k)=2 ! no need to divide by area
        kd = denom_ijts(kx)
        if(kd.gt.0) then
          ijtype(k)=3  ! ratio; set denominator aij2
          aij2(:,:,k) = taijs(:,:,kd)*
     &       ( real(idacc(ia_ijts(kx)),kind=8)/
     &         real(idacc(ia_ijts(kd)),kind=8) )
        endif

      end do

#ifdef TRACERS_COSMO
      if (n_Be7.gt.0 .and. n_Be10.gt.0 .and. n_Pb210.gt.0) then
C**** Be10/Be7
        k=k+1
        ijtype(k) = 3 !perform a ratio
        name(k) = "be10be7_ij"
        lname(k) = "surface ratio Be10 to Be7"
        units(k) = " "
        irange(k) = ir_0_18
        iacc(k) = ia_srf
        iord(k) = 1
        aij1(:,:,k) = taijn(:,:,tij_surf,n_Be10) !set as numerator
        aij2(:,:,k) = taijn(:,:,tij_surf,n_Be7) !set as denom
        scale(k) = 1.
C**** Be7/Pb210
        k=k+1
        ijtype(k) = 3 !perform a ratio
        name(k) = "be7pb210_ij"
        lname(k) = "surface ratio Be7 to Pb210"
        units(k) = "mBq/mBq "
        irange(k) = ir_0_180
        iacc(k) = ia_srf
        iord(k) = 2
        scale(k) = 1.d0
        aij1(:,:,k) = taijn(:,:,tij_surf,n_Be7) !numerator
C*** scale by (Be7decay/mm_Be7)/(Pb210decay/mm_Pb210) to convert to mBq
        aij1(:,:,k)=aij1(:,:,k)*trdecay(n_Be7)*tr_mm(n_Pb210)
     *             /trdecay(n_Pb210)/tr_mm(n_Be7)
        aij1(:,:,k)=scale(k)*aij1(:,:,k)
        aij2(:,:,k) = taijn(:,:,tij_surf,n_Pb210) !denominator
      end if
#endif

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-d18O)
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** precipitation
        k=k+1
        ijtype(k) = 3
        name(k) = "prec_ij_dex"
        lname(k) = "Deuterium excess in precip"
        units(k) = "per mil"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 1
        aij1(:,:,k) = 1d3*(taijn(:,:,tij_prec,n2)/trw0(n2)-
     *                8.*taijn(:,:,tij_prec,n1)/trw0(n1)+
     *                7.*taijn(:,:,tij_prec,n_water))
        aij2(:,:,k) = taijn(:,:,tij_prec,n_water)
        scale(k) = 1.
C**** ground concentration
        k=k+1
        ijtype(k) = 3
        name(k) = "grnd_ij_dex"
        lname(k) = "Deuterium excess at Ground"
        units(k) = "per mil"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 2
        aij1(:,:,k) = 1d3*(taijn(:,:,tij_grnd,n2)/trw0(n2)-
     *                8.*taijn(:,:,tij_grnd,n1)/trw0(n1)+
     *                7.*taijn(:,:,tij_grnd,n_water))
        aij2(:,:,k) = taijn(:,:,tij_grnd,n_water)
        scale(k) = 1.
      end if

C****
C**** Calculations of D17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
C**** precipitation
        k=k+1
        ijtype(k) = 3
        name(k) = "prec_ij_D17O"
        lname(k) = "D17O excess in precip"
        units(k) = "per meg"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 1
        do j=1,jm
          do i=1,im
            if (taijn(i,j,tij_prec,n_water).gt.0) then
              aij1(i,j,k) = 1d6*taijn(i,j,tij_prec,n_water)*
     *             (log(taijn(i,j,tij_prec,n2)/trw0(n2))-
     *             0.529d0*log(taijn(i,j,tij_prec,n1)/trw0(n1))-
     *             0.471d0*log(taijn(i,j,tij_prec,n_water)))
              aij2(i,j,k) = taijn(i,j,tij_prec,n_water)
            else
              aij1(i,j,k)=0.
              aij2(i,j,k)=0.
            end if
          end do
        end do
        scale(k) = 1.
C**** ground concentration
        k=k+1
        ijtype(k) = 3
        name(k) = "grnd_ij_D17O"
        lname(k) = "D17O excess at Ground"
        units(k) = "per meg"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 2
        do j=1,jm
          do i=1,im
            if (taijn(i,j,tij_grnd,n_water).gt.0) then
              aij1(i,j,k) = 1d6*taijn(i,j,tij_grnd,n_water)*
     *             (log(taijn(i,j,tij_grnd,n2)/trw0(n2))-
     *             0.529d0*log(taijn(i,j,tij_grnd,n1)/trw0(n1))-
     *             0.471d0*log(taijn(i,j,tij_grnd,n_water)))
              aij2(i,j,k) = taijn(i,j,tij_grnd,n_water)
            else
              aij1(i,j,k)=0.
              aij2(i,j,k)=0.
            end if
          end do
        end do
        scale(k) = 1.
      end if
#endif

      nmaplets = k
      Qk(k+1:ktmax)=.false.

      xlb=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
C****
      DAYS=(Itime-Itime0)/FLOAT(nday)
C**** Collect the appropriate weight-arrays in WT_IJ
      wt_ij(:,:,:) = 1.

C**** Print out 6-map pages
      do n=1,nmaplets
        if (mod(n-1,6) .eq. 0) then
c**** print header lines
          write (6,'("1",a)') xlabel
          write (6,902) jyear0,amon0,jdate0,jhour0,
     *      year,amon,date,hour,itime,days
        end if
        kcolmn = 1 + mod(n-1,3)
        if (kcolmn .eq. 1) line=' '
c**** Find, then display the appropriate array
        if (Iord(n) .gt. 0 .and. Qk(n)) then
          call ijt_mapk (ijtype(n),aij1(1,1,n),aij2(1,1,n),smap,smapj,gm
     *         ,jgrid,scale(n),iacc(n),irange(n),name(n),lname(n)
     *         ,units(n))
          title=trim(lname(n))//' ('//trim(units(n))//')'
          call maptxt(smap,smapj,gm,irange(n),title,line,kcolmn,nlines)
c assuming igrid=jgrid for now
          if(qdiag) call pout_ij(title//xlb,name(n),lname(n),units(n),
     *                            smap,smapj,gm,jgrid,jgrid)
          Qk(n) = .false.
        end if

c**** copy virtual half-page to paper if appropriate
        if (kcolmn.eq.3 .or. n.eq.nmaplets) then
          do k=1,nlines
            write (6,'(a133)') line(k)
          end do
        end if
      end do

      if (.not.qdiag) then
        deallocate(aij1,aij2)
        RETURN
      endif

C**** produce binary files of remaining fields if appropriate
      do n=1,ktmax
        if (Qk(n)) then
          call ijt_mapk (ijtype(n),aij1(1,1,n),aij2(1,1,n),smap,smapj,gm
     *         ,jgrid,scale(n),iacc(n),irange(n),name(n),lname(n)
     *         ,units(n))
          title=trim(lname(n))//' ('//trim(units(n))//')'
          call pout_ij(title//xlb,name(n),lname(n),units(n),smap,smapj,
     *         gm,jgrid,jgrid) ! assuming igrid=jgrid for now
        end if
      end do
      if(qdiag) call close_ij

      deallocate(aij1,aij2)
      RETURN
C****
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      END SUBROUTINE DIAGIJt


      SUBROUTINE DIAGIJLt
!@sum  DIAGIJLt produces 3D lat-lon fields as maplets (6/page) or full-page
!@+    digital maps, and binary (netcdf etc) files (if qdiag=true)
!@auth Jean Lerner (adapted from work of G. Russell,R. Ruedy)
!@ESMF This routine should only be called from a serial region.
!@     It is NOT parallelized.
      USE RESOLUTION, only : im,jm,lm
      use model_com, only: modelEclock
      USE MODEL_COM, only: jhour0,jdate0,amon,amon0
     *     ,jyear0,nday,itime,itime0,xlabel,lrunid,idacc
      use OldTracer_mod, only: trw0, src_dist_index
      USE TRACER_COM, only: ntm, n_water
#ifdef TRACERS_SPECIAL_O18
      use tracer_com, only: n_h2o18, n_hdo, n_h2o17
#endif

      USE TRDIAG_COM, only : taijln, taijls, sname_ijlt, lname_ijlt,
     *     units_ijlt, sname_ijt, lname_ijt, units_ijt, scale_ijt,
     *     ir_ijlt, ia_ijlt, scale_ijlt, ktaijl
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
     &     ,to_per_mil
#endif
      USE DIAG_SERIAL, only : MAPTXT, scale_ijlmap
      USE CONSTANT, only : teeny
      USE DIAG_COM, only : ijkgridc,ctr_ml,ir_m45_130,
     *                     kdiag,ir_log2,ia_src,wt_ij,
     *                     qdiag
      USE MDIAG_COM, only: acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE ATM_COM, only : pmidl00
      use filemanager
      IMPLICIT NONE

      CHARACTER XLB*24,titlel(lm)*80
      CHARACTER*11 CPRESS(LM)
      CHARACTER*3 CLEV(LM)
      INTEGER i,j,l,kxlb,k

!TODO fix kludge - NTM replaced with MAXNTM to allow ktmax to be a parameter
      integer, parameter :: MAXNTM=1000
      integer, parameter :: ktmax = maxntm+ktaijl
#ifdef TRACERS_SPECIAL_O18
     *     + 2        ! include dexcess + D17O diags
#endif

!@var Qk: if Qk(k)=.true. field k still has to be processed
      logical, dimension (ktmax) :: Qk
!@var Iord: Index array, fields are processed in order Iord(k), k=1,2,..
!@+     only important for fields 1->nmaplets which appear in printout
!@+     Iord(k)=0 indicates that a blank space replaces a maplet
      INTEGER nmaplets
      INTEGER, DIMENSION(ktmax) :: ijtype,Iord,irange,iacc
      CHARACTER(len=lname_strlen), DIMENSION(ktmax) :: lname
      CHARACTER(len=sname_strlen), DIMENSION(ktmax) :: name
      CHARACTER(len=units_strlen), DIMENSION(ktmax) :: units
      integer, dimension(ktmax) :: lgrid
      REAL*8, DIMENSION(ktmax) :: scale
      REAL*8, DIMENSION(:,:,:,:), allocatable :: aijl1,aijl2
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(IM,JM,LM) :: SMAPIJL
      REAL*8, DIMENSION(JM) :: SMAPJ
      REAL*8, DIMENSION(JM,LM) :: SMAPJL
      REAL*8, DIMENSION(LM) :: SMAPL
!@var LINE virtual half page (with room for overstrikes)
      CHARACTER*133 LINE(53)
      INTEGER ::  kx,N,kcolmn,nlines,jgrid,n1,n2,nn
      REAL*8 :: DAYS,gm
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      if (kdiag(8).ge.1) return

      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
      DAYS=(Itime-Itime0)/FLOAT(nday)
      DO L=1,LM
        WRITE(CPRESS(L),'(F8.3,A3)') pmidl00(l),' mb'
        WRITE(CLEV(L),'(I3)') l
      END DO

      allocate(aijl1(IM,JM,LM,ktmax),aijl2(IM,JM,LM,ktmax))

C**** Initialise
      do k=1,ktmax
        lgrid(k) = ctr_ml
        lname(k) = 'no output'
      enddo
      Qk = .true.
      k = 0

C**** Fill in the undefined pole box duplicates
      do i=2,im
        taijln(i,1,:,:) = taijln(1,1,:,:)
        taijln(i,jm,:,:) = taijln(1,jm,:,:)
        taijls(i,1,:,:) = taijls(1,1,:,:)
        taijls(i,jm,:,:) = taijls(1,jm,:,:)
      end do

C**** Fill in maplet indices for tracer concentrations
      do n=1,NTM
        if (src_dist_index(n)/=0) cycle
        k = k+1
        iord(k) = n
        ijtype(k) = 1
        name(k) = sname_ijt(n)
        lname(k) = lname_ijt(n)
        units(k) = units_ijt(n)
        irange(k) = ir_log2
        iacc(k) = ia_src
        scale(k) = scale_ijt(n)
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0) then
           ijtype(k) = 3
           scale(k) = 1.
        end if
#endif
        do l=1,lm
          aijl1(:,:,l,k) = taijln(:,:,l,n)
          aijl2(:,:,l,k) = 1.
#ifdef TRACERS_WATER
          if (to_per_mil(n).gt.0) then
            aijl1(:,:,l,k)=1d3*(taijln(:,:,l,n)-taijln(:,:,l,n_water)
     *           *trw0(n))
            aijl2(:,:,l,k)=taijln(:,:,l,n_water)*trw0(n)
          end if
#endif
        end do
      end do

C**** Fill in maplet indices for 3D tracer specials
      do kx=1,ktaijl
        if (index(lname_ijlt(kx),'unused').gt.0) cycle
        k = k+1
        iord(k) = kx
        ijtype(k) = 2
        name(k) = sname_ijlt(kx)
        lname(k) = lname_ijlt(kx)
        units(k) = units_ijlt(kx)
        irange(k) = ir_ijlt(kx)
        iacc(k) = ia_ijlt(kx)
        scale(k) = scale_ijlt(kx)
        do l=1,lm
          aijl1(:,:,l,k) = taijls(:,:,l,kx)
          aijl2(:,:,l,k) = 1.
        end do
      end do

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-d18O)
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** water vapour
        k=k+1
        ijtype(k) = 3
        lname(k) = "Deuterium excess water vapour"
        name(k) = "wvap_ij_dex"
        units(k) = "per mil"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 1  !???l+2
        scale(k) = 1.
        do l=1,lm
          aijl1(:,:,l,k) = 1d3*(taijln(:,:,l,n2)/trw0(n2)-
     *         8.*taijln(:,:,l,n1)/trw0(n1)+
     *         7.*taijln(:,:,l,n_water))
          aijl2(:,:,l,k) = taijln(:,:,l,n_water)
        end do
      end if

C****
C**** Calculations of D17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
C**** water vapour
        k=k+1
        ijtype(k) = 3
        lname(k) = "D17O excess water vapour"
        name(k) = "wvap_ij_D17O"
        units(k) = "per meg"
        irange(k) = ir_m45_130
        iacc(k) = ia_src
        iord(k) = 2   !??? l+2
        scale(k) = 1.
        do l=1,lm
          do j=1,jm
            do i=1,im
              if (taijln(i,j,l,n_water).gt.0) then
                aijl1(i,j,l,k) = 1d6*taijln(i,j,l,n_water)*
     *               (log(taijln(i,j,l,n2)/trw0(n2))-
     *               0.529d0*log(taijln(i,j,l,n1)/trw0(n1))-
     *               0.471d0*log(taijln(i,j,l,n_water)))
                aijl2(i,j,l,k) = taijln(i,j,l,n_water)
              else
                aijl1(i,j,l,k)=0.
                aijl2(i,j,l,k)=0.
              end if
            end do
          end do
        end do
      end if
#endif

      if (k.gt.ktmax) print*,"K < KTMAX in diagijlt",k,ktmax
      nmaplets = k
      Qk(k+1:ktmax)=.false.
C**** Collect the appropriate weight-arrays in WT_IJ
      wt_ij(:,:,:) = 1.

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      if (qdiag) call open_ijl(trim(acc_period)//'.ijlt'/
     *     /XLABEL(1:LRUNID),im,jm,lm,
     &     nmaplets,name,lname,units,lgrid)

C**** Print out 6-map pages
      do n=1,nmaplets*lm
        nn=1+(n-1)/lm
        if (mod(n-1,6) .eq. 0) then
c**** print header lines
          write (6,'("1",a)') xlabel
          write (6,902) jyear0,amon0,jdate0,jhour0,
     *      year,amon,date,hour,itime,days
        end if
        kcolmn = 1 + mod(n-1,3)
        if (kcolmn .eq. 1) line=' '
c**** Find, then display the appropriate array
        if (Iord(nn) .gt. 0 .and. Qk(nn)) then
          l=1+mod(n-1,lm)
          call ijt_mapk (ijtype(nn),aijl1(1,1,l,nn),aijl2(1,1,l,nn),smap
     *         ,smapj,gm,jgrid,scale(nn),iacc(nn),irange(nn),name(nn)
     *         ,lname(nn),units(nn))
          titlel(l)=trim(lname(nn))//' Level '//clev(l)//' ('/
     *         /trim(units(nn))//')'
          call maptxt(smap,smapj,gm,irange(nn),titlel(l),line,kcolmn
     *         ,nlines)
        end if
c**** copy virtual half-page to paper if appropriate
        if (kcolmn.eq.3 .or. n.eq.nmaplets*lm) then
          do k=1,nlines
            write (6,'(a133)') line(k)
          end do
        end if

C**** for every diag, output all levels at once
        if (l.eq.lm) then
           call scale_ijlmap(ijtype(nn),aijl1(1,1,1,nn),aijl2(1,1,1,nn),
     *       scale(nn),1,idacc(iacc(nn)),idacc(iacc(nn)),smapijl,smapjl)
          smapl=0               ! tmp
          if (qdiag) call pout_ijl(titlel,name(nn),lname(nn)
     *         ,units(nn),smapijl,smapjl,smapl,ijkgridc)
        end if
      end do

      if(qdiag) call close_ijl

      deallocate(aijl1,aijl2)
      RETURN
C****
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      END SUBROUTINE DIAGIJLt

      subroutine IJt_MAPk(nmap,aij1,aij2,smap,smapj,gm,jgrid
     *     ,scale,iacc,irange,name,lname,units)
!@sum ijt_MAPk returns the map data and related terms for the k-th field
!@+   for tracers and tracer sources/sinks
      USE CONSTANT, only: teeny
      USE RESOLUTION, only : im,jm
      USE MODEL_COM, only:idacc
      USE GEOM, only: dxyp
      USE TRACER_COM, only: ntm
      USE DIAG_COM
      USE DIAG_SERIAL, only : IJ_avg
      USE MDIAG_COM, only:
     &     sname_strlen,units_strlen,lname_strlen
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap,aij1,aij2
      REAL*8, DIMENSION(JM) :: smapj
      integer i,j,iwt,jgrid,irange,nmap,iacc
      character(len=sname_strlen) :: name
      character(len=units_strlen) :: units
      character(len=lname_strlen) :: lname
      real*8 :: gm,nh,sh, byiacc,scale
!@var isumz,isumg = 1 or 2 if zon,glob sums or means are appropriate
      integer isumz,isumg,k1

C****
C**** Extract useful local domain parameters from "grid"
C****
      isumz = 2 ; isumg = 2  !  default: in most cases MEANS are needed
      iwt = 1 ; jgrid = 1

      adenom = 1.
      byiacc = 1./(idacc(iacc)+teeny)
c**** tracer amounts (divide by area) and Sources and sinks
      if (nmap.eq.1) then
        do j=1,JM
        do i=1,im
          anum(i,j)=aij1(i,j)*byiacc*scale/dxyp(j)
        end do
        end do
c**** tracer sums and means (no division by area)
      else if (nmap.eq.2) then
        do j=1,JM
        do i=1,im
          anum(i,j)=aij1(i,j)*byiacc*scale
        end do
        end do
c**** ratios (i.e. per mil diags)
      else if (nmap.eq.3) then
        do j=1,JM
        do i=1,im
          anum(i,j)=aij1(i,j)*scale
          adenom(i,j)=aij2(i,j)
        end do
        end do
C**** PROBLEM
      else  ! should not happen
        write (6,*) 'no ij map type defined for ijt_index',name
        call stop_model(
     &       'ijt_mapk: undefined extra ij_field for tracers',255)
      end if

c**** Find final field and zonal, global, and hemispheric means
  100 continue
      call ij_avg (anum,adenom,wt_ij(1,1,iwt),jgrid,isumz,isumg, ! in
     *             smap,smapj,gm,nh,sh)                    ! out
      return
      end subroutine ijt_mapk

#endif /* cubed sphere skipping prt routines */

#else  /* TRACERS_ON */

      subroutine tracers_are_off
      return
      end subroutine tracers_are_off

#endif /* TRACERS_ON */
