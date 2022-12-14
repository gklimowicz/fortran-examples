      program convert_rsf_files
C**** convert rsf files from model II' (B399) format to modelE format
C**** compile with: gmake conv_rsf.o
C**** f90 -o conv_rsf conv_rsf.o *.o -O2 -64 -mips4 -static \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE CONSTANT, only : lhm,shi,by3
      USE MODEL_COM, only : im,jm,lm,wm,u,v,t,p,q,xlabel,iowrite
     *     ,iowrite_mon,focean,nday,itime,itimei,itimee,itime0,iyear1
      USE SOMTQ_COM
      USE GHY_COM, only : snowe,tearth,wearth,aiearth,snoage,wbare,wvege
     *     ,htbare,htvege,snowbv,ngm,evap_max_ij,fr_sat_ij,qg_ij
      use veg_com, only : cint,qfol,cnc_ij
      USE RAD_COM, only : rqt,lm_req
      USE CLOUDS_COM, only : ttold,qtold,svlhx,rhsav,cldsav
      USE DIAG_COM, only : keynr,tsfrez
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
     *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,
     *     ustar_pbl,egcm,w2gcm,tgvavg,qgavg
      USE STATIC_OCEAN, only : tocean,z1o,sss0
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi,pond_melt,flag_dsws
      USE SEAICE, only : ace1i,xsi,ac2oim,ssi0,tfrez
      USE LANDICE_COM, only : tlandi,snowli
      USE LAKES_COM, only : flake
      USE FILEMANAGER
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60              ,clabel*156
      INTEGER IARGC,iu_AIC,I,J,L,N,ioerr,iu_TOPO   ,jc(100)
      REAL*8 TAUX,X                                ,rc(161)
      REAL*8 MSI1,tfo
      INTEGER ItimeX
!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5
      integer :: conv_index(im,jm,2)
      character*120 cindex_file, gic_file, argv
      integer iu_GIC, iu_CI, i0, j0
      logical :: extend_gh = .false., dump_gic = .false.

      IF (IARGC().lt.2) THEN
        PRINT*,"Convert rsf files from old format to new"
        PRINT*,"To create a restart file for the standard land mask:"
        PRINT*,"conv_rsf in_file rsf_file [-g gic_file] [-e ind_file]"
        PRINT*,"To create a restart file with the extended land"
        PRINT*,"surface data (which can be used with any land mask):"
        PRINT*,"use '-e conv_indices_file' option"
        PRINT*,"To create a GIC file use '-g gic_file' option"
        STOP
      END IF

      CALL GETARG(1,infile)
      CALL GETARG(2,outfile)

      i = 3
      do while (i <= IARGC())
        call getarg(i,argv); i = i+1
        if ( argv == "-e" ) then
          extend_gh = .true.
          call getarg(i,cindex_file); i = i+1
          cycle
        endif
        if ( argv == "-g" ) then
          dump_gic = .true.
          call getarg(i,gic_file); i = i+1
          cycle
        endif
      enddo


      if ( extend_gh ) then
        call openunit (trim(cindex_file), iu_CI, .true. , .true. )
        read(iu_CI) conv_index
        call closeunit (iu_CI)
      endif

      call openunit (trim(infile), iu_AIC, .true. , .true. )

      READ (iu_AIC,ERR=800,END=810) TAUX,JC,CLABEL,RC,KEYNR,
     *     U,V,T,P,Q,
     2     ((TOCEAN(1,I,J),I=1,IM),J=1,JM),RSI,MSI,
     *     (((TOCEAN(L,I,J),I=1,IM),J=1,JM),L=2,3),
     *     SNOWI,SNOWE,
     *     ((HSI(1,I,J),I=1,IM),J=1,JM),TEARTH,WEARTH,AIEARTH,
     *     ((HSI(2,I,J),I=1,IM),J=1,JM),((X,I=1,IM),J=1,JM),
     *     (((SNOAGE(L,I,J),I=1,IM),J=1,JM),L=1,3),SNOWLI,
     *     (((TLANDI(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     (((HSI(L,I,J),I=1,IM),J=1,JM),L=3,4),
     *     (((wbare(L,I,J),I=1,IM),J=1,JM),L=1,NGM),
     *     (((wvege(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((htbare(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((htvege(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((snowbv(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg,  ustar_pbl
     *           ,uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,(((TTOLD(L,I,J)
     *     ,I=1,IM),J=1,JM),L=1,LM),(((QTOLD(L,I,J),I=1,IM),J=1,JM),L=1
     *     ,LM),(((SVLHX(L,I,J),I=1,IM),J=1,JM),L=1,LM),(((RHSAV(L,I,J)
     *     ,I=1,IM),J=1,JM),L=1,LM),WM,(((CLDSAV(L,I,J),I=1,IM),J=1,JM)
     *     ,L=1,LM),((((TMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9)
     *     ,((((QMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9),(((RQT(L,I
     *     ,J),I=1,IM),J=1,JM),L=1,LM_REQ)
      call closeunit (iu_AIC)

      xlabel=clabel(1:132)
      ItimeX=NINT(TAUX)
      NDAY = 24
      IYEAR1 = JC(41)
      itimei = itimex
      itimee = itimex
      itime0 = itimex
      print*,iyear1,ItimeX,xlabel

C**** read in FLAKE/FOCEAN data
      call openunit ("TOPO", iu_TOPO, .true., .true.)
      CALL READT (iu_TOPO,0,IM*JM,FOCEAN,1) ! ocean fraction
      CALL READT (iu_TOPO,0,IM*JM,FLAKE ,1) ! Lake fraction
      call closeunit (iu_TOPO)

C**** convert sea ice temperatures into enthalpy
C**** and initialize sea ice salinity to 3.2 ppt (0 in snow & lake ice).
      DO J=1,JM
        DO I=1,IM
          IF (RSI(I,J).gt.0) THEN
            MSI1=SNOWI(I,J)+ACE1I
            IF (FOCEAN(I,J).gt.0) THEN
              TFO=tfrez(sss0)   ! use default salinity
              SSI(3:4,I,J)=SSI0 * XSI(3:4)*MSI(I,J)
              IF (ACE1I*XSI(1).gt.SNOWI(I,J)*XSI(2)) THEN
                SSI(1,I,J)=SSI0 * (ACE1I-(ACE1I+SNOWI(I,J))* XSI(2))
                SSI(2,I,J)=SSI0 * (ACE1I+SNOWI(I,J))* XSI(2)
              ELSE
                SSI(1,I,J)=0.
                SSI(2,I,J)=SSI0 * ACE1I
              END IF
              HSI(1:2,I,J)=(SHI*MIN(HSI(1:2,I,J),TFO)-LHM)*XSI(1:2)*MSI1
     *             +LHM*SSI(1:2,I,J)
              HSI(3:4,I,J)=(SHI*MIN(HSI(3:4,I,J),TFO)-LHM)*XSI(3:4)
     *             *MSI(I,J)+LHM*SSI(3:4,I,J)
            ELSE
              HSI(1:2,I,J)=(SHI*MIN(HSI(1:2,I,J),0d0)-LHM)*XSI(1:2)*MSI1
              HSI(3:4,I,J)=(SHI*MIN(HSI(3:4,I,J),0d0)-LHM)*XSI(3:4)
     *             *MSI(I,J)
              SSI(1:4,I,J) = 0
            END IF
          ELSE
            MSI(I,J)=AC2OIM
            SNOWI(I,J)=0.
            IF (FOCEAN(I,J).gt.0) THEN
              TFO=tfrez(sss0)   ! use default salinity
              HSI(1:2,I,J) = (SHI*TFO-LHM*(1.-SSI0))*XSI(1:2)*ACE1I
              HSI(3:4,I,J) = (SHI*TFO-LHM*(1.-SSI0))*XSI(3:4)*AC2OIM
              SSI(1:2,I,J)=SSI0 * ACE1I  * XSI(1:2)
              SSI(3:4,I,J)=SSI0 * AC2OIM * XSI(3:4)
            ELSE
              HSI(1:2,I,J) = -LHM*XSI(1:2)*ACE1I
              HSI(3:4,I,J) = -LHM*XSI(3:4)*AC2OIM
              SSI(1:4,I,J) = 0.
            END IF
          END IF
          POND_MELT(I,J)=0.
          FLAG_DSWS(I,J)=.FALSE.
        END DO
      END DO


c     initialize the 3-d turbulent kinetic energy to be used in
c     the subroutine diffus.
      do j=1,jm
        do i=1,im
          do l=1,lm
            egcm(l,i,j)=egcm_init_max/(float(l)**2)
            w2gcm(l,i,j)=egcm(l,i,j)*by3
          end do
        end do
      end do

C**** set TSFREZ to zero (proper initalisation in init_diag)
      TSFREZ(:,:,:)=0.

C**** set default values for evaporation limiting arrays
      evap_max_ij(:,:) = 1.d0
      fr_sat_ij(:,:) = 1.d0
      qg_ij(:,:) = 0.d0

C**** extending ground hydrology
      if ( extend_gh ) then
        print *,"Extending ground hydrology data"
        do j=1,jm
          do i=1,im
            if ( conv_index(i,j,1) .ne. 0 ) then
              i0 = conv_index(i,j,1)
              j0 = conv_index(i,j,2)
              snowe(i,j)=snowe(i0,j0)
              tearth(i,j)=tearth(i0,j0)
              wearth(i,j)=wearth(i0,j0)
              aiearth(i,j)=aiearth(i0,j0)
              wbare(:,i,j) = wbare(:,i0,j0)
              wvege(:,i,j) = wvege(:,i0,j0)
              htbare(:,i,j)= htbare(:,i0,j0)
              htvege(:,i,j)= htvege(:,i0,j0)
              snowbv(:,i,j)= snowbv(:,i0,j0)
            endif
          end do
        end do
      endif

C**** set default foliage values
      Qfol=3.D-6
      Cint=0.0127D0
      cnc_ij=0.d0

      call io_rsf(trim(outfile),ItimeX,iowrite_mon,ioerr)
      print*,ioerr
      print*,"New rsf file written out to ",trim(outfile)

      if ( dump_gic ) then
        call openunit (trim(gic_file),iu_GIC, .true.,.false.)
        ioerr=-1
        call io_ocean  (iu_GIC,iowrite,ioerr)
        call io_seaice (iu_GIC,iowrite,ioerr)
        call io_earth  (iu_GIC,iowrite,ioerr)
        call io_soils  (iu_GIC,iowrite,ioerr)
        call io_landice(iu_GIC,iowrite,ioerr)
        call closeunit(iu_GIC)
        print *,ioerr
        print*,"New GIC file written out to ",trim(gic_file)
      endif

      stop
 800  print*,"Error reading in AIC file "
      stop
 810  print*,"Error reading in AIC file - end reached"
      stop
      end
