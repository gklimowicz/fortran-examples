      program convert_rsf_files
C**** convert rsf files from model II' (B399) format to modelE format
C**** compile with: gmake conv_rsf.o
C**** f90 -o conv_rsf conv_rsf.o *.o -O2 -64 -mips4 -static \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE CONSTANT, only : lhm,shi
      USE MODEL_COM, only : im,jm,lm  !  ,wm,u,v,t,p,q,xlabel
     *     ,iowrite_mon,focean,nday,itime,itimei,itimee,itime0,iyear1
     *     ,ioread,iowrite,ioreadnt
!     USE SOMTQ_COM
      USE GHY_COM, only : snowe,tearth,wearth,aiearth,w_ij
     *     ,ht_ij,snowbv,ngm,evap_max_ij,fr_sat_ij,qg_ij
!     USE RAD_COM, only : rqt,lm_req
!     USE CLOUDS_COM, only : ttold,qtold,svlhx,rhsav,cldsav
!     USE DIAG_COM, only : keynr,tsfrez
!     USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
!    *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,
!    *     ustar_pbl,egcm,tgvavg,qgavg
!     USE STATIC_OCEAN, only : tocean,z1o,sss0
!     USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi,pond_melt,flag_dsws
!     USE SEAICE, only : ace1i,xsi,ac2oim,ssi0,tfrez
!     USE LANDICE_COM, only : tlandi,snowli
!     USE LAKES_COM, only : flake
      USE FILEMANAGER
      USE DOMAIN_DECOMP_1D, ONLY : init_app,init_grid,grid,finish_app
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60              ,clabel*156
      INTEGER IARGC,iu_GIC,I,J,L,N,ioerr,iu_TOPO   ,jc(100)
      REAL*8 TAUX,X                                ,rc(161)
      REAL*8 MSI1,tfo
      INTEGER ItimeX
!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5
      integer :: conv_index(im,jm,2)
      character*120 cindex_file
      integer iu_CI, i0, j0
      logical :: extend_gh = .false.
      character*1, allocatable :: dummy(:), dummy2(:)

      ! hack to read irrelevant records for different resolutions
      if ( im==72 ) then
        allocate( dummy(106064), dummy2(331280) )
      else if ( im==144 ) then
        allocate( dummy(414800), dummy2(1296080) )
      else
        call stop_model("This resolution is not supported",255)
      endif

      IF (IARGC().lt.2) THEN
        PRINT*,"Convert GIC file to the extended format"
        PRINT*,"To create a restart file for the standard land mask:"
        PRINT*,"conv_rsf filename output_file"
        PRINT*,"To create a restart file with the extended land"
        PRINT*,"surface data (which can be used with any land mask):"
        PRINT*,"ext_gic filename output_file conv_indices_file"
        STOP
      END IF

        call init_app()
        call init_grid(grid,im,jm,lm)
        call alloc_drv()

      IF (IARGC() >= 3 ) extend_gh = .true.

      CALL GETARG(1,infile)
      CALL GETARG(2,outfile)

      if ( extend_gh ) then
        CALL GETARG(3,cindex_file)
        call openunit (trim(cindex_file), iu_CI, .true. , .true. )
        read(iu_CI) conv_index
        call closeunit (iu_CI)
      endif

c??   iu_GIC=9
c??   OPEN(iu_GIC,FILE=trim(infile),FORM="UNFORMATTED",STATUS="OLD")
      call openunit (trim(infile), iu_GIC, .true. , .true. )
      !call io_rsf(iu_GIC,ItimeX,ioread,ioerr)
        ioerr=-1
        read(iu_GIC) dummy ! ignore first line (ocean ic done in init_OCEAN)
        !call io_seaice (iu_GIC,ioreadnt,ioerr)
        read(iu_GIC) dummy2
        call io_earth  (iu_GIC,ioreadnt,ioerr)
        call io_soils  (iu_GIC,ioreadnt,ioerr)
        call io_landice(iu_GIC,ioreadnt,ioerr)
        if (ioerr.eq.1) then
          WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
          call stop_model("INPUT: GIC READ IN ERROR",255)
        end if
      call closeunit (iu_GIC)
      print *, ioerr
      print *, 'Read input file'


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
              w_ij(:,:,i,j) = w_ij(:,:,i0,j0)
              ht_ij(:,:,i,j)= ht_ij(:,:,i0,j0)
              snowbv(:,i,j)= snowbv(:,i0,j0)
            endif
          end do
        end do
      endif


c??   OPEN(iu_GIC,FILE=trim(outfile),
c??  *     FORM="UNFORMATTED",STATUS="UNKNOWN")
      call openunit (trim(outfile),iu_GIC, .true.,.false.)

      !call io_rsf(iu_GIC,ItimeX,iowrite_mon,ioerr)
      write(iu_GIC) dummy
      !call io_seaice (iu_GIC,iowrite,ioerr)
      write(iu_GIC) dummy2
      call io_earth  (iu_GIC,iowrite,ioerr)
      call io_soils  (iu_GIC,iowrite,ioerr)
      call io_landice(iu_GIC,iowrite,ioerr)

      close (iu_GIC)
      print*,ioerr
      print*,"New rsf file written out to ",trim(outfile)
      call finish_app()
      stop
 800  print*,"Error reading in file"
 810  print*,"Error reading in file"
      stop
      end
