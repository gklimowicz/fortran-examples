      program convert_river_direction_files
!@sum RDconvert switch between grid-independent and ascii versions
!@+   of the river direction files
!@+   probably not safe on multiple processors (so use NPES=1!)
!@auth Gavin Schmidt
      USE CONSTANT, only : undef
      USE MODEL_COM, only : im,jm,lm,focean
      USE DOMAIN_DECOMP_1D, ONLY : init_app,init_grid
      USE DOMAIN_DECOMP_ATM, ONLY : grid,readt_parallel
     *     ,WRITE_PARALLEL,get,am_i_root
      USE LAKES, only : KDIREC,KD911,IFLOW,JFLOW,IFL911,JFL911
      USE LAKES_COM, only : NAMERVR
      USE FILEMANAGER, only : openunit,closeunit,nameunit
      USE GEOM, only : lon2d_dg,lat2d_dg,geom_b,imaxj,lonlat_to_ij
     *     ,lat_to_j
      IMPLICIT NONE
      character*80 filein,arg
      logical ascii, have_south_pole, have_north_pole
!@var I,J,I72 loop variables
!@var iwrap true if I direction is periodic and has no halo
      logical :: iwrap=.true.
      INTEGER I,J,I72,INM,KD,KDE,jmax_fill,jmin_fill,nrvr
      INTEGER, PARAMETER :: nrvrmx=42
      INTEGER I_0,I_1,J_0,J_1,J_0H,J_1H,I_0H,I_1H
      INTEGER iu_RVR  !@var iu_RVR unit number for river direction file
      INTEGER iu_TOPO
      integer get_dir
      CHARACTER TITLEI*80, TITLE2*80,CDIREC(IM,JM)*1, TITLE*80
      character(len=300) :: out_line
      REAL*4, dimension(im,jm) :: down_lat,down_lon,down_lat_911
     *     ,down_lon_911
      REAL*8, allocatable, dimension(:,:) :: down_lat_loc
     *     ,down_lon_loc,down_lat_911_loc,down_lon_911_loc
      REAL*4, dimension(nrvrmx) :: lat_rvr,lon_rvr
      INTEGER, dimension(nrvrmx) :: irvr,jrvr
      INTEGER, dimension(2) :: ij
      REAL*8, dimension(2) :: ll
      LOGICAL, DIMENSION(IM,jM) :: NODIR
      real*8 deg
C**** Usage: RDconvert -b RD_file.bin   => ascii version
C****        RDconvert -a RD_file.ascii => binary version

      call getarg(2,filein)
      call getarg(1,arg)
      if (arg .eq. "-a") then ! convert from ascii to bin
        ascii=.true.
      else ! from bin to ascii
        ascii=.false.
      end if

      call init_app()
      call init_grid(grid, im, jm, lm, CREATE_CAP=.true.)
      call alloc_drv()
      call GEOM_B

      CALL GET(GRID, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO= J_0H, J_STOP_HALO= J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),FOCEAN,1)
      call closeunit(iu_TOPO)

      if (ascii) then
C****
C**** Read in CDIREC: Number = octant direction
C****                 upper case = river mouth
C****                 lower case = emergency octant direction 
      call openunit(filein,iu_RVR,.false.,.true.)
      READ  (iu_RVR,'(A80)') TITLEI
      WRITE (out_line,*) 'River Direction file read: ',TITLEI
      CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
      READ  (iu_RVR,910)
      DO I72=1,1+(IM-1)/72
        DO J=JM, 1, -1
          READ  (iu_RVR,911) (CDIREC(I,J),I=72*(I72-1)+1,MIN(IM,I72*72))
        END DO
      END DO

C**** read in named rivers (if any)
      READ (iu_RVR,*,END=10)
      READ (iu_RVR,'(A80)',END=10) TITLE2
      READ (iu_RVR,*,END=10)
      IF (TITLE2.eq."Named River Mouths:") THEN
        DO I=1,NRVRMX,5
          READ(iu_RVR,'(5(A8,1X))',END=10) NAMERVR(I:MIN(NRVRMX,I+4))
        END DO
      END IF
 10   call closeunit (iu_RVR)
c      write(*,*) NAMERVR

C**** Create integral direction array KDIREC/KD911 from CDIREC

      ! Use unusual loop bounds to fill KDIREC/KD911 in halo
      DO J=MAX(1,J_0-1),MIN(JM,J_1+1)
      DO I=I_0, I_1
C**** KD: -16 = blank, 0-8 directions, 9 recirculation, 
C****      >9 named rivers, >48 emerg. dir
        KD= ICHAR(CDIREC(I,J)) - 48
        KDE=KD
        NODIR(I,J) = (KD.eq.-16) .or. (KD.gt.9 .and. KD.lt.48)
        IF (KDE.gt.48) THEN
           KDE=KDE-48  ! emergency direction for no outlet boxes
           KD=0        ! normally no outlet
        END IF

C**** Default direction is down (if ocean box), or no outlet (if not)
C**** Also ensure that all ocean boxes are done properly
        IF ((KD.lt.0 .or. KD.gt.9)) THEN ! .or. FOCEAN(I,J).eq.1.) THEN
          KDIREC(I,J)=0
          KD911(I,J) =0
        ELSE
          KDIREC(I,J) = KD
          KD911(I,J)  = KDE
        END IF
C**** Check for specified river mouths
        IF (KD.GE.17 .AND. KD.LE.42) THEN
          IF (FOCEAN(I,J).le.0) THEN
            WRITE(6,*)
     *       "Warning: Named river outlet must be in ocean",i
     *           ,j,FOCEAN(I,J),CDIREC(I,J)
          END IF
        END IF
      END DO
      END DO

      write(*,*) "IM, JM=",IM,JM
      INM=0
      DO J=1,JM
      DO I=1,IM
C**** KD: -16 = blank, 0-8 directions, 9 recirculation, >9 named rivers
        KD= ICHAR(CDIREC(I,J)) - 48
C**** Check for specified river mouths
        IF (KD.GE.17 .AND. KD.LE.42) THEN
          INM=INM+1        
          IF (CDIREC(I,J).ne.NAMERVR(INM)(1:1)) THEN
            WRITE(6,*)
     *           "Warning: Named river in RVR does not correspond"
     *           //" with letter in direction file. Please check"
            WRITE(6,*) "INM, CDIREC, NAMERVR = ",INM,CDIREC(I,J)
     *           ," ",NAMERVR(INM)
            NAMERVR(INM)=CDIREC(I,J)  ! set default
          END IF
          lat_rvr(inm)=lat2d_dg(I,J)
          lon_rvr(inm)=lon2d_dg(I,J)
          print*,inm,NAMERVR(INM),lat_rvr(inm),lon_rvr(inm)
        END IF
      END DO
      END DO
      NRVR=INM
C****
C**** From each box calculate the downstream river box
C****
      ! odd bounds to fill IFLOW and JFLOW in halo
      if(have_south_pole) then
        jmin_fill=2
      else
        jmin_fill=J_0H
      endif
      if(have_north_pole) then
        jmax_fill=JM-1
      else
        jmax_fill=J_1H
      endif
      DO J=JMIN_FILL,JMAX_FILL
        DO I=I_0H,I_1H
c should put in a check here whether we are at a nonexistent SW/NW/SE/NE
c halo corner of a cubed sphere face - see later whether necessary here.
          SELECT CASE (KDIREC(I,J))
          CASE (0)
            IFLOW(I,J) = I
            JFLOW(I,J) = J
          CASE (1)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J+1
            IF(I.eq.IM .and. iwrap)  IFLOW(I,J) = 1
          CASE (2)
            IFLOW(I,J) = I
            JFLOW(I,J) = J+1
          CASE (3)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J+1
            IF(I.eq.1 .and. iwrap)  IFLOW(I,J) = IM
          CASE (4)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J
            IF(I.eq.1 .and. iwrap)  IFLOW(I,J) = IM
          CASE (5)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J-1
            IF(I.eq.1 .and. iwrap)  IFLOW(I,J) = IM
          CASE (6)
            IFLOW(I,J) = I
            JFLOW(I,J) = J-1
          CASE (7)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J-1
            IF(I.eq.IM .and. iwrap)  IFLOW(I,J) = 1
          CASE (8)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J
            IF(I.eq.IM .and. iwrap)  IFLOW(I,J) = 1
          CASE (9) 
            IFLOW(I,J) = -9999
            JFLOW(I,J) = -9999
          END SELECT

C****
          SELECT CASE (KD911(I,J))   ! emergency directions
          CASE (1)
            IFL911(I,J) = I+1
            JFL911(I,J) = J+1
            IF(I.eq.IM .and. iwrap)  IFL911(I,J) = 1
          CASE (2)
            IFL911(I,J) = I
            JFL911(I,J) = J+1
          CASE (3)
            IFL911(I,J) = I-1
            JFL911(I,J) = J+1
            IF(I.eq.1 .and. iwrap)  IFL911(I,J) = IM
          CASE (4)
            IFL911(I,J) = I-1
            JFL911(I,J) = J
            IF(I.eq.1 .and. iwrap)  IFL911(I,J) = IM
          CASE (5)
            IFL911(I,J) = I-1
            JFL911(I,J) = J-1
            IF(I.eq.1 .and. iwrap)  IFL911(I,J) = IM
          CASE (6)
            IFL911(I,J) = I
            JFL911(I,J) = J-1
          CASE (7)
            IFL911(I,J) = I+1
            JFL911(I,J) = J-1
            IF(I.eq.IM .and. iwrap)  IFL911(I,J) = 1
          CASE (8)
            IFL911(I,J) = I+1
            JFL911(I,J) = J
            IF(I.eq.IM .and. iwrap)  IFL911(I,J) = 1
          END SELECT
        END DO
      END DO

C**** Poles are special cases
      IF (HAVE_SOUTH_POLE) Then
        DO I=1,IM
          IF(KDIREC(I,1).eq.2)  THEN
            IFLOW(1,1) = I
            JFLOW(1,1) = 2
          END IF
        END DO
        IFLOW(2:IM,1)=IFLOW(1,1)
        JFLOW(2:IM,1)=JFLOW(1,1)
        
        DO I=1,IM
          IF(KD911(I,1).eq.2)  THEN
            IFL911(1,1) = I
            JFL911(1,1) = 2
          END IF
        END DO
        IFL911(2:IM,1)=IFL911(1,1)
        JFL911(2:IM,1)=JFL911(1,1)
      END IF

      IF (HAVE_NORTH_POLE) Then
        DO I=1,IM
          IF(KDIREC(I,JM).eq.6)  THEN
            IFLOW(1,JM) = I
            JFLOW(1,JM) = JM-1
          END IF
        END DO
        IFLOW(2:IM,JM)=IFLOW(1,JM)
        JFLOW(2:IM,JM)=JFLOW(1,JM)
        
        DO I=1,IM
          IF(KD911(I,JM).eq.6)  THEN
            IFL911(1,JM) = I
            JFL911(1,JM) = JM-1
          END IF
        END DO
        IFL911(2:IM,JM)=IFL911(1,JM)
        JFL911(2:IM,JM)=JFL911(1,JM)
      END IF

C**** Next to pole too
      IF (J_0H.LE.2 .and. J_1H.GE.2) THEN
        DO I=1,IM
          IF(KDIREC(I,2).eq.6)  THEN
            IFLOW(I,2) = 1
            JFLOW(I,2) = 1
          END IF
          IF(KD911(I,2).eq.6)  THEN
            IFL911(I,2) = 1
            JFL911(I,2) = 1
          END IF
        END DO
      END IF
      IF (J_0H.LE.JM-1 .and. J_1H.GE.JM-1) THEN
        DO I=1,IM
          IF(KDIREC(I,JM-1).eq.2)  THEN
            IFLOW(I,JM-1) = 1
            JFLOW(I,JM-1) = JM
          END IF
          IF(KD911(I,JM-1).eq.6)  THEN
            IFL911(I,JM-1) = 1
            JFL911(I,JM-1) = JM
          END IF
        END DO
      END IF


      do j=1,jm
         do i=1,im
C**** calculate lat/lon pairs for downstream box.
            if (.not. nodir(i,j)) then
               if (iflow(i,j) .eq. -9999) then
                  down_lat(i,j)=-9999.
                  down_lon(i,j)=-9999.
               else
                  down_lat(i,j)=lat2d_dg(iflow(i,j),jflow(i,j))
                  down_lon(i,j)=lon2d_dg(iflow(i,j),jflow(i,j))
               endif
            else
               down_lat(i,j)=undef
               down_lon(i,j)=undef
            end if
            if (.not. nodir(i,j) .and. ifl911(i,j).gt.0) then
               down_lat_911(i,j)=lat2d_dg(ifl911(i,j),jfl911(i,j))
               down_lon_911(i,j)=lon2d_dg(ifl911(i,j),jfl911(i,j))
            else
               down_lat_911(i,j)=undef
               down_lon_911(i,j)=undef
            end if
         end do
      end do

C**** output binary file
      if (am_i_root()) then
      call openunit(trim(filein)//".bin",iu_RVR,.true.,.false.)
      write(iu_RVR) titlei
      write(iu_RVR) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *     ,lon_rvr(1:nrvr)
      title="Latitude of downstream river direction box"
      write(iu_RVR) title,down_lat
      title="Longitude of downstream river direction box"
      write(iu_RVR) title,down_lon
      title="Latitude of emergency downstream river direction box"
      write(iu_RVR) title,down_lat_911
      title="Longitude of emergency downstream river direction box"
      write(iu_RVR) title,down_lon_911

      call closeunit(iu_RVR)
      end if
      else  ! binary ==> ascii

C**** read in binary file:

        allocate(down_lat_loc(I_0H:I_1H,J_0H:J_1H),
     *       down_lon_loc(I_0H:I_1H,J_0H:J_1H),
     *       down_lat_911_loc(I_0H:I_1H,J_0H:J_1H),
     *       down_lon_911_loc(I_0H:I_1H,J_0H:J_1H) )

      call openunit(trim(filein),iu_RVR,.true.,.true.)
      read(iu_RVR) titlei
      read(iu_RVR) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *     ,lon_rvr(1:nrvr)
      CALL READT_PARALLEL(grid,iu_RVR,NAMEUNIT(iu_RVR),down_lat_loc,1)
      CALL READT_PARALLEL(grid,iu_RVR,NAMEUNIT(iu_RVR),down_lon_loc,1)
      CALL READT_PARALLEL(grid,iu_RVR,NAMEUNIT(iu_RVR),down_lat_911_loc
     *     ,1)
      CALL READT_PARALLEL(grid,iu_RVR,NAMEUNIT(iu_RVR),down_lon_911_loc
     *     ,1)
      call closeunit(iu_RVR)

      IFLOW=0 ; JFLOW=0
      IFL911=0 ; JFL911=0

C**** calculate I,J points corresponding to lat,lon
      do J=J_0,J_1
        do I=I_0,I_1
           if (down_lon_loc(i,j).ne.undef .and. 
     &          down_lon_loc(i,j).ne. -9999) then
              ll(1)=down_lon_loc(i,j)
              ll(2)=down_lat_loc(i,j)
              call lonlat_to_ij(ll,ij) 
              IFLOW(I,J)=ij(1) ; JFLOW(I,J)=ij(2)
              ll(1)=down_lon_911_loc(i,j)
              ll(2)=down_lat_911_loc(i,j)
              call lonlat_to_ij(ll,ij)
              IFL911(I,J)=ij(1) ; JFL911(I,J)=ij(2)
           elseif (down_lon_loc(i,j) == -9999.) then
              IFLOW(i,j)=-9999
              JFLOW(i,j)=-9999
              ll(1)=down_lon_911_loc(i,j)
              ll(2)=down_lat_911_loc(i,j)
              call lonlat_to_ij(ll,ij)
              IFL911(I,J)=ij(1) ; JFL911(I,J)=ij(2)
           end if
         END DO
      END DO

      do inm=1,nrvr
        ll(1)=lon_rvr(inm)
        ll(2)=lat_rvr(inm)
        call lonlat_to_ij(ll,ij)
        IRVR(INM)=ij(1) ; JRVR(INM)=ij(2)
      end do

C**** calculate river directions from IFLOW,JFLOW

      KDIREC=-16
      KD911=-16
      DO J=J_0,J_1
        DO I=I_0,I_1
          IF (IFLOW(I,J).gt.0) KDIREC(I,J)=get_dir(I,J,IFLOW(I,J)
     *         ,JFLOW(I,J),IM,JM)
          IF (IFL911(I,J).gt.0) KD911(I,J)=get_dir(I,J,IFL911(I,J)
     *         ,JFL911(I,J),IM,JM)
          IF (IFLOW(I,J) == -9999) KDIREC(I,J)=9
        END DO
      END DO

C**** output to file
      call openunit(trim(filein)//".asc",iu_RVR,.false.,.false.)
      write(iu_RVR,'(A)') trim(titlei)
      write(iu_RVR,*) 
      DO J=J_0,J_1
        DO I=I_0,I_1
          CDIREC(I,J)=CHAR(KDIREC(I,J)+48)
          IF (KDIREC(I,J).eq.0. .and. KD911(I,J).gt.0) THEN
            CDIREC(I,J)=CHAR(KD911(I,J)+96)
          END IF
          DO INM=1,NRVR
            IF (I.eq.IRVR(INM) .and. J.eq.JRVR(INM)) THEN
              CDIREC(I,J)=NAMERVR(INM)(1:1)
              EXIT
            END IF
          END DO
        END DO
      END DO
      DO I72=1,1+(IM-1)/72
        DO J=JM, 1, -1
          write (iu_RVR,'(72A)') (CDIREC(I,J),I=72*(I72-1)+1,MIN(IM,I72
     *         *72))
        END DO
      END DO
      write(iu_RVR,*)
      write(iu_RVR,'(A)') trim(title2)
      write(iu_RVR,'(A)')
     *     "(read in with '(5(A8,X))', first letter matches"
     *     //" character in direction array)"
      do inm=1,nrvr,5
        WRITE(iu_RVR,'(5(A8,1X))') NAMERVR(INM:MIN(NRVR,INM+4))
      end do

      call closeunit(iu_RVR)

      end if

      if (am_i_root()) then
        call write_river_dir_to_giss_format
      end if
      stop
C****
 910  FORMAT (A72)
 911  FORMAT (72A1)

      contains

      subroutine write_river_dir_to_giss_format
C
C This subroutine is intended to make easy to navigate 
C at the ascii file. The output is binary file that can 
C be see by ijprt. THIS FILE IS ONLY FOR HUMAN VIEW.
C You can see the indexes i,j and 
C geographical coordinates each box and value of river
C direction:
C                    2(b)
C               3(c) ^  1(a)
C                  \ | /
C            4(d)<-- .-->8(h)    0 means desert, 9 means recirculation
C                  / | \
C               5(c) v  7(g)
C                    6(f)
C
C For ijprt values a-h  converted (-1)-(-8),
C mouths of the rivers  A-Z ==> 901-926.
C So if you see number 903 it means "C" (third letter in alphabet),
C 926 means Z(last 26th letter) ==> mouth Zambezi river. 
C 999 means that box 100% cover ocean.
C
C
      character(len=1), dimension(45) :: cvalue = (/      
     1    ' ','0'                          
     2  , '1','2','3','4','5','6','7','8','9'  
     3  , 'a','b','c','d','e','f','g','h'  
     4  , 'A','B','C','D','E','F','G','H'  
     5  , 'I','J','K','L','M','N','O','P'  
     6  , 'Q','R','S','T','U','V','W','X'  
     7  , 'Y','Z'
     8  /)
      real*4,           dimension(45) :: rvalue = (/       
     1     999. , 0.                           
     2  ,   1.,   2.,  3.,  4.,  5.,  6.,  7.,  8., 9.  
     3  ,  -1.,  -2., -3., -4., -5., -6., -7., -8.  
     4  ,  901.,902.,903.,904.,905.,906.,907.,908.  
     5  ,  909.,910.,911.,912.,913.,914.,915.,916.  
     6  ,  917.,918.,919.,920.,921.,922.,923.,924.  
     7  ,  925.,926.
     8  /)
      real*4, dimension(im,jm) ::  PDIREC
      integer                  :: i,j,k
      character(len=80)        :: ijprt_file, title
      character(len= 1)        :: symb
      logical                  :: found

      PDIREC(:,:) = 999.

      do i = 1, im
      do j = 1, jm
         symb = CDIREC(I,J)
         found = .false.
         do k = 1, 45
           if( symb == cvalue(k) ) then
             PDIREC(I,J) = rvalue(k) 
             found = .true.
             exit 
           end if
         end do
         if( .not. found ) then
             write(6,*) ' File River Direction has wrong symbol=', symb
             write(6,*) ' Good values are 0-9, a-h, A-Z'
             STOP
         end if
      end do
      end do

      ijprt_file="ijprt_river_dir.bin"
      write(6,*)
      write(6,*) "You might want to look at file:"
      write(6,'(2x,A)') trim(ijprt_file)
      write(6,*) "by GISS utility: ijprt"

      call openunit(trim(ijprt_file),iu_RVR,.true.,.false.)
      title="999->Ocean, 0-9 -> Octan dir, " //
     a      "(-1)-(-8) -> a-h emergency dir, " //
     b      " 901-926 -> mouths(A-Z)"
      write(iu_RVR) title, PDIREC
      call closeunit(iu_RVR)

      end subroutine write_river_dir_to_giss_format
      end

