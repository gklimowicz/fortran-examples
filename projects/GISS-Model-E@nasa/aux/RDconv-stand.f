      program convert_river_direction_files
!@sum RDconvert switch between binary and ascii versions
!@+   of the river direction files
!@auth Gavin Schmidt. Standalone version D. Gueyffier
      IMPLICIT NONE
      integer, parameter :: im=144,jm=90
      real*8,parameter :: undef=-1.d30
      character*80 filein,arg,fin
      logical ascii
!@var I,J,I72 loop variables
!@var iwrap true if I direction is periodic and has no halo
      logical :: iwrap=.true.
      INTEGER I,J,I72,INM,KD,KDE,nrvr
      INTEGER, PARAMETER :: nrvrmx=42
      integer get_dir
      CHARACTER TITLEI*80, TITLE2*80,TITLE*80,zfile*80
      character*1, :: CDIREC(im,jm)
      character(len=300) :: out_line
      REAL*4, dimension(im,jm) :: down_lat,down_lon,
     *     down_lat_911,down_lon_911,FOCEAN
      real*8 :: lat_dg(JM),lon_dg(IM)
      integer, dimension(im,jm) :: KDIREC,KD911,IFLOW,JFLOW,
     *     IFL911,JFL911
      CHARACTER*8, DIMENSION(NRVRMX) :: NAMERVR
      real*8 :: dlat_dg,dlon_dg,fjeq
      REAL*4, dimension(nrvrmx) :: lat_rvr,lon_rvr
      INTEGER, dimension(nrvrmx) :: irvr,jrvr
      LOGICAL, dimension(im,jm) :: NODIR
C**** Usage: RDconvert -b RD_file.bin   => ascii version
C****        RDconvert -a RD_file.ascii => binary version

      call getarg(2,filein)
      call getarg(1,arg)
      if (arg .eq. "-a") then ! convert from ascii to bin
        ascii=.true.
      else ! from bin to ascii
        ascii=.false.
      end if

      zfile="Z144X90N_nocasp.1"
      zfile=trim(zfile)
      open(10,FILE=zfile,FORM='unformatted', STATUS='unknown')
      read(10) title,FOCEAN
      close(10)
      write(*,*) "ZFILE:",title

      dlon_dg = 360./dble(IM)
      dlat_dg = 180./dble(JM)  
      lat_dg(1)  = -90.
      lat_dg(JM) = 90.
      FJEQ=.5*(1+JM)
      do J=2,JM-1
        lat_dg(J)  = dlat_dg*(J-FJEQ)
      enddo
      lon_dg(1) = -180.+360./(2.*FLOAT(IM))
      do i=2,IM
        lon_dg(i) = lon_dg(i-1)+360./FLOAT(IM)
      enddo

      if (ascii) then
C****
C**** Read in CDIREC: Number = octant direction
C****                 upper case = river mouth
C****                 lower case = emergency octant direction 
         open(unit=20,file=filein,form='formatted',status='old',
     &     action='read')

         READ  (20,'(A80)') TITLEI
         WRITE (6,*) 'River Direction file read: ',TITLEI
         READ  (20,910)
         DO I72=1,1+(IM-1)/72
           DO J=JM, 1, -1
             READ (20,911) (CDIREC(I,J),I=72*(I72-1)+1,MIN(IM,I72*72))
           END DO
         END DO

C**** read in named rivers (if any)
         READ (20,*,END=10)
         READ (20,'(A80)',END=10) TITLE2
         READ (20,*,END=10)
         IF (TITLE2.eq."Named River Mouths:") THEN
            DO I=1,NRVRMX,5
               READ(20,'(5(A8,1X))',END=10) NAMERVR(I:MIN(NRVRMX,I+4))
            END DO
         END IF
 10      close(20)
         write(*,*) NAMERVR

C**** Create integral direction array KDIREC/KD911 from CDIREC

      DO J=1,JM
      DO I=1,IM
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
          lat_rvr(inm)=lat_dg(J)
          lon_rvr(inm)=lon_dg(I)
          print*,inm,NAMERVR(INM),lat_rvr(inm),lon_rvr(inm)
        END IF
      END DO
      END DO
      NRVR=INM
C****
C**** From each box calculate the downstream river box
C****

      DO J=2,JM-1
        DO I=1,IM
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
      

C**** Next to pole too

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


      do j=1,jm
         do i=1,im
C**** calculate lat/lon pairs for downstream box.
            if (.not. nodir(i,j)) then
               if (iflow(i,j) .eq. -9999) then
                  down_lat(i,j)=-9999.
                  down_lon(i,j)=-9999.
               else
                  down_lat(i,j)=lat_dg(jflow(i,j))
                  down_lon(i,j)=lon_dg(iflow(i,j))
               endif
            else
               down_lat(i,j)=undef
               down_lon(i,j)=undef
            end if
            if (.not. nodir(i,j) .and. ifl911(i,j).gt.0) then
               down_lat_911(i,j)=lat_dg(jfl911(i,j))
               down_lon_911(i,j)=lon_dg(ifl911(i,j))
            else
               down_lat_911(i,j)=undef
               down_lon_911(i,j)=undef
            end if
         end do
      end do

C**** output binary file
      fin=(trim(filein))//".bin"
      open(30,FILE=fin,
     &     FORM='unformatted', STATUS='unknown')
      write(30) titlei
      write(30) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *     ,lon_rvr(1:nrvr)
      title="Latitude of downstream river direction box"
      write(30) title,down_lat
      title="Longitude of downstream river direction box"
      write(30) title,down_lon
      title="Latitude of emergency downstream river direction box"
      write(30) title,down_lat_911
      title="Longitude of emergency downstream river direction box"
      write(30) title,down_lon_911
      close(30)

      else                      ! binary ==> ascii
C**** read in binary file:
         fin=trim(filein)
         open(40,FILE=fin,FORM='unformatted', STATUS='unknown')
         read(40) titlei
         read(40) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *        ,lon_rvr(1:nrvr)
         read(40) title,down_lat
         read(40) title,down_lon
         read(40) title,down_lat_911
         read(40) title,down_lon_911
         close(40)

      IFLOW=0 ; JFLOW=0
      IFL911=0 ; JFL911=0

C**** calculate I,J points corresponding to lat,lon
      do J=1,JM
        do I=1,IM
           if (down_lon(i,j).ne.undef .and. 
     &          down_lon(i,j).ne. -9999) then
             
              IFLOW(I,J)= nint( .5*(im + 1) + 
     &             down_lon(i,j)/dlon_dg )
              JFLOW(I,J) = nint( .5*(jm + 1) + 
     &             down_lat(i,j)/dlat_dg )
c              write(*,*) I,J,IFLOW(I,J),JFLOW(I,J)
              if(JFLOW(I,J)>jm) JFLOW(I,J)=jm
              if(JFLOW(I,J)<1 ) JFLOW(I,J)=1

              IFL911(I,J)= nint( .5*(im + 1) + 
     &             down_lon_911(i,j)/dlon_dg )
              JFL911(I,J) = nint( .5*(jm + 1) + 
     &             down_lat_911(i,j)/dlat_dg )
              if(JFL911(I,J)>jm) JFL911(I,J)=jm
              if(JFL911(I,J)<1 ) JFL911(I,J)=1

           elseif (down_lon(i,j) == -9999.) then
              IFLOW(i,j)=-9999
              JFLOW(i,j)=-9999

              IFL911(I,J)= nint( .5*(im + 1) + 
     &             down_lon_911(i,j)/dlon_dg )
              JFL911(I,J) = nint( .5*(jm + 1) + 
     &             down_lat_911(i,j)/dlat_dg )
              if(JFL911(I,J)>jm) JFL911(I,J)=jm
              if(JFL911(I,J)<1 ) JFL911(I,J)=1

           end if
         END DO
      END DO

      do i=1,nrvr
              IRVR(i)= nint( .5*(im + 1) + 
     &             lon_rvr(i)/dlon_dg )
              JRVR(i) = nint( .5*(jm + 1) + 
     &             lat_rvr(i)/dlat_dg )
              if(JRVR(i)>jm) JRVR(i)=jm
              if(JRVR(i)<1 ) JRVR(i)=1
      end do

C**** calculate river directions from IFLOW,JFLOW

      KDIREC=-16
      KD911=-16
      DO J=1,JM
        DO I=1,IM
          IF (IFLOW(I,J).gt.0) KDIREC(I,J)=get_dir(I,J,IFLOW(I,J)
     *         ,JFLOW(I,J),IM,JM)
          IF (IFL911(I,J).gt.0) KD911(I,J)=get_dir(I,J,IFL911(I,J)
     *         ,JFL911(I,J),IM,JM)
          IF (IFLOW(I,J) == -9999) KDIREC(I,J)=9
        END DO
      END DO

C**** output to file
      fin=(trim(filein))//".asc"
      open(unit=50,file=fin,
     &     form='formatted',status='unknown')
      
      write(50,'(A)') trim(titlei)
      write(50,*) 
      DO J=1,JM
        DO I=1,IM
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
          write (50,'(72A)') (CDIREC(I,J),I=72*(I72-1)+1,MIN(IM,I72
     *         *72))
        END DO
      END DO
      write(50,*)
      write(50,'(A)') trim(title2)
      write(50,'(A)')
     *     "(read in with '(5(A8,X))', first letter matches"
     *     //" character in direction array)"
      do inm=1,nrvr,5
        WRITE(50,'(5(A8,1X))') NAMERVR(INM:MIN(NRVR,INM+4))
      end do

      close(50)

      endif

      call write_river_dir_to_giss_format
      
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
      character*80 :: fin
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

      fin=trim(ijprt_file)
      open(60,file=fin,
     &     FORM='unformatted', STATUS='unknown')
      title="999->Ocean, 0-9 -> Octan dir, " //
     a      "(-1)-(-8) -> a-h emergency dir, " //
     b      " 901-926 -> mouths(A-Z)"
      write(60) title, PDIREC
      close(60)

      end subroutine write_river_dir_to_giss_format
      end

      integer function get_dir(I,J,ID,JD,IM,JM)
!@sum get_dir derives the locally orientated river direction
      integer, intent(in) :: I,J,ID,JD,IM,JM
      integer ::  DI,DJ


      DI=I-ID
      IF (DI.eq.IM-1) DI=-1
      IF (DI.eq.1-IM) DI=1
      DJ=J-JD
      get_dir=-99
c      write(*,*) "DI, DJ=",DI,DJ
      if (DI.eq.-1 .and. DJ.eq.-1) then
        get_dir=1
      elseif (DI.eq.-1 .and. DJ.eq.0) then
        get_dir=8
      elseif (DI.eq.-1 .and. DJ.eq.1) then
        get_dir=7
      elseif (DI.eq.0 .and. DJ.eq.1) then
        get_dir=6
      elseif (DI.eq.0 .and. DJ.eq.0) then
        get_dir=0
      elseif (DI.eq.0 .and. DJ.eq.-1) then
        get_dir=2
      elseif (DI.eq.1 .and. DJ.eq.-1) then
        get_dir=3
      elseif (DI.eq.1 .and. DJ.eq.0) then
        get_dir=4
      elseif (DI.eq.1 .and. DJ.eq.1) then
        get_dir=5
      end if
      if (J.eq.JM) then         ! north pole
         if (DI.eq.0) then
            get_dir=6
         else
            get_dir=8
         end if
      elseif (J.eq.1) then      ! south pole
        if (DI.eq.0) then
           get_dir=2
        else
           get_dir=8
        end if
      elseif (J.eq.JM-1) then
         if (JD.eq.JM) get_dir=2
      elseif (J.eq.2) then
         if (JD.eq.1) get_dir=6
      end if
      if (get_dir.eq.-99) then
        print*,"get_dir error",i,j,id,jd
        get_dir=0
      end if
         
      return
      end function get_dir
