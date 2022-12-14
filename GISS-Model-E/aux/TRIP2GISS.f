
      program trip2giss
!@sum converts river directions from TRIP format to GISS format
!@+   see http://hydro.iis.u-tokyo.ac.jp/~taikan/TRIPDATA/TRIPDATA.html
!@+   1->2
!@+   8->3
!@+   7->4
!@+   6->5
!@+   5->6
!@+   4->7
!@+   3->8
!@+   2->1
!@+   9->8 if land cell
!@+   9->undef (SPACE) if fractional ocean cell
!@+   compile using ifort -convert big_endian TRIP2GISS.f -o trip2giss
      integer, parameter :: im=360,jm =180
      character*80 filein,titlei,title,title2,name
      real*4 :: FOCEAN(IM,JM)
      character*1 :: CDIREC(im,jm) 
      integer :: iu_RVR, I72, kd

      iu_TOPO=30

      name="Z1X1N"
      open(iu_TOPO,FILE=name,FORM='unformatted', STATUS='old')
      read(iu_TOPO) title,FOCEAN
      close(iu_TOPO)

      WRITE(*,*) FOCEAN

      iu_RVR=20
      filein="RD_modelE_O.RVR"
      open(iu_RVR,FILE=filein,FORM='formatted', STATUS='old')
      read(iu_RVR,'(A80)') TITLEI
      write (*,*) 'River Direction file read: ',TITLEI

      READ  (iu_RVR,910)
      DO I72=1,1+(IM-1)/72
        DO J=JM, 1, -1
          READ  (iu_RVR,911) (CDIREC(I,J),I=72*(I72-1)+1,MIN(IM,I72*72))
        END DO
      END DO

 910  FORMAT (A72)
 911  FORMAT (72A1)

      READ (iu_RVR,*,END=10)
      READ (iu_RVR,'(A80)',END=10) TITLE2
      READ (iu_RVR,*,END=10)

 10   close (iu_RVR)

      do j=1,JM
         do i=1,IM
            KD= ICHAR(CDIREC(I,J)) - 48
c            write(*,*) kd
            if (kd .eq. 1) CDIREC(i,j)=ACHAR(50)
            if (kd .eq. 8) CDIREC(i,j)=ACHAR(51)
            if (kd .eq. 7) CDIREC(i,j)=ACHAR(52)
            if (kd .eq. 6) CDIREC(i,j)=ACHAR(53)
            if (kd .eq. 5) CDIREC(i,j)=ACHAR(54)
            if (kd .eq. 4) CDIREC(i,j)=ACHAR(55)
            if (kd .eq. 3) CDIREC(i,j)=ACHAR(56)
            if (kd .eq. 2) CDIREC(i,j)=ACHAR(49)
            if (kd .eq. 9) then
               if (FOCEAN(i,j) .gt. 1.e-6 .or. FOCEAN(i+1,j) .gt. 1.e-6  
     &              .or. FOCEAN(i+1,j+1) .gt. 1.e-6 
     &              .or. FOCEAN(i,j+1) .gt. 1.e-6  
     &              .or. FOCEAN(i-1,j+1) .gt. 1.e-6 
     &              .or. FOCEAN(i-1,j) .gt. 1.e-6  
     &              .or. FOCEAN(i-1,j-1) .gt. 1.e-6 
     &              .or. FOCEAN(i,j-1) .gt. 1.e-6
     &              .or. FOCEAN(i+1,j-1) .gt. 1.e-6   ) then  
                  CDIREC(i,j)=ACHAR(32)  !SPACE
               else
                  CDIREC(i,j)=ACHAR(48)
               endif
               write(*,*) "found 9"
            endif             
         enddo
      enddo

      name=trim(filein)//".GISS"
      open(iu_RVR,FILE=name,FORM='formatted', STATUS='unknown')
      write(iu_RVR,'(A)') trim(titlei)
      write(iu_RVR,*) 

      DO I72=1,1+(IM-1)/72
         DO J=JM,1,-1
            write (iu_RVR,'(72A)') (CDIREC(I,J),I=72*(I72-1)+1,
     &           MIN(IM,I72*72))
        END DO
      END DO
      write(iu_RVR,*)
      write(iu_RVR,'(A)') trim(title2)
      close(iu_RVR)

      end
