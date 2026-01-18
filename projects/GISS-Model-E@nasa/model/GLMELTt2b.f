      program glmeltt2b
      implicit none
      integer, parameter :: IM=72,JM=46
      integer ncols
      character*1 :: char(256)
      character*80 :: fname,oname
      CHARACTER*1, DIMENSION(IM,JM) :: CGLM   ! global array
      integer :: iu_GL,Irow,i,j,len,status
      real*4, allocatable :: RGLM(:,:)

      iu_GL=20
      fname="GLMELT_4X5.OCN"

      open(unit=iu_GL,FILE=fname,FORM ="FORMATTED",
     &     status ="UNKNOWN")

      write(*,*) "usage: converts txt file to binary. 
     &     Title (1st line) MUST have same number of characters as 
     &     all other lines in the file (for example 72)"
      do 
      read(iu_GL,'(256A1)',size=ncols,advance='no',eor=800) char
      write(*,*) char
      enddo
  800 write(*,*) "ncols",ncols

      READ  (iu_GL,*)

C**** assumes a fix width slab - will need adjusting for CS
      DO Irow=1,1+(IM-1)/ncols
        DO J=JM,1,-1
          READ (iu_GL,200)
     *         (CGLM(I,J),I=ncols*(Irow-1)+1,MIN(IM,Irow*ncols))
        END DO
      END DO
  200 FORMAT(<ncols>A1)
C****

      allocate(RGLM(IM,JM))

      do i=1,IM
         do j=1,JM
            if (CGLM(i,j) .eq. "1") then
               RGLM(i,j)=1.0
            else
               RGLM(i,j)=0
            endif
         enddo
      enddo

c      write(*,*) "RGLM=",RGLM
      close(iu_GL)

      oname=trim(fname)//".bin"
      write(*,*) oname
      write(*,*) "size iglm",size(RGLM)
      open(unit=3200, FILE=oname,FORM='unformatted',STATUS='unknown')


      write(3200) char(1:80),RGLM
      close(unit=3200)
      deallocate(RGLM)

      end program glmeltt2b
