      program vegtypet2b
!@sum converting txt file vegtype.global and LAI files to binary. 
!@+ We adopt a format which is compatible with 
!@+ the binary format used elsewhere in the model:
!@+ do K=1:NVEGTYPE
!@+ read(iu_data) TITLE_VEG_TYPE(K), IUSE_glob(1:IM,1:JM,K)
!@! enddo 
!@  note that the upper bound is now NVEGTYPE. 
!@+ The main change is that now the file contains many entries where IUSE_glob=0.
!@+ We use an increasing ordering of the vegetation types (do K=1:NVEGTYPE). This 
!@+ provides a global definition of the ordering instead of the local definition used
!@+ in the native txt files, i.e
!@+ NATIVE:  
!@+    cell i=22, j=11,  7  0  24  25  52  53  8  54  73 263 100 413  12  12 127
!@+ INCREASING ORDER: 
!@+    cell i=22, j=11,  7  0  8  24  25  52  53  54  73  12 263 100 413  12 127
!@+   
!@+ When reading the new binary file the arrays IREG_glob and ILAND_GLOB can be deduced 
!@+ simply by checking which entries of IUSE_glob are nonzero
!@+ 
!@+ The new binary file vegtype.global is used in a separate routine (regridinput) to regrid 
!@+ the fractions IUSE_glob to the cubed-sphere in such a way that sum(IUSE(I0,J0,:),3) = 1000 
!@+ for any (I0,J0) on the cubed-sphere surface.
!@+ 
!@auth Denis Gueyffier
!@+ 
c      integer, parameter :: IM=72, JM=46
      implicit none
      integer, parameter :: IM=720, JM=360
      INTEGER, PARAMETER :: NPOLY   = 20,
     &                      NTYPE   = 16,
     &                      NVEGTYPE= 74
      INTEGER, PARAMETER :: NWAT=6
      integer, dimension(IM,JM) :: IREG_glob
      integer, dimension(IM,JM,NTYPE) :: ILAND_glob,IUSE_glob
      real*4, dimension(IM,JM,NVEGTYPE) :: FUSE_glob,XLAI_TEMP1,
     &     XLAI_TEMP2,FUSE_TEMP1,FUSE_TEMP2
      real*8, dimension(im,jm,ntype) :: XLAI_glob
      real*4, dimension(im,jm,ntype) :: XLAI_glob4
      real*4, dimension(im,jm,nvegtype) :: XLAI_out
      real*4, dimension(im,jm) :: ones, FOCEAN

      INTEGER :: I,J,K,L,iu_data,iveg,imonth,index,in,jn
      character*80 :: fname,oname,TITLE_VEG_TYPE,TITLE
      character(len=2) :: vegtype,c2month 
      character(len=1) :: c1month 
      real*4, parameter :: skip = -9.999
      real*4 :: nd

      write(6,*) 'READING land types and fractions'
      write(6,*) 'and converting to binary'
      fname='vegtype.global'

      open(20,FILE=fname,FORM="FORMATTED",status="UNKNOWN")

 100  READ(20,'(20I4)',end=110) I,J,IREG_glob(I,J),
     &     (ILAND_glob(I,J,K),K=1,IREG_glob(I,J)),
     &     (IUSE_glob(I,J,K),K=1,IREG_glob(I,J))
      GO TO 100
 110  CONTINUE
      close(20)

      FUSE_glob(:,:,:)= 0.

c***  Warning : ILAND_GLOB takes values between 0 and 73, and 0 isn't a friendly bound 
c***  for fortran arrays
      do I=1,IM
         do J=1,JM
            do K=1,IREG_glob(I,J)
               FUSE_glob(I,J,ILAND_GLOB(I,J,K)+1)=
     &              IUSE_glob(I,J,K)/1000.
            enddo
         enddo
      enddo
      
      open(unit=3100, FILE='ZHXH00',FORM='unformatted',
     &     STATUS='unknown')
       
      read(3100) title,FOCEAN
      close(3100)

        do I=1,IM
           do J=1,JM
              do K=1,NVEGTYPE
                 if (FOCEAN(I,J) .gt. 0.
     &                .and. FUSE_glob(I,J,K) .eq. 0.) then
                    FUSE_temp1(I,J,K)=skip
                 else
                    FUSE_temp1(I,J,K)=FUSE_glob(I,J,K)
                 endif
              enddo
           enddo
        enddo
        write(*,*) "THERE"
        do k=1,nvegtype
           call fillin(FUSE_TEMP2(:,:,K),
     &          FUSE_TEMP1(:,:,K),im,jm,skip,.false.)
        enddo

c***  Make sure that fractions at a particular cell add up to 1
      do J=1,JM
      do I=1,IM
         if (FUSE_TEMP2(I,J,54) .le. 0.001 
     &        .and. FUSE_TEMP2(I,J,71) .ge. .9995 ) then
            write(*,*) "tundra < 0.001 landice",i,j,
     &        FUSE_TEMP2(I,J,54),FUSE_TEMP2(I,J,71) 
            FUSE_TEMP2(I,J,71)=FUSE_TEMP2(I,J,71)+FUSE_TEMP2(I,J,54)
     &        -0.001
            FUSE_TEMP2(I,J,54)=0.001
            write(*,*) "after",FUSE_TEMP2(I,J,54),FUSE_TEMP2(I,J,71)
         endif
         if ( FUSE_TEMP2(I,J,1) .le. 0.) FUSE_TEMP2(I,J,1)=0.
         if (abs(sum(FUSE_TEMP2(I,J,:))-1.0) .gt. 0.001) then 
c            write(6,*) "Warning veg. fractions at cell:",i,j,
c     &        "do not add up to 1. Difference is more than 1 per mil"
            FUSE_TEMP2(I,J,:)=FUSE_TEMP2(I,J,:)/sum(FUSE_TEMP2(I,J,:))
         endif
      enddo
      enddo

      oname=trim(fname)//".bin"
      write(*,*) oname

      open(unit=3200, FILE=oname,FORM='unformatted',STATUS='unknown')

      do iveg=1,NVEGTYPE
         write(vegtype,'(i2)') iveg
         TITLE_VEG_TYPE="veg type="//vegtype
         write(3200) TITLE_VEG_TYPE,FUSE_TEMP2(:,:,iveg)
      enddo
 
      close(unit=3200)

      open(unit=3200, FILE='sumfuse',FORM='unformatted',
     &     STATUS='unknown')
      TITLE="sum fuse"
      write(3200) TITLE,SUM(FUSE_TEMP2,3)

      close(unit=3200)

C     Read lai:
      do imonth=1,12
         if (imonth .lt. 10) then
            write(c1month,'(i1)') imonth
            c2month="0"//c1month
         else
            write(c2month,'(i2)') imonth
         endif
         fname='lai'//c2month//'.global'
         write(*,*) fname
         open(20,FILE=fname,FORM="FORMATTED",status="UNKNOWN")

         XLAI_glob(:,:,:)=0.

 10      READ(20,"(3I3,20F5.1)",END=20) I,J,INDEX,
     &        (XLAI_glob(I,J,K),K=1,INDEX)
         GOTO 10
 20     close(20)
        
        XLAI_glob4=XLAI_glob

        XLAI_out(:,:,:)=0.
        
        do I=1,IM
           do J=1,JM
              do K=1,IREG_glob(I,J)
                 XLAI_out(I,J,ILAND_GLOB(I,J,K)+1)=XLAI_GLOB4(I,J,K)
              enddo
           enddo
        enddo

        do I=1,IM
           do J=1,JM
              do K=1,NVEGTYPE
                 if (FOCEAN(I,J) .gt. 0. 
     &                .and. XLAI_out(I,J,K) .eq. 0.) then
                    XLAI_temp1(I,J,K)=skip
                 else
                    XLAI_temp1(I,J,K)=XLAI_out(I,J,K)
                 endif
              enddo
           enddo
        enddo
        
        do k=1,nvegtype
           call fillin(XLAI_TEMP2(:,:,K),
     &          XLAI_TEMP1(:,:,K),im,jm,skip,.false.)
        enddo 


        oname=trim(fname)//".bin"
        write(*,*) oname
        open(unit=3200, FILE=oname,FORM='unformatted',
     &       STATUS='unknown')
        
        do iveg=1,NVEGTYPE
           write(vegtype,'(i2)') iveg
           TITLE_VEG_TYPE="LAI corresponding to veg type="//vegtype
           write(3200) TITLE_VEG_TYPE,XLAI_temp2(:,:,iveg)
        enddo
        
        close(3200)

      enddo

      open(unit=3500, FILE="sumlai",FORM='unformatted',
     &     STATUS='unknown')
      title="sum veg"
      write(3500) title,sum(XLAI_temp2,3)
      close(3500)

      end program vegtypet2b


      subroutine fillin (fil,dat,im,jm,skip,vrb)
C**** nearest neighbor fill of 'skip's for all points where pe>0
      real fil(im,jm),dat(im,jm)
      logical vrb

      write(*,*) im,jm
      do j=1,jm
        do i=1,im
           fil(i,j)=dat(i,j)
           if(dat(i,j) .eq. skip) then
             call getnn(im,jm,dat,skip,i,j,in,jn,nd)
             fil(i,j)=dat(in,jn)
             if(vrb) write(*,*) i,j,' filled from',in,jn,dat(in,jn),nd
           end if
        end do
      end do

      return
      end

      subroutine getnn(im,jm,dat,skip,i0,j0,in,jn,nd)
C**** dat(i0,j0)=skip; find (i,j) with dat(i,j).ne.skip 'closest'
C**** to point (i0,j0) - 'distance' = (j-j0)^2 + min(|i-i0|,im-|i-i0|)
C**** taking into account that lat is more important than lon: (in,jn)
      real dat(im,jm)

C****    the distance function
      ndist(ix,jx,iy,jy)=(jx-jy)**2+ 
!     &     (ix-iy)**2
     &     min(abs(ix-iy),im-abs(ix-iy))

      nd=jm*jm+im       ! larger than ndist for any pair of points
      do 20 id=1,im/2                   ! first try j=j0
      if(id.ge.nd) go to 40
      do 10 m=-1,1,2
      i=i0+m*id
      if(i.gt.im) i=i-im
      if(i.lt.1) i=i+im
      if (dat(i,j0).eq.skip) go to 10
        nd0=id  !  =ndist(i0,j0,i,j0)
        if (nd0.lt.nd) then
          in=i
          jn=j0
          nd=nd0
          go to 40
        end if
   10 continue
   20 continue

   40 do 100 jd=1,jm-1
      if (jd*jd.ge.nd) return
      do 60 m=-1,1,2
      j=j0+m*jd
      if(j.gt.jm) go to 60
      if(j.le.0) go to 60
      do 50 i=1,im
      if (dat(i,j).eq.skip) go to 50
        nd0=ndist(i0,j0,i,j)
        if (nd0.lt.nd) then
          in=i
          jn=j
          nd=nd0
        end if
   50 continue
   60 continue
  100 continue

      return
      end
