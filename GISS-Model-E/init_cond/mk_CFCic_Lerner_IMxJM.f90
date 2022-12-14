!!!! mk_CFCic   7/2009     J.Lerner
!!!! Create first layer initial distribution for Lerner/Rind CFC
!!!!
      implicit none
      REAL*8, ALLOCATABLE, DIMENSION(:) :: LAT_DG
      REAL*4, ALLOCATABLE, DIMENSION(:,:) :: CFCic
      real*8 dlat_dg,fjeq,xlat,steps
      integer i,j,im,jm,k,jS,jN,nargs
      character*80 title,fileout
      character*4 cim,cjm

      NARGS = IARGC()
      IF(NARGS.lt.2)  then
        write(*,*) 'USAGE: mk_CFCic im jm'
        write(*,*) 'EXAMPLE: mk_CFCic 144 90'
        stop
      end if
      CALL GETARG (1,cim)
      read (cim,*) im
      CALL GETARG (2,cjm)
      read (cjm,*) jm
      cim = adjustl(cim)
      cjm = adjustl(cjm)
      fileout = 'CFCic_Lerner_'//trim(cim)//'x'//trim(cjm)
      write(*,*) fileout
      title = 'Initial conditions for CFC. IM, JM= '//cim//cjm

      allocate (lat_dg(jm),CFCic(im,jm))

!**** LATITUDES (degrees)
      DLAT_DG=180./REAL(JM)                   ! even spacing (default)
      IF (JM.eq.46) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole for 4x5
!c    IF (JM.eq.24) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole, orig 8x10
      IF (JM.eq.24) DLAT_DG=180./REAL(JM-1.5) ! 1/4 box at pole, 'real' 8x10

      FJEQ=.5*(1+JM)
      LAT_DG(1)=-90.
      LAT_DG(JM)=90.
      DO J=2,JM-1
        LAT_DG(J)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO

      DO J=1,JM
        if (lat_dg(j).ge.-20.) then
          jS = j
        exit
        end if
      enddo
      DO J=JM,1,-1
        if (lat_dg(j).le.20.) then
          jN = j
        exit
        end if
      enddo
      steps = jN-jS+1
      write(*,*) jn,js,steps

      DO I=1,im
      DO J=1,js-1
          CFCic(i,j) = 220.d-12*136.5/29.029
        enddo
        DO J=jN+1,jm
          CFCic(i,j) = 235.d-12*136.5/29.029
        enddo
        DO J=js,jn
          XLAT = (J-jS+.5)/steps
          if(i.eq.1) write(*,*) j,xlat
          CFCic(i,j) = (220.d-12 + XLAT*15.d-12)*136.5/29.029
        enddo
      enddo

      do j=1,jm
        write(*,*) j,lat_dg(j),CFCic(1,j)
      end do

      open (20,file=fileout,form='unformatted')
      write(20) title,CFCic
      stop
      end

!**** original code from modelII
!     DO 410 J=1,18                                                     
! 410 CFC11(J) = 220.E-12*136.5/29.029                                  
!     DO 420 J=29,46                                                    
! 420 CFC11(J) = 235.E-12*136.5/29.029                                  
!     DO 430 J=19,28                                                    
!     XLAT = (J-18.5)/10.                                               
! 430 CFC11(J) = (220.E-12 + XLAT*15.E-12)*136.5/29.029 


