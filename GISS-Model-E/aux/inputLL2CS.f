      program regrid_input
!@sum Input file regridding routines. Uses x_2grids derived type
!@    to define source and target grids.
!@auth Denis Gueyffier
      use regrid
      implicit none
      integer, parameter :: imt=90,jmt=90
      type (x_2grids) :: xll2cs_TOPO,
     &     xll2cs_OSST,xll2cs_GIC,xll2cs_SICE,
     &     xll2cs_SOIL,xll2cs_GLMELT,
     &     xll2cs_VEGFRAC,xll2cs_LAI,xll2cs_AIC,
     &     xll2cs_VEG,xll2cs_GIC2,
     &     xll2cs_4x5,xll2cs_2x2h,
     &     xll2cs_1x1q,xll2cs_1x1,xll2cs_halfdeg

      integer ::  ims_4x5, jms_4x5,
     &     ims_2x2h, jms_2x2h,
     &     ims_1x1q, jms_1x1q,
     &     ims_1x1, jms_1x1,
     &     ims_halfdeg, jms_halfdeg
      integer, parameter :: ntilessource=1,ntilestarget=6

      ims_4x5=72
      jms_4x5=46
      ims_2x2h=144
      jms_2x2h=90
      ims_1x1q=288
      jms_1x1q=180
      ims_1x1=360
      jms_1x1=180
      ims_halfdeg=720
      jms_halfdeg=360

      call init_regrid(xll2cs_4x5,ims_4x5,jms_4x5,
     &     ntilessource,imt,jmt,ntilestarget)
      call init_regrid(xll2cs_2x2h,ims_2x2h,jms_2x2h,
     &     ntilessource,imt,jmt,ntilestarget)
      call init_regrid(xll2cs_1x1q,ims_1x1q,jms_1x1q,
     &     ntilessource,imt,jmt,ntilestarget)
      call init_regrid(xll2cs_1x1,ims_1x1,jms_1x1,
     &     ntilessource,imt,jmt,ntilestarget)
      call init_regrid(xll2cs_halfdeg,ims_halfdeg,jms_halfdeg,
     &     ntilessource,imt,jmt,ntilestarget)

c     288x180 for SICE, OSST, TOPO, VEG
      xll2cs_TOPO=xll2cs_1x1q   
      xll2cs_OSST=xll2cs_1x1q
      xll2cs_SICE=xll2cs_1x1q
      xll2cs_VEG=xll2cs_1x1q

c     360x180 for GIC, AIC, GLMELT
      xll2cs_GIC=xll2cs_1x1
      xll2cs_AIC=xll2cs_1x1
      xll2cs_GLMELT=xll2cs_1x1 

c     144x90 for GIC2
      xll2cs_GIC2=xll2cs_2x2h

c      call regridTOPO(xll2cs_TOPO)
c      call regridOSST(xll2cs_OSST)
ccc      call testOSST()
c      call regridSICE(xll2cs_SICE)
c      call regridVEG(xll2cs_VEG)
c      call regridSOIL(xll2cs_halfdeg)
c      call regridSOIL(xll2cs_1x1)
c      call regridGLMELT(xll2cs_GLMELT)
c      call regridVEGFRAC_LAI(xll2cs_halfdeg)
c      call regridVEGFRAC_LAI(xll2cs_4x5)
c      call regridGIC(xll2cs_GIC,xll2cs_GIC2)
c      call regridAIC(xll2cs_AIC)
      call regridRDSCAL(xll2cs_halfdeg)

      end program regrid_input
c*

      subroutine regridTOPO(x2grids)
      use regrid
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE(nrecmax),name,ncfile
      type (x_2grids), intent(inout) :: x2grids
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      real*8, allocatable :: tsource(:,:,:,:)
      real*4, allocatable :: ttargr4(:,:,:,:)
      real*4, allocatable :: tsourc4(:,:,:,:)
      integer :: iu_TOPO,i,j,k,irec,iunit,imt,jmt,ntt,ims,jms,nts,
     &     maxrec,status,vid,fid,ir

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget
      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource


      write(*,*) "imt jmt ntt",imt,jmt,ntt

      allocate( ttargglob(imt,jmt,ntt,nrecmax),
     &     ttargr4(imt,jmt,ntt,nrecmax),
     &     ones(imt,jmt,ntt), 
     &     tsource(ims,jms,nts,nrecmax),
     &     tsourc4(ims,jms,nts,nrecmax)  )
      
      iu_TOPO=20
      name="Z1QX1N"

      write(*,*) name
      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "iu_TOPO",iu_TOPO
c      call read_regrid_4D_1R(x2grids,iu_TOPO,TITLE,ttargglob,maxrec)

      tsource(:,:,:,:)=0.
      ttargglob(:,:,:,:)=0.

      write(*,*) "ims, jms, nts=",ims,jms,nts
            
      irec=1

      do
         read(unit=iu_TOPO,END=30) TITLE(irec), tsourc4(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= tsourc4(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      write(*,*) "maxrec",maxrec
      
      do ir=1,maxrec
         call do_regrid(x2grids,tsource(:,:,:,ir),
     &        ttargglob(:,:,:,ir))
        write(*,*) "TITLE",TITLE(ir)
      enddo

c
c     CONSISTENCY CHECKS: 1) FOCEAN+FLAKE+FGRND+FGICE=1 
c                         2) IF FOCEAN(i,j) > 0 set FGRND=FGRND+FLAKE, FLAKE=0

      ones(:,:,:)=ttargglob(:,:,:,1)
     &     +ttargglob(:,:,:,2)
     &     +ttargglob(:,:,:,3)
     &     +ttargglob(:,:,:,4)
      
      if (any( abs(ones(1:imt,1:jmt,1:ntt) - 1.0) .gt. 1.d-6 ) ) then 
         write(*,*) "WARNING FOCEAN+FLAKE+FGRND+FGICE=1 BROKEN"
      endif

      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if ( (ttargglob(i,j,k,1) .gt. 0.) .and.
     &              (ttargglob(i,j,k,2) .gt. 0. ) ) then
                  write(*,*) "FGRND=FGRND+FLAKE, FLAKE=0"
                  ttargglob(i,j,k,3)=ttargglob(i,j,k,2)
     &                 +ttargglob(i,j,k,3)
                  ttargglob(i,j,k,2)=0.0
               endif
            enddo
         enddo
      enddo

      close(iu_TOPO)

      name="Z_CS90"

      write(*,*) name

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='unknown')

      ttargr4=ttargglob

      do ir=1,maxrec
         write(unit=iu_TOPO) TITLE(ir), ttargr4(:,:,:,ir)
         write(*,*) TITLE(ir)
      enddo

      close(iu_TOPO)


c      ncfile="topo6tiles.nc"

c      status = nf_open(trim(ncfile),nf_write,fid)
c      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN FILE"
c      status = nf_inq_varid(fid,'zatmo',vid)
c      write(*,*) NF_STRERROR(status)
c      status = nf_put_var_double(fid,vid,ttargglob(:,:,:,1))
c      write(*,*) "STATUS",NF_STRERROR(status),"<<"
c      status = nf_close(fid)

c      ncfile="topo.nc"

c      status = nf_open(trim(ncfile),nf_write,fid)
c      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN FILE"
c      status = nf_inq_varid(fid,'zatmo',vid)
c      write(*,*) NF_STRERROR(status)
c      status = nf_put_var_double(fid,vid,ttargglob(:,:,2,1))
c      write(*,*) "STATUS",NF_STRERROR(status),"<<"
c      status = nf_close(fid)

      write(*,*) "end topo"

      deallocate(ttargglob,ones,ttargr4,tsource)
   
      end subroutine regridTOPO 
c*

      subroutine regridOSST(x2grids)
c     Reto has extracted a 288x180 data set from Hadley data
c     this is compatible with focean mask in Gary's Z288x180 topo file
      use regrid
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 :: TITLE,name
      real*4 :: FOCEAN(x2grids%imsource,x2grids%jmsource)
      real*4 OSTmean(x2grids%imsource,x2grids%jmsource),
     &     OSTend(x2grids%imsource,x2grids%jmsource)
      integer iu_OSST,iu_TOPO,i,j,k
      real*8 :: missing

      if (x2grids%imtarget .eq. 32) then
         if (x2grids%imsource .eq. 144 .and. 
     &        x2grids%jmsource .eq. 90) then
c*    Read ocean fraction on input grid
      iu_TOPO=19
      name="Z144X90N_nocasp.1"
      open(iu_TOPO,FILE=name,FORM='unformatted', STATUS='old')
      read(iu_TOPO) title,FOCEAN
      close(iu_TOPO)

      iu_OSST=20
      name="OST_144x90.1876-1885avg.HadISST1.1"
c
      write(*,*) name
      missing=-9.999999171244748e33

      endif 
      endif

      if (x2grids%imtarget .eq. 90) then
         if (x2grids%imsource .eq. 288 .and. 
     &        x2grids%jmsource .eq. 180) then
c*    Read ocean fraction on input grid
      iu_TOPO=19
      name="Z288X180N"
      open(iu_TOPO,FILE=name,FORM='unformatted', STATUS='old')

      read(iu_TOPO) title,FOCEAN
      close(iu_TOPO)

      iu_OSST=20
      name="OST_288x180.1975-1984avg.HadISST1.1"

      write(*,*) name
      missing=-9.999999171244748e33

      endif 
      endif

      open(iu_OSST, FILE=name,FORM='unformatted', STATUS='old')

      call regrid_mask_OSST(x2grids,name,iu_OSST,FOCEAN,missing)

      close(iu_OSST)

      end subroutine regridOSST
c


      subroutine regridSICE(x2grids)
c     Reto has extracted a 288x180 data set from Hadley data
c     this is compatible with focean mask in Gary's Z288x180 topo file

      use regrid
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE,name,outunformat, TITLE2(nrecmax)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tin(:,:,:),tout(:,:,:)
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
      integer iu_SICE,i,j,k,ims,jms,nts,imt,jmt,ntt,iuout,
     &    maxrec,ir

      iu_SICE=19

      if (x2grids%imtarget .eq. 32) then
         if (x2grids%imsource .eq. 144 .and. 
     &        x2grids%jmsource .eq. 90) then
            name="SICE4X5.B.1876-85avg.Hadl1.1"
         endif
      endif   
      if (x2grids%imtarget .eq. 90) then
         if (x2grids%imsource .eq. 288 .and. 
     &        x2grids%jmsource .eq. 180) then
            name="SICE_288x180.1975-1984avg.HadISST1.1"    
         endif
      endif
      write(*,*) name

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),
     &     tin(ims,jms,nts),
     &     tout(imt,jmt,ntt))

      tsource(:,:,:)=0.0
      
      if (ims .eq. 72 .and. jms .eq. 46) then
      name="SICE4X5.B.1876-85avg.Hadl1.1"
      else if (ims .eq. 288 .and. jms .eq. 180) then
      name="SICE_288x180.1975-1984avg.HadISST1.1"
      endif
      write(*,*) name

      open (iu_SICE, FILE=name,FORM='unformatted', STATUS='old')

      read(unit=iu_SICE) TITLE, tin
      tsource=tin


      call do_regrid(x2grids,tsource(:,:,:),ttargglob)
      tout(:,:,:)=ttargglob(:,:,:)
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      write(unit=iuout) TITLE,tout(:,:,:)
      write(*,*) TITLE
      
      open( 33, FILE="sice1",
     &     FORM='unformatted', STATUS="UNKNOWN")

      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt),tbig(imt,jmt,2,ntt) )

      tsource1(:,:,:,:)=0.0
      tsource2(:,:,:,:)=0.0

      call read_recs_2R(tsource1,tsource2,iu_SICE,TITLE2,
     &     maxrec,ims,jms,nts)

      write(*,*) "maxrec",maxrec

      do ir=1,maxrec
         call do_regrid(x2grids,tsource1(:,:,:,ir),ttargglob1)
         tout1(:,:,:)=ttargglob1(:,:,:)
         call do_regrid(x2grids,tsource2(:,:,:,ir),ttargglob2)
         tout2(:,:,:)=ttargglob2(:,:,:)
         tbig(:,:,1,:)=tout1(:,:,:)
         tbig(:,:,2,:)=tout2(:,:,:)
         write(unit=iuout) TITLE2(ir),tbig
         write(33) TITLE2(ir),tout1(:,:,:)
         write(*,*) "TITLE",TITLE2(ir)
      enddo

      close(iuout)

      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2,
     &     tin,tout,tbig)

      close(iu_SICE)
      
      end subroutine regridSICE

      subroutine regridVEG(x2grids)
c     for 1x1 resolution: Jeff uses VEG=V360X180_no_crops.rep
c     It is identical to V144X90_no_crops.ext (144X90 data 
c     was just transfered to 360X180 grid without any change)
      use regrid
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_VEG
	
      iu_VEG=20
      name=" V288X180_no_crops.ext"

      open(iu_VEG,FILE=name,FORM='unformatted', STATUS='old')
      
      call read_regrid_write_veg(x2grids,name,iu_VEG)
           
      close(iu_VEG)
      
      end subroutine regridVEG
c*


      subroutine read_regrid_write_veg(x2grids,name,iuin)
      use regrid
      implicit none
      real*8,parameter :: pi = 3.1415926535897932d0 
      real*8,parameter :: twopi = 2d0*pi         
      real*8,parameter :: radian = pi/180d0      
      real*8,parameter :: radius = 6371000.
      type(x_2grids), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),arrsum(:,:,:),
     &     tout(:,:,:,:)
      real*4, allocatable :: t4(:,:,:),data(:,:,:)
      real*8 :: alpha
      real*8 :: totalfrac(10)
      real*8, allocatable :: axyp(:,:),
     &     fgroundLL(:,:),fgroundCS(:,:,:)
      real*4, allocatable :: fgroundLL4(:,:),fgroundCS4(:,:,:)
      character*80 :: TITLE(10),titleZ
      character*80 :: outunformat,outfile
      integer :: ir,ims,jms,nts,imt,jmt,ntt,i,j,k,l,iveg
      real*8 ::   FJEQ,DLAT_DG,DLAT,DLON,SINV,SINVm1,vsum
      real*8, allocatable :: dxyp(:)


      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt

      allocate (
     &     tsource(ims,jms,nts),data(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),tout(imt,jmt,10,ntt),
     &     t4(imt,jmt,ntt),arrsum(imt,jmt,ntt),
     &     dxyp(jms),axyp(imt,jmt),
     &     fgroundLL(ims,jms),fgroundCS(imt,jmt,ntt),
     &     fgroundLL4(ims,jms),fgroundCS4(imt,jmt,ntt)
     &     )
    
      arrsum=0.
            
      outunformat=trim(name)//".CS"
      write(*,*) outunformat
      
      iuout=40
      outfile="Z_CS90"
      open( iuout, FILE=outfile,
     &     FORM='unformatted', STATUS="UNKNOWN")
      write(*,*) outfile
      read(iuout) titleZ,fgroundCS4
      read(iuout) titleZ,fgroundCS4
      read(iuout) titleZ,fgroundCS4
      fgroundCS=fgroundcs4
      write(*,*) titleZ
      close(iuout)
      
      iuout=50
      outfile="Z288X180N"
      open( iuout, FILE=outfile,
     &     FORM='unformatted', STATUS="UNKNOWN")
      write(*,*) outfile
      read(iuout) titleZ,fgroundLL4
      read(iuout) titleZ,fgroundLL4
      read(iuout) titleZ,fgroundLL4
      fgroundll=fgroundll4
      write(*,*) titleZ
      close(iuout)
      
c     computing area of latlon cells 
      FJEQ=.5*(1+jms)
      DLAT_DG=180./REAL(jms)    ! even spacing (default)
      DLAT=DLAT_DG*radian
      DLON=TWOPI/REAL(ims)
      
      SINV    = Sin (DLAT*(1+.5-FJEQ))
      DXYP(1) = RADIUS*RADIUS*DLON*(SINV+1)
      
      SINVm1  = Sin (DLAT*(jms-.5-FJEQ))
      DXYP(jms)= RADIUS*RADIUS*DLON*(1-SINVm1)
      
      DO J=2,jms-1
         SINVm1  = Sin (DLAT*(J-.5-FJEQ))
         SINV    = Sin (DLAT*(J+.5-FJEQ))
         DXYP(J) = RADIUS*RADIUS*DLON*(SINV-SINVm1)
c         write(*,*) "dxyp=",DXYP(J)
      END DO
c     end area
      
      do ir=1,10
         read(unit=iuin) TITLE(ir),data
         write(*,*) "TITLE, ir",TITLE(ir),ir
         tsource(1:ims,1:jms,1)=data(1:ims,1:jms,1)*
     &         fgroundll(1:ims,1:jms)
        
  
         totalfrac(ir)=0.
         do j=1,jms
            do i=1,ims
               totalfrac(ir)=totalfrac(ir)
     &              +tsource(i,j,1)*dxyp(j)
c     &                 *fgroundll(i,j)
            enddo
         enddo
         WRITE(*,*) "ll tot surface vetype ",ir," =",totalfrac(ir)
         
         call do_regrid(x2grids,tsource,ttargglob)
         tout(:,:,ir,:)=ttargglob(:,:,:)
      enddo

c     divide by sum of veg. fractions
       do k=1,ntt
         do j=1,jmt
            do i=1,imt
               vsum=(tout(i,j,1,k)+
     &                 tout(i,j,2,k)+tout(i,j,3,k)+
     &                 tout(i,j,4,k)+tout(i,j,5,k)+
     &                 tout(i,j,6,k)+tout(i,j,7,k)+
     &                 tout(i,j,8,k)+tout(i,j,9,k)+
     &                 tout(i,j,10,k)  ) 
               if (vsum .gt. 0.) then
                 tout(i,j,1:10,k)=tout(i,j,1:10,k)/vsum
               else
                  if (fgroundCS(i,j,k) .gt. 0.)  
     &           write(*,*)  "WARNING!!!!!!!!",i,j,k,fgroundCS(i,j,k)
               endif
              
            enddo
         enddo
       enddo

      iuout=iuout+1
      
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      
c     check bare soil fraction if v(1)+v(10) < 0.01 set v(1)=v(10)=0, 
c     if v(1)+v(10) > 0.99, set v(1)+v(10)=1, v(2)=...=v(9)=0
      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if (tout(i,j,1,k)+tout(i,j,10,k) .le. 0.01 .and.
     &              tout(i,j,1,k)+tout(i,j,10,k) .gt. 0. ) then
c                  write(*,*) " v(1)+v(10) < 0.01 at",i,j,k
                  alpha=1./(1.-(tout(i,j,1,k)+tout(i,j,10,k)) )
                  tout(i,j,1,k)  = 0.
                  tout(i,j,10,k) = 0.
                  do l=2,9
                     tout(i,j,l,k)=alpha*tout(i,j,l,k)
                  enddo
               endif
               if (tout(i,j,1,k)+tout(i,j,10,k) .ge. 0.99) then
c                  write(*,*) " v(1)+v(10) > 0.95 at",i,j,k
                  alpha=1./(1.- (
     &                 tout(i,j,2,k)+tout(i,j,3,k)+
     &                 tout(i,j,4,k)+tout(i,j,5,k)+
     &                 tout(i,j,6,k)+tout(i,j,7,k)+
     &                 tout(i,j,8,k)+tout(i,j,9,k) ) )
                  tout(i,j,1,k)  = alpha*tout(i,j,1,k)
                  tout(i,j,10,k) = alpha*tout(i,j,10,k) 
                  do l=2,9
                     tout(i,j,l,k) = 0.
                  enddo
               endif
            enddo
         enddo
      enddo
      
c     fix all fractions
      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               do iveg=1,10
                  if (tout(i,j,iveg,k) .le. 0.01 .and.
     &                 tout(i,j,iveg,k) .gt. 0. ) then
c                     write(*,*) " >>>v(",iveg,") < 0.01 at",i,j,k
                     alpha=1./( 1.- tout(i,j,iveg,k) )
                     do l=1,10
                        tout(i,j,l,k)=alpha*tout(i,j,l,k)
                     enddo
                     tout(i,j,iveg,k)  = 0.
                  endif
               enddo
            enddo
         enddo
      enddo
c     
      
      do ir=1,10
         arrsum(:,:,:)=arrsum(:,:,:)+tout(:,:,ir,:) 
         t4(:,:,:)=tout(:,:,ir,:) 
         write(unit=iuout) TITLE(ir),t4
      enddo
      
      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if (abs( arrsum(i,j,k) - 1.0) .gt. 0.0001) then
c                  write(*,*) "arrsum(",i,",",j,",",k,")=",
c     &                 arrsum(i,j,k)
               endif
            enddo
         enddo
      enddo
      
      call areacs(axyp,imt,jmt)

      do iveg=1,10
         totalfrac(iveg)=0.
         do k=1,ntt
            do j=1,jmt
               do i=1,imt
                  totalfrac(iveg)=totalfrac(iveg)
     &                 +tout(i,j,iveg,k)*axyp(i,j)
     &                    *fgroundCS(i,j,k)
               enddo
            enddo
         enddo
         WRITE(*,*) "TOTAL surface covered by veg. of type ",
     &        iveg," =", totalfrac(iveg)
         
      enddo
      
      close(iuout) 
      
      deallocate(tsource,ttargglob,tout,t4,data,arrsum,axyp,
     &     dxyp,fgroundCS,fgroundLL)
      
      end subroutine read_regrid_write_veg
c*

      
      subroutine regridSOIL(x2grids)
c
c     for 1x1 resolution: Jeff uses SOIL=S360X180_0098M.rep
c
      use regrid
      implicit none
      type (x_2grids), intent(in) :: x2grids
      real*4, allocatable :: dz(:,:,:),ftext(:,:,:,:),
     &     ftextk(:,:,:,:),sl(:,:)
      real*4, allocatable :: dzout(:,:,:,:),ftextout(:,:,:,:,:),
     &     ftextkout(:,:,:,:,:),slout(:,:,:),bigarrout(:,:,:,:)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLE,name,outunformat
      integer iu_SOIL,iuout,ims,jms,nts,imt,jmt,ntt,i,j,k,l,inds

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

      iu_SOIL=20
c      name="S360X180_0098M.rep"
      name="SOIL720X360_Reynolds"
      
      open(iu_SOIL,FILE=name,FORM='unformatted', STATUS='old')
      
      allocate (dz(ims,jms,6),ftext(ims,jms,5,6),
     &     ftextk(ims,jms,5,6),sl(ims,jms),
     &     dzout(imt,jmt,6,ntt),ftextout(imt,jmt,5,6,ntt),
     &     ftextkout(imt,jmt,5,6,ntt),
     &     slout(imt,jmt,ntt),bigarrout(imt,jmt,67,ntt) )

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt) )

      read(iu_SOIL) dz,ftext,ftextk,sl

      close(iu_SOIL)
            
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat

c      write(*,*) dz
      iuout=iu_SOIL+1
      
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      
      do k=1,6
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=dz(i,j,k)
            enddo
         enddo
         call do_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         dzout(:,:,k,:)=ttargglob(:,:,:)
      enddo

      do k=1,6
         do l=1,5
            tsource(:,:,1)=ftext(:,:,l,k)
            call do_regrid(x2grids,tsource,ttargglob)
            ftextout(:,:,l,k,:)=ttargglob(:,:,:)
         enddo
      enddo


      do k=1,6
         do l=1,5
            tsource(:,:,1)=ftextk(:,:,l,k)
            call do_regrid(x2grids,tsource,ttargglob)
            ftextkout(:,:,l,k,:)=ttargglob(:,:,:)
         enddo
      enddo

      tsource(:,:,1)=sl(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      slout(:,:,:)=ttargglob(:,:,:)
      
      do k=1,6
         bigarrout(:,:,k,:)=dzout(:,:,k,:)
      enddo

      inds=6
      do k=1,6
         do l=1,5
            bigarrout(:,:,inds+l+5*(k-1),:)=ftextout(:,:,l,k,:)
         enddo
      enddo

      inds=inds+30
      do k=1,6
         do l=1,5
            bigarrout(:,:,inds+l+5*(k-1),:)=ftextkout(:,:,l,k,:)
         enddo
      enddo

      inds=inds+30
      bigarrout(:,:,inds+1,:)=slout(:,:,:)

      write(unit=iuout) bigarrout
         
      close(iuout)
      

      end subroutine regridSOIL
c*
      
      

      subroutine regridGLMELT(x2grids)
c     Gavin has made a 1x1 file ascii file GLMELT_360X180.OCN
c     convert it to binary using GLMELTt2b.f
      use regrid
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_GLMELT,imt,jmt

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      
 
      name="GLMELT_360X180.OCN.bin"

      iu_GLMELT=20
      
      open (iu_GLMELT,FILE=name,FORM='unformatted', STATUS='old')

      call read_regrid_write_GLMELT(x2grids,name,iu_GLMELT)

      close(iu_GLMELT)
      
      end subroutine regridGLMELT
c*



      subroutine read_regrid_write_GLMELT(x2grids,name,iuin)
      use regrid
      implicit none
      type(x_2grids), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,itile,i,j

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:,:)=0.0

      
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)

      write(*,*) "maxrec",maxrec

      outunformat=trim(name)//".CS"

      write(*,*) outunformat

      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      
      do ir=1,maxrec
         call do_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         do itile=1,ntt
           do i=1,imt
             do j=1,jmt
             if (ttargglob(i,j,itile) .le. 0.) then
               tout(i,j,itile)=0.
             else
               tout(i,j,itile)=1.
             endif
            enddo
           enddo
          enddo
      
          write(unit=iuout) TITLE(ir),tout(:,:,:)
          write(*,*) "TITLE",TITLE(ir)
      
      enddo

      close(iuout)

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_GLMELT
c*



      subroutine regridVEGFRAC_LAI(x2grids)
c
c     regriding vegetation fraction and LAI used by tracers code
c
      use regrid
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE,name,oname,TITLEFILE
      type (x_2grids), intent(inout) :: x2grids
      real*8, allocatable :: vfraccs(:,:,:,:),ones(:,:,:),
     &     laiglob(:,:,:,:),real8ll(:,:,:),vfracll(:,:,:),
     &     laics(:,:,:)
      real*4, allocatable :: real4ll(:,:)
      integer :: iu_VEGFRAC,iu_VEGFRACCS,i,j,k,irec,iunit,imt,jmt,
     &     ntt,maxrec,status,vid,fid,ir,imlon,jmlat
      integer :: iu_LAI,iu_LAICS,imonth,ios
      character(len=2) :: c2month 
      character(len=1) :: c1month 

      imlon = x2grids%imsource
      jmlat = x2grids%jmsource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "imt jmt ntt",imt,jmt,ntt

      allocate(
     &     vfracll(imlon,jmlat,nrecmax),
     &     vfraccs(imt,jmt,ntt,nrecmax),
     &     real8ll(imlon,jmlat,1),
     &     laics(imt,jmt,ntt),
     &     real4ll(imlon,jmlat),
     &     ones(imt,jmt,ntt) )

      vfraccs(:,:,:,:) = 0.0
     
      iu_VEGFRAC=20
      name="vegtype.global.bin"
      write(*,*) name
      open( iu_VEGFRAC, FILE=name,FORM='unformatted', STATUS='old')
      oname=trim(name)//".CS"
      write(*,*) oname
      open( iu_VEGFRACCS, FILE=oname,FORM='unformatted', 
     &     STATUS='unknown')

      ir = 0
      do
        read(iu_VEGFRAC,iostat=ios) title,real4ll
        if(ios.ne.0) exit
        ir = ir + 1
        vfracll(:,:,ir) = real4ll
        real8ll(:,:,1) = real4ll
        call do_regrid(x2grids,real8ll,vfraccs(:,:,:,ir))
        write(iu_VEGFRACCS) TITLE,real(vfraccs(:,:,:,ir),kind=4)
      enddo
      maxrec = ir
      close(iu_VEGFRAC)
      close(iu_VEGFRACCS)
      
c***  
c***  CONSISTENCY CHECKS: 1) sum(fractions)=1 
c***  

      ones=sum(vfraccs(1:imt,1:jmt,1:ntt,1:maxrec),4)
      
      if (any( abs(ones(1:imt,1:jmt,1:ntt) - 1.0) .gt. 1.d-3 ) ) then 
         write(*,*) "WARNING FRACTIONS DO NOT ADD UP TO 1 
     &        EVERYWHERE ON CUBED SPHERE"
         write(*,*) "ones=",ones
      endif

      iu_LAI=21
      do imonth=1,12
         write(c2month,'(i2.2)') imonth
         name='lai'//c2month//'.global.bin'
         write(*,*) name
         oname=trim(name)//".CS"
         open( iu_LAI, FILE=name,FORM='unformatted', STATUS='old')
         open( iu_LAICS, FILE=oname,FORM='unformatted', 
     &        STATUS='unknown')
         do ir=1,maxrec
           read(iu_LAI) TITLE, real4ll
           real8ll(:,:,1) = real4ll*vfracll(:,:,ir) ! multiply lai by weight
           call do_regrid(x2grids,real8ll,laics)
           where(vfraccs(:,:,:,ir).gt.0.)
             laics = laics / vfraccs(:,:,:,ir)  ! divide by regridded weight
           end where
           write(iu_LAICS) TITLE, real(laics,kind=4)
         enddo
         close(iu_LAI)
         close(iu_LAICS)
      enddo

      deallocate(vfracll,vfraccs,laics,ones,real4ll,real8ll)
   
      end subroutine regridVEGFRAC_LAI
c*


      subroutine regridGIC(x2grids,x2grids2)
c
c     We use Jeff's 1x1 file GIC=GIC.360X180.DEC01.1.rep
c     and replace the ice initial conditions
c  SICE02         R8 F(im,jm),H(4,im,jm),snw,msi,ssi(4),pond_melt,L flag_dsws
c     by Larissa's E42F40oQ32 run in DEC1971.iceE44F40oQ32 
      use regrid
      use ncio, only : defvar, write_data

      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids,x2grids2

c*    read
      real*8, allocatable :: Tocn(:,:,:),MixLD(:,:)        ! OCN01
      real*8, allocatable :: F0(:,:),H0(:,:,:),snw0(:,:),msi0(:,:),
     &     ssi0(:,:,:),pond_melt0(:,:)
      logical, allocatable :: flag_dsws0(:,:)               ! SICE02
      real*8, allocatable :: F(:,:),H(:,:,:),snw(:,:),msi(:,:),
     &     ssi(:,:,:),pond_melt(:,:)
      logical, allocatable :: flag_dsws(:,:)               ! SICE02
      real*8, allocatable :: snowe(:,:),Te(:,:),WTRe(:,:),ICEe(:,:),
     &     SNOage(:,:,:),evmax(:,:),fsat(:,:),gq(:,:)      ! EARTH01
      real*8, allocatable :: W(:,:,:,:),HT(:,:,:,:),SNWbv(:,:,:) ! SOILS03
      real*8, allocatable :: SNOW(:,:),T(:,:,:),
     &     MDWN(:,:),EDWN(:,:)
      real*8 ::  ACCPDA, ACCPDG,  EACCPDA, EACCPDG           !GLAI

c*    write
      real*8, allocatable :: Tocn_out(:,:,:,:),MixLD_out(:,:,:)   ! OCN01
      real*8, allocatable :: F_out(:,:,:),H_out(:,:,:,:),
     &     snw_out(:,:,:),msi_out(:,:,:),ssi_out(:,:,:,:),
     &     pond_melt_out(:,:,:), flag_dsws_out(:,:,:)                ! SICE02
      real*8, allocatable :: snowe_out(:,:,:),Te_out(:,:,:),
     &     WTRe_out(:,:,:), ICEe_out(:,:,:),SNOage_out(:,:,:,:),
     &     evmax_out(:,:,:),fsat_out(:,:,:),gq_out(:,:,:)         ! EARTH01
      real*8, allocatable :: W_out(:,:,:,:,:),
     &     HT_out(:,:,:,:,:),SNWbv_out(:,:,:,:)                   ! SOILS02
      real*8, allocatable :: SNOW_out(:,:,:),T_out(:,:,:,:),
     &     MDWN_out(:,:,:),EDWN_out(:,:,:)        !GLAI
      real*8, allocatable :: tsource(:,:,:),tsource2(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLEOCN01,TITLESICE02,TITLEEARTH01,TITLESOILS03,
     &     TITLEGLAIC01,name,outnc,name2
      integer :: iu_GIC,iuout,ims,jms,nts,imt,jmt,ntt,ims2,jms2,nts2,
     &     iu_ICE
      integer :: i,j,k,l,m,fid,status,ntiles,im,jm,d2,d3
     

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget
      ims2=x2grids2%imsource
      jms2=x2grids2%jmsource
      nts2=x2grids2%ntilessource
      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt
      write(*,*) "ims2,jms2,nts2",ims2,jms2,nts2

      name="GIC.360X180.DEC01.1.rep"
      name2="1DEC1971.iceE44F40oQ32"
      iu_GIC=20
      iu_ICE=30

      allocate (Tocn(3,ims,jms),MixLD(ims,jms),
     &  F0(ims,jms),H0(4,ims,jms),snw0(ims,jms),msi0(ims,jms),
     &  ssi0(4,ims,jms),pond_melt0(ims,jms),flag_dsws0(ims,jms),
     &  F(ims2,jms2),H(4,ims2,jms2),snw(ims2,jms2),msi(ims2,jms2),
     &  ssi(4,ims2,jms2),pond_melt(ims2,jms2),flag_dsws(ims2,jms2),
     &  snowe(ims,jms),Te(ims,jms),WTRe(ims,jms),ICEe(ims,jms), 
     &  SNOage(3,ims,jms),evmax(ims,jms),fsat(ims,jms),gq(ims,jms),
     &  W(7,3,ims,jms),HT(7,3,ims,jms),
     &  SNWbv(2,ims,jms),
     &  SNOW(ims,jms),T(2,ims,jms),
     &  MDWN(ims,jms),EDWN(ims,jms))

      allocate (Tocn_out(3,imt,jmt,ntt),MixLD_out(imt,jmt,ntt),
     &     F_out(imt,jmt,ntt),H_out(4,imt,jmt,ntt),
     &     snw_out(imt,jmt,ntt),msi_out(imt,jmt,ntt),
     &     ssi_out(4,imt,jmt,ntt),pond_melt_out(imt,jmt,ntt),
     &     flag_dsws_out(imt,jmt,ntt),
     &     snowe_out(imt,jmt,ntt),Te_out(imt,jmt,ntt),
     &     WTRe_out(imt,jmt,ntt),ICEe_out(imt,jmt,ntt), 
     &     SNOage_out(3,imt,jmt,ntt),evmax_out(imt,jmt,ntt),
     &     fsat_out(imt,jmt,ntt),gq_out(imt,jmt,ntt),
     &     W_out(7,3,imt,jmt,ntt),
     &     HT_out(7,3,imt,jmt,ntt),
     &     SNWbv_out(3,imt,jmt,ntt),
     &     SNOW_out(imt,jmt,ntt),T_out(2,imt,jmt,ntt),
     &     MDWN_out(imt,jmt,ntt),EDWN_out(imt,jmt,ntt))

      allocate (tsource(ims,jms,nts),
     &     tsource2(ims2,jms2,nts2),
     &     ttargglob(imt,jmt,ntt) )


      open(iu_GIC,FILE=name,FORM='unformatted', STATUS='old')

      read(iu_GIC) TITLEOCN01, Tocn,MixLD
      write(*,*) TITLEOCN01
      read(iu_GIC) 
     &  TITLESICE02,F0,H0,snw0,msi0,ssi0,pond_melt0,flag_dsws0
      write(*,*) TITLESICE02
      read(iu_GIC) TITLEEARTH01, snowe,Te,WTRe, ICEe, SNOage,evmax,
     &     fsat,gq
      write(*,*) TITLEEARTH01
      read(iu_GIC) TITLESOILS03, W,HT,SNWbv
      write(*,*) TITLESOILS03
      read(iu_GIC) TITLEGLAIC01, SNOW,T,MDWN,EDWN,
     &     ACCPDA,ACCPDG,EACCPDA,EACCPDG
      write(*,*) TITLEGLAIC01

      close(iu_GIC)

      open(iu_ICE,FILE=name2,FORM='unformatted', STATUS='old')
      read(iu_ICE) TITLESICE02,F,H,snw,msi,ssi,pond_melt,flag_dsws
      write(*,*) TITLESICE02
      close(iu_ICE)
            
      do k=1,3
         tsource(:,:,1)=Tocn(k,:,:)
         call do_regrid(x2grids,tsource,ttargglob)
         Tocn_out(k,:,:,:)=ttargglob(:,:,:)
      enddo
     
      tsource(:,:,1)=MixLD(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      MixLD_out(:,:,:)=ttargglob(:,:,:)

      tsource2(:,:,1)=F(:,:)

      call do_regrid(x2grids2,tsource2,ttargglob)
      F_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         tsource2(:,:,1)=H(k,:,:)
         call do_regrid(x2grids2,tsource2,ttargglob)
         H_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource2(:,:,1)=snw(:,:)
      call do_regrid(x2grids2,tsource2,ttargglob)
      snw_out(:,:,:)=ttargglob(:,:,:)

      tsource2(:,:,1)=msi(:,:)
      call do_regrid(x2grids2,tsource2,ttargglob)
      msi_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         tsource2(:,:,1)=ssi(k,:,:)
         call do_regrid(x2grids2,tsource2,ttargglob)
         ssi_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource2(:,:,1)=pond_melt(:,:)
      call do_regrid(x2grids2,tsource2,ttargglob)
      pond_melt_out(:,:,:)=ttargglob(:,:,:)

      tsource2(:,:,1) = 0.d0

      do j=1,jms2
         do i=1,ims2
            if (flag_dsws(i,j)) then
               tsource2(i,j,1)=1.d0
            else
               tsource2(i,j,1)=0.d0
            endif
         enddo
      enddo

      call do_regrid(x2grids2,tsource2,ttargglob)


c***  Compatibility
      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if (ttargglob(i,j,k) .ge. 0.5d0) then
                  flag_dsws_out(i,j,k)=1.
               else
                  flag_dsws_out(i,j,k)=0.
               endif
            enddo
         enddo
      enddo

      write(*,*) "HERE"

      tsource(:,:,1)=snowe(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      snowe_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=Te(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      Te_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=WTRe(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      WTRe_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=ICEe(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      ICEe_out(:,:,:)=ttargglob(:,:,:)

      do k=1,3
         tsource(:,:,1)=SNOage(k,:,:)
         call do_regrid(x2grids,tsource,ttargglob)
         SNOage_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=evmax(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      evmax_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=fsat(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      fsat_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=gq(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      gq_out(:,:,:)=ttargglob(:,:,:)
      
      
      do k=1,7
         do m=1,3
            tsource(:,:,1)=W(k,m,:,:)
            call do_regrid(x2grids,tsource,ttargglob)
            W_out(k,m,:,:,:)=ttargglob(:,:,:)
         enddo
      enddo

      do k=1,7
         do m=1,3
            tsource(:,:,1)=HT(k,m,:,:)
            call do_regrid(x2grids,tsource,ttargglob)
            HT_out(k,m,:,:,:)=ttargglob(:,:,:)
         enddo
      enddo

      do k=1,2
         tsource(:,:,1)=SNWbv(k,:,:)
         call do_regrid(x2grids,tsource,ttargglob)
         SNWbv_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      SNWbv_out(3,:,:,:) = 0.

      tsource(:,:,1)=SNOW(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      SNOW_out(:,:,:)=ttargglob(:,:,:)

      do k=1,2
         tsource(:,:,1)=T(k,:,:)
         call do_regrid(x2grids,tsource,ttargglob)
         T_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=MDWN(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      MDWN_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=EDWN(:,:)
      call do_regrid(x2grids,tsource,ttargglob)
      EDWN_out(:,:,:)=ttargglob(:,:,:)
      

c***  Write Netcdf file
#ifdef TRACERS_WATER
      write(*,*) "STOP TRACERS WATER NOT IMPLEMENTED IN regridinput"
      stop
#endif      
      

      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc
      
      status = nf_create(outnc,nf_clobber,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      
    
c***  Define OCN variables
      call defvar(fid,imt,jmt,ntt,Tocn_out,'tocean(d3,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,MixLD_out,'z1o(im,jm,tile)'
     &     ,.true.)
c***  Define SICE variables
      call defvar(fid,imt,jmt,ntt,F_out,'rsi(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,H_out,'hsi(lmi,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,snw_out,'snowi(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,msi_out,'msi(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,ssi_out,'ssi(lmi,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,pond_melt_out,
     &     'pond_melt(im,jm,tile)',.true.)
      call defvar(fid,imt,jmt,ntt,flag_dsws_out,
     &     'flag_dsws(im,jm,tile)',.true.)
c***  Define EARTH variables
      call defvar(fid,imt,jmt,ntt,snowe_out,'snowe(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,Te_out,'tearth(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,WTRe_out,'wearth(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,ICEe_out,'aiearth(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,SNOage_out,'snoage(d3,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,evmax_out,
     &     'evap_max_ij(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,fsat_out,'fr_sat_ij(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,gq_out,'qg_ij(im,jm,tile)'
     &     ,.true.)
c***  Define SOIL variables     ! this is the old SOIL02 version, implement the new version
      call defvar(fid,imt,jmt,ntt,W_out,
     &     'w_ij(zero_to_ngm,ls_nfrac,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,HT_out,
     &     'ht_ij(zero_to_ngm,ls_nfrac,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,SNWbv_out,
     &     'snowbv(ls_nfrac,im,jm,tile)'
     &     ,.true.)
c***  Define GLAIC variables
      call defvar(fid,imt,jmt,ntt,SNOW_out,'snowli(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,T_out,'tlandi(d2,im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,MDWN_out,'mdwnimp(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,EDWN_out,'edwnimp(im,jm,tile)'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,ACCPDA,'accpda'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,ACCPDG,'accpdg'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,EACCPDA,'eaccpda'
     &     ,.true.)
      call defvar(fid,imt,jmt,ntt,EACCPDG,'eaccpdg'
     &     ,.true.)

      status = nf_enddef(fid)
      if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"

      write(*,*) "finished defvar"

c***  Write OCN variables
      call write_data(fid,'tocean',Tocn_out)
      call write_data(fid,'z1o',MixLD_out)
      write(*,*) "finished writing OCN variables"
c***  Write SICE variables
      call write_data(fid,'rsi',F_out)
      call write_data(fid,'hsi',H_out)
      call write_data(fid,'snowi',snw_out)
      call write_data(fid,'msi',msi_out)
      call write_data(fid,'ssi',ssi_out)
      call write_data(fid,'pond_melt',pond_melt_out)
      call write_data(fid,'flag_dsws',flag_dsws_out)
      write(*,*) "finished writing SICE variables"
c***  Write EARTH variables
      call write_data(fid,'snowe',snowe_out)
      call write_data(fid,'tearth',Te_out)
      call write_data(fid,'wearth',WTRe_out)
      call write_data(fid,'aiearth',ICEe_out)
      call write_data(fid,'snoage',SNOage_out)
      call write_data(fid,'evap_max_ij',evmax_out)
      call write_data(fid,'fr_sat_ij',fsat_out)
      call write_data(fid,'qg_ij',gq_out)
      write(*,*) "finished writing EARTH variables"
c***  Write SOIL variables     ! this is the old SOIL02 version, implement the new version
      call write_data(fid,'w_ij',W_out)
      call write_data(fid,'ht_ij',HT_out)
      call write_data(fid,'snowbv',SNWbv_out)
      write(*,*) "finished writing SOIL variables"
c***  Write GLAIC variables
      call write_data(fid,'snowli',SNOW_out)
      call write_data(fid,'tlandi',T_out)
      call write_data(fid,'mdwnimp',MDWN_out)
      call write_data(fid,'edwnimp',EDWN_out)
      call write_data(fid,'accpda',ACCPDA)
      call write_data(fid,'accpdg',ACCPDG)
      call write_data(fid,'eaccpda',EACCPDA)
      call write_data(fid,'eaccpdg',EACCPDG)
      write(*,*) "finished writing GLAIC variables"

      deallocate (Tocn_out,MixLD_out,F_out,H_out,
     &     snw_out,msi_out,ssi_out,pond_melt_out,
     &     flag_dsws_out,snowe_out,Te_out,
     &     WTRe_out,ICEe_out,SNOage_out,evmax_out,
     &     fsat_out,gq_out,W_out,
     &     HT_out,SNWbv_out,
     &     SNOW_out,T_out,MDWN_out,EDWN_out )

      status = nf_close(fid)
      deallocate (Tocn,MixLD,
     &     F0,H0,snw0,msi0,
     &     ssi0,pond_melt0,flag_dsws0)

      deallocate (F,H,snw,msi,
     &     ssi,pond_melt)

      deallocate(flag_dsws)

      deallocate (snowe,Te,WTRe,ICEe, 
     &     SNOage,evmax,fsat,gq,
     &     W,HT,SNWbv,SNOW,T,MDWN,EDWN )

      deallocate (tsource,tsource2,ttargglob)     
      
      write(*,*) "end regrid GIC"

      end subroutine regridGIC
c*



      subroutine regridAIC(x2grids)
c
c     for 1x1 resolution : Jeff uses AIC=AIC.RES_X40.D771201N.rep
c
      use regrid
      use ncio, only : defvar,write_data
      implicit none
      include 'netcdf.inc'
      type(x_2grids), intent(in) :: x2grids
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:)
      real*4, allocatable :: ts4(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),tcopy(:,:,:),
     &     tcopy2(:,:,:),tcopy3(:,:,:),tcopy4(:,:,:)
      character*80 :: TITLE,outunformat,outnc,name
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,
     &     status,fid,vid
      integer iu_AIC

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt r8",
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts),
     &     ts4(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),tcopy(imt,jmt,ntt),
     &     tcopy2(imt,jmt,ntt), 
     &     tcopy3(imt,jmt,ntt), 
     &     tcopy4(imt,jmt,ntt))

      tsource(:,:,:)=0.0

      iu_AIC=20

      if (ims .eq. 72 .and. jms .eq. 46) then
         name="AIC.RES_M20A.D771201"
      elseif (ims .eq. 360 .and. jms .eq. 180) then
         name="AIC.RES_X40.D771201N.rep"
      endif

      open(iu_AIC,FILE=name,FORM='unformatted', STATUS='old')

      outunformat=trim(name)//".CS"
      write(*,*) outunformat

      iuout=21
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      irec=1

      write(*,*) "irec=",irec
      do
         read(unit=iu_AIC,END=30) TITLE, ts4
         tsource=ts4
         write(*,*) "TITLE, irec",TITLE,irec
         call do_regrid(x2grids,tsource,ttargglob)

         write(unit=iuout) TITLE,real(ttargglob,KIND=4)
         irec=irec+1
      enddo
   
 30   continue

      maxrec=irec-1

      write(*,*) "maxrec",maxrec

      close(iuout) 
      close(iu_AIC)


c      write(*,*) "here w r4"
c
c      outnc=trim(name)//"-CS.nc"
c      write(*,*) outnc
c
c         write(*,*) "TCOPY=",tcopy
c         status = nf_create(outnc,nf_clobber,fid)
c         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
c
c      call defvar(fid,imt,jmt,ntt,tcopy,'press(im,jm,tile)')
c      call defvar(fid,imt,jmt,ntt,tcopy2,'surf_temp(im,jm,tile)')
c      call defvar(fid,imt,jmt,ntt,tcopy3,'u974(im,jm,tile)')
c      call defvar(fid,imt,jmt,ntt,tcopy4,'v974(im,jm,tile)')
c
c     status = nf_enddef(fid)
c     if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
c
c      call write_data(fid,'press',tcopy)
c      call write_data(fid,'surf_temp',tcopy2)
c      call write_data(fid,'u974',tcopy3)
c      call write_data(fid,'v974',tcopy4)
c
c      status = nf_close(fid)

      deallocate(tsource,ts4,ttargglob,tcopy,tcopy2,
     &     tcopy3,tcopy4)
 
      end subroutine regridAIC
c*


      subroutine read_recs_1R(tsource,iuin,TITLE,maxrec,im1,jm1,ntl)
      use regrid
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*4 :: data(im1,jm1,ntl,nrecmax)
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      write(*,*) "iuin",iuin

      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), data(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= data(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)


      end subroutine read_recs_1R
c*

      subroutine read_recs_1R_r4_r8(tsource,iuin,TITLE,maxrec,
     *     im1,jm1,ntl)
      use regrid
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      real*4 :: ts4(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec


      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), ts4(:,:,:,irec)
         tsource(:,:,:,irec)=ts4(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)


      end subroutine read_recs_1R_r4_r8
c*


      subroutine read_recs_1R_r8(tsource,iuin,TITLE,maxrec,im1,jm1,
     *     ntl)
      use regrid
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      
      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), tsource(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_1R_r8
c*


      subroutine read_recs_2R(tsource1,tsource2,iuin,
     &     TITLE,maxrec,im1,jm1,ntl)
      use regrid
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*4 :: data1(im1,jm1,ntl,nrecmax),
     &     data2(im1,jm1,ntl,nrecmax)
      real*8, intent(inout) :: tsource1(im1,jm1,ntl,nrecmax),
     &     tsource2(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), data1(:,:,:,irec),
     &        data2(:,:,:,irec)

         tsource1(:,:,:,irec)= data1(:,:,:,irec)
         tsource2(:,:,:,irec)= data2(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_2R
c*


      subroutine read_regrid_write_4D_1R(x2grids,name,iuin)
      use regrid
      implicit none
      type(x_2grids), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:)
      real*4, allocatable :: tsourc4(:,:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt

      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt),tsourc4(ims,jms,nts,nrecmax))
      tsource(:,:,:,:)=0.0

      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec),tsourc4(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)=tsourc4(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat

      iuout=iuin+1

      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call do_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         tout=ttargglob
         write(*,*) "test=",tout(58,13,3)
         write(unit=iuout) TITLE(ir),tout
         write(*,*) "TITLE",TITLE(ir)
      enddo

      close(iuout) 

      deallocate(tsource,ttargglob,tout,tsourc4)

      end subroutine read_regrid_write_4D_1R
c*


      subroutine read_regrid_write_4D_1R_r8(x2grids,name,iuin)
      use regrid
      use ncio, only : defvar,write_data
      implicit none
      include 'netcdf.inc'
      type(x_2grids), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),tcopy(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat,outnc
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,
     &     status,fid,vid

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt r8",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),tcopy(imt,jmt,ntt) )
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R_r4_r8(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call do_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         if (ir .eq. 1) tcopy=ttargglob
         write(unit=iuout) TITLE(ir),ttargglob(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
      enddo

      close(iuout) 

      write(*,*) "here w r8"

      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc

      write(*,*) "TCOPY=",tcopy
      status = nf_create(outnc,nf_clobber,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      
      call defvar(fid,imt,jmt,ntt,tcopy,'press(im,jm,tile)')

      status = nf_enddef(fid)
      if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      
      call write_data(fid,'press',tcopy)

      deallocate(tsource,ttargglob,tcopy)

      end subroutine read_regrid_write_4D_1R_r8
c*



      subroutine read_regrid_write_4D_1R_rmax(x2grids,name,iuin,rmax)
      use regrid
      implicit none
      type(x_2grids), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer, intent(in):: rmax
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),arrsum(:,:,:)
      real*4, allocatable :: tout(:,:,:),data(:,:,:,:)
      character*80, allocatable :: TITLE(:)
      character*80 :: outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),arrsum(ims,jmt,ntt),
     &     tout(imt,jmt,ntt),
     &     data(ims,jms,nts,rmax),
     &     TITLE(rmax))
      tsource(:,:,:,:)=0.0
      data(:,:,:,:)=0.0
      
      do irec=1,rmax
         read(unit=iuin) TITLE(irec), data(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= data(:,:,:,irec)
      enddo
      maxrec=rmax
            
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      arrsum(:,:,:)=0.

      do ir=1,maxrec
         call do_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         arrsum(:,:,:)=arrsum(:,:,:)+ttargglob(:,:,:)
         tout(:,:,:)=ttargglob(:,:,:)
         write(unit=iuout) TITLE(ir),tout(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
      enddo

      write(*,*) "SUM ARRAY=",arrsum

      close(iuout) 

      deallocate(tsource,ttargglob,arrsum,tout,data,TITLE)

      end subroutine read_regrid_write_4D_1R_rmax
c*


      subroutine read_regrid_4D_1R(x2grids,iuin,TITLE,
     &     ttargglob,maxrec)
      use regrid
      implicit none
      type(x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iuin
      integer, intent(inout) :: maxrec
      real*8, allocatable :: tsource(:,:,:,:)
      real*8 :: ttargglob(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget,nrecmax)
      character*80, intent(inout) :: TITLE(nrecmax)
      integer :: ir,ims,jms,nts

      
      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource

      allocate (tsource(ims,jms,nts,nrecmax))
      tsource(:,:,:,:)=0.
      ttargglob(:,:,:,:)=0.

      write(*,*) "ims, jms, nts=",ims,jms,nts
            
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      do ir=1,maxrec
         call do_regrid(x2grids,tsource(:,:,:,ir),
     &        ttargglob(:,:,:,ir))
c        write(*,*) "TITLE",TITLE(ir)
      enddo

      deallocate(tsource)

      end subroutine read_regrid_4D_1R
c*



      subroutine read_regrid_write_4D_2R(x2grids,name,iuin)
      use regrid
      type(x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iuin
      character*80 name
      integer :: iuout
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
      real*4, allocatable :: data1(:,:,:,:),data2(:,:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     data1(ims,jms,nts,nrecmax),
     &     data2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt),
     &     tbig(imt,jmt,2,ntt) )

      tsource1(:,:,:,:)=0.0
      tsource2(:,:,:,:)=0.0

      irec=1
      do
         read(unit=iuin,END=30) TITLE(irec), data1(:,:,:,irec),
     &        data2(:,:,:,irec)
         write(*,*) "tile",TITLE(irec)
         tsource1(:,:,:,irec)= data1(:,:,:,irec)
         tsource2(:,:,:,irec)= data2(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1
      close(iuin)

      write(*,*) "maxrec",maxrec
      write(*,*) "TITLE RECS",TITLE(:)

      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=21
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      open( 33, FILE="ost1",
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call do_regrid(x2grids,tsource1(:,:,:,ir),ttargglob1)
         tout1(:,:,:)=ttargglob1(:,:,:)
         call do_regrid(x2grids,tsource2(:,:,:,ir),ttargglob2)
         tout2(:,:,:)=ttargglob2(:,:,:)
         tbig(:,:,1,:)=tout1
         tbig(:,:,2,:)=tout2
         write(unit=iuout) TITLE(ir),tbig
         write(unit=33) TITLE(ir),tout1
         write(*,*) "TITLE",TITLE(ir)
      enddo
      
      close(iuout) 
      close(33)
      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2,
     &     tbig,data1,data2)

      end subroutine read_regrid_write_4D_2R
c*


      subroutine regrid_mask_OSST(x2grids,name,iuin,
     &    FOCEAN_LL,missing)
c
c     regrid OSST using FOCEAN mask
c
      use regrid
      implicit none
      type(x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iuin
      character*80, intent(in) :: name
      integer :: iuout,ims,jms,nts,imt,jmt,ntt
      real*4,intent(in) :: FOCEAN_LL(x2grids%imsource,
     &     x2grids%jmsource,x2grids%ntilessource)
      real*8, allocatable :: foc8(:,:,:)
      real*4, allocatable :: FOCEAN_CS(:,:,:)
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*4, allocatable :: data1(:,:,:,:),data2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*8, allocatable :: tt1(:,:,:),tt2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
      real*8, intent(in) :: missing
      character*80 :: TITLE(nrecmax),outunformat,T2,filefcs
      integer :: ir,irec,maxrec,i,j,k

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt

      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     data1(ims,jms,nts,nrecmax),
     &     data2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tt1(0:imt+1,0:jmt+1,ntt),
     &     tt2(0:imt+1,0:jmt+1,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt),
     &     tbig(imt,jmt,2,ntt),
     &     FOCEAN_CS(imt,jmt,ntt),
     &     foc8(ims,jms,nts)
     &     )

      irec=1    
      do
         read(unit=iuin,END=30) TITLE(irec), data1(:,:,:,irec),
     &        data2(:,:,:,irec)

         tsource1(:,:,:,irec)= data1(:,:,:,irec)
         tsource2(:,:,:,irec)= data2(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      write(*,*) "maxrec",maxrec
      write(*,*) "TITLE RECS",TITLE(:)


      foc8=FOCEAN_LL

      outunformat=trim(name)//".CS"
      write(*,*) outunformat
      if (imt .eq. 32) then
      filefcs="Z_CS32"
      else if (imt .eq. 90) then
      filefcs="Z_CS90"
      endif
      open(50, FILE=filefcs,
     &     FORM='unformatted', STATUS="UNKNOWN")
      read(50) T2,FOCEAN_CS
      write(*,*) T2
      close(50)


      write(*,*) outunformat
      iuout=21
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      open( 33, FILE="ost1",
     &     FORM='unformatted', STATUS="UNKNOWN")
      
      do ir=1,maxrec
         call do_regrid_wt(x2grids,foc8,missing,
     &        tsource1(:,:,:,ir),ttargglob1)
         call do_regrid_wt(x2grids,foc8,missing,
     &        tsource2(:,:,:,ir),ttargglob2)

         tt1(1:imt,1:jmt,:)=ttargglob1(:,:,:)
         tt2(1:imt,1:jmt,:)=ttargglob2(:,:,:)

         tt1(0,1:jmt,:)=ttargglob1(1,:,:)
         tt1(imt+1,1:jmt,:)=ttargglob1(imt,:,:)
         tt1(1:imt,0,:)=ttargglob1(:,1,:)
         tt1(1:imt,jmt+1,:)=ttargglob1(:,jmt,:)

         tt2(0,1:jmt,:)=ttargglob2(1,:,:)
         tt2(imt+1,1:jmt,:)=ttargglob2(imt,:,:)
         tt2(1:imt,0,:)=ttargglob2(:,1,:)
         tt2(1:imt,jmt+1,:)=ttargglob2(:,jmt,:)

c*       fix inconsistencies between input and target FOCEAN masks
         do k=1,6
           do j=1,jmt
             do i=1,imt
                if (FOCEAN_CS(i,j,k) .gt. 0) then
                   if ( tt1(i,j,k) .eq. missing) then
                      if (tt1(i+1,j,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i+1,j,k)
                      elseif (tt1(i+1,j+1,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i+1,j+1,k)
                      elseif (tt1(i,j+1,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i,j+1,k)
                      elseif (tt1(i-1,j+1,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i-1,j+1,k) 
                      elseif (tt1(i-1,j,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i-1,j,k)
                      elseif (tt1(i-1,j-1,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i-1,j-1,k)
                      elseif (tt1(i,j-1,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i,j-1,k)
                      elseif (tt1(i+1,j-1,k) .gt. missing) then
                         tt1(i,j,k)= tt1(i+1,j-1,k)
                      else
                         write(*,*) "PROBLEM WITH INCONSISTENT FOCEAN"
                      endif
                      write(*,*) "fixed->",focean_cs(i,j,k),i,j,k,
     &                     tt1(i,j,k)
                  endif
                  if ( tt2(i,j,k) .eq. missing) then
                      if (tt2(i+1,j,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i+1,j,k)
                      elseif (tt2(i+1,j+1,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i+1,j+1,k)
                      elseif (tt2(i,j+1,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i,j+1,k)
                      elseif (tt2(i-1,j+1,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i-1,j+1,k) 
                      elseif (tt2(i-1,j,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i-1,j,k)
                      elseif (tt2(i-1,j-1,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i-1,j-1,k)
                      elseif (tt2(i,j-1,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i,j-1,k)
                      elseif (tt2(i+1,j-1,k) .gt. missing) then
                         tt2(i,j,k)= tt2(i+1,j-1,k)
                      else
                         write(*,*) "PROBLEM WITH INCONSISTENT FOCEAN"
                      endif
                      write(*,*) "fixed->",focean_cs(i,j,k),i,j,k,
     &                     tt2(i,j,k)
                   endif
                else                   
                   tt1(i,j,k)=missing
                   tt2(i,j,k)=missing
                   write(*,*) "fixed pure land cell",i,j,k,
     &                  focean_CS(i,j,k)
                endif
             enddo
          enddo
       enddo
                      

       write(*,*) "aft loop incons"
       tbig(:,:,1,:)=tt1(1:imt,1:jmt,:)
       tbig(:,:,2,:)=tt2(1:imt,1:jmt,:)
       write(unit=iuout) TITLE(ir),tbig
       write(unit=33) TITLE(ir),tbig(:,:,1,:)
       write(*,*) "TITLE",TITLE(ir)
      enddo
      
      close(iuout)
      close(33)
      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tt1,tt2,
     &     tbig,data1,data2,FOCEAN_CS,foc8)

      end subroutine regrid_mask_OSST


      subroutine testOSST
      character*80 :: name,T2,TITLE(12)
      real*4 :: FOCEAN(32,32,6)
      real*4 :: tbig(32,32,2,6)
      real*4, parameter :: missing = -9.999999171244748e33
      integer iu_OSST,ir,failed
      
      iu_OSST=20
      name="OSST_CS32"
      write(*,*) name
      open(iu_OSST, FILE=name,FORM='unformatted', STATUS='old')
      
      open(21,FILE="Z_CS32",FORM='unformatted', STATUS='old')
      read(21) T2,FOCEAN
      close(21)
      
      do ir=1,12
         failed=0
         read(unit=iu_OSST) TITLE(ir),tbig
         write(*,*) TITLE(ir)
         do k=1,6
            do j=1,32
               do i=1,32
                  if (tbig(i,j,1,k) .eq. missing) then
c     write(*,*) "mask cell",i,j,k
                     if (FOCEAN(i,j,k) .gt. 0) then
                        failed=1
                     endif
                  endif
                  if (tbig(i,j,2,k) .eq. missing) then
c     write(*,*) "mask cell",i,j,k
                     if (FOCEAN(i,j,k) .gt. 0) then
                        failed=1
                     endif
                  endif
               enddo
            enddo
         enddo
         if (failed.ne.0) then
            write(*,*) "record ",ir," failed test"
         else
            write(*,*) "record ",ir," succeeded test"
         endif
      enddo
      
      close(iu_OSST)

      end subroutine testOSST

      module area
      contains
      real*8 function aint(x,y)
c calculates the area integral from the center of a cube face
c to the point x,y.
      implicit none
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: x,y
c      real*8 :: aint
      real*8 :: tangx,tangy
      tangx = tan(g*x)
      tangy = tan(g*y)
      aint = atan(2.*tangx*tangy/sqrt(1.+2.*(tangx*tangx+tangy*tangy)))
      return
      end function aint
      end module area

      subroutine areacs(axyp,im,jm)
      use area, only : aint
      IMPLICIT NONE
      SAVE
      integer :: im,jm
      real*8 :: x,y,xp1,yp1
      real*8, intent(out) :: axyp(im,jm)
      real*8,parameter :: radius = 6371000.
      integer :: i,j

      do j=1,jm
         do i=1,jm
            x = -1d0 + 2d0*(i-1)/im
            y = -1d0 + 2d0*(j-1)/im
            
            xp1= -1d0 + 2d0*i/im
            yp1= -1d0 + 2d0*j/im
            
            axyp(i,j) = aint(xp1,yp1)-aint(x,yp1)
     &           -aint(xp1,y)+aint(x,y)
            axyp(i,j) = axyp(i,j)*radius*radius
         enddo
      enddo

      end subroutine areacs


      subroutine regridRDSCAL(x2grids)
c     
c     regridding scalar quantities used to derive river directions 
c     using simple sampling scheme (non conservative) 
c
      use regrid
      implicit none
      type(x_2grids), intent(in) :: x2grids
      character*80 :: TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,
     &     TITLE6,TITLE7,name
      real*8,allocatable :: ttargglob(:,:,:)
      real*8,allocatable :: ts(:,:,:),ts2(:,:,:),lon_dg(:),lat_dg(:)
      real*4,allocatable :: ttargupst4(:,:,:)
      real*4,allocatable :: ttargdist4(:,:,:)
      real*4,allocatable :: tsourc4(:,:,:),tlondg(:,:,:),tlatdg(:,:,:)
      integer, allocatable :: bassId(:,:,:),tbassId(:,:,:)
      integer:: iu_RDSCAL,i,j,k,iunit,imt,jmt,ntt,ims,jms,nts,
     &     status,nijt,ierr,ilon,jlat,tile
      real*8 :: dlat_dg,dlon_dg,x,y,area
      integer :: bId_max,imax,jmax
      real*8,parameter :: pi = 3.1415926535897932d0 
      real*8,parameter :: twopi = 2d0*pi         
      real*8, parameter :: shiftwest = twopi/36.
      real*8,parameter :: radian = pi/180d0      
      REAL*8 :: fjeq



      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "imt jmt ntt ims jms nts",imt,jmt,ntt,ims,jms,nts

      allocate( 
     &     ttargglob(imt,jmt,ntt),
     &     ttargupst4(imt,jmt,ntt),
     &     ttargdist4(imt,jmt,ntt),
     &     tbassId(imt,jmt,ntt),
     &     tlondg(imt,jmt,ntt),
     &     tlatdg(imt,jmt,ntt),
     &     tsourc4(ims,jms,nts),
     &     ts(ims,jms,nts),
     &     ts2(ims,jms,nts),
     &     bassId(ims,jms,nts),
     &     lon_dg(ims),lat_dg(jms)
     &     )
      


      FJEQ=.5*(1+jms)
      DLAT_DG=180./REAL(jms)
      LON_DG(1) = -180.+360./(2.*FLOAT(ims))
      DO I=2,ims
        LON_DG(I) = LON_DG(I-1)+360./FLOAT(ims)
      END DO
      LAT_DG(1)=-90.
      LAT_DG(jms)=90.
      DO J=2,jms-1
        LAT_DG(J)=DLAT_DG*(J-FJEQ)    ! primary (tracer) latitudes
      END DO

c*

      iu_RDSCAL=20
      name="STN-30p.bin"

      write(*,*) name
      open( iu_RDSCAL, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "ims, jms, nts=",ims,jms,nts
      read(iu_RDSCAL) TITLE1
      write(*,*) TITLE1
      read(iu_RDSCAL) TITLE2
      write(*,*) TITLE2
      read(iu_RDSCAL) TITLE3,bassId(:,:,:)
      write(*,*) TITLE3
      read(iu_RDSCAL) TITLE4,tsourc4(:,:,:)
      write(*,*) TITLE4
      ts(:,:,:)=tsourc4(:,:,:)
      read(iu_RDSCAL) TITLE5,tsourc4(:,:,:)
      write(*,*) TITLE5
      ts2(:,:,:)=tsourc4(:,:,:)
      do k=1,ntt
      do j=1,jmt
      do i=1,imt
c*    for each cell, choose id of bassin with greatest area
      call maxarea_bassId(x2grids,bassId,i,j,k,area,bId_max
     &     ,ims,jms,nts,imax,jmax)
      tbassId(i,j,k)=bId_max

      if (i .ge. 31 .and. i .le. 33 .and. j .ge. 14 .and.
     &     j .le. 16 .and. k == 3) then
         write(*,*) i,j,k,bId_max,ts(imax,jmax,1)
      endif

c*    choose distance to the ocean of latlon cell with bass Id = bId_max 
c*    and with largest exchange grid area 
      ttargdist4(i,j,k)=ts(imax,jmax,1)
      ttargupst4(i,j,k)=ts2(imax,jmax,1)
      tlondg(i,j,k)=lon_dg(imax)
      tlatdg(i,j,k)=lat_dg(jmax)
      enddo
      enddo
      enddo
      close (iu_RDSCAL)

      name="STN_C90"
c      name="STN_CS32"

      TITLE6="longitude of corresponding latlon grid cell"
      TITLE7="latitude of corresponding latlon grid cell"
      open( iu_RDSCAL, FILE=name,FORM='unformatted', STATUS='unknown')
      write(iu_RDSCAL) TITLE1
      write(iu_RDSCAL) TITLE2
      write(iu_RDSCAL) TITLE3,tbassId(:,:,:)
      write(iu_RDSCAL) TITLE4,ttargdist4(:,:,:)
      write(iu_RDSCAL) TITLE5,ttargupst4(:,:,:)
      write(iu_RDSCAL) TITLE6,tlondg(:,:,:)
      write(iu_RDSCAL) TITLE7,tlatdg(:,:,:)
      close(iu_RDSCAL)
      

      deallocate(ttargglob,ts,ts2,ttargupst4,ttargdist4,bassId,
     &             tbassId,lon_dg,lat_dg)

      end subroutine regridRDSCAL
c*

      subroutine csxy2ll(x,y,tile,lon,lat)
c converts x,y,tile to lon,lat (radians)
c This routine places the center of tile 1 at the IDL.
      implicit none
      real*8 :: x,y ! input
      integer :: tile ! input
      real*8 :: lon,lat ! output
      real*8,parameter :: pi = 3.1415926535897932d0 
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      real*8 :: gx,gy,tmpx,tmpy,cosgx,tangx,tangy,coslon
      gx = g*x
      gy = g*y
      if(tile.eq.4 .or. tile.eq.5) then ! 90 deg rotation
        tmpx = gx
        tmpy = gy
        gx = +tmpy
        gy = -tmpx
      elseif(tile.eq.6) then ! tile 6 = tile 3 flipped around the axis x+y=0
        tmpx = gx
        tmpy = gy
        gx = -tmpy
        gy = -tmpx
      endif
      if(tile.eq.3 .or. tile.eq.6) then
        tangx = tan(gx)
        tangy = tan(gy)
        lat = atan(1d0/sqrt(2d0*(tangx**2 +tangy**2)+1d-40))
        lon = atan2(tangy,tangx)
        if(tile.eq.6) lat = -lat
      else
        cosgx = cos(gx)
        coslon = cosgx/sqrt(2d0-cosgx*cosgx)
        lat = atan(coslon*sqrt(2d0)*tan(gy))
        lon = sign(acos(coslon),gx)
! add longitude offset to tiles 1,2,4,5. integer arithmetic.
        lon = lon + (mod(tile,3)-1)*pi/2. -pi*(1-(tile/3))
        if(lon.lt.-pi) lon=lon+2.*pi
      endif
      return
      end subroutine csxy2ll

      subroutine maxarea_bassId(x2grids,bassId,ics,jcs,tile,
     &     area_max,bId_max,ims,jms,nts,imax,jmax)
!@sum  for each cell choose bassin Id with greatest area
!@auth Denis Gueyffier
      use regrid
      implicit none
      type (x_2grids), intent(in) :: x2grids
      
      integer, parameter :: bassMax=6119
      real*8 :: area_max, area_cumul(bassMax+1),aream,bIdM
      integer :: n,ics,jcs,tile,itile,icc,jcc,il,jl,
     &     bId,ims,jms,nts,imax,jmax
      integer :: bId_max(1)
      integer :: bassId(ims,jms,nts)
c      real*8 :: area(x2grids%imtarget,x2grids%jmtarget,
c     &     x2grids%ntilestarget)

      area_max = 0.d0
      area_cumul = 0.d0
c      area=0.0d0

c      write(*,*) "ics,jcs,tile=",ics,jcs,tile

      do n=1,x2grids%xgrid%ncells

         itile=x2grids%xgrid%tile(n)
         icc=x2grids%xgrid%ijcub(1,n)
         jcc=x2grids%xgrid%ijcub(2,n)
         il=x2grids%xgrid%ijlatlon(1,n)
         jl=x2grids%xgrid%ijlatlon(2,n)
         
c         area(icc,jcc,itile) = area(icc,jcc,itile)
c     &        + x2grids%xgrid%xgrid_area(n)

         if ( (itile == tile)    .and. (icc == ics) 
     &        .and. (jcc == jcs) ) then
c     Sum up area for each separate bassin in the cell.
c     bassin index shifted by 1 to make room for bass Id = -9999   
            if (bassId(il,jl,1) == -9999) then
               area_cumul(1)=
     &                 area_cumul(1)
     &              + x2grids%xgrid%xgrid_area(n)
            endif
            do bId=1,bassMax
               if (bassId(il,jl,1) == bId) then
               area_cumul(bId+1)=
     &                 area_cumul(bId+1)
     &              + x2grids%xgrid%xgrid_area(n)
               endif
            enddo
         endif
      enddo

      area_max = maxval(area_cumul(:))
      bId_max  = maxloc(area_cumul(:))-1
      if (bId_max(1) == 0 ) bId_max(1)=-9999

c     test area
c      if (sum(area_cumul(:))-area(ics,jcs,tile) .gt. 1.e-10) then
c         write(*,*) "warning area",ics,jcs,tile,bId_max
c         stop
c      endif


      bIdM=bId_max(1)   ! id of bassin with max area

c     Now find x-grid cell belonging to bassin id = bIdM
c     with the largest area

      aream=0.0d0

      do n=1,x2grids%xgrid%ncells

         itile=x2grids%xgrid%tile(n)
         icc=x2grids%xgrid%ijcub(1,n)
         jcc=x2grids%xgrid%ijcub(2,n)
         il=x2grids%xgrid%ijlatlon(1,n)
         jl=x2grids%xgrid%ijlatlon(2,n)
         
         if ( (itile == tile)    
     &        .and. (icc == ics) 
     &        .and. (jcc == jcs) 
     &        .and. (bassId(il,jl,1) == bIdM) 
     &        ) then

            if (x2grids%xgrid%xgrid_area(n) .ge. aream) then
               aream=x2grids%xgrid%xgrid_area(n)
               imax=il
               jmax=jl
            endif
         endif
      enddo

c      write(*,*) ics,jcs,tile,aream/sum(area_cumul(:))


      end subroutine 
