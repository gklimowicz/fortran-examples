#define GCM_DRV

      module ent_prescribed_drv_geo
!@sum ent_prescr_veg_geo   Module of routines to read geographic cover, max LAI, canopy heights, and calculate plant densities and carbon pools.  Also, routine to read geographic initial LAI.
!@auth N.Y.Kiang

      use ent_const
      use ent_types
      use ent_pfts
      
      implicit none
      private
      save

      public IJMULT,IJDIV,IJADD,IJSUB
      public init_entvegdata_geo
      public init_ent_laimax_geo

      !* Maybe get this from GCM
!      real*8, parameter :: undef=-1.d30 !no value




      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES                 *
      !***************************************************

      contains

!------------------------------------------------------------------------


      !---- Utility functions -----------------------------------

      logical function NONUM(num)
      real*8 :: num
      if ((num.le.undef).or.(isNaN(num))) then
         NONUM = .true.
      else
         NONUM = .false.
      endif
      end function NONUM

      function NONUM4(num) Result(NONUMres)
      real*4 :: num
      logical :: NONUMres
      if ((num.le.undef).or.(isNaN(num))) then
         NONUMres = .true.
      else
         NONUMres = .false.
      endif
      end function NONUM4


      function CHECKNUM(num) Result(CHECKres)
!@sum Return 0.0 if num is undef or NaN
      real*8 :: num
      real*8 :: CHECKres
      if (NONUM(num)) then
         CHECKres = 0.d0
      else
         CHECKres = num
      endif
      end function CHECKNUM


      function CHECKNUM4(num) Result(CHECKres)
!@sum Return 0.0 if num is undef or NaN
      real*4 :: num
      real*4 :: CHECKres
      if (NONUM4(num)) then
         CHECKres = 0.0
      else
         CHECKres = num
      endif
      end function CHECKNUM4

      function DIV0(num, denom) Result(DIVres)
!@sum Returns 0 if divide by zero or divide by undef.
      real*8 :: num, denom, DIVres
      real*8 :: n,d

      n = CHECKNUM(num)
      d = CHECKNUM(denom)
      if (d.eq.0.d0) then
         DIVres = 0.0
      else
         DIVres = n/d
      endif

      end function DIV0


      function IJDIV(I0,I1,J0,J1,num, denom) Result(DIVres)
!@sum IJ array division, return zero if divide-by-zero or -undef.
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: num(I0:I1,J0:J1),denom(I0:I1,J0:J1)
      real*8 :: DIVres(I0:I1,J0:J1)
      integer :: i,j

      do i=I0,I1
         do j=J0,J1
            DIVres(i,j) = DIV0(num(i,j),denom(i,j))
         enddo
      enddo
      end function IJDIV

      function IJMULT(I0,I1,J0,J1,fac1, fac2) Result(MULTres)
!@sum IJ array multiplication, return zero if multiply by undef.

      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*8 :: MULTres(I0:I1,J0:J1)
      integer :: i,j
      real*8 :: m1, m2

      do i=I0,I1
         do j=J0,J1
            m1 = CHECKNUM(fac1(i,j))
            m2 = CHECKNUM(fac2(i,j))
            MULTres(i,j) = m1 * m2
         enddo
      enddo
      end function IJMULT

      function IJSUB(I0,I1,J0,J1,fac1, fac2) Result(SUBres)
!@sum IJ array subtraction, with check for undef.      
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*8 :: SUBres(I0:I1,J0:J1)
      integer :: i,j
      real*8 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM(fac1(i,j))
            a2 = CHECKNUM(fac2(i,j))
            SUBres(i,j) = a1 - a2
         enddo
      enddo
      end function IJSUB

      function IJADD(I0,I1,J0,J1,fac1, fac2) Result(ADDres)
!@sum IJ array addition, with check for undef.
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*8 :: ADDres(I0:I1,J0:J1)
      !---------
      integer :: i,j
      real*8 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM(fac1(i,j))
            a2 = CHECKNUM(fac2(i,j))
            ADDres(i,j) = a1 + a2
         enddo
      enddo
      end function IJADD


      function IJADD4(I0,I1,J0,J1,fac1, fac2) Result(ADDres)
!@sum IJ array addition, with check for undef. real*4
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*4 :: fac1(I0:I1,J0:J1),fac2(I0:I1,J0:J1)
      real*4 :: ADDres(I0:I1,J0:J1)
      !---------
      integer :: i,j
      real*8 :: a1, a2

      do i=I0,I1
         do j=J0,J1
            a1 = CHECKNUM4(fac1(i,j))
            a2 = CHECKNUM4(fac2(i,j))
            ADDres(i,j) = a1 + a2
         enddo
      enddo
      end function IJADD4

!      real*8 DIMENSION(I0:I1,J0:J1) function IJzero(datij)
!      real*8 :: datij(I0:I1,J0:J1)
!      !---------
!      integer :: i,j
!
!      do i=I0,I1
!         do j=J0,J1
!            IJzero(i,j) = CHECKNUM(dat(i,j))
!         enddo
!      enddo
!
!      end function IJzero


      !-----------ROUTINES FOR READING INPUT FILES---------------!

      subroutine init_ent_laimax_geo(IM,JM,I0,I1,J0,J1,laimax)
!@sum Read the maximum annual LAI for vegetation structure from allometry.
!@auth N.Y.Kiang
      use FILEMANAGER, only : openunit,closeunit,nameunit !for VEG_PROGNOSTIC
      use ent_const,only : N_COVERTYPES
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8 :: laimax(N_COVERTYPES,I0:I1,J0:J1) 
      !----------
      character*80 :: title
      real*4 :: buf(IM,JM)
      integer :: iu_LAImax
      integer :: k !@var cover type

      call openunit("LAImax",iu_LAImax,.true.,.true.)

      do k=1,N_COVERTYPES
        read(iu_LAImax) title , buf
        laimax(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      end do
      !print *,"laimax", laimax(:,I0:I1,J0:J1) !#DEBUG

      call closeunit(iu_LAImax)

      !* Return lai for each vegetation type.
      end subroutine init_ent_laimax_geo


!-----------------------------------------------------------------

      subroutine init_entvegdata_geo( IM,JM,I0,I1,J0,J1
     i     ,laidata,hdata,laimaxgeo
     o     ,popdata,dbhdata
     o     ,craddata,cpooldata)  !,rootprofdata)--##LATER MUST SCALE w/HEIGHT.
!     &     ,soil_color,soil_texture,
!     &     Tpooldata)

!@sum get_entvegdata_geo Called by main driver to initialize
!@+   geographic vegetation structure.
!@+   Calls internal routines to read in LAI max data 
!@+   and calculate vegetation population & cpools
!@+   using allometry functions.
!@+   Outputs diagnostic files of plant density and biomass on startup.

      use ent_const
      use ent_pfts
      use allometryfn, only : dbh2Cdead,height2dbh, Cfol_fn,Csw_fn
     &     ,hw_fract,Crown_rad_allom, Crown_rad_max_from_density
     &     ,allom_plant_cpools,init_Clab

      implicit none
      integer, intent(in) :: im, jm, i0, i1, j0, j1
      !Vegetation structure arrays input
      real*8,intent(in),optional :: laidata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in) :: laimaxgeo(N_COVERTYPES,I0:I1,J0:J1) 

      !Vegetation structure arrays calculated
      real*8,intent(out) :: dbhdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) !gC
      real*8,intent(out) :: popdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: craddata(N_COVERTYPES,I0:I1,J0:J1)
!      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)

      !-----Local---------------------------------------

      character*80 :: TITLE, TITLES(N_COVERTYPES)
      integer :: skip

      !Arrays for calculating veg structure.
      real*4 :: L1(I0:I1,J0:J1)
      real*4 :: L2(I0:I1,J0:J1)
      real*4 :: L3(I0:I1,J0:J1)
      real*4 :: L4(I0:I1,J0:J1)
      real*8 :: LAmax(N_COVERTYPES,I0:I1,J0:J1)!Leaf area max (m2) per plant
      real*8 :: LA(N_COVERTYPES,I0:I1,J0:J1)!Leaf area actual (m2) per plant

      integer :: i,j,k,p

      call flush(6)
      !******* CALCULATE DENSITY AND BIOMASS POOLS *****************
!      call get_ent_laimax_geo(IM,JM,I0,I1,J0,J1,laimaxgeo) !Read file

      !* DBH, cpools, LA/plant,population density, crown radius
      dbhdata(:,:,:) = 0.d0
      cpooldata(:,:,:,:) = 0.d0 
      popdata(:,:,:) = 0.d0
      craddata = 0.d0
      do p=1,N_PFT
         print *,"p",p
         do i=I0,I1
            do j=J0,J1
!               print *,'hdata(p,i,j)',p,i,j,hdata(p,i,j)  !##NK
!               print *,'laidata(p,i,j)',p,i,j,laidata(p,i,j)
!               print *,'laimaxgeo(p,i,j)',p,i,j,laimaxgeo(p,i,j)
            if ((hdata(p,i,j).gt.0.d0).and.
     &            (laimaxgeo(p,i,j).gt.0.d0)) then
               dbhdata(p,i,j) = height2dbh(p,hdata(p,i,j))
!		print *,'height2dbh',p,i,j,hdata(p,i,j),laimaxgeo(p,i,j)
!               if (pfpar(p)%woody) then
!                  if ((hdata(p,i,j).gt.0.d0).and.
!     &                 (dbhdata(p,i,j).eq.0.d0)) then
!     &                 print *, 'WARNING: h>0,dbh=0' , p,i,j
!     &                 hdata(p,i,j),dbhdata(p,i,j)
!                  endif
!               endif

!               !cpooldata(p,HW,i,j) = Cstem_fn(p,dbhdata(p,i,j))
!               cpooldata(p,HW,i,j)=dbh2Cdead(p,dbhdata(p,i,j))*hw_fract
!		print *,'dbh2Cdead'
!               cpooldata(p,SW,i,j) = 
!!     &             Csw_fn(p,dbhdata(p,i,j),hdata(p,i,j),laidata(p,i,j))
!     &               Csw_fn(p,dbhdata(p,i,j),hdata(p,i,j)
!     &               ,max(0.5d0*laimaxgeo(p,i,j),laidata(p,i,j)))
               call flush(6)
               
               !LAmax(p,i,j) = dbh2Cfol(p,dbhdata(p,i,j)) * pfpar(p)%sla
               LAmax(p,i,j) = 0.001d0 * !gC to kgC for sla in m2/kgC
     &              Cfol_fn(p,dbhdata(p,i,j),hdata(p,i,j))*pfpar(p)%sla
               popdata(p,i,j) = DIV0(laimaxgeo(p,i,j),LAmax(p,i,j))
!               ! Get LA actual per plant
!!               if (.not.present(laidata)) then
!!                  call calc_leafarea_plant_hemi(
!!     &                 IM,JM,i,j,p,LAmax(p,i,j),LA(p,i,j))
!!               else
!                  LA(p,i,j) = DIV0(laidata(p,i,j), popdata(p,i,j))
!		  print *,'LA',p,i,j
!!               endif
!               cpooldata(p,FOL,i,j) = LA(p,i,j)/pfpar(p)%sla*1.d3
!                print *, 'FOL',p,i,j               
!               !cpooldata(p,CR,i,j) = Ccroot_fn(p,Cstem(p,i,j))
!               cpooldata(p,CR,i,j) = cpooldata(p,HW,i,j)*
!     &              (1-hw_fract)/hw_fract
!                print *, 'CR',p,i,j
!               !cpooldata(p,CR,i,j) = Cfroot_fn(p,Cfol(p,i,j))
!               cpooldata(p,FR,i,j) = cpooldata(p,FOL,i,j)


               call allom_plant_cpools(p
     &              ,laidata(p,i,j),hdata(p,i,j)
     &              ,dbhdata(p,i,j),popdata(p,i,j)
     &              ,cpooldata(p,:,i,j))

               call init_Clab(p,dbhdata(p,i,j)
     &              ,hdata(p,i,j),cpooldata(p,LABILE,i,j)) !LABILE

               craddata(p,i,j) = min(Crown_rad_allom(p,hdata(p,i,j)),
     &              Crown_rad_max_from_density(popdata(p,i,j)))
               call flush(6)
            endif
            enddo
         enddo
      enddo
      write(*,*) "Calculated biomass pools, LA, density, crad, lai"
      !write(1231) cpooldata(7,HW,:,:)

      !Return for ent_cell_set:
      !vegdata=coverdata, popdata, laidata,
      !hdata, dbhdata, craddata, cpooldata

      !These get from ent_prescribed_drv.f
      !nmdata, rootprofdata,
      !soil_color, albedodata, soil_texture, Tpool_ini
      !do_soilinit, do_phenology_activegrowth

      write(*,*) "Done calculating Ent GVSD arrays."

      ! Write a diagnostics file.  DEBUGGING CHECKS.
      !do k=1,N_COVERTYPES
      !   TITLE = trim(Ent_title(k))//' density (#/m2)'
      !   write(9950) TITLE, popdata(k,:,:)
      !enddo
 !!     L1(:,:) = 0.d0
 !!     L2(:,:) = 0.d0
 !!     L3(:,:) = 0.d0
 !!     L4(:,:) = 0.d0
 !!     open(9950,file='Ent16dens.ij',form='unformatted',status='unknown')
 !!     open(9951,file='Ent16cpools.ij',
 !!    &     form='unformatted',status='unknown')
 !!     do p=1,N_COVERTYPES
 !!           TITLE = trim(Ent_title(p)(1:50))//' popdata (#/m2)'
 !!           L2 = popdata(p,:,:)
 !!           write(9951) TITLE, L2
 !!           L3(:,:) = 0.0
 !!        do k = 1,N_BPOOLS
 !!           L2(:,:) = IJMULT(I0,I1,J0,J1
 !!    &           ,cpooldata(p,k,:,:),popdata(p,:,:))*1.d-3 !gC to kg-C
 !!           TITLE = trim(Ent_title(p)(1:50))//'  '//Ent_cpool_title(k)//
 !!    &           '   '//'Biomass (kgC/m2)'
 !!           write(9950) TITLE, L2(:,:)
 !!           TITLE = trim(Ent_title(p)(1:50))//'  '//
 !!    &           Ent_cpool_title(k)//'  (kgC/plant)'
 !!           L4(:,:) =cpooldata(p,k,:,:)*1.d-3
 !!           write(9951) TITLE, L4
 !!           L3(:,:) = IJADD4(I0,I1,J0,J1,L3(:,:),L2(:,:))
 !!           L1(:,:) = IJADD4(I0,I1,J0,J1,L1(:,:),L2(:,:))
 !!        enddo
 !!        TITLE = trim(Ent_title(p)(1:50))//'  Biomass  '//'(kgC/m2)'
 !!        write(9950) TITLE, L3
 !!     enddo
 !!     TITLE = 'Total Biomass   (kgC/m2)'
 !!     write(9950) TITLE, L1(:,:)
 !!     write(*,*) 'Wrote calculated density and Cpools to file.'
 !!     close(9950)
 !!     close(9951)

      end subroutine init_entvegdata_geo

      end module ent_prescribed_drv_geo


