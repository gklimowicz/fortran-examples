#include "rundeck_opts.h"
      module bio_inicond_mod
      
      implicit none

      contains

      subroutine bio_inicond_read(filename, array, depth,ogrid)

      use pario, only : par_open,par_close
     &     ,read_data,read_dist_data,get_dimlens
      USE DOMAIN_DECOMP_1D, only: DIST_GRID


      implicit none
 
      type(DIST_GRID), intent(in), optional :: ogrid
      character(len=*), intent(in) :: filename

      real*8, dimension(:,:,:), allocatable, intent(out) :: array
      real*8, dimension(:), allocatable, intent(out), optional :: depth
      integer :: fid, dlens(7), ndims

      fid=par_open(ogrid, filename, 'read')
      call get_dimlens(ogrid, fid, 'array', ndims, dlens)
      if ((dlens(1)/=ogrid%im_world).or.(dlens(2)/=ogrid%jm_world))
     &   call stop_model('dimension mismatch: '//trim(filename), 255)
      allocate(array(ogrid%i_strt:ogrid%i_stop,
     &                  ogrid%j_strt:ogrid%j_stop, dlens(3)))
      if (present(depth)) then
        allocate(depth(dlens(3)))
        call read_data(ogrid, fid, 'depth', depth, bcast_all=.true.)
      endif
      call read_dist_data(ogrid, fid, 'array', array)
      call par_close(ogrid, fid)
      end subroutine bio_inicond_read

!NOT FOR HYCOM: lmm passed to subroutine
      subroutine bio_inicond(filename,fldo,kdm,im,ogrid,ip,lmm)

      use obio_com, only: ze 
      USE DOMAIN_DECOMP_1D, only :DIST_GRID


      implicit none
  
      type(DIST_GRID), intent(in) :: ogrid
      integer, intent(in) :: kdm,im, 
     &                  lmm(im,ogrid%j_strt_halo:ogrid%j_stop_halo)      
      real, intent(in) :: ip(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &                      ogrid%j_strt_halo:ogrid%j_stop_halo)
      character(len=*), intent(in) :: filename
      real, dimension(ogrid%i_strt:ogrid%i_stop,
     &    ogrid%j_strt:ogrid%j_stop, kdm), intent(out) ::  fldo

      real*8, dimension(:,:,:), allocatable :: array
      real*8, dimension(:), allocatable :: depth, nodc_d
      real*8, dimension(:), allocatable :: dummy
      logical :: regrid
      integer :: i, j, k, kmax, nodc_kmax

      interface
        Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ, missing, foo)
        integer :: km,lm
        Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM), RZ(LM), missing
        logical foo
        end Subroutine VLKtoLZ
      end interface

      call bio_inicond_read(filename, array, depth,ogrid)
      regrid=kdm/=size(depth)
      if (.not.regrid) regrid=all(abs(depth-
     &                      ze(ogrid%i_strt, ogrid%j_strt, :))<1d0)
      if (regrid) then
#ifdef OBIO_ON_GISSocean
        allocate(dummy(size(fldo, 3)))
        do i=ogrid%i_strt,ogrid%i_stop
          do j=ogrid%j_strt,ogrid%j_stop
            if (ip(i, j)==0) cycle
            call vlktolz(size(depth), lmm(i, j), depth, ze(i, j, :),
     &           array(i, j, :), fldo(i, j, :), dummy, -999999.,.false.)
          end do
        end do
#else
        fldo=-9999.d0
        allocate(nodc_d(size(depth)+1))
        do j=ogrid%j_strt,ogrid%j_stop
        do i=ogrid%i_strt,ogrid%i_stop
          if (ip(i,j)==0) cycle
          kmax=1
          do k=1,kdm
            if (ze(i, j, k) .gt. ze(i, j, k-1)) kmax=k+1
          enddo
          do k=1,size(depth)
            if (depth(k) .le. ze(i, j, min(20,kmax))) then
              nodc_d(k)=depth(k)
              nodc_kmax=k
            endif
          enddo
          nodc_d(nodc_kmax+1)=ze(i, j, min(20,kmax))
          call remap1d_plm(array(i,j,1:nodc_kmax),nodc_d,nodc_kmax,
     .             fldo(i,j,1:kdm),ze(i, j, :),kdm,.false.,i,j)
        enddo
        enddo
#endif
      else
        fldo=array
      endif
      end subroutine bio_inicond

      end module bio_inicond_mod


c ----------------------------------------------------------------
      subroutine obio_init(kdm,dtsrc,ogrid,dlatm)
c --- biological/light setup
c ----------------------------------------------------------------
c 
      USE FILEMANAGER, only: openunit,closeunit,file_exists
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, DIST_GRID
      USE timestream_mod, only : init_stream
      USE obio_dim
      USE obio_incom
      USE bio_inicond_mod, only: bio_inicond_read
      USE obio_forc, only : atmFe,stream_atmFe,alk,surfN
      USE obio_com, only : npst,npnd,WtoQ,obio_ws,P_tend,D_tend
     .                    ,C_tend,wsdet,gro,obio_deltath,obio_deltat
     .                    ,sday
#ifdef OBIO_RUNOFF
!    .                    ,rnitrmflo_loc
     .                    ,rnitrconc_loc
     .                    ,rdicconc_loc
     .                    ,rdocconc_loc
     .                    ,rsiliconc_loc
     .                    ,rironconc_loc
     .                    ,rpocconc_loc
     .                    ,ralkconc_loc
#endif


#ifdef STANDALONE_OCEAN
      USE obio_forc, only: Eda,Esa
#endif
      USE pario
      USE model_com, only: modelEclock
      USE ocalbedo_mod, only: lam, ocalbedo_init=>init

      implicit none  

      integer,intent(in) :: kdm
      real, intent(in) :: dtsrc,dlatm
      type(DIST_GRID), intent(in) :: ogrid
    

      integer i,j,k
      integer jyear, jday
      integer iu_bio
      integer nt,nl
      integer imon,ihr,nrec,ichan
      integer lambda,ic
      integer icd,ntr,ich,ih,iu_fac
      integer fid
      real saw,sbw,sac,sbc
      real*4  facirr4(nh,nch,5,ncd)

      real planck,c,hc,oavo,rlamm,rlam450,Sdom,rlam,hcoavo
     .    ,rnn,rbot,pi
     .    ,dummy

      character*50 title
!     character*50 cfle
      character cacbc*11,cabw*10
      character*80 filename,fn,filename1,filename2

      data cacbc,cabw /'acbc25b.dat','abw25b.dat'/

c 
      if (AM_I_ROOT()) print*, 'Ocean Biology setup starts'

! time steps
      obio_deltath = dtsrc/3600.d0  !time step in hours
      obio_deltat = obio_deltath*3600.d0    !time step in s  !July 2016
      call modelEclock%get(year=jyear, dayOfYear=jday)
      if (AM_I_ROOT()) 
     . print*, 'Ocean Biology time step =',obio_deltat

c  Read in constants, light data
c  Computes constants over entire run of model, reads in required
c  data files, and otherwise obtains one-time-only information
c  necessary for the run.

c  Degrees to radians conversion
      pi = dacos(-1.0D0)
      pi2 = pi*2.0
      rad = 180.0D0/pi

#ifdef STANDALONE_OCEAN
      if (AM_I_ROOT()) then
      print*, '    '
      print*, 'reading OASIM data.....'
      print*, '    '
      endif

      !reading Eda
      filename1='oasimdirect1'
      filename2='oasimdirect2'
      call obio_edaesa_g(filename1,filename2)
#endif

      do nt = 1,nchl
       rkn(nt) = 0.0
       rks(nt) = 0.0
       rkf(nt) = 0.0
      enddo
c
c  Phytoplankton group parameters
      do nt = 1,nchl
       obio_wsd(nt)    = 0.0
       obio_wss(nt)    = 0.0
      enddo
      do nt = 1,nchl
       rmumax(nt) = 0.0
        rik(1,nt) = 0.0
        rik(2,nt) = 0.0
        rik(3,nt) = 0.0
      enddo
      Pdeep(1) = 32.0    !bottom BC for nitrate
      Pdeep(2) = 0.1     !ammonium
      Pdeep(3) = 60.0    !silica
      Pdeep(4) = 0.6     !iron from Archer and Johnson 2000
      do nt = nnut+1,ntyp
       Pdeep(nt) = 0.0    !chl and herbivores
      enddo
      do nt = 1,ndet
       detdeep(nt) = 0.0 !detritus
      enddo
      cardeep(1) = 0.0   !DOC
      cardeep(2) = 2330.0  !DIC uM(C) from Goyet et al (2000)
c
c  Carbon:chl ratios for different adaptation states
      cchl(1) = 25.0
      cchl(2) = 50.0
      cchl(3) = 80.0
c      cchl(1) = 20.0
c      cchl(2) = 60.0
c      cchl(3) = 100.0
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)
      ro2c_DET = -150.0/106.0  !@PL O2:C ratio respiration,degradation,remineralization (uM:uM), (Anderson, 1995)
      ro2c_NH4 = -118.0/106.0  !@PL O2:N ratio ammonium production (uM:uM), (Dunne et al., 2013)
      ro2c_NO3 = -150.0/106.0   !@PL O2:N ratio nitrate production (uM:uM), (Anderson, 1995)

!change: March 15, 2010
       bn = cchl(2)/cnratio         !N:chl ratio (uM/ugl)
       bf = cchl(2)/cfratio         !Fe:chl ratio (nM/ugl)
       cchlratio = cchl(2)          !C:chl ratio (ugl/ugl)
       mgchltouMC = cchlratio/uMtomgm3
c
!!#if NCHL_DEFINED > 0
      if (nchl > 0) then
c  Diatoms
      nt = 1
!change: March 15, 2010
!     rmumax(nt) = 1.50       !u max in /day at 20C
      rmumax(nt) = 2.00       !u max in /day at 20C
#ifdef newpp_May10_2018
      rmumax(nt) = 2.50       !u max in /day at 20C
#endif
#ifdef newpp_May10_2018a
      rmumax(nt) = 2.75       !u max in /day at 20C
#endif
#ifdef newpp_May10_2018b
      rmumax(nt) = 3.00       !u max in /day at 20C
#endif
#ifdef newpp_May10_2018c
      rmumax(nt) = 3.00       !u max in /day at 20C
#endif

#ifdef OBIO_ON_GISSocean
      obio_wsd(nt)    = 0.75  !sinking rate in m/day
#else
      obio_wsd(nt)    = 0.50  !sinking rate in m/day   !!change Oct27,2008
#endif

      rik(1,nt)  = 90.0       !low light-adapted Ik (<50 uE/m2/s)
      rik(2,nt)  = 93.0       !medium light-adapted Ik (50-200 uE/m2/s)
      rik(3,nt)  = 184.0      !high light adapted Ik (>200 uE/m2/s)
      rkn(nt) = 1.0           !M-M half-sat constant for nitrogen
      rks(nt) = 0.2           !M-M half-sat constant for silica
      rkf(nt) = 0.12          !M-M half-sat constant for iron

      endif
!!#endif
!!#if NCHL_DEFINED > 1
      if (nchl > 1) then
c  Chlorophytes
      nt = 2
      rmumax(nt) = rmumax(nt-1)*0.840
      obio_wsd(nt)    = 0.25
      rik(1,nt)  = rik(1,nt-1)*1.077
      rik(2,nt)  = rik(2,nt-1)*0.935
      rik(3,nt)  = rik(3,nt-1)*0.781
      rkn(nt) = rkn(nt-1)*0.75
      rkn(nt) = rkn(nt-1)*0.667   !1/3 distance bet. cocco and dia
      rkf(nt) = rkf(nt-1)*0.835   !midway between cocco's and diatoms
      rkf(nt) = rkf(nt-1)*0.779   !1/3 distance bet. cocco and dia

      endif
!!#endif
!!#if NCHL_DEFINED > 2
      if (nchl > 2) then
c  Cyanobacteria
      nt = 3
      rmumax(nt) = rmumax(nt-2)*0.670
      obio_wsd(nt)    = 0.0085
      rik(1,nt)  = rik(1,nt-2)*0.723
      rik(2,nt)  = rik(2,nt-2)*0.710
      rik(3,nt)  = rik(3,nt-2)*0.256
      rkn(nt) = rkn(nt-2)*0.50
      rkf(nt) = rkf(nt-2)*0.67  !equals cocco

      endif
!!#endif
!!#if NCHL_DEFINED > 3
      if (nchl > 3) then
c  Coccolithophores
      nt = 4
       rmumax(nt) = rmumax(nt-3)*0.755   !E. huxleyi only
c      rmumax(nt) = rmumax(nt-3)*0.781   !E. huxleyi only (no Sunda/Hunts)
      obio_wsd(nt)    = 0.82
      obio_wsd(nt)    = 0.648
      rik(1,nt)  = rik(1,nt-3)*0.623
      rik(2,nt)  = rik(2,nt-3)*0.766
      rik(3,nt)  = rik(3,nt-3)*0.899
      rkn(nt) = rkn(nt-3)*0.5
      rkf(nt) = rkf(nt-3)*0.67

      endif
!!#endif
!!#if NCHL_DEFINED > 4
      if (nchl > 4) then
c  Dinoflagellates
      nt = 5
      rmumax(nt) = rmumax(nt-4)*0.335
      obio_wsd(nt)    = 0.0
      rik(1,nt)  = rik(1,nt-4)*1.321
      rik(2,nt)  = rik(2,nt-4)*1.381
      rik(3,nt)  = rik(3,nt-4)*1.463
      rkn(nt) = rkn(nt-4)*1.0
      rkf(nt) = rkf(nt-4)*0.67

      endif
!!#endif
      do nt = 1,nchl
       obio_wsh(nt) = obio_wsd(nt)/24.0  !convert to m/hr
       obio_wss(nt) = obio_wsd(nt)/sday  !convert to m/s     !July 2016
      enddo

c  Detrital sinking rates m/h  -> m/s   July 2016
      !default
!change: March 10, 2010
!     wsdeth(1) = 30.0/24.0     !nitrogen
!     wsdeth(1) = 20.0/24.0     !nitrogen
      wsdeth(1) = 20.0/sday     !nitrogen

!     wsdeth(2) = 50.0/24.0     !silica
      wsdeth(2) = 50.0/sday     !silica

!     wsdeth(3) = 20.0/24.0     !iron
!change June 1, 2010
!     wsdeth(3) =  5.0/24.0     !iron
      wsdeth(3) =  5.0/sday     !iron
!endofchange
c
c  Detrital remineralization rates /s
!change: March 10, 2010
!     remin(1) = 0.010/24.0            !nitrogen
!     remin(1) = 0.020/24.0            !nitrogen
      remin(1) = 0.020/sday            !nitrogen
#ifdef increaseNremin
!     remin(1) = 0.5/24.0            !nitrogen
      remin(1) = 0.5/sday            !nitrogen
#endif
#ifdef increaseNremin2
!AR5 preprocessor option
!     remin(1) = 0.1/24.0            !nitrogen
      remin(1) = 0.1/sday            !nitrogen
#endif
#ifdef newpp_May10_2018
      remin(1) = 0.025/sday            !nitrogen
#endif
#ifdef newpp_May10_2018a
      remin(1) = 0.0275/sday            !nitrogen
#endif
#ifdef newpp_May10_2018b
      remin(1) = 0.03/sday            !nitrogen
#endif
#ifdef newpp_May10_2018c
      remin(1) = 0.1/sday            !nitrogen
#endif
#ifdef increaseNremin3
!     remin(1) = 0.3/24.0            !nitrogen
      remin(1) = 0.3/sday            !nitrogen
#endif
#ifdef increaseNremin4
!     remin(1) = 0.5/24.0            !nitrogen
      remin(1) = 0.5/sday            !nitrogen
#endif
#ifdef increaseNremin5
!     remin(1) = 0.5/24.0            !nitrogen
      remin(1) = 0.9/sday            !nitrogen
#endif
#ifdef increaseNremin6
      remin(1) = 1.5/sday            !nitrogen
#endif


!     remin(2) = 0.0001/24.0           !silica
      remin(2) = 0.0001/sday           !silica
#ifdef increaseSremin
!AR5 preprocessor option
!     remin(2) = 0.002/24.0           !silica
      remin(2) = 0.002/sday           !silica
#endif
!     remin(3) = 0.020/24.0            !iron
!change June 1, 2010
!     remin(3) = 0.50/24.0            !iron
      remin(3) = 0.50/sday            !iron
!endofchange
#ifdef increaseIremin
!AR5 preprocessor option
!     remin(3) = 0.70/24.0            !iron
      remin(3) = 0.70/sday            !iron
#endif
#ifdef decreaseIremin
      remin(3) = 0.20/sday            !iron
#endif
#ifdef TRACERS_degC
      tdegC = 0.96/sday                ! transfer rate constant for degradable carbon
#endif

!     fescavrate(1) = 2.74E-5/24.0      !low fe scavenging rate/s
      fescavrate(1) = 2.74E-5/sday      !low fe scavenging rate/s     !July 2016
      fescavrate(2) = 50.0*fescavrate(1) !high fe scavenging rate/s -> rate/s
c
c  (originally done inside lidata subroutine of obio_daysetrad)
c  Reads in radiative transfer data: specifically
c  water data (seawater absorption and total scattering coefficients,
c  and chl-specific absorption and total scattering data for
c  several phytoplankton groups).  PAR (350-700) begins at index 3,
c  and ends at index 17.
c     
c  Water data files
!     cfle = cabw                       
!     open(4,file='/explore/nobackup/aromanou/2.0deg/'//cfle
!    .      ,status='old',form='formatted')


c  Phytoplankton group chl-specific absorption and total scattering
c  data.  Chl-specific absorption data is normalized to 440 nm; convert
c  here to actual ac*(440)
!     cfle = cacbc
!     open(4,file='/explore/nobackup/aromanou/2.0deg/'//cfle
!    .      ,status='old',form='formatted')
      call ocalbedo_init
      call openunit('cfle2',iu_bio)
      do ic = 1,6
       read(iu_bio,'(a50)')title
      enddo
      do nt = 1,nchl
       read(iu_bio,'(a50)')title
       do nl = 1,19
        read(iu_bio,30)lambda,sac,sbc
        ac(nt,nl) = sac
        bc(nt,nl) = sbc
       enddo
       do nl = 20,nlt
        ac(nt,nl) = 0.0
        bc(nt,nl) = 0.0
       enddo
      enddo
      call closeunit(iu_bio)
 30   format(i4,2f10.4)

!ifst part from daysetrad.f
c      h = 6.6256E-34   !Plancks constant J sec
       planck = 6.6256E-34   !Plancks constant J sec
       c = 2.998E8      !speed of light m/sec
c      hc = 1.0/(h*c)
       hc = 1.0/(planck*c)
       oavo = 1.0/6.023E23   ! 1/Avogadros number
       hcoavo = hc*oavo
       do nl = npst,npnd
        rlamm = float(lam(nl))*1.0E-9  !lambda in m
        WtoQ(nl) = rlamm*hcoavo        !Watts to quanta conversion
       enddo
       !CDOM absorption exponent
       rlam450 = 450.0
       Sdom = 0.014
       do nl = 1,nlt
        if (lam(nl) .eq. 450)nl450 = nl
        rlam = float(lam(nl))
        excdom(nl) = exp(-Sdom*(rlam-rlam450))
       enddo
       if (nl450.eq.0) stop 'obio_init: nl450=0'

!ifst part from edeu.f
       bbw = 0.5            !backscattering to forward scattering ratio
       rmus = 1.0/0.83      !avg cosine diffuse down
       Dmax = 500.0         !depth at which Ed = 0

       rnn = 1.341
       rmuu = 1.0/0.4             !avg cosine diffuse up
       rbot = 0.0                 !bottom reflectance
       rd = 1.5   !these are taken from Ackleson, et al. 1994 (JGR)
       ru = 3.0

c  Read in factors to compute average irradiance
! (this part originally done inside obio_edeu)
      if (AM_I_ROOT()) then
      print*, '    '
      print*,'Reading factors for mean irradiance at depth...'
      print*,'nh,nch,ncd=',nh,nch,ncd
      endif

      call openunit('facirr',iu_fac)
      do icd=1,ncd
       do ntr=1,5
        do ich=1,nch
         do  ih=1,nh
          read(iu_fac,*)facirr4(ih,ich,ntr,icd)
           facirr(ih,ich,ntr,icd)=1.D0*facirr4(ih,ich,ntr,icd)
         enddo
        enddo
       enddo
      enddo
      call closeunit(iu_fac)

!ifst part from ptend.f
       do k=1,kdm

         do nt=1,nchl
          obio_ws(k,nt) = 0.0
          gro(k,nt) = 0.0
         enddo
 
         do nt=1,ndet
          D_tend(k,nt) = 0.0
          wsdet(k,nt) = 0.0
         enddo
 
         do nt=1,ncar
          C_tend(k,nt) = 0.0
         enddo
       enddo
       do nt=1,nchl
        obio_ws(kdm+1,nt)=0.0
       enddo
       do nt=1,ndet
        wsdet(kdm+1,nt) = 0.0
       enddo

#ifdef exp_wsdiat
! exponential profile coefficients for diatoms 
         adiat_exp = 0.01   !diatoms  a coef
         bdiat_exp = 5.0    !diatoms  b coef
#endif

#ifdef exp_wsdet
! exponential profile coefficients for detritus
         adet_exp(1) = 2.0    !nitrogen a coef
#ifdef w_PL_1
         adet_exp(1) = adet_exp(1) + 0.915
#endif

#ifdef w_PL_2
         adet_exp(1) = adet_exp(1) + 2.0*0.915
#endif

#ifdef w_PL_3
         adet_exp(1) = adet_exp(1) + 5.0*0.915
#endif
         bdet_exp(1) = 3.0    !nitrogen b coef
         adet_exp(2) = 3.0    !silica a coef
         bdet_exp(2) = 6.0    !silica b coef
         adet_exp(3) = 1.0    !iron a coef
         bdet_exp(3) = 2.0    !iron b coef
#endif
 
!read in atmospheric iron deposition (this will also be changed later...)
      if (AM_I_ROOT()) then
      print*, '    '
      print*, 'reading iron data.....'
      print*, '    '
      endif

      if(file_exists('ironflux')) then ! read netcdf format flux
        fid = par_open(ogrid,'ironflux','read')
        call read_dist_data(ogrid,fid,'ironflux',atmFe)
        call par_close(ogrid,fid)
      else
       call init_stream(ogrid,stream_atmFe,'atmFe_inicond_new',
     &  'array', 0d0, 1d9,"linm2m",jyear,jday)
#ifdef Relax2SurfN
        allocate(surfn(ogrid%i_strt:ogrid%i_stop,
     &                       ogrid%j_strt:ogrid%j_stop))
        call bio_surfN('nitrates_inicond',surfn,ogrid,dlatm)
#endif
      endif ! netcdf iron or not

#ifdef TRACERS_Alkalinity
! Alkalinity will be read in from obio_bioinit
! don't do anything here
#else
!read in alkalinity annual mean file
      if (ALK_CLIM.eq.1) then      !read from climatology
        call init_alk(alk,kdm,ogrid)
      else      !set to zero, obio_carbon sets alk=tabar*sal/sal_mean
        alk = 0.
      endif
#endif

#ifdef OBIO_RUNOFF
! read in nutrient concentrations, already regridded to model grid
        if (AM_I_ROOT()) then
        print*, '    '
        print*, 'reading nutrient runoff data.....'
        print*, '    '
        endif
!       filename='rnitr_mflo'
        filename='rnitr_conc'
        fid=par_open(ogrid,filename,'read')
!       call read_dist_data(ogrid,fid,'din',rnitrmflo_loc)
        call read_dist_data(ogrid,fid,'din',rnitrconc_loc)
        call par_close(ogrid,fid)
        filename='rdic_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'dic',rdicconc_loc)
        call par_close(ogrid,fid)
        write(*,*)'reading dic from',filename
        filename='rdoc_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'doc',rdocconc_loc)
        call par_close(ogrid,fid)
        filename='rsili_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'sil',rsiliconc_loc)
        call par_close(ogrid,fid)
        filename='riron_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'fe',rironconc_loc)
        call par_close(ogrid,fid)
        filename='rpoc_conc'
        fid=par_open(ogrid,filename,'read')
        call read_dist_data(ogrid,fid,'poc',rpocconc_loc)
        call par_close(ogrid,fid)
!       filename='ralk_conc'
!       fid=par_open(ogrid,filename,'read')
!       call read_dist_data(ogrid,fid,'alk',ralkconc_loc)
!       call par_close(ogrid,fid)
#endif

! printout some key information
      if (AM_I_ROOT()) then
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'           INITIALIZATION                         '

      write(*,'(a,i5)') 'OBIO - NUMBER OF TRACERS=',ntrac
      write(*,'(a,f10.2)') 'OBIO - dtsrc=',dtsrc

      write(*,*)'ALK_CLIM = ', ALK_CLIM
      if (ALK_CLIM.eq.0) write(*,*) 'ALKALINITY, from SALINITY'
      if (ALK_CLIM.eq.1) write(*,*) 'ALKLNTY, GLODAP annmean'
      if (ALK_CLIM.eq.2) write(*,*) 'ALKALINITY prognostic'

#ifdef pCO2_ONLINE
      print*, 'PCO2 is computed online and not through lookup table'
#else
      print*, 'PCO2 is computed through lookup table'
#endif
      write(*,'(a,4e12.4)')'obio_init, sinking rates for chl (per s): ',
     .    obio_wss(1),obio_wss(2),obio_wss(3),obio_wss(4)
      write(*,'(a,3e12.4)')'obio_init, settl rates for detr (per s): ',
     .      wsdeth(1),  wsdeth(2),  wsdeth(3)

       write(*,'(a,3(f8.6,1x))'), 'OBIO remin rates (per day)=',
     . remin(1)*3600.*24.,remin(2)*3600.*24.,remin(3)*3600.*24.

#ifdef OBIO_RUNOFF
       write(*,*) 'obio-river turned on'
#endif
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      endif

      return
      end subroutine obio_init

c------------------------------------------------------------------------------
#ifdef STANDALONE_OCEAN
      subroutine obio_edaesa_g(filename1,filename2)
!read in eda and esa
!read in a field and convert to ocean grid (using Gary Russel's routine) 

      USE FILEMANAGER, only: openunit,closeunit
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,unpack_data
      USE OCEANR_DIM, only : ogrid

      USE OCEANRES, only : imo,jmo,lmo
      USE OCEAN, only : oDLATM=>DLATM,LMOM=>LMM,ZOE=>ZE,FOCEAN
      USE obio_forc, only: Eda_glob, Esa_glob, Eda,Esa

      implicit none


      integer, parameter :: igrd=360,jgrd=180,kgrd=33
      integer, parameter :: igrd2=288
      integer, parameter :: nmo=12,nhr=12
      integer i,j,k,l,n,lm
      integer iu_file,lgth
      real data1(igrd,jgrd)
      real data2(igrd,jgrd)
      real data_mask(igrd,jgrd)

      integer imon,ihr

      logical vrbos

      character*80 filename1,filename2

!--------------------------------------------------------------
      if ( AM_I_ROOT() ) then
      lgth=len_trim(filename1)
      print*, 'obio-init: reading from file...',filename1(1:lgth)
      call openunit(filename1,iu_file,.false.,.true.)

      do imon=1,nmo; do ihr=1,nhr; do k=1,kgrd;
          do i=1,igrd2; do j=1,jgrd
          read(iu_file,'(e12.4)')Eda_glob(i,j,k,ihr,imon)
          enddo; enddo    ! i,j-loop
      enddo; enddo;enddo    ! k,imon
      print*, 'completed reading ',filename1(1:lgth)
      call closeunit(iu_file)
      lgth=len_trim(filename2)
      print*, 'obio-init: reading from file...',filename2(1:lgth)
      call openunit(filename2,iu_file,.false.,.true.)

      do imon=1,nmo; do ihr=1,nhr; do k=1,kgrd;
          do i=1,igrd2; do j=1,jgrd
          read(iu_file,'(e12.4)')Esa_glob(i,j,k,ihr,imon)
          enddo; enddo    ! i,j-loop
      enddo; enddo;enddo    ! k,imon
      call closeunit(iu_file)
      endif   !AM_I_ROOT
!--------------------------------------------------------------

      call unpack_data(ogrid, Eda_glob, Eda)
      call unpack_data(ogrid, Esa_glob, Esa)

      end subroutine obio_edaesa_g
#endif /*  STANDALONE_OCEAN */

c------------------------------------------------------------------------------

      subroutine bio_surfn(filename,fldo,ogrid,dlatm)
 
      use bio_inicond_mod, only: bio_inicond_read
      use dictionary_mod, only: sync_param
      USE DOMAIN_DECOMP_1D, only : DIST_GRID


      implicit none

      type(DIST_GRID), intent(in) :: ogrid
      character(len=*), intent(in) :: filename
      real, intent(out) :: fldo(ogrid%I_STRT:ogrid%I_STOP,
     .          ogrid%J_STRT:ogrid%J_STOP)
      real, intent(in) :: dlatm      

      real*8, dimension(:, :, :), allocatable :: array

      call bio_inicond_read(filename, array)
      fldo=array(:, :, 1)
      end subroutine bio_surfn


c------------------------------------------------------------------------------


