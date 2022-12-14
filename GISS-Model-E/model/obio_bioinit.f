#include "rundeck_opts.h"

!NOT FOR HYCOM: mo, DXYPO, trmo,n_abioDIC,im,jm, number traces and lmm
!passed to subroutine
      subroutine obio_bioinit(kdm,n_abioDIC,DXYPO,ogrid,
     &                    MO,trmo,num_tracers,
     &                   im,jm,rlon2D,rlat2D,ip,lmm)

 
!note: 
!obio_bioinit is called only for a cold start and reads in INITIAL conditions and interpolates them
!to the ocean grid. Such fields are nitrates,silicate,dic
!obio_init  is called for every start of the run and reads in BOUNDARY conditions and interpolates 
!them to ocean grid. such fields are iron,alkalinity,chlorophyl

! based on /g6/aromanou/Watson_new/BioInit/rstbio.F

!  Makes initialization data files for biological variables.
!  This is for the global model.  To subset, use the routines in
!  /u2/gregg/bio/biodat/subreg.
!  Uses NOAA 2001 atlas for NO3 and SiO2 distributions.
!  Includes initial iron distributions.
 
      USE FILEMANAGER, only: openunit,closeunit


      USE CONSTANT, only: tf

      USE obio_dim
      USE obio_incom
      USE obio_forc, only: avgq
      USE obio_com, only: gcmax, tracer

!@PL 
      USE OFLUXES, only: oAPRESS,ocnatm

      use ocean, only : t3d,s3d,r3d 
!@PL
 
      use obio_com, only: ze
      use bio_inicond_mod, only : bio_inicond
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT ,DIST_GRID

      implicit none

      integer, intent(in) :: num_tracers,kdm,n_abioDIC,
     &                       im,jm
      type(DIST_GRID), intent(in) :: ogrid
      integer, intent(in) :: lmm(im,ogrid%j_strt_halo:ogrid%j_stop_halo)
      real, intent(in) ::DXYPO(jm),
     &                   MO(im,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm),
     &                   ip(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &                      ogrid%j_strt_halo:ogrid%j_stop_halo),
     &                   rlon2D(im,jm) ,rlat2D(im,jm)
      real, intent(out) :: trmo(im,
     &         ogrid%j_strt_halo:ogrid%j_stop_halo,kdm,num_tracers)


      real, parameter, dimension(0:13) :: fer_values=
     &  (/ 0E0,
     &     3.5E-3,      !Antarctic
     &     20.0E-3,     !South Indian
     &     4.5E-3,      !South Pacific
     &     20.0E-3,     !South Atlantic
     &     22.5E-3,     !Equatorial Indian
     &     4.0E-3,      !Equatorial Pacific
     &     20.0E-3,     !Equatorial Atlantic
     &     25.0E-3,     !North Indian
     &     6.0E-3,      !North Central Pacific
     &     20.0E-3,     !North Central Atlantic
     &     6.0E-3,      !North Pacific
     &     10.0E-3,     !North Atlantic
     &     20.0E-3 /)   !Mediterranean

      integer i,j,k,l

!@PL sum vars for abiotic O2 initialization
!@PL var O2SAT0 reference O2 value at equilibrium with atmosphere at sea level pressure = 1 atm (mmol/kg)
!@PL var xo2  mol fraction of O2 in dry air - constant value
!@PL var pH2O  air humidity (atm)
!@PL var Ksol  solubility constant (mmol/kg/atm)
!@PL var O2sat O2 at equilibrium with atmosphere 

      real*8 :: ta2,ta3,ta4,ta5,ksol,xo2,ps,ps2,SLP,
     &        tk100,O2sat,pH2O
      integer nir(nrg),nt,I_0,I_1,J_0,J_1

      integer, ALLOCATABLE, DIMENSION(:,:)   :: ir
      real,  ALLOCATABLE, DIMENSION(:,:,:) :: fer,dic
#ifdef TRACERS_Ocean_O2
          real,  ALLOCATABLE, DIMENSION(:,:,:) :: o2_in
          real*8,  ALLOCATABLE, DIMENSION(:,:,:) :: ta,O2SAT0 
#endif

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP

      ALLOCATE(ir(i_0:i_1,j_0:j_1))
      ALLOCATE(fer(i_0:i_1,j_0:j_1,kdm))
      allocate(dic(i_0:i_1,j_0:j_1,kdm))
#ifdef TRACERS_Ocean_O2
          allocate(o2_in(i_0:i_1,j_0:j_1,kdm))
          allocate(ta(i_0:i_1,j_0:j_1,kdm))
          allocate(O2SAT0(i_0:i_1,j_0:j_1,kdm))
#endif

      tracer(:,:,:,1:ntyp)=0.d0
      Fer(:,:,:) = 0.d0
      dic(:,:,:) = 0.d0

!NOT FOR HYCOM: lmm passed to subroutine
      call bio_inicond('nitrates_inicond',tracer(:,:,:,1),
     &                 kdm,im,ogrid,ip,lmm)
      call bio_inicond('silicate_inicond',tracer(:,:,:,3),
     &                 kdm,im,ogrid,ip,lmm)
#ifdef TRACERS_Ocean_O2
      call bio_inicond('oxygen_inicond',
     &     o2_in(:,:,:),kdm,im,ogrid,ip,lmm) !@PL o2_frac is frac. saturation from a climatology (e.g., WOA2013v2)
                                               !@PL o2_in=o2 concentraiton (mmol/kg) when using GLODAP O2 for initilization
#endif


#ifdef obio_TRANSIENTRUNS
! in the transient runs keep the dic read in from RSF/AIC file
      dic = tracer(:,:,:,15)
#else
! otherwise take from rundeck
!NOT FOR HYCOM: lmm passed to subroutine
      call bio_inicond('dic_inicond',dic,kdm,im,ogrid,ip,lmm)
#endif

#ifdef TRACERS_Alkalinity
!check the value of ntrac
      call init_alk(tracer(:,:,:,ntyp+ndet+ncar+nalk),
     &     kdm,ogrid,im,ip,lmm)
#endif

!     /archive/u/aromanou/Watson_new/BioInit/iron_ron_4x5.asc
!     /archive/u/aromanou/Watson_new/BioInit/CHL_WG_4x5


!!these rno3 and so2_init fields are not correct. There are void points due to
!mismatch of the noaa grid and the hycom grid. To fill in have to do
!the interpolation in matlab (furtuna).
!at the same time use dps and interpolate to layer depths from the model

      do k=1,kdm
       do j=ogrid%j_strt,ogrid%j_stop
        do i=ogrid%i_strt,ogrid%i_stop
          if(tracer(i,j,k,1).le.0.)tracer(i,j,k,1)=0.085d0
          if(tracer(i,j,k,3).le.0.)tracer(i,j,k,3)=0.297d0
          if (dic(i,j,k).le.0.0) dic(i,j,k)=1837.0
          dic(i,j,k)=max(dic(i,j,k),1837.0)   !set minimum =1837
        enddo
       enddo
      enddo

      if (AM_I_ROOT()) 
     .  write(*,'(a,2e12.4)')'BIO: bioinit: dic min-max=',
     .       minval(dic),maxval(dic)

c  Obtain region indicators
c     write(6,*)'calling fndreg...'
      call fndreg(ir,ogrid,im,jm,rlon2D,rlat2D,ip)
 
c  Define Fe:NO3 ratios by region, according to Fung et al. (2000)
c  GBC.  Conversion produces nM Fe, since NO3 is as uM

      do j=ogrid%j_strt,ogrid%j_stop
        do i=ogrid%i_strt,ogrid%i_stop
          fer(i,j,:)=fer_values(ir(i,j))
        enddo  !i-loop
      enddo  !j-loop
 
c  Create arrays 
      if (AM_I_ROOT())
     .  write(6,*)'Creating bio restart data for ',ntyp,' arrays and'
     . ,kdm,'  layers...'

      do j=ogrid%j_strt,ogrid%j_stop
      do i=ogrid%i_strt,ogrid%i_stop
        if (ip(i,j)==0) cycle

        do k=1,kdm
          !Nitrate
          !read earlier from file

          !Ammonium
          tracer(i,j,k,2) = 0.5
          !!!if (ze(i,j,k).gt.4000.d0)tracer(i,j,k,2) = Pdeep(2)

          !Silica
          !read earlier from file
          !!!if (ze(i,j,k).gt. 4000.d0)tracer(i,j,k,2) = Pdeep(3)

          !Iron
          tracer(i,j,k,4) = Fer(i,j,k)*tracer(i,j,k,1)  !Fung et al. 2000
!          if (ir(nw) .eq. 3)then
!           P(i,j,k,4) = 0.04*float(k-1) + 0.2
!           P(i,j,k,4) = 0.06*float(k-1) + 0.2
!           P(i,j,k,4) = 0.08*float(k-1) + 0.2
!           P(i,j,k,4) = min(P(i,j,k,4),0.65)
!           P(i,j,k,4) = min(P(i,j,k,4),0.75)
!          endif
          if (ir(i,j) .eq. 1)then
           tracer(i,j,k,4) = Fer(i,j,k)*0.5*tracer(i,j,k,1)
          endif
          tracer(i,j,k,4) = max(tracer(i,j,k,4),0.01)
          !!!if (ze(i,j,k).gt. 4000.)tracer(i,j,k,4) = Pdeep(4)

          !Herbivores
          do nt = nnut+1,ntyp-nzoo
           tracer(i,j,k,nt) = 0.05
          enddo
          do nt = ntyp-nzoo+1,ntyp
           tracer(i,j,k,nt) = 0.05  !in chl units mg/m3
!          tracer(i,j,k,nt) = 0.05*50.0  !in C units mg/m3
          enddo

         enddo
      end do
      end do


c  Detritus (set to 0 for start up)
      if (AM_I_ROOT()) write(6,*)'Detritus...'
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

      do j=j_0,j_1
       do i=i_0,i_1
         do k=1,kdm
           if (ze(i,j,k)>ze(i,j,k-1)) then
           !only detritus components
             tracer(i,j,k,ntyp+1) = tracer(i,j,k,1)*0.25*cnratio !as carbon
             tracer(i,j,k,ntyp+2) = tracer(i,j,k,3)*0.1
             tracer(i,j,k,ntyp+3) = tracer(i,j,k,4)*0.25
             tracer(i,j,k,ntyp+1) = 0.0
             tracer(i,j,k,ntyp+2) = 0.0
             tracer(i,j,k,ntyp+3) = 0.0
!@PL nondegradable carbon
#ifdef TRACERS_degC
             tracer(i,j,k,ndimndegC) = 0.0
#endif
          endif
         enddo
        enddo
      enddo



c  Carbon (set to 0 for start up)
c   DIC is derived from GLODAP.  Using mean H from exp601,
c   mean DIC for these values is computed.  Surface DIC is taken
c   as the mean for 020m deeper than the mixed layer, converted from
c   uM/kg to uM
      if (AM_I_ROOT()) write(6,*)'Carbon...'
c    conversion from uM to mg/m3
      do j=j_0,j_1
       do i=i_0,i_1
         do k=1,kdm
          tracer(i,j,k,ntyp+ndet+1) = 0.0
          tracer(i,j,k,ntyp+ndet+2) = 0.0
!@PL
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
          tracer(i,j,k,ntyp+ndet+ncar+nalk+no2) = 0.d0
#endif
#ifdef TRACERS_abio_O2
          tracer(i,j,k,ntyp+ndet+ncar+nalk+no2+nabo2) = 0.d0
#endif
          ta(i,j,k) = 0.d0
          O2SAT0(i,j,k) = 0.d0
#endif
!@Pl
         enddo
       enddo
      enddo

      !only carbon components
      do j=j_0,j_1
       do i=i_0,i_1
         if (ip(i,j)==0) cycle
         do k = 1,kdm
          tracer(i,j,k,ntyp+ndet+2) = dic(i,j,k)
     .       * 1024.5d0 * 0.001                               ! convert micromole/kg to mili-mol/m3
!initialize abioDIC
!@PL changed 1024.5 to r3d(k,i,j)
      if (n_abioDIC.ne.0) 
     .    trmo(i,j,k,n_abioDIC) = tracer(i,j,k,ntyp+ndet+2)  ! mili-mol/m3
     .                          * 1.d-06 *
     .                          12.d0*MO(I,J,K)*DXYPO(J)/1024.5d0
#ifdef TRACERS_Ocean_O2
!@PL          !oxygen   
!@PL sum: initialize abiotic O2 with values at equlibrium with atmosphere at all depths, biotic O2 wit GLODAPv2
!@PL auth Paul Lerner
!@PL var O2SAT0 reference O2 value at equilibrium with atmosphere at sea level pressure = 1 atm
!@PL var xo2  mol fraction of O2 in dry air - constant value
!@PL var pH2O  air humidity (atm)
!@PL var Ksol  solubility constant (mmol/kg/atm)
!@PL var O2sat O2 at equilibrium with atmosphere 
!@PL Ts is placeholder !
c      ta(i,j,k) = log((298.15d0-t3d(k,i,j))/(t3d(k,i,j) + tf)) !@PL scaled T from Garcia and Gordon, 1992
c      ta2 = ta(i,j,k)*ta(i,j,k)
c      ta3 = ta2*ta(i,j,k)
c      ta4 = ta3*ta(i,j,k)
c      ta5 = ta4*ta(i,j,k)
c      ps = s3d(k,i,j)*1000.d0
c      ps2 = ps*ps
c      tk100 = 100.d0/(t3d(1,i,j) + tf) !@PL denominator is in K

c      SLP = ((oAPRESS(i,j)/100.d0)+stdslp) !@PL oAPRESS is pressure anomoly in Pa, stdslp in hPa

c          O2SAT0(i,j,k) =(1.d-03)*exp(5.80818d0 + 3.20684d0*ta(i,j,k)  +
c     .     4.11890d0*ta2 + 4.93845d0*ta2 +
c     .     4.93845d0*ta3 + 1.01567d0*ta4 +1.41575d0*ta5 -
c     .     ps * (0.00701211d0 - 0.00725958d0*ta(i,j,k) +
c     .     0.00793334d0*ta2 - 0.000554491d0*ta3) -0.000000132412d0*ps2) !@PL units: mmol/kg

c        xo2 = 2.0946d-1 !@PL mol fraction O2 in dry air

c#ifdef pH2O 
c         pH2O = ocnatm%QSAVG(i,j)*SLP/1013.25d0
c#else
c         pH2O = exp(24.4543d0 - 67.4509d0*(tk100) -
c     .           4.8489d0*log(1.d0/tk100)
c     .        - (0.000544d0 * s3d(1,i,j)*1000.d0)) !@PL water vapor pres used in correction term, units: atm
c#endif
c        Ksol = O2SAT0(i,j,k) / (xO2*((stdslp/1013.25d0)-pH2O)) !@PL units: mmol/kg/atm
        call init_abo2(t3d(k,i,j),s3d(k,i,j),oAPRESS(i,j)
     &           ,ocnatm%QSAVG(i,j),O2sat)   !@PL mmol/kg 

#ifdef TRACERS_bio_O2
        tracer(i,j,k,ntyp+ndet+ncar+nalk+no2) = o2_in(i,j,k)
#endif
#ifdef TRACERS_abio_O2
        tracer(i,j,k,ntyp+ndet+ncar+nalk+no2+nabo2)=O2sat ! initial conditions for abiotic O2 are saturation
#endif
!@PLtest
!@PL         print * , '#biotracers = ',ntrac
!@PLtest
  
#endif


         enddo
c         car(i,j,k,1) = 3.0  !from Bissett et al 1999 (uM(C))
c         car(i,j,k,1) = 0.0  !from Walsh et al 1999
        enddo
      enddo



c  Light saturation data
      avgq = 0.0

      do j=j_0,j_1
        do i=i_0,i_1
          if (ip(i,j)==0) cycle
          do k=1,kdm
            avgq(i,j,k) = 25.0
          enddo
        enddo
      enddo

c  Coccolithophore max growth rate
      do j=j_0,j_1
       do k=1,kdm
        do i=i_0,i_1
          gcmax(i,j,k) = 0.0
        enddo
       enddo
      enddo
 
      return
      end subroutine obio_bioinit

      subroutine init_alk(alk,kdm,ogrid,im,ip,lmm)
      use obio_com, only : ze
      use bio_inicond_mod, only: bio_inicond
      USE DOMAIN_DECOMP_1D, only : DIST_GRID
      
      implicit none
      
      integer, intent(in) ::  kdm,im
      type(DIST_GRID), intent(in) :: ogrid
      integer, intent(in) :: lmm(im,ogrid%j_strt_halo:ogrid%j_stop_halo)
      real, intent(in) :: ip(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &                      ogrid%j_strt_halo:ogrid%j_stop_halo)
      real, intent(out) :: alk(ogrid%i_strt:ogrid%i_stop,
     &                         ogrid%j_strt:ogrid%j_stop,kdm)
  
      integer :: i, j, k

      call bio_inicond('alk_inicond',alk,kdm,im,ogrid,ip,lmm)
      do k=1,kdm
      do j=ogrid%j_strt,ogrid%j_stop
      do i=ogrid%i_strt,ogrid%i_stop
        if (alk(i,j,k).lt.0.) then
          if (ze(i, j, k).le.150.) alk(i,j,k)=2172.      !init neg might be under ice,
          if (ze(i, j, k).gt.150. .and. ze(i, j, k).lt.1200.) !meq/m^3
     &                           alk(i,j,k)=2200.
          if (ze(i, j, k).ge.1200.) alk(i,j,k)=2300.
        endif
      enddo
      enddo
      enddo
      return
      end subroutine init_alk

!@PL initialize o2 for abiotic O2 at saturation

      subroutine init_abo2(t,s,oap,pH2Oin,O2sat)
      USE CONSTANT, only: tf
      USE obio_incom


      implicit none

      real, intent(in) ::  t,s,oap,pH2Oin
      real, intent(out) :: O2sat

      real :: ta,ta2,ta3,ta4,ta5,ps,ps2,tk100,SLP
     &        ,xo2,O2SAT0,pH2O,Ksol

      !@PL sum: initialize abiotic O2 with values at equlibrium with atmosphere at all depths, biotic O2 wit GLODAPv2
!@PL auth Paul Lerner
!@PL var O2SAT0 reference O2 value at equilibrium with atmosphere at sea level pressure = 1 atm
!@PL var xo2  mol fraction of O2 in dry air - constant value
!@PL var pH2O  air humidity (atm)
!@PL var Ksol  solubility constant (mmol/kg/atm)
!@PL var O2sat O2 at equilibrium with atmosphere 
!@PL t is temperature (deg C), s is salinity (psu) !
      ta = log((298.15d0-t)/(t + tf)) !@PL scaled T from Garcia and Gordon, 1992
      ta2 = ta*ta
      ta3 = ta2*ta
      ta4 = ta3*ta
      ta5 = ta4*ta
      ps = s
      ps2 = ps*ps
      tk100 = 100.d0/(t + tf) !@PL denominator is in K

      SLP = ((oap/100.d0)+stdslp) !@PL oap is pressure anomoly in Pa, stdslp in hPa

          O2SAT0 =(1.d-03)*exp(5.80818d0 + 3.20684d0*ta  +
     .     4.11890d0*ta2 + 4.93845d0*ta2 +
     .     4.93845d0*ta3 + 1.01567d0*ta4 +1.41575d0*ta5 -
     .     ps * (0.00701211d0 - 0.00725958d0*ta +
     .     0.00793334d0*ta2 - 0.000554491d0*ta3) -0.000000132412d0*ps2) !@PL units: mmol/kg

        xo2 = 2.0946d-1 !@PL mol fraction O2 in dry air

#ifdef pH2O 
         pH2O = pH2Oin*SLP/1013.25d0
#else
         pH2O = exp(24.4543d0 - 67.4509d0*(tk100) -
     .           4.8489d0*log(1.d0/tk100)
     .        - (0.000544d0 * s*1000.d0)) !@PL water vapor pres/ used in correction term, units: atm
#endif
        Ksol = O2SAT0 / (xo2*((stdslp/1013.25d0)-pH2O)) !@PL units: mmol/kg/atm
        O2sat = Ksol * ((SLP/1013.25d0) - pH2O)*xo2   !@PL mmol/kg 

      end subroutine init_abo2
c------------------------------------------------------------------------------
      subroutine fndreg(ir,ogrid,im,jm,rlon2D,rlat2D,ip)

c  Finds nwater indices corresponding to significant regions,
c  and defines arrays.  Variables representative of the regions will 
c  later be kept in these array indicators.
c  Regions are defined as follows:
c        1 -- Antarctic
c        2 -- South Indian
c        3 -- South Pacific
c        4 -- South Atlantic
c        5 -- Equatorial Indian
c        6 -- Equatorial Pacific
c        7 -- Equatorial Atlantic
c        8 -- North Indian
c        9 -- North Central Pacific
c       10 -- North Central Atlantic
c       11 -- North Pacific
c       12 -- North Atlantic
c       13 -- Mediterranean/Black Seas
 

      USE obio_dim
      use obio_com, only: ze
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,DIST_GRID

      implicit none

      integer, intent(in) :: im,jm
      type(DIST_GRID), intent(in) :: ogrid
      integer, intent(out) :: ir(ogrid%i_strt:ogrid%i_stop,
     &                           ogrid%j_strt:ogrid%j_stop)
      real, intent(in) :: rlon2D(im,jm),rlat2D(im,jm),
     &                    ip(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &                      ogrid%j_strt_halo:ogrid%j_stop_halo)

      integer i,j,l
      integer iant,isin,ispc,isat,iein,iepc,ieat,incp
     .       ,inca,inat,imed,inin,inpc,nr,ntot

      real rlat,rlon

      integer nir(nrg)

      real antlat,rnpolat
      data antlat,rnpolat /-40.0, 40.0/


c  Set up indicators
      iant = 1   !antarcic region
      isin = 2   !south indian ocean
      ispc = 3   !south pacific
      isat = 4   !south atlantic
      iein = 5   !equatorial indian ocean
      iepc = 6   !equatorial pacific ocean
      ieat = 7   !equatorial atlantic ocean
      inin = 8   !north indian ocean
      incp = 9   !north-central pacific
      inca = 10  !north-central atlantic
      inpc = 11  !north pacific
      inat = 12  !north atlantic
      imed = 13  !mediterranean/black sea
 
c  Initialize region indicator array
      ir = 0
      nir = 0
 
c  Find nwater values corresponding to regions
       do j=ogrid%j_strt,ogrid%j_stop
       do i=ogrid%i_strt,ogrid%i_stop

        rlon=rlon2D(i,j)
        rlat=rlat2D(i,j)
 
         if (AM_I_ROOT()) then
         if (i.eq.1) then
           write(*,*) 'rlon and rlat',i,j,rlon, rlat
         endif
         endif

        if (rlon .gt. 180)rlon = rlon-360.0

        if (ip(i,j)==0) cycle
        if (ze(i,j,1)<=ze(i,j,0)) cycle

c   Antarctic region
        if (rlat .le. antlat)then
         ir(i,j) = iant
         nir(ir(i,j)) = nir(ir(i,j))+1
        endif

c   South Indian region
        if (rlat .le. -30.0 .and. rlat .gt. antlat)then
         if (rlon .gt. 20.0 .and. rlon .lt. 150.0)then
          ir(i,j) = isin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -10.0 .and. rlat .gt. -30.0)then
         if (rlon .gt. 20.0 .and. rlon .lt. 142.5)then
          ir(i,j) = isin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   South Pacific region
        if (rlat .le. -10.0 .and. rlat .gt. antlat)then
         if (rlon .le. -70.0 .and. rlon .ge. -180.0)then
          ir(i,j) = ispc
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
         if (rlon .ge. 150.0 .and. rlon .lt. 180.0)then
          ir(i,j) = ispc
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -10.0 .and. rlat .gt. -30.0)then
         if (rlon .ge. 142.5 .and. rlon .lt. 180.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = ispc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   South Atlantic region
        if (rlat .le. -10.0 .and. rlat .gt. antlat)then
         if (rlon .gt. -70.0 .and. rlon .le. 20.0)then
          ir(i,j) = isat
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Equatorial Indian Ocean region
        if (rlat .le. -8.0 .and. rlat .gt. -10.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 142.5)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -6.0 .and. rlat .gt. -8.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 108.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iein
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif
        if (rlat .le. -4.0 .and. rlat .gt. -6.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 105.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -2.0 .and. rlat .gt. -4.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 103.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 2.0 .and. rlat .gt. -2.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 104.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 4.0 .and. rlat .gt. 2.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 103.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 7.0 .and. rlat .gt. 4.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 102.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .lt. 10.0 .and. rlat .gt. 7.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 99.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Equatorial Pacific region
        if (rlat .lt. 10.0 .and. rlat .gt. -10.0)then
         if (rlon .gt. 90.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iepc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
         if (rlon .gt. -180.0 .and. rlon .lt. -70.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iepc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   Equatorial Atlantic region
        if (rlat .lt. 10.0 .and. rlat .gt. -10.0)then
         if (rlon .lt. 20.0 .and. rlon .ge. -74.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = ieat
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   North Indian Ocean
        if (rlat .le. 30.0 .and. rlat .ge. 10.0)then
         if (rlon .gt. 20.0 .and. rlon .lt. 99.0)then
          ir(i,j) = inin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   North Central Pacific region
        if (rlat .ge. 10.0 .and. rlat .le. rnpolat)then
         if (rlon .gt. 99.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
         if (rlon .le. -100.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .ge. 10.0 .and. rlat .lt. 18.0)then
         if (rlon .gt. -100.0 .and. rlon .le. -90.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .ge. 10.0 .and. rlat .lt. 14.0)then
         if (rlon .gt. -90.0 .and. rlon .le. -84.5)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Mediterranean/Black Seas region
        if (rlat .ge. 30.0 .and. rlat .le. 43.0)then
         if (rlon .ge. -5.5 .and. rlon .lt. 60.0)then
          ir(i,j) = imed
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .gt. 43.0 .and. rlat .le. 48.0)then
         if (rlon .ge. 0.0 .and. rlon .lt. 70.0)then
          ir(i,j) = imed
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   North Central Atlantic region
        if (rlat .ge. 7.0 .and. rlat .le. 40.0)then !pickup Carib.Sea
         if (ir(i,j) .eq. 0)then
          ir(i,j) = inca
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c  North Pacific and Atlantic
        if (rlat .gt. rnpolat)then
c  North Pacific
         if (rlon .lt. -105.0)then
          ir(i,j) = inpc
          nir(ir(i,j)) = nir(ir(i,j))+1
         else if (rlon .gt. 120.0)then
          ir(i,j) = inpc
          nir(ir(i,j)) = nir(ir(i,j))+1
         else
c  North Atlantic
          if (ir(i,j) .ne. imed)then
           ir(i,j) = inat
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

       end do
       end do
 
c  Set nir to minimum 1 value to prevent error in division
      do nr = 1,nrg
       nir(nr) = max(nir(nr),1)
      enddo
 
c  Total up points for check
      if (AM_I_ROOT()) then
      ntot = 0
      do nr = 1,nrg
       ntot = ntot + nir(nr)
       write(6,*)'Region, no. points = ',nr,nir(nr)
      enddo
      write(6,*)'Total ocean points = ',ntot
      endif
 
      return
      end subroutine fndreg
c------------------------------------------------------------------------------

      subroutine remap1d_plm(yold,xold,kold,ynew,xnew,knew,vrbos,i,j)
c
c --- consider two stepwise constant functions -yold,ynew- whose
c --- discontinuities are at abscissa values -xold,xnew- respectively.
c --- treat -ynew- as unknown. solve for -ynew- under the condition that
c --- the integral over y*dx is preserved (integration based on PLM).
c
      implicit none
      integer,intent(IN) :: kold,knew
      integer,intent(IN) :: i,j                !  current location in horiz.grid
      real,intent(IN)    :: yold(kold),xold(kold+1),xnew(knew+1)
      real,intent(OUT)   :: ynew(knew)
      logical,intent(IN) :: vrbos        !  if true, print diagnostics
      integer k,ko,n
      real colin,clout,slope,wgta,wgtb,wgtc,yinteg,
     .     ylft(kold),yrgt(kold),xlo,xhi,ra,rb,ya,yb,q,plmslp,
     .     yrka,ylk,yrk,ylkb
      external plmslp
      logical at_top
      real,parameter    :: onemu=1.e-6, acurcy=1.e-6, flag=-999.
      integer,parameter :: iter=2
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- old profile:     x         y',
     .  (k,xold(k),yold(k),k=1,kold),kold+1,xold(kold+1)
 101  format (2i5,a/(i30,f12.1,f10.2))
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      colin=0.
      clout=0.
      do 3 k=1,kold
 3    colin=colin+yold(k)*(xold(k+1)-xold(k))
c
c --- replace each flat segment of stairstep curve by
c --- a slanting segment, using PLM-type limiters.
c
      ylft(   1)=yold(   1)
      yrgt(   1)=yold(   1)
      ylft(kold)=yold(kold)
      yrgt(kold)=yold(kold)
      do 6 n=1,iter                        !  iterate to optimize limiters
      do 2 k=2,kold-1
      if (n.eq.1) then
        yrka=yold(k-1)
        ylk=yold(k)
        yrk=yold(k)
        ylkb=yold(k+1)
      else
        yrka=yrgt(k-1)
        ylk=ylft(k)
        yrk=yrgt(k)
        ylkb=ylft(k+1)
      end if
      wgta=max(onemu,xold(k  )-xold(k-1))
      wgtb=max(onemu,xold(k+1)-xold(k  ))
      wgtc=max(onemu,xold(k+2)-xold(k+1))
      if (k.eq.     1) wgta=onemu
      if (k.eq.kold-1) wgtc=onemu
      slope=plmslp((wgtb*yrka+wgta*ylk)/(wgtb+wgta),
     .     yold(k),(wgtb*ylkb+wgtc*yrk)/(wgtb+wgtc))
      ylft(k)=yold(k)-slope
 2    yrgt(k)=yold(k)+slope
      if (vrbos) print '(8x,a,12x,a,5x,a/(i3,f9.2,5x,2f9.2))',
     .  'y','ylft','yrgt',(ko,yold(ko),ylft(ko),yrgt(ko),ko=1,kold)
 6    continue
      if (vrbos) print '(a/(10f8.2))',
     .  '  target values:',(ynew(ko),ko=1,knew)
c
c --- y in k-th interval now varies from ylft at xold(k) to yrgt at xold(k+1).
c --- find ynew(k) by requiring
c --- that the integral over y*dx from xnew(k) to xnew(k+1) be preserved.
c
      at_top=.true.
      do 4 k=1,knew
      yinteg=0.
      xlo=xnew(k  )
      xhi=xnew(k+1)
ccc      if (vrbos) print '(a,2f9.3)','xlo,xhi =',xlo,xhi
      if (xhi.gt.xlo) then
        at_top=.false.
        do 5 ko=1,kold
        if (xold(ko).ge.xhi) go to 1
c --- integrate over sloping portions of y(x) curve:
        ra=max(xlo,min(xhi,xold(ko  )))
        rb=max(xlo,min(xhi,xold(ko+1)))
cnat    ya=ylft(k)
cnat    yb=yrgt(k)
        ya=ylft(ko)
        yb=yrgt(ko)
        wgta=flag
        wgtb=flag
        if (xold(ko+1).ne.xold(ko)) then
          if (ra.ge.xold(ko).and.ra.le.xold(ko+1)) then
            wgta=(xold(ko+1)-ra)/(xold(ko+1)-xold(ko))
            ya=ylft(ko)*wgta+yrgt(ko)*(1.-wgta)
          end if
          if (rb.ge.xold(ko).and.rb.le.xold(ko+1)) then
            wgtb=(xold(ko+1)-rb)/(xold(ko+1)-xold(ko))
            yb=ylft(ko)*wgtb+yrgt(ko)*(1.-wgtb)
          end if
        end if
        yinteg=yinteg+.5*(ya+yb)*(rb-ra)
ccc        if (vrbos) print '(2i4,4f9.3,3f11.1)',
ccc     .    k,ko,ra,rb,wgta,wgtb,ya,yb,yinteg
 5      continue
        yinteg=yinteg+yb*(xhi-rb)
ccc        if (vrbos) print '(2i4,4f9.3,3f11.1)',
ccc     .    k,0,rb,xhi,wgta,wgtb,yb,yb,yinteg
 1      ynew(k)=yinteg/(xhi-xlo)
      else if (at_top) then
        ynew(k)=yold(   1)
      else                              !  at end
        ynew(k)=yold(kold)
      end if
ccc      if (vrbos) print '(a,f11.1)','ynew =',ynew(k)
      clout=clout+ynew(k)*(xnew(k+1)-xnew(k))
 4    continue
c
      if (abs(clout-colin).gt.acurcy*10.*xold(kold+1))
     .  write (*,100) i,j,' remap1d - column intgl.error',
     .   colin,clout,(clout-colin)/colin
 100  format (2i5,a,2es14.6,es9.1)
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- new profile:     x         y',
     .  (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)

      return
      end
c
c
      function plmslp(ylft,ymid,yrgt)
c
c --- get slope at point 'ymid' for piecewise linear interpolation
      if (ymid.le.min(ylft,yrgt) .or.
     .    ymid.ge.max(ylft,yrgt)) then
        plmslp=0.
      else if ((yrgt-ylft)*(ylft+yrgt-2.*ymid).gt.0.) then
        plmslp=ymid-ylft
      else
        plmslp=yrgt-ymid
      end if
      return
      end

      subroutine remap1d_pcm(yold,xold,kold,ynew,xnew,knew,vrbos,i,j)
c
c --- consider two stepwise constant functions -yold,ynew- whose
c --- discontinuities are at abscissa values -xold,xnew- respectively.
c --- treat -ynew- as unknown. solve for -ynew- under the condition that
c --- the integral over y*dx is preserved (integration based on PCM).
c
      implicit none
      integer,intent(IN) :: kold,knew
      integer,intent(IN) :: i,j                !  current location in horiz.grid
      real,intent(IN)    :: yold(kold),xold(kold+1),xnew(knew+1)
      real,intent(OUT)   :: ynew(knew)
      logical,intent(IN) :: vrbos        !  if true, print diagnostics
      integer k,ko
      real colin,clout,colmx,yinteg,xlo,xhi,xa,xb
      logical at_top
      real,parameter    :: acurcy=1.e-6
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- old profile:     x         y',
     .  (k,xold(k),yold(k),k=1,kold),kold+1,xold(kold+1)
 101  format (2i5,a/(i30,f12.1,f10.3))
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      colin=0.
      clout=0.
      colmx=0.
      do 3 k=1,kold
      colmx=max(colmx,abs(yold(k)))
 3    colin=colin+yold(k)*(xold(k+1)-xold(k))
c
c --- find ynew(k) by requiring that the integral over y*dx
c --- from xnew(k) to xnew(k+1) be preserved.
c
      at_top=.true.
      do 4 k=1,knew
      yinteg=0.
      xlo=xnew(k  )
      xhi=xnew(k+1)
      if (xhi.gt.xlo) then
        at_top=.false.
        xb=xlo
        do 5 ko=1,kold
        xa=xb
        xb=min(xhi,max(xlo,xold(ko+1)))
        yinteg=yinteg+yold(ko)*(xb-xa)
        if (xa.ge.xhi) go to 1
 5      continue
        xa=xb
        xb=xhi
        yinteg=yinteg+yold(kold)*(xb-xa)
 1      ynew(k)=yinteg/(xhi-xlo)
        clout=clout+yinteg
      else if (at_top) then
        ynew(k)=yold(   1)
      else                              !  at end
        ynew(k)=yold(kold)
      end if
 4    continue
c
      if (abs(clout-colin).gt.acurcy*colmx*xold(kold+1))
     .  write (*,100) i,j,' remap1d - column intgl.error',
     .   colin,clout,(clout-colin)/colin
 100  format (2i5,a,2es14.6,es9.1)
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- new profile:     x         y',
     .  (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)
      return
      end
