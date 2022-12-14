c   -----------------------------------------------------------------------------
      module hycom_arrays
      implicit none
      real, allocatable ::
     . u(:,:,:),v(:,:,:)		! velocity components
     .,dp(:,:,:),dpold(:,:,:)		! layer thickness
     .,dpu(:,:,:),dpv(:,:,:)		! layer thickness at u,v points
     .,p(:,:,:)				! interface pressure
     .,pu(:,:,:),pv(:,:,:)		! interface pres. at u,v points
     .,latij(:,:,:),lonij(:,:,:)	! latitude/longitude
     .,corio(:,:)			! coriolis parameter
     .,potvor(:,:)			! potential vorticity
     .,temp(:,:,:)			! temperature
     .,saln(:,:,:)			! salinity
     .,th3d(:,:,:)			! potential density
     .,thstar(:,:,:)			! virtual potential density
     .,wgtkap(:,:)			! scale factor of ref2
     .,psikk(:,:)			! initl.montg.pot. in bottom layer
     .,thkk(:,:)			! initl.thstar in bottom layer
     .,dpmixl(:,:,:)			! mixed layer depth
     .,srfhgt(:,:)			! sea surface height
c
      real, allocatable ::
     . montg(:,:,:)			! montgomery potential
     .,defor1(:,:),defor2(:,:)		! deformation components
     .,ubavg(:,:),vbavg(:,:)	        ! barotropic velocity
     .,pbavg(:,:)			! barotropic pressure
     .,ubrhs(:,:),vbrhs(:,:)		! rhs of barotropic u,v eqns.
     .,utotm(:,:),vtotm(:,:)		! total (barotrop.+baroclin.)...
     .,utotn(:,:),vtotn(:,:)		! ...velocities at 2 time levels
     .,uflux(:,:),vflux(:,:)		! horizontal mass fluxes
     .,uflux1(:,:),vflux1(:,:)		! more mass fluxes
     .,uflux2(:,:),vflux2(:,:)		! more mass fluxes
     .,uflux3(:,:),vflux3(:,:)		! more mass fluxes
     .,uflx(:,:,:),vflx(:,:,:)		! more mass fluxes
     .,bolusu(:,:,:),bolusv(:,:,:)	! thickness (bolus) fluxes
c
      real, allocatable ::
     .   uav(:,:,:),  vav(:,:,:)
     .,dpuav(:,:,:),dpvav(:,:,:)
     .,temav(:,:,:),salav(:,:,:)
     .,th3av(:,:,:), dpav(:,:,:)
     .,ubavav(:,:),vbavav(:,:)
     .,pbavav(:,:),sfhtav(:,:)
     .,uflxav(:,:,:),vflxav(:,:,:)
     .,diaflx(:,:,:)                    ! time integral of diapyc.flux
     .,salflav(:,:),brineav(:,:)
     .,eminpav(:,:),surflav(:,:)
     .,tauxav(:,:),tauyav(:,:)
     .,ufxcum(:,:,:),vfxcum(:,:,:)
     .,dpinit(:,:,:),dpmxav(:,:),oiceav(:,:)
c
      real, allocatable ::
     . util1(:,:),util2(:,:)		! arrays for temporary storage
     .,util3(:,:),util4(:,:)		! arrays for temporary storage
c
     .,scpx(:,:),scpy(:,:)		! mesh size at p pts in x,y dir.
     .,scux(:,:),scuy(:,:)		! mesh size at u pts in x,y dir.
     .,scvx(:,:),scvy(:,:)		! mesh size at v pts in x,y dir.
     .,scqx(:,:),scqy(:,:)		! mesh size at q pts in x,y dir.
     .,scu2(:,:),scv2(:,:)		! grid box size at u,v pts
     .,scp2(:,:),scq2(:,:)		! grid box size at p,q pts
     .,scuxi(:,:),scvyi(:,:)		! inverses of scux,scvy
     .,scp2i(:,:),scq2i(:,:)		! inverses of scp2,scq2
c
     .,pgfx(:,:),pgfy(:,:)		! horiz. presssure gradient
     .,gradx(:,:),grady(:,:)		! horiz. presssure gradient
     .,depthu(:,:),depthv(:,:)		! bottom pres. at u,v points
     .,pvtrop(:,:)			! pot.vort. of barotropic flow
     .,depths(:,:)			! water depth
     .,drag(:,:)			! bottom drag
     .,glue(:,:)			! regional viscosity enhancement
     .,dampu(:,:),dampv(:,:)		! coastal wave damping coeff.
c
      real, allocatable ::
     . uja(:,:),ujb(:,:)		! velocities at lateral
     .,via(:,:),vib(:,:)		!          neighbor points
     .,pbot(:,:)			! bottom pressure at t=0
     .,tracer(:,:,:,:)			! inert tracer (optional)
     .,diadff(:,:,:)			! effective diapycnal diffusivity
     .,tprime(:,:)			! temp.change due to surflx
     .,sgain(:,:)			! salin.changes from diapyc.mix.
     .,surflx(:,:)			! surface thermal energy flux
     .,salflx(:,:)			! surface salinity flux
     .,thkice(:,:)			! grid-cell avg. ice thknss (cm)
     .,covice(:,:)			! ice coverage (rel.units)
c    .,temice(:,:)			! ice surf.temp.
c    .,odhsi(:,:)			! heat borrowed from frozen
     .,odmsi(:,:)			! newly formed ice
     .,omlhc(:,:)
     .,dmfz(:,:)			! ice mass due to freezing
     .,thmix(:,:),tmix(:,:),smix(:,:),umix(:,:),vmix(:,:)
c
     .,ylo(:,:),southfl(:,:),xlo(:,:),eastfl(:,:)

      integer, allocatable, dimension (:,:) ::
     .  klist         !k-index of layer below mixl'r
     .,ijlist         !global ij index
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing functions
c --- windf       include wind stress in forcing functions
c --- relax       activate lateral boundary nudging
c --- trcout      advect tracer and save results in history/restart file
c --- dotrcr      perform column physics operations on tracer array(s)
c
      real, allocatable ::
     . taux(:,:)                          !  wind stress in x direction
     .,tauy(:,:)                          !  wind stress in y direction
c    .,wndspd(:,:,:)                      !  wind speed (tke source)
c    .,airtmp(:,:,:)                      !  pseudo air temperature
c    .,vapmix(:,:,:)                      !  atmosph. vapor mixing ratio
c    .,oprec(:,:)                         !  precipitation
c    .,oevap(:,:)                         !  evaportation
     .,oemnp(:,:)                         !  e - p 
     .,oflxa2o(:,:),oice(:,:)
     .,ustar(:,:)                         ! surface friction velocity
     .,ustarb(:,:)                        ! bottom friction velocity
     .,osalt(:,:)                         ! saltflux from SI(kg/m*m)
     .,freshw(:,:)                        !  river & glacier runoff
     .,diafor(:,:)                        !  imposed diapycnal forcing
c
      contains

      subroutine alloc_hycom_arrays
      USE HYCOM_DIMEN, only : I_0H,I_1H,J_0H,J_1H,kdm,ntrcr

c     print *,"alloc_hycom_arrays: halos", I_0H,I_1H,J_0H,J_1H

      allocate( 
     . u(I_0H:I_1H,J_0H:J_1H,2*kdm),v(I_0H:I_1H,J_0H:J_1H,2*kdm) 
     .,dp(I_0H:I_1H,J_0H:J_1H,2*kdm),dpold(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,dpu(I_0H:I_1H,J_0H:J_1H,2*kdm),dpv(I_0H:I_1H,J_0H:J_1H,2*kdm) 
     .,p(I_0H:I_1H,J_0H:J_1H,kdm+1) 
     .,pu(I_0H:I_1H,J_0H:J_1H,kdm+1),pv(I_0H:I_1H,J_0H:J_1H,kdm+1) 
     .,latij(I_0H:I_1H,J_0H:J_1H,4),lonij(I_0H:I_1H,J_0H:J_1H,4) 
     .,corio(I_0H:I_1H,J_0H:J_1H) 
     .,potvor(I_0H:I_1H,J_0H:J_1H) 
     .,temp(I_0H:I_1H,J_0H:J_1H,2*kdm) 
     .,saln(I_0H:I_1H,J_0H:J_1H,2*kdm) 
     .,th3d(I_0H:I_1H,J_0H:J_1H,2*kdm) 
     .,thstar(I_0H:I_1H,J_0H:J_1H,2*kdm) 
     .,wgtkap(I_0H:I_1H,J_0H:J_1H) 
     .,psikk(I_0H:I_1H,J_0H:J_1H) 
     .,thkk(I_0H:I_1H,J_0H:J_1H) 
     .,dpmixl(I_0H:I_1H,J_0H:J_1H,2) 
     .,srfhgt(I_0H:I_1H,J_0H:J_1H) ) 
c 
      allocate( 
     . montg(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,defor1(I_0H:I_1H,J_0H:J_1H),defor2(I_0H:I_1H,J_0H:J_1H) 
     .,ubavg(I_0H:I_1H,J_0H:J_1H),vbavg(I_0H:I_1H,J_0H:J_1H) 
     .,pbavg(I_0H:I_1H,J_0H:J_1H) 
     .,ubrhs(I_0H:I_1H,J_0H:J_1H),vbrhs(I_0H:I_1H,J_0H:J_1H) 
     .,utotm(I_0H:I_1H,J_0H:J_1H),vtotm(I_0H:I_1H,J_0H:J_1H) 
     .,utotn(I_0H:I_1H,J_0H:J_1H),vtotn(I_0H:I_1H,J_0H:J_1H) 
     .,uflux(I_0H:I_1H,J_0H:J_1H),vflux(I_0H:I_1H,J_0H:J_1H) 
     .,uflux1(I_0H:I_1H,J_0H:J_1H),vflux1(I_0H:I_1H,J_0H:J_1H) 
     .,uflux2(I_0H:I_1H,J_0H:J_1H),vflux2(I_0H:I_1H,J_0H:J_1H) 
     .,uflux3(I_0H:I_1H,J_0H:J_1H),vflux3(I_0H:I_1H,J_0H:J_1H) 
     .,uflx(I_0H:I_1H,J_0H:J_1H,kdm),vflx(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,bolusu(I_0H:I_1H,J_0H:J_1H,kdm),bolusv(I_0H:I_1H,J_0H:J_1H,kdm) ) 
c 
      allocate( 
     .   uav(I_0H:I_1H,J_0H:J_1H,kdm),  vav(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,dpuav(I_0H:I_1H,J_0H:J_1H,kdm),dpvav(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,temav(I_0H:I_1H,J_0H:J_1H,kdm),salav(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,th3av(I_0H:I_1H,J_0H:J_1H,kdm), dpav(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,ubavav(I_0H:I_1H,J_0H:J_1H),vbavav(I_0H:I_1H,J_0H:J_1H)
     .,pbavav(I_0H:I_1H,J_0H:J_1H),sfhtav(I_0H:I_1H,J_0H:J_1H) 
     .,uflxav(I_0H:I_1H,J_0H:J_1H,kdm),vflxav(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,diaflx(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,salflav(I_0H:I_1H,J_0H:J_1H),brineav(I_0H:I_1H,J_0H:J_1H)
     .,eminpav(I_0H:I_1H,J_0H:J_1H) 
     .,surflav(I_0H:I_1H,J_0H:J_1H) 
     .,tauxav(I_0H:I_1H,J_0H:J_1H),tauyav(I_0H:I_1H,J_0H:J_1H) 
     .,ufxcum(I_0H:I_1H,J_0H:J_1H,kdm),vfxcum(I_0H:I_1H,J_0H:J_1H,kdm)
     .,dpinit(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,dpmxav(I_0H:I_1H,J_0H:J_1H),oiceav(I_0H:I_1H,J_0H:J_1H) ) 
c 
      allocate( 
     . util1(I_0H:I_1H,J_0H:J_1H),util2(I_0H:I_1H,J_0H:J_1H) 
     .,util3(I_0H:I_1H,J_0H:J_1H),util4(I_0H:I_1H,J_0H:J_1H)
c 
     .,scpx(I_0H:I_1H,J_0H:J_1H),scpy(I_0H:I_1H,J_0H:J_1H) 
     .,scux(I_0H:I_1H,J_0H:J_1H),scuy(I_0H:I_1H,J_0H:J_1H) 
     .,scvx(I_0H:I_1H,J_0H:J_1H),scvy(I_0H:I_1H,J_0H:J_1H) 
     .,scqx(I_0H:I_1H,J_0H:J_1H),scqy(I_0H:I_1H,J_0H:J_1H) 
     .,scu2(I_0H:I_1H,J_0H:J_1H),scv2(I_0H:I_1H,J_0H:J_1H) 
     .,scp2(I_0H:I_1H,J_0H:J_1H),scq2(I_0H:I_1H,J_0H:J_1H) 
     .,scuxi(I_0H:I_1H,J_0H:J_1H),scvyi(I_0H:I_1H,J_0H:J_1H) 
     .,scp2i(I_0H:I_1H,J_0H:J_1H),scq2i(I_0H:I_1H,J_0H:J_1H)
c  
     .,pgfx(I_0H:I_1H,J_0H:J_1H),pgfy(I_0H:I_1H,J_0H:J_1H) 
     .,gradx(I_0H:I_1H,J_0H:J_1H),grady(I_0H:I_1H,J_0H:J_1H) 
     .,depthu(I_0H:I_1H,J_0H:J_1H),depthv(I_0H:I_1H,J_0H:J_1H) 
     .,pvtrop(I_0H:I_1H,J_0H:J_1H) 
     .,depths(I_0H:I_1H,J_0H:J_1H) 
     .,drag(I_0H:I_1H,J_0H:J_1H) 
     .,glue(I_0H:I_1H,J_0H:J_1H) 
     .,dampu(I_0H:I_1H,J_0H:J_1H),dampv(I_0H:I_1H,J_0H:J_1H) ) 
c
       allocate(  
     . uja(I_0H:I_1H,J_0H:J_1H),ujb(I_0H:I_1H,J_0H:J_1H) 
     .,via(I_0H:I_1H,J_0H:J_1H),vib(I_0H:I_1H,J_0H:J_1H) 
     .,pbot(I_0H:I_1H,J_0H:J_1H) 
     .,tracer(I_0H:I_1H,J_0H:J_1H,2*kdm,ntrcr) 
     .,diadff(I_0H:I_1H,J_0H:J_1H,kdm) 
     .,tprime(I_0H:I_1H,J_0H:J_1H) 
     .,sgain(I_0H:I_1H,kdm) 
     .,surflx(I_0H:I_1H,J_0H:J_1H) 
     .,salflx(I_0H:I_1H,J_0H:J_1H) 
     .,thkice(I_0H:I_1H,J_0H:J_1H) 
     .,covice(I_0H:I_1H,J_0H:J_1H) 
c    .,temice(I_0H:I_1H,J_0H:J_1H) 
c    .,odhsi(I_0H:I_1H,J_0H:J_1H) 
     .,odmsi(I_0H:I_1H,J_0H:J_1H) 
     .,omlhc(I_0H:I_1H,J_0H:J_1H) 
     .,dmfz(I_0H:I_1H,J_0H:J_1H)
     .,thmix(I_0H:I_1H,J_0H:J_1H)
     .,tmix(I_0H:I_1H,J_0H:J_1H)
     .,smix(I_0H:I_1H,J_0H:J_1H)
     .,umix(I_0H:I_1H,J_0H:J_1H)
     .,vmix(I_0H:I_1H,J_0H:J_1H) 
c 
     .,ylo    (I_0H:I_1H,J_0H:J_1H) 
     .,xlo    (I_0H:I_1H,J_0H:J_1H) 
     .,southfl(I_0H:I_1H,J_0H:J_1H) 
     .,eastfl (I_0H:I_1H,J_0H:J_1H) )

      allocate( klist(I_0H:I_1H,J_0H:J_1H) 
     .        ,ijlist(I_0H:I_1H,J_0H:J_1H) )
c  
      allocate(  
     . taux(I_0H:I_1H,J_0H:J_1H) 
     .,tauy(I_0H:I_1H,J_0H:J_1H) 
c    .,wndspd(I_0H:I_1H,J_0H:J_1H,4) 
c    .,airtmp(I_0H:I_1H,J_0H:J_1H,4) 
c    .,vapmix(I_0H:I_1H,J_0H:J_1H,4) 
c    .,oprec(I_0H:I_1H,J_0H:J_1H) 
c    .,oevap(I_0H:I_1H,J_0H:J_1H) 
     .,oemnp(I_0H:I_1H,J_0H:J_1H) 
     .,oflxa2o(I_0H:I_1H,J_0H:J_1H),oice(I_0H:I_1H,J_0H:J_1H) 
     .,ustar(I_0H:I_1H,J_0H:J_1H) 
     .,ustarb(I_0H:I_1H,J_0H:J_1H) 
     .,osalt(I_0H:I_1H,J_0H:J_1H)
c 
     .,freshw(I_0H:I_1H,J_0H:J_1H) 
     .,diafor(I_0H:I_1H,J_0H:J_1H) ) 
c 

      u = 0
      v = 0
      dp = 0
      dpold = 0
      dpu = 0
      dpv = 0
      p = 0
      pu = 0
      pv = 0
      latij = 0
      lonij = 0
      corio = 0
      potvor = 0
      temp = 0
      saln = 0
      th3d = 0
      thstar = 0
      wgtkap = 0
      psikk = 0
      thkk = 0
      dpmixl = 1.    ! TNL: avoid NaN for the first step
      srfhgt = 0
      montg = 0
      defor1 = 0
      defor2 = 0
      ubavg = 0
      vbavg = 0
      pbavg = 0
      ubrhs = 0
      vbrhs = 0
      utotm = 0
      vtotm = 0
      utotn = 0
      vtotn = 0
      uflux = 0
      vflux = 0
      uflux1 = 0
      vflux1 = 0
      uflux2 = 0
      vflux2 = 0
      uflux3 = 0
      vflux3 = 0
      uflx = 0
      vflx = 0
      bolusu = 0
      bolusv = 0
      uav = 0
      vav = 0
      dpuav = 0
      dpvav = 0
      temav = 0
      salav = 0
      th3av = 0
      dpav = 0
      ubavav = 0
      vbavav = 0
      pbavav = 0
      sfhtav = 0
      uflxav = 0
      vflxav = 0
      diaflx = 0
      salflav = 0
      brineav = 0
      eminpav = 0
      surflav = 0
      tauxav = 0
      tauyav = 0
      ufxcum = 0
      vfxcum = 0
      dpinit = 0
      dpmxav = 0
      oiceav = 0
      util1 = 0
      util2 = 0
      util3 = 0
      util4 = 0
      scpx = 0
      scpy = 0
      scux = 0
      scuy = 0
      scvx = 0
      scvy = 0
      scqx = 0
      scqy = 0
      scu2 = 0
      scv2 = 0
      scp2 = 0
      scq2 = 0
      scuxi = 0
      scvyi = 0
      scp2i = 0
      scq2i = 0
      pgfx = 0
      pgfy = 0
      gradx = 0
      grady = 0
      depthu = 0
      depthv = 0
      pvtrop = 0
      depths = 0
      drag = 0
      glue = 0
      dampu = 0
      dampv = 0
      uja = 0
      ujb = 0
      via = 0
      vib = 0
      pbot = 0
      tracer = 0
      diadff = 0
      tprime = 0
      sgain = 0
      surflx = 0
      salflx = 0
      odmsi = 0
      omlhc = 0
      dmfz = 0
      taux = 0
      tauy = 0
      oemnp = 0
      oflxa2o = 0
      oice = 0
      ustar = 0
      ustarb = 0
      osalt = 0
      freshw = 0
      diafor = 0
      ijlist = 0
      klist = 0

      end subroutine alloc_hycom_arrays

      end module hycom_arrays
c
c> Revision history:
c>
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c-----------------------------------------------------------------------------
