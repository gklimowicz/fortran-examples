#include "rundeck_opts.h"
      module HYCOM_SCALARS
      USE HYCOM_DIM_GLOB

      implicit none

      private

c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing functions
c --- windf       include wind stress in forcing functions
c --- relax       activate lateral boundary nudging
c --- trcout      advect tracer and save results in history/restart file
c --- dotrcr      perform column physics operations on tracer array(s)
c
      logical, public:: diagno,thermo,windf,relax,trcout,dotrcr
c
      real, public :: time,time0,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3,
     . area,avgbot,ocnvol,watcum=0.,empcum=0.,slfcum=0.,brncum=0.,
     . sala2o,tavini,tmean0,smean0,tmean1,smean1,
     . zonarea(3)=0.,zonwat(3)=0.,zonemp(3)=0.,zonsfl(3)=0.,
     . zonbrn(3)=0.,zonsqi(3)=0.,zonrun(3)=0.,zondm(3)=0.
c
      integer, public ::  nstep,nstep0,nstepi,lstep,l0,l1,l2,l3,ls0,ls1
     .             ,ls2,ls3,oddev
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --  'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'trcfrq' = number of time steps between tracer transport calculations
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- slip = +1  for free-slip boundary cond., slip = -1  for non-slip cond.
c --- 'diagfq' = number of days between model diagnostics (incl.output)
c --- 'ntracr' = number of time steps between tracer transport
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'thkmin' = minimum mixed-layer thickness
c --- 'thkbot' = thickness of bottom boundary layer
c --- 'botmin' = minimum topo depth
c --- 'ekman'  = thickness of ekman layer
c --- 'sigjmp' = minimum density jump at mixed-layer bottom
c ---' salmin' = minimum salinity allowed in an isopycnic layer
c --- 'acurcy' = permissible roundoff error in column integral calc.
c --- 'nhr   ' = coupling freq. in hours
c
      dimension theta(kdm),salmin(kdm),dplist(kdm)
      real, dimension(kdm) :: pr1d
      real, public ::
     &     theta,thbase,baclin,batrop,veldff,temdff,viscos,
     &     vertmx,h1,slip,cbar,diagfq,wuv1,wuv2,wts1,wts2,dplist,
     &     acurcy, wbaro,thkmin,thkbot,botmin,ekman,sigjmp,salmin
      public pr1d, init_pr1d
c
      integer, public ::       trcfrq,ntracr,nhr,mixfrq
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'g'      = gravity acceleration
c --- 'cd'     = drag coefficient
c --- 'ct'     = thermal transfer coefficient
c --- 'epsil'  = small nonzero number used to prevent division by zero
c
      real, public ::
     &     tenm,onem,tencm,onecm,onemm,onemu,g,csubp,spcifh,cd,ct,
     &     airdns,evaplh,thref,epsil,huge,radian,pi
c
      character*60, public ::
     &     flnmdep,flnmrsi,flnmrso,flnmarc,flnmfor,flnmovt
     &            ,flnmini,flnmriv,flnmbas,flnmdia,flnmlat
     &            ,flnminp,flnmint,flnmins
     &            ,flnmcoso,flnmcosa,flnma2o,flnmo2a,flnmcellsz

c --- opening the bering strait requires information exchange across a
c --- 'u' face represented in 2 different locations in the tri-pole grid.
c --- the 2 locations are the northern tip of the bering 'inlet' (truncated
c --- bering channel) on the pacific side and the southern (i.e., upper)
c --- tip of the bering inlet in the panam pole patch
c
c --- ipacn,jpac:  grid point north of inlet head on pacific side
c --- ipacs,jpac:  grid point south of inlet head on pacific side
c --- iatln,jatl:  grid point north of inlet head on arctic ocean side
c --- iatls,jatl:  grid point south of inlet head on arctic ocean side
c
c --- thus, the pairs [(ipacn,jpac),(iatln,jatl)],[(ipacs,jpac),(iatls,jatl)]
c --- refer to identical grid cells in physical space.
c
      logical, public, parameter :: beropn=.true.   !true if bering strait open
     .                             ,kappa =.false.  !true to include thermobaricity
#ifdef HYCOM2deg
      integer, public, parameter :: ipacn=67,ipacs=68,jpac= 95
      integer, public, parameter :: iatln= 2,iatls= 1,jatl=156
#endif
#ifdef HYCOM1degRefined
      integer, public, parameter :: ipacn=137,ipacs=138,jpac=189
      integer, public, parameter :: iatln= 2,iatls= 1,jatl=312
#endif
#ifdef HYCOM1degUnrefined
      integer, public, parameter :: ipacn=137,ipacs=138,jpac=189
      integer, public, parameter :: iatln= 2,iatls= 1,jatl=312
#endif
c-----------------------------------------------------------------------------
c
c --- layer densities (theta units):
c
#ifdef HYCOM26layers
      data theta/
     . 24.35,26.07,27.50,28.67,29.61,30.35,30.92,31.35,31.67,31.90
     .,32.06,32.17,32.25,32.31,32.36,32.40,32.43,32.46,32.49,32.52
     .,32.54,32.56,32.58,32.60,32.62,32.64/     ! 26-sig1
      data dplist/
     .     5.,  7.,  9., 11., 13., 15., 17., 19., 22., 26.,
     .    31., 37., 45., 55., 67., 81., 98.,118.,141.,168.,
     .   199.,234.,274.,319.,369.,425./         ! total 2805m
#endif

#ifdef HYCOM32layers
      data theta/
     . 20.82,21.64,22.46,23.28,24.10,24.92,25.74,26.55,27.34,28.10,
     . 28.82,29.49,30.10,30.64,31.10,31.48,31.78,32.01,32.18,32.30,
     . 32.38,32.43,32.46,32.48,32.50,32.52,32.54,32.56,32.58,32.60,
     . 32.62,32.64 /                            ! 32-sig1a
      data dplist/
     .   2.0,   3.0,  4.0,  5.0,  6.0,  8.0, 10.0, 12.0, 15.0, 18.0,
     .  21.0,  25.0, 29.0, 33.0, 38.0, 43.0, 48.0, 54.0, 60.0, 66.0,
     .  73.0,  81.0, 90.0,101.0,114.0,129.0,146.0,166.0,189.0,215.0,
     . 244.0, 276. /
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'diagfq' = number of days between model diagnostics (incl.output)
c --- 'equatn' = the i index of the equator
      real, public :: equatn
#ifdef HYCOM2deg
      data baclin,batrop/3600.,100./,diagfq/365./          ! 2deg full global
      data equatn/122./
#endif
#ifdef HYCOM1degRefined
      data baclin,batrop/1800., 50./,diagfq/365./          ! 1deg full global
      data equatn/243./
#endif
#ifdef HYCOM1degUnrefined
      data baclin,batrop/1800., 50./,diagfq/365./          ! 1deg full global
      data equatn/229./
#endif
c
c --- 'thkdff' = diffusion velocity (m/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (m/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (m/s) for temp/salin. mixing
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'vertmx' = scale velocity for vertical momentum mixing (m/s)
      data veldff/.1/,temdff/.02/,viscos/0.3/,vertmx/0./
c
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- 'thkmin' = minimum mixed-layer thickness (m)
c --- 'acurcy' = permissible roundoff error in column integral calc.
      data h1/98060./,thkmin/5./,acurcy/1.e-11/,botmin/30./
c
c --- slip=+1  for free-slip boundary cond., slip=-1  for non-slip cond.
c --- 'cbar'   = rms flow speed (m/s) for linear bottom friction law
c --- 'thkbot' = thickness of bottom boundary layer (m)
c --- 'ekman'  = thickness of ekman layer (m)
c --- 'sigjmp' = minimum density jump at mixed-layer bottom (theta units)
      data slip/-1./,cbar/0.1/,thkbot/10./,ekman/30./,sigjmp/.01/
c
c --- weights for time smoothing
ccc   data wuv1,wuv2/.5,.25/
      data wuv1,wuv2/.75,.125/
ccc   data wts1,wts2/.875,.0625/
ccc   data wts1,wts2/.9375,.03125/
ccc   data wts1,wts2/.96875,.015625/
      data wts1,wts2/.984375,.0078125/
      data wbaro/.125/
c
c --- layer thicknesses in units of pressure (kg/m/sec^2):
      data tenm, onem, tencm, onecm, onemm, onemu
     .   /98060.,9806.,980.6, 98.06, 9.806,.0098/
      data radian/57.2957795/,pi/3.1415926536/
c
c --- 'g'      = gravitational acceleration (m/s^2)
c --- 'csubp'  = specific heat of air at constant pressure (j/kg/deg)
c --- 'spcifh' = specific heat of sea water (j/kg/deg)
c --- 'airdns' = air density at sea level (kg/m^3)
c --- 'evaplh' = latent heat of evaporation (j/kg)
c --- 'thref'  = reference value of specific volume (m^3/kg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
      data g/9.806/,csubp/1005.7/,spcifh/4185./
      data airdns/1.2/,evaplh/2.47e6/,thref/1.e-3/,epsil/1.e-11/
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- thermo      use thermodynamic forcing functions
c --- windf       use wind stress forcing function
c --- relax       activate lateral boundary nudging
c
#if (defined TRACERS_HYCOM_Ventilation) \
 || (defined TRACERS_OceanBiology) \
 || (defined TRACERS_AGE_OCEAN) \
 || (defined TRACERS_OCEAN_WATER_MASSES) \
 || (defined TRACERS_ZEBRA)
      data thermo/.true./, windf/.true./,relax/.false./,trcout/.true./
#else
      data thermo/.true./, windf/.true./,relax/.false./,trcout/.false./
#endif
c
c
c --- use 'huge' to initialize array portions that the code should never access
      data huge/1.e33/
!!! temporary value for debug (for checksums to make sense)
!!!      data huge/0.e0/
      data nhr/1/                        ! couple every nhr hours
      data oddev/-1/
c
c --- i/o file names
c
c     flnmdep = name/location of basin depth array
c     flnmint = name/location of initial -t- field
c     flnmins = name/location of initial -s- field
c     flnminp = name/location of initial -p- field
c     flnmfor = name/location of the forcing functions
c     flnmdia = name/location of the diapycnal flux forcing fields
c     flnmriv = name/location of freshwater (river runoff) fields
c     flnmbas = name/location of basin mask file used in overtn diagno
c     flnmrsi = location (pathname) of restart file (input)
c     flnmrso = location (pathname) of restart file (output)
c     flnmarc = location (pathname) of archive files
c     flnmovt = location (pathname) of ovtn.xxxxxx files
c     flnmlat = location (pathname) of lat/lon at vorticity points
c     flnmscp2= name/location of scp2
c
      data flnmlat    /'latlonij'/
      data flnmdep    /'hycomtopo'/
      data flnmint    /'temp_ini'/
      data flnmins    /'salt_ini'/
      data flnminp    /'pout_ini'/
      data flnmbas    /'ibasin'/
      data flnma2o    /'wgt_a2o'/
      data flnmo2a    /'wgt_o2a'/
      data flnmcoso   /'cososino'/
      data flnmovt    /'./'/
      data flnmcellsz /'hycom_cellsz'/

c --- grid point where detailed diagnostics are desired:
      integer, public :: itest=-1, jtest=-1    !overwritten by values in rundeck
c
c --- ocean mixed layer schemes
      integer, public :: iocnmx=2              !overwritten by value in rundeck
c
c --- initial jerlov water type (1 to 5; 0 to use KPAR)
      integer, public :: jerlv0=1              !overwritten by values in rundeck
c
c --- brntop/brnbot:top/bottom of depth interval over which to distribute brine
      real, public :: brntop=50., brnbot=200.  !overwritten by values in rundeck
c
c --- ocnmx_factor_s/ocnmx_factor_t:factor to reduce difs/dift in mxkprf.f
      real, public :: ocnmx_factor_s=.1, ocnmx_factor_t=.1
c
c --- 'diapyn' = diapycnal diffusivity times buoyancy freq. (m^2/s^2)
c --- 'diapyc' = diapycnal diffusivity (m^2/s)
      real, public :: diapyn=2.e-7, diapyc=.2e-4 !overwritten by values in rundeck
c
c --- 'thkdff' = diffusion velocity (m/s) for thickness diffusion
      real, public :: thkdff = 0.05             !overwritten by values in rundeck
c
      real, public :: stdsal=34.7, h_glb_cum(2)=0., s_glb_cum(2)=0.
c --- choices for bolus velocity (interface smoothing); overwritten by values in rundeck
c --- 1 = true, 0 = false
      integer, public :: bolus_biharm_constant=0
      integer, public :: bolus_laplc_constant =1
      integer, public :: bolus_laplc_exponential=0

      contains

      subroutine init_pr1d
      integer :: k
#ifdef HYCOM26layers
! 200 m step from 0 to 5000m:
      pr1d(:)=(/(5000.*float(k-1)/float(kdm-1)* onem,k=1,kdm)/) ! isobaric depth levels
#endif
#ifdef HYCOM32layers
! 160 m step from 0 to 4960 meters (31 layers) plus last layer from 4960 until bottom
      pr1d(:)=(/(4960.*float(k-1)/float(kdm-1)* onem,k=1,kdm)/) ! isobaric depth levels
#endif
      end subroutine init_pr1d

      end module HYCOM_SCALARS
