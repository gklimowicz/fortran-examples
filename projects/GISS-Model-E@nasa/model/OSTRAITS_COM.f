#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif

      module straits
      use oceanres,  only: lmo
      use seaice, only: lmi
      implicit none
      save

      INTEGER :: NMST=0  !@var NMST no. of ocean straits

!@var WIST,ZST width, nominal depth of strait (m)
!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      real*8, allocatable :: wist(:),zst(:),xst(:,:),yst(:,:)

!@var IST,JST i,j coordinates of ends of straits
!@var LMST no. of levels in strait
      integer, allocatable :: lmst(:),ist(:,:),jst(:,:)

!@var name_st Names of straits
      character(len=20), allocatable :: name_st(:)

!@var DIST distance along strait (m)
!@var DISTPG distance between centre points of adjoining ocean boxes (m)
!@var MMST mass of water in strait (kg)
!@var MUST mass flux of water in strait (kg/s)
!@var G0MST,GXMST,GZMST pot. enthalpy of water in strait (+ moments) (J)
!@var S0MST,SXMST,SZMST salinity of water in strait (+ moments) (kg)
!@var RSIST Sea ice fraction in strait
!@var RSIXST Center of sea ice in strait (m)
!@var MSIST Mass of ice within strait (kg)
!@var HSIST Enthalpy of ice within strait (J)
!@var SSIST Salinity of ice within strait (kg)
!@param USIFAC ratio of strait sea ice velocity to current
!@var TRMST,TXMST,TZMST tracer amount in strait (+ moments) (kg)
!@var TRSIST tracer amount in with strait (kg)
#ifdef OCN_GISS_TURB
!@var OTKEST turbulent kinetic energy in strait (m/s)^2
#endif

      Real*8 ::
     *    USIFAC = .1d0  !  ratio of strait sea ice velocity to current

      real*8, dimension(:), allocatable ::
     &       DIST, !  fixed length of strait (m)
     &     DISTPG  !  fixed distance between ocean cell centers
!     &      RSIST, !  horizontal sea ice cover in strait
!     &     RSIXST  !  center of sea ice cover in strait

      real*8, dimension(:,:), allocatable ::
     &   MMST, !  mass of water in strait (kg)
     &   MUST, !  mass flux of water in strait (kg/s)
     &  G0MST, !  pot. enthalpy of water in strait (J)
     &  GXMST, !  west-east gradient of pot. enthalpy (J)
     &  GZMST, !  down. vert. gradient of pot. enthalpy (J)
     &  S0MST, !  mass of salt in strait (kg)
     &  SXMST, !  west-east gradient of salt (kg)
     &  SZMST  !  down. vert. gradient of salt (kg)
!     &  MSIST, !  2 mass layers of sea ice in strait (kg)
!     &  HSIST, !  LMI layers of heat content in strait (J)
!     &  SSIST  !  LMI layers of salt in strait (kg)

#ifdef OCN_GISS_TURB
      real*8, dimension(:,:), allocatable ::
     &  OTKEST!  turbulent kinetic energy in strait (m/s)^2
#endif

!@var QTYE workspace holding values of QTY at endpoints of straits
      integer, dimension(:,:), allocatable :: LMMe
      real*8, dimension(:,:), allocatable :: OPRESE,HOCEANE
      real*8, dimension(:,:,:), allocatable ::
     &     MOE, G0ME,GXME,GYME,GZME, S0ME,SXME,SYME,SZME

!@var kn2: Where kn2>0, the index pair k,n [k=1,2;n=1,nmst] corresponds
!@+        to the same i,j endpoint as index pair kn2(1:2,k,n).
!@+        Currently, an i,j endpoint can be shared by only 2 straits.
      integer, dimension(:,:,:), allocatable :: kn2

#ifdef TRACERS_OCEAN
!@var TRMST,TXMST,TZMST tracer amount in strait (+ moments) (kg)
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::
     &     TRMST,TXMST,TZMST !(LMO,NMST,NTM)
!@var TRME,TXME,TYME,TZME tracers at the endpoints of straits
      REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE ::
     &     TRME,TXME,TYME,TZME !(2,NMST,LMO,NTM)
#endif
!#ifdef TRACERS_WATER
!      Real*8, ALLOCATABLE :: TRSIST(:,:,:) !(NTM_ATM,LMI,NMST)
!#endif


      end module straits

      subroutine alloc_straits
      use straits
      use filemanager, only : openunit,closeunit,file_exists
      implicit none
      integer :: iu,ios,n
      character(len=20) :: name
      integer :: ij1(2),ij2(2),lm,iter,ier
      real*8 :: width,depth,xy1(2),xy2(2)
      namelist/strait/ name,ij1,ij2,lm,width,depth,xy1,xy2

      if(file_exists('OSTRAITS')) then
        call openunit('OSTRAITS',iu,.false.,.true.)
      else
        iu = -999
      endif

      do iter=1,2 ! first iteration counts straits, second fills data
        nmst = 0
        if(iu.ne.-999) then
          do
            n = nmst+1
            depth = -1d30
            read(iu,nml=strait,iostat=ios)
            if(ios.ne.0) exit
            nmst = nmst + 1
            if(iter.eq.2) then
              name_st(n) = name
              wist(n) = width
              zst(n) = depth
              lmst(n) = lm
              ist(n,1) = ij1(1)
              jst(n,1) = ij1(2)
              ist(n,2) = ij2(1)
              jst(n,2) = ij2(2)
              xst(n,1) = xy1(1)
              yst(n,1) = xy1(2)
              xst(n,2) = xy2(1)
              yst(n,2) = xy2(2)
            endif
          enddo
          rewind(iu)            ! use rewind_parallel?
        endif
        if(iter.eq.1) then      ! allocate once we know how many straits exist
          allocate(wist(nmst),xst(nmst,2),yst(nmst,2))
          allocate(lmst(nmst),ist(nmst,2),jst(nmst,2))
          allocate(name_st(nmst),zst(nmst))
        endif
      enddo
      if(iu.ne.-999) then
        call closeunit(iu)
      endif

      allocate(
     &     MMST(LMO,NMST),
     &     MUST(LMO,NMST),
     &     G0MST(LMO,NMST),
     &     GXMST(LMO,NMST),
     &     GZMST(LMO,NMST),
     &     S0MST(LMO,NMST),
     &     SXMST(LMO,NMST),
     &     SZMST(LMO,NMST),
     &     DIST(NMST),
     &     DISTPG(NMST),
!     &     RSIST(NMST),
!     &     RSIXST(NMST),
!     &     MSIST(2,NMST),
!     &     HSIST(LMI,NMST),
!     &     SSIST(LMI,NMST),
#ifdef OCN_GISS_TURB
     &     OTKEST(LMO,NMST),
#endif
     &     OPRESE(2,nmst),
     &     HOCEANE(2,nmst),
     &     LMMe(2,nmst),      
     &     MOE(2,NMST,LMO),
     &     G0ME(2,NMST,LMO),
     &     GXME(2,NMST,LMO),
     &     GYME(2,NMST,LMO),
     &     GZME(2,NMST,LMO),
     &     S0ME(2,NMST,LMO),
     &     SXME(2,NMST,LMO),
     &     SYME(2,NMST,LMO),
     &     SZME(2,NMST,LMO),
     &     kn2(2,2,nmst),
     &     stat=ier)


      end subroutine alloc_straits
