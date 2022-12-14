#include "rundeck_opts.h" 
module TRACERS_SEASALT

!===============================================================================
implicit none
!===============================================================================
!@dbparam tune_ss1, tune_ss2 factors to tune seasalt sources
real*8 :: tune_ss1=1.d0, tune_ss2=1.d0
#ifdef TRACERS_AEROSOLS_OCEAN
!@var OC_SS_enrich_fact OCocean enrichment factor of seasalt1
      real*8, ALLOCATABLE, DIMENSION(:,:) :: OC_SS_enrich_fact
#endif  /* TRACERS_AEROSOLS_OCEAN */
!===============================================================================

contains

!===============================================================================
subroutine alloc_seasalt_sources(grid)

use domain_decomp_atm, only : dist_grid, getDomainBounds
implicit none

type (dist_grid), intent(in) :: grid
integer :: J_1H, J_0H, I_1H, I_0H

call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
I_0H = grid%I_STRT_HALO
I_1H = grid%I_STOP_HALO

#ifdef TRACERS_AEROSOLS_OCEAN
allocate(OC_SS_enrich_fact(I_0H:I_1H,J_0H:J_1H))
#endif  /* TRACERS_AEROSOLS_OCEAN */

end subroutine alloc_seasalt_sources
!===============================================================================
subroutine read_seasalt_sources(swind,itype,ibin,i,j,ss,tr)
!@sum determines wind-speed dependent oceanic seasalt source
!@auth Dorothy Koch
! want kg seasalt/m2/s, for now in 2 size bins
use TimeConstants_mod, only: SECONDS_PER_DAY
USE GEOM, only: axyp
use model_com, only: modelEclock
use lakes_com, only: flake
#ifdef TRACERS_TOMAS
USE TOMAS_EMIS, only : scalesizeSalt
#endif
!USE FLUXES, only: gtemp !Jaegle
use Dictionary_mod, only: sync_param

implicit none

REAL*8 erate,swind_cap
integer jread,nb
integer, INTENT(IN)::itype,ibin,i,j
REAL*8, INTENT(IN)::swind
REAL*8, INTENT(OUT)::ss
character*8, intent(in) :: tr

ss=0.
erate=0.d0
if (flake(i,j) /= 0.d0) return ! if there are lakes, there is no ocean

  if (itype.eq.1) then
! Monahan 1971, bubble source, important for small (<10um) particles
!  swind_cap=swind !modelE
!  if (swind.gt.10.d0) swind_cap=10.d0 !modelE
!  erate= 1.373d0 * swind_cap**(3.41d0) !modelE
!  erate=1.373d0*swind**3.41d0 !Monahan
!  erate=(swind/10.d0)**2.5d0 !Lewis and Schwartz
!  erate=(0.3d0+0.1d0*gtemp(1,1,i,j)-0.0076d0*gtemp(1,1,i,j)**2+&
!         0.00021d0*gtemp(1,1,i,j)**3)*1.373d0*swind**3.41d0 !Jaegle
  erate=1.373d0*swind**3.41d0 !Gong
! units are kg salt/m2/s
#ifdef TRACERS_TOMAS 
  ss=tune_ss1*erate*scalesizeSalt(ibin)
#else
  if (ibin.eq.1) then ! submicron (0.1 < r_d < 1.)
!    ss=tune_ss1*erate*2.11d-14 !modelE
!    ss=tune_ss1*erate*1.750298d-14 !Monahan
!    ss=tune_ss1*erate*3.431418d-11 !Lewis and Schwartz
!    ss=tune_ss1*erate*1.336825d-14 !Jaegle
    ss=tune_ss1*erate*1.336825d-14 !Gong
#ifdef TRACERS_AEROSOLS_OCEAN
    if (trim(tr).eq.'OCocean') then
      ss=ss*OC_SS_enrich_fact(i,j)
    else
      ss=ss*(1.d0-OC_SS_enrich_fact(i,j))
    endif
#endif  /* TRACERS_AEROSOLS_OCEAN */
  else ! supermicron (1. < r_d < 4.)
!    ss=tune_ss2*erate*7.78d-14 !modelE
!    ss=tune_ss2*erate*7.097854d-14 !Monahan
!    ss=tune_ss2*erate*1.650674d-9 !Lewis and Schwartz
!    ss=tune_ss2*erate*8.676763d-14 !Jaegle
    ss=tune_ss2*erate*8.676763d-14 !Gong
  endif
#endif  /* TRACERS_TOMAS */
endif

end subroutine read_seasalt_sources
!===============================================================================
#ifdef TRACERS_AEROSOLS_OCEAN
subroutine read_seawifs_chla(imon)
!@sum read_seawifs_chla Reads in SeaWiFS chlorophyll-a concentration
!@auth Kostas Tsigaridis

use filemanager, only: openunit,closeunit,nameunit
use domain_decomp_atm, only: grid,get,readt_parallel
implicit none

real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                  grid%j_strt_halo:grid%j_stop_halo) :: seawifs_chla
integer :: iu, imon
integer :: j,i,j_0,j_1,i_0,i_1

call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
call getDomainBounds(grid, i_strt=i_0, i_stop=i_1)

! read chla SeaWiFS data
call openunit('SeaWiFS_chla',iu,.true.,.true.)
call readt_parallel(grid,iu,nameunit(iu),seawifs_chla,imon)
call closeunit(iu)

! SS enrichment factor of OC (Vignati et al., AE, 2009, 670)
! if chla is more than 1.43 mg m-3 the enrichment factor remains
! constant, 76%.
do j=j_0,j_1; do i=i_0,i_1
  OC_SS_enrich_fact(i,j)=min(seawifs_chla(i,j),1.43d0)*0.435d0+0.13805d0
enddo ; enddo

end subroutine read_seawifs_chla
#endif  /* TRACERS_AEROSOLS_OCEAN */
!===============================================================================

end module TRACERS_SEASALT
