#include "hycom_mpi_hacks.h"
      subroutine dpudpv(mmnn)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
      USE HYCOM_DIM_GLOB, only : kk,jj,isu,ifu,ilu,isv,ifv,ilv
     &     ,jchunk
      USE HYCOM_ARRAYS_GLOB
      implicit none
c
c --- ----------------------------------
c --- define layer depth at  u,v  points
c --- ----------------------------------
c
      integer :: i,k,j,l, ja
c
      integer mmnn,kmn
c
      do k=1,kk
        call cpy_p(p(1,1,k+1))
      end do
c
      do j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 155 k=1,kk
      kmn=k+mmnn
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
 154  dpu(i,j,kmn)=max(0.,
     .           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     .           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
 155  dpv(i,j,kmn)=max(0.,
     .           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,ja ,k+1)))-
     .           min(depthv(i,j),.5*(p(i,j,k  )+p(i,ja ,k  ))))
c
      end do
      return
      end
c
c
c> Revision history:
c>
c> Sep. 2000 - (dpu,dpv,depthu,depthv,p) no longer passed as arguments
c> Mar. 2006 - added bering strait exchange logic
c--------------------------------------------------------------------------

      subroutine pardpudpv(mmnn)
c
c --- parallel version
c --- version 2.8 -- cyclic and noncyclic b.c. combined
      USE HYCOM_DIM, only : kk,jj,isu,ifu,ilu,isv,ifv,ilv
     &     ,jchunk,ogrid,I_0H,I_1H,J_0H,J_1H,J_0,J_1
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, ONLY: HALO_UPDATE, SOUTH
      implicit none
c
c --- ----------------------------------
c --- define layer depth at  u,v  points
c --- ----------------------------------
c
      integer :: i,k,j,l, ja
c
      integer mmnn,kmn
c
      do k=1,kk
        call cpy_p_par(p(I_0H,J_0H,k+1))
      end do

      CALL HALO_UPDATE(ogrid,p, FROM=SOUTH)
c
      do j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 155 k=1,kk
      kmn=k+mmnn
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
 154  dpu(i,j,kmn)=max(0.,
     .           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     .           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
 155  dpv(i,j,kmn)=max(0.,
     .           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,ja ,k+1)))-
     .           min(depthv(i,j),.5*(p(i,j,k  )+p(i,ja ,k  ))))
c
      end do
      return
      end
c
c
c> Revision history:
c>
c> Sep. 2000 - (dpu,dpv,depthu,depthv,p) no longer passed as arguments
c> Mar. 2006 - added bering strait exchange logic
