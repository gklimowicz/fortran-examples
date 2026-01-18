      subroutine dpthuv
c
c --- define water depth (bottom pressure) at  u,v  points and barotp.pot.vort.
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
      USE HYCOM_DIM_GLOB, only : jj,ifu,isu,ilu,isv,ifv,ilv,isq,ifq,ilq
     &     ,jchunk
      USE HYCOM_SCALARS, only : onem
      USE HYCOM_ARRAYS_GLOB
      implicit none
c
      integer i,j,l,ja,jb
c
      real uvdep,a,b,damp
      data damp/5.e-6/		!  inverse time scale for coastal wave damping
c
c --- function for determining depth at u,v points
      uvdep(a,b)=min(a,b)
c
      call cpy_p(pbot)
c
      do j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
c
      do 151 l=1,isu(j)
      do 151 i=ifu(j,l),ilu(j,l)
      depthu(i,j)=uvdep(pbot(i,j),pbot(i-1,j))
      dampu(i,j)=damp*exp(-depthu(i,j)/(50.*onem))
      pvtrop(i,j  )=corio(i,j  )*2./(pbot(i,j)+pbot(i-1,j))
 151  pvtrop(i,jb )=corio(i,jb )*2./(pbot(i,j)+pbot(i-1,j))
c
      do 152 l=1,isv(j)
      do 152 i=ifv(j,l),ilv(j,l)
      depthv(i,j)=uvdep(pbot(i,j),pbot(i,ja ))
      dampv(i,j)=damp*exp(-depthv(i,j)/(50.*onem))
      pvtrop(i  ,j)=corio(i  ,j)*2./(pbot(i,j)+pbot(i,ja ))
 152  pvtrop(i+1,j)=corio(i+1,j)*2./(pbot(i,j)+pbot(i,ja ))
      end do
c
      do j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 153 l=1,isq(j)
      do 153 i=ifq(j,l),ilq(j,l)
 153  pvtrop(i,j)=corio(i,j)*4./(pbot(i,j  )+pbot(i-1,j  )
     .                          +pbot(i,ja )+pbot(i-1,ja ))
      end do
c
      return
      end
c
c
c> Revision history:
c>
c> May  2000 - changed i/j loop nesting to j/i
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Sep. 2000 - made sure loop 153 is executed after loops 151/152 are finished
c> Oct. 2000 - added calc. of damping factors simulating coastal wave breaking 
c> Mar. 2006 - added bering strait exchange logic
