!#include "hycom_mpi_hacks.h"
      subroutine advfct(iord,fld,u,v,scal,scali,dt,fco,fc)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
c --- this is advem.f with MPDATA replaced by 2nd/4th order FCT
c --- if iord=1, scheme reduces to simple donor cell scheme
c
c  fld    - transported mixing ratio, e.g., salinity or temperature
c  u,v    - mass fluxes satisfying continuity equation
c  scal   - spatial increments (squared)
c  scali  - inverse of scal
c  dt     - temporal increment
c  fco,fc - depth of the layer at previous and new time step
c
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : itest,jtest,onemu

      implicit none
      integer i,j,l,ia,ib,ja,jb
c
      real fld(idm,jdm),u(idm,jdm),v(idm,jdm),scal(idm,jdm),
     .     scali(idm,jdm),fco(idm,jdm),fc(idm,jdm)
      real fmx(idm,jdm),fmn(idm,jdm),flp(idm,jdm),fln(idm,jdm),
     .     flx(idm,jdm),fly(idm,jdm),uan(idm,jdm),van(idm,jdm),
     .     flxdiv(idm,jdm),clipj(jdm),vlumj(jdm)
      real dt,q,clip,vlume,amount,bfore,after,epsil
      integer iord,jaa
      logical wrap,recovr
      data recovr/.true./
c
      parameter (epsil=1.e-11)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional code for checking conservation properties
ccc      bfore=0.
ccc      do 14 j=1,jj
ccc      do 14 l=1,isp(j)
ccc      do 14 i=ifp(j,l),ilp(j,l)
ccc 14   bfore=bfore+fld(i,j)*fco(i,j)*scal(i,j)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- compute low-order and antidiffusive (high-minus-low order) fluxes
c
      call cpy_p(fld)
c
      do 3 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      jaa=mod(j-3+jj,jj)+1
c
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      if (u(i,j).ge.0.) then
        q=fld(i-1,j)
      else
        q=fld(i  ,j)
      end if
      flx(i,j)=u(i,j)*q
      q=fld(i,j)+fld(i-1,j)				!  2nd order
      if (ip(i+1,j)+iu(i-1,j).eq.2)
     .  q=1.125*q-.125*(fld(i+1,j)+fld(i-2,j))		!  4th order
 2    uan(i,j)=.5*q*u(i,j)-flx(i,j)
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
      if (v(i,j).ge.0.) then
        q=fld(i,ja )
      else
        q=fld(i,j  )
      end if
      fly(i,j)=v(i,j)*q
      q=fld(i,ja )+fld(i,j)				!  2nd order
      if (ip(i,jb )+iv(i,ja).eq.2)
     .  q=1.125*q-.125*(fld(i,jb )+fld(i,jaa))		!  4th order
 3    van(i,j)=.5*q*v(i,j)-fly(i,j)
c
      do 22 j=1,jj
      do 22 l=1,isp(j)
      flx(ifp(j,l)  ,j)=0.
      flx(ilp(j,l)+1,j)=0.
      uan(ifp(j,l)  ,j)=0.
      uan(ilp(j,l)+1,j)=0.
 22   continue
c
      do 33 i=1,ii1
      wrap=jfv(i,1).eq.1	! true if j=1 and j=jj are both water points
      do 33 l=1,jsp(i)
      j=jfp(i,l)
      if (j.gt.1 .or. .not.wrap) then
        fly(i,j)=0.
        van(i,j)=0.
      end if
      j=mod(jlp(i,l),jj)+1
      if (j.gt.1 .or. .not.wrap) then
        fly(i,j)=0.
        van(i,j)=0.
      end if
 33   continue
c
      do 11 j=1,jj
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      ja=mod(j-2+jj,jj)+1
      if (ip(i,ja).eq.0) ja=j
      jb=mod(j     ,jj)+1
      if (ip(i,jb).eq.0) jb=j
      fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb))
 11   fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb))
c
cdiag i=itest
cdiag j=jtest
cdiag ja=mod(j+2+jj,jj)+1
cdiag jb=mod(j     ,jj)+1
cdiag write (*,101) 'advem (1)',i,j,fld(i-1,j),flx(i,j),fld(i,ja )
cdiag.,fly(i,j),fld(i,j),fly(i,jb ),fld(i,jb ),flx(i+1,j),fld(i+1,j)
  101 format(a,2i5,f20.3/es39.2/f21.3,es9.2,f9.3,es9.2,f9.3/
     .  es39.2/f39.3)
c
      do 61 j=1,jj
      jb=mod(j     ,jj)+1
      vlumj(j)=0.
      clipj(j)=0.
      do 61 l=1,isp(j)
      do 61 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*dt
     .   *scali(i,j)
      q=fld(i,j)*fco(i,j)-flxdiv(i,j)
      amount=max(fmn(i,j)*fc(i,j),min(q,fmx(i,j)*fc(i,j)))
      if (recovr) then
        vlumj(j)=vlumj(j)+scal(i,j)*fc(i,j)
        clipj(j)=clipj(j)+(q-amount)*scal(i,j)
      end if
 61   fld(i,j)=(fld(i,j)*onemu+amount)/(onemu+fc(i,j))
c
      if (iord.le.1) go to 100
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes
c
      do 12 j=1,jj
      jb=mod(j     ,jj)+1
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      flp(i,j)=(fmx(i,j)-fld(i,j))*fc(i,j)
     ./((max(0.,uan(i,j))-min(0.,uan(i+1,j))
     .  +max(0.,van(i,j))-min(0.,van(i,jb ))+epsil)*dt*scali(i,j))
c
 12   fln(i,j)=(fmn(i,j)-fld(i,j))*fc(i,j)
     ./((min(0.,uan(i,j))-max(0.,uan(i+1,j))
     .  +min(0.,van(i,j))-max(0.,van(i,jb ))-epsil)*dt*scali(i,j))
c
c---- limit antidiffusive fluxes
c
      call cpy_p(flp)
      call cpy_p(fln)
c
      do 8 j=1,jj
      ja=mod(j-2+jj,jj)+1
c
      do 7 l=1,isu(j)
      do 7 i=ifu(j,l),ilu(j,l)
      if (uan(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i-1,j))
      else
        clip=min(1.,fln(i,j),flp(i-1,j))
      end if
 7    flx(i,j)=uan(i,j)*clip
c
      do 8 l=1,isv(j)
      do 8 i=ifv(j,l),ilv(j,l)
      if (van(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i,ja ))
      else
        clip=min(1.,fln(i,j),flp(i,ja ))
      end if
 8    fly(i,j)=van(i,j)*clip
c
cdiag i=itest
cdiag j=jtest
cdiag ja=mod(j+2+jj,jj)+1
cdiag jb=mod(j     ,jj)+1
cdiag write (*,101) 'advem (2)',i,j,fld(i-1,j),flx(i,j),fld(i,ja )
cdiag.,fly(i,j),fld(i,j),fly(i,jb ),fld(i,jb ),flx(i+1,j),fld(i+1,j)
c
      do 62 j=1,jj
      jb=mod(j     ,jj)+1
      do 62 l=1,isp(j)
      do 62 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*dt
     .   *scali(i,j)
      q=fld(i,j)*fc(i,j)-flxdiv(i,j)
      amount=max(fmn(i,j)*fc(i,j),min(q,fmx(i,j)*fc(i,j)))
      if (recovr) clipj(j)=clipj(j)+(q-amount)*scal(i,j)
 62   fld(i,j)=(fld(i,j)*onemu+amount)/(onemu+fc(i,j))
c
  100 continue
c
c --- revover 'clipped' amount and return to field
c
      if (recovr) then
        vlume=0.
        clip=0.
        do 19 j=1,jj
        vlume=vlume+vlumj(j)
 19     clip=clip+clipj(j)
c
        if (vlume.ne.0.) then
          clip=clip/vlume
cdiag     write (*,'(a,es11.3)') 'tracer drift in advem:',-clip
          do 13 j=1,jj
          do 13 l=1,isp(j)
          do 13 i=ifp(j,l),ilp(j,l)
 13       fld(i,j)=fld(i,j)+clip
        end if
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional code for checking conservation properties
ccc      after=0.
ccc      do 15 j=1,jj
ccc      do 15 l=1,isp(j)
ccc      do 15 i=ifp(j,l),ilp(j,l)
ccc 15   after=after+fld(i,j)*fc(i,j)*scal(i,j)
ccc      write (*,'(a,3es14.6,e11.1)') 'advem conservation:',
ccc     .  bfore,after,after-bfore,(after-bfore)/bfore
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
c
c
c> Revision history:
c>
c> Mar. 2000 - removed 'cushn' and added logic to assure global conservation
c> Apr. 2000 - conversion to SI units
c> Apr. 2000 - changed i/j loop nesting to j/i
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Sep. 2000 - fixed cyclicity problem in loop 33
c> Nov. 2004 - switched from mpdata to fct
c> Mar. 2006 - added bering strait exchange logic
