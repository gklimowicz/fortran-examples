#include "hycom_mpi_hacks.h"
      subroutine advem(iord,fld,u,v,scal,scali,dt,fco,fcn)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
      USE HYCOM_DIM
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT,HALO_UPDATE,NORTH,SOUTH,
     &                             GLOBALSUM
      USE HYCOM_SCALARS, only: itest,jtest,onecm,onemu

      implicit none
c
c combined monotone scheme, for details see section 3.3 (eqs. 34 to 37)
c in smolarkiewicz and clark, 1986, j.comput.phys.,67,no 2, p. 396-438
c and smolarkiewicz and grabowski, 1989, j.comput.phys.
c  fld     - transported mixing ratio, e.g., salinity or temperature
c  u,v     - mass fluxes satisfying continuity equation
c  scal    - grid cell size
c  scali   - inverse of scal
c  dt      - temporal increment
c  fco,fcn - depth of the layer at previous and new time step
c
      real,intent(INOUT) :: fld(idm,J_0H:J_1H)
      real,intent(IN)    :: u(idm,J_0H:J_1H),v(idm,J_0H:J_1H),
     .     scal(idm,J_0H:J_1H),scali(idm,J_0H:J_1H),
     .     fco(idm,J_0H:J_1H),fcn(idm,J_0H:J_1H),dt
      real fmx(idm,J_0H:J_1H),fmn(idm,J_0H:J_1H),
     .     flp(idm,J_0H:J_1H),fln(idm,J_0H:J_1H),
     .     flx(idm,J_0H:J_1H),fly(idm,J_0H:J_1H)
      real u1(idm,J_0H:J_1H),v1(idm,J_0H:J_1H),
     .     flxdiv(idm,J_0H:J_1H),vlumj(J_0H:J_1H),
     .     bforj(J_0H:J_1H),aftrj(J_0H:J_1H),
     .     thko(idm,J_0H:J_1H),thkn(idm,J_0H:J_1H)

      real q,clip,vlume,clipped,unclipd,bfore,after,div,offset
      integer iord,i,j,l,n,ia,ib,ja,jb,ip1,im1,jp1,jm1
      logical,parameter :: recovr=.true.
c
c --- if iord=1, scheme reduces to simple donor cell scheme.
c
      if (recovr) then
        do 14 j=J_0,J_1
        bforj(j)=0.
        do 14 l=1,isp(j)
        do 14 i=ifp(j,l),ilp(j,l)
 14     bforj(j)=bforj(j)+fld(i,j)*fco(i,j)*scal(i,j)
      end if

c --- generate auxiliary layer depths consistent with mass flux divergence
      do 12 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      div=(u(i+1,j)-u(i,j)+v(i,jb )-v(i,j))*scali(i,j)*dt
      thkn(i,j)=.5*(fco(i,j)+fcn(i,j)-div)
      thko(i,j)=.5*(fco(i,j)+fcn(i,j)+div)
      offset=min(0.,thko(i,j),thkn(i,j))
      thkn(i,j)=thkn(i,j)-offset
 12   thko(i,j)=thko(i,j)-offset
c
c --- compute low-order and part of antidiffusive fluxes
c
      call cpy_p_par(fld)

      CALL HALO_UPDATE(ogrid,fld, FROM=SOUTH+NORTH)
c
      do 11 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do i=1,ii
        fmx(i,j)=0.
        fmn(i,j)=0.
        flp(i,j)=0.
        fln(i,j)=0.
        flx(i,j)=0.
        fly(i,j)=0.
        u1 (i,j)=0.
        v1 (i,j)=0.
      end do

      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      u1(i,j)=.5*abs(u(i,j))*(fld(i,j)-fld(i-1,j))
      if (u(i,j).ge.0.) then
        q=fld(i-1,j)
      else
        q=fld(i  ,j)
      end if
    2 flx(i,j)=u(i,j)*q
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
      v1(i,j)=.5*abs(v(i,j))*(fld(i,j)-fld(i,ja ))
      if (v(i,j).ge.0.) then
        q=fld(i,ja )
      else
        q=fld(i,j  )
      end if
    3 fly(i,j)=v(i,j)*q
c
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      ja = PERIODIC_INDEX(j-1, jj)
      if (ip(i,ja).eq.0) ja=j
      jb = PERIODIC_INDEX(j+1, jj)
      if (ip(i,jb).eq.0) jb=j
      fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb))
   11 fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb))
c
cdiag i=itest
cdiag j=jtest
cdiag write (*,'("advem (1)",2i5,f22.3/es39.2/f21.3,es9.2,f9.3,
cdiag.es9.2,f9.3/es39.2/f39.3)') i,j,fld(i-1,j),u(i,j),fld(i,j-1),
cdiag.v(i,j),fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
      do 61 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      vlumj(j)=0.
      do 61 l=1,isp(j)
      do 61 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*dt
     .   *scali(i,j)
      unclipd=fld(i,j)*thko(i,j)-flxdiv(i,j)
      clipped=max(fmn(i,j)*thkn(i,j),min(unclipd,fmx(i,j)*thkn(i,j)))
      if (recovr) vlumj(j)=vlumj(j)+scal(i,j)*fcn(i,j)
   61 fld(i,j)=(fld(i,j)*onemu+clipped)/(onemu+thkn(i,j))
c
      if (iord.le.1) go to 100

      CALL HALO_UPDATE(ogrid,flxdiv, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,thko,   FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,thkn,   FROM=SOUTH)
c
c --- finish computation of antidiffusive fluxes
c
      do 8 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 7 l=1,isu(j)
      do 7 i=ifu(j,l),ilu(j,l)
    7 flx(i,j)=u1(i,j)-u(i,j)*(flxdiv(i,j)+flxdiv(i-1,j))
     .   /(thko(i,j)+thko(i-1,j)+thkn(i,j)+thkn(i-1,j)+onemu)
c
      do 8 l=1,isv(j)
      do 8 i=ifv(j,l),ilv(j,l)
    8 fly(i,j)=v1(i,j)-v(i,j)*(flxdiv(i,j)+flxdiv(i,ja ))
     .   /(thko(i,j)+thko(i,ja )+thkn(i,j)+thkn(i,ja )+onemu)
c
c---- limit antidiffusive fluxes
c
      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
      do 16 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 16 l=1,isp(j)
      do 16 i=ifp(j,l),ilp(j,l)
      flp(i,j)=(fmx(i,j)-fld(i,j))*thkn(i,j)*scal(i,j)/( (onemu
     .  -min(0.,flx(i+1,j))+max(0.,flx(i,j))
     .  -min(0.,fly(i,jb ))+max(0.,fly(i,j)) )*dt)
      fln(i,j)=(fld(i,j)-fmn(i,j))*thkn(i,j)*scal(i,j)/( (onemu
     .  +max(0.,flx(i+1,j))-min(0.,flx(i,j))
     .  +max(0.,fly(i,jb ))-min(0.,fly(i,j)) )*dt)
   16 continue
c
      call cpy_p_par(flp)
      call cpy_p_par(fln)

      CALL HALO_UPDATE(ogrid,fln, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,flp, FROM=SOUTH)
c
      do 18 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 17 l=1,isu(j)
      do 17 i=ifu(j,l),ilu(j,l)
      flx(i,j)=max(0.,flx(i,j))*min(1.,flp(i,j),fln(i-1,j))
     .        +min(0.,flx(i,j))*min(1.,flp(i-1,j),fln(i,j))
   17 continue
c
      do 18 l=1,isv(j)
      do 18 i=ifv(j,l),ilv(j,l)
      fly(i,j)=max(0.,fly(i,j))*min(1.,flp(i,j),fln(i,ja ))
     .        +min(0.,fly(i,j))*min(1.,flp(i,ja ),fln(i,j))
   18 continue

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
cdiag i=itest
cdiag j=jtest
cdiag write (*,'("advem (2)",2i5,f22.3/es39.2/f21.3,es9.2,f9.3,
cdiag.es9.2,f9.3/es39.2/f39.3)') i,j,fld(i-1,j),u(i,j),fld(i,ja ),
cdiag.v(i,j),fld(i,j),v(i,jb ),fld(i,jb ),u(i+1,j),fld(i+1,j)
c
      do 62 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 62 l=1,isp(j)
      do 62 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*dt
     .   *scali(i,j)
      unclipd=fld(i,j)*thkn(i,j)-flxdiv(i,j)
      clipped=max(fmn(i,j)*thkn(i,j),min(unclipd,fmx(i,j)*thkn(i,j)))
   62 fld(i,j)=(fld(i,j)*onemu+clipped)/(onemu+thkn(i,j))
c
  100 continue
c
c --- revover 'clipped' amount and return to field
c
      if (recovr) then
        call GLOBALSUM(ogrid,vlumj,vlume, all=.true.)
        if (vlume.gt.0.) then
          do 15 j=J_0,J_1
          aftrj(j)=0.
          do 15 l=1,isp(j)
          do 15 i=ifp(j,l),ilp(j,l)
 15       aftrj(j)=aftrj(j)+fld(i,j)*fcn(i,j)*scal(i,j)
          call GLOBALSUM(ogrid,bforj,bfore, all=.true.)
          call GLOBALSUM(ogrid,aftrj,after, all=.true.)
          clip=(bfore-after)/vlume
          do 13 j=J_0,J_1
          aftrj(j)=0.
          do 13 l=1,isp(j)
          do 13 i=ifp(j,l),ilp(j,l)
          fld(i,j)=fld(i,j)+clip
 13       aftrj(j)=aftrj(j)+fld(i,j)*fcn(i,j)*scal(i,j)
        end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- optional:
cc      call GLOBALSUM(ogrid,aftrj,after, all=.true.)
cc
cc      if (AM_I_ROOT())
cc   .   write (*,'(a,2f17.11,es11.2)') 'advem conservation:',
cc   .   bfore/vlume,after/vlume,(after-bfore)/vlume
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      end if			! recovr
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
c> Jul  2017 - new 'recovr' logic to assure exact global conservation

