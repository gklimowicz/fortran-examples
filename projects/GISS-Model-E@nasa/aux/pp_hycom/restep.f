      subroutine restep(uold,vold,temold,salold,trcold,dnsold,pold,
     .                  unew,vnew,temnew,salnew,trcnew,dnsnew,pnew,
     .                  theta,kold,knew)
c
c --- convert an idm x jdm array of stairstep (i.e., piecewise constant)
c --- u/v/T/S/dens/tracer profiles into an idm x jdm array of stairstep
c --- profiles constrained to have prescribed density ('theta') steps.
c
c --- input  variables: uold,vold,temold,salold,trcold,dnsold,pold,kold,theta
c --- output variables: unew,vnew,temnew,salnew,trcnew,dnsnew,pnew,knew
c
      use hycom_dimen
      use const_proc
      implicit none
c
      integer kold,knew,ko,kmx,i1,j1,ja,jb
c
      real uold(idm,jdm,kold),vold(idm,jdm,kold),temold(idm,jdm,kold),
     .     salold(idm,jdm,kold),dnsold(idm,jdm,kold),
     .     trcold(idm,jdm,kold),pold(idm,jdm,kold+1)
c
      real unew(idm,jdm,knew),vnew(idm,jdm,knew),temnew(idm,jdm,knew),
     .     salnew(idm,jdm,knew),dnsnew(idm,jdm,knew),
     .     trcnew(idm,jdm,knew),pnew(idm,jdm,knew+1),theta(knew)
c
      parameter (kmx=40)
      real cloutu(idm),cloutv(idm),cloutt(idm),clouts(idm),cloutr(idm),
     .     colinu(idm),colinv(idm),colint(idm),colins(idm),colinr(idm),
     .     pinteg,uinteg,vinteg,tinteg,sinteg,qinteg,phi,plo,pa,pb,siga,
     .     sigb,q,acurcy,oldsig(idm,kmx),puold,pvold,punew,pvnew
      logical at_top(idm),abort
      data acurcy/1.e-2/
css   data acurcy/1.e-9/
      data abort/.false./
c
c --- statement functions for computing interface pressure at u/v points:
c
ccc      puold(i,i1,j,k)=min(.5*(pold(i,j,k   )+pold(i1,j,k   )),
ccc     .                        pold(i,j,kold),pold(i1,j,kold))
ccc      punew(i,i1,j,k)=min(.5*(pnew(i,j,k   )+pnew(i1,j,k   )),
ccc     .                        pnew(i,j,knew),pnew(i1,j,knew))
ccc      pvold(i,j,j1,k)=min(.5*(pold(i,j,k   )+pold(i,j1,k   )),
ccc     .                        pold(i,j,kold),pold(i,j1,kold))
ccc      pvnew(i,j,j1,k)=min(.5*(pnew(i,j,k   )+pnew(i,j1,k   )),
ccc     .                        pnew(i,j,knew),pnew(i,j1,knew))
      puold(i,i1,j,k)=.5*(pold(i,j,k)+pold(i1,j,k))
      punew(i,i1,j,k)=.5*(pnew(i,j,k)+pnew(i1,j,k))
      pvold(i,j,j1,k)=.5*(pold(i,j,k)+pold(i,j1,k))
      pvnew(i,j,j1,k)=.5*(pnew(i,j,k)+pnew(i,j1,k))
c
      if (max(kold,knew).gt.kmx) then
        write (lp,'(a,i4)') 'restep error -- increase kmx to',
     .    max(kold,knew)
        abort=.true.
      end if
c
      do 1 j=1,jdm
c
      do 1 l=1,isp(j)
c
      do 2 k=1,knew
      do 2 i=ifp(j,l),ilp(j,l)
 2    dnsnew(i,j,k)=theta(k)
c
c --- remove density inversions from input profile
      do 29 i=ifp(j,l),ilp(j,l)
 29   oldsig(i,1)=dnsold(i,j,1)
      do 30 k=2,kold
      do 30 i=ifp(j,l),ilp(j,l)
 30   oldsig(i,k)=max(oldsig(i,k-1),dnsold(i,j,k))
c
      do 3 i=ifp(j,l),ilp(j,l)
 101  format (2i5,a/(30x,i3,f9.3,f10.2,2f9.3))
 102  format (2i5,a/(30x,i3,2(0p,f10.2,1p,e10.2)))
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,101) itest,jtest,
     .  '  restep -- old profile:   dens      dpth     temp     saln',
     .  (k,oldsig(i,k),pold(i,j,k+1),temold(i,j,k),salold(i,j,k),
     .  k=1,kold)
        ja=mod(j-2+jdm,jdm)+1
        write (lp,102) itest,jtest,
     .  '  restep -- old profile:    dpthu     u         dpthv     v',
     .  (k,puold(i,i-1,j,k+1),uold(i,j,k),
     .     pvold(i,j,ja ,k+1),vold(i,j,k),k=1,kold)
      end if
      dnsnew(i,j,   1)=min(oldsig(i,1),oldsig(i,kold),theta(   1))
      dnsnew(i,j,knew)=max(oldsig(i,1),oldsig(i,kold),theta(knew))
      pnew(i,j,     1)=pold(i,j,     1)
      pnew(i,j,knew+1)=pold(i,j,kold+1)
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      cloutr(i)=0.
      colinr(i)=oldsig(i,kold)*(pold(i,j,kold+1)-pold(i,j,kold))
      colint(i)=temold(i,j,kold)*(pold(i,j,kold+1)-pold(i,j,kold))
 3    colins(i)=salold(i,j,kold)*(pold(i,j,kold+1)-pold(i,j,kold))
c
      do 9 k=1,kold-1
      do 9 i=ifp(j,l),ilp(j,l)
      colinr(i)=colinr(i)+oldsig(i,k)*(pold(i,j,k+1)-pold(i,j,k))
      colint(i)=colint(i)+temold(i,j,k)*(pold(i,j,k+1)-pold(i,j,k))
 9    colins(i)=colins(i)+salold(i,j,k)*(pold(i,j,k+1)-pold(i,j,k))
c
c --- find interface depth pnew(k+1) separating layers k and k+1 by requiring 
c --- that integral over p*d(sigma) from dnsnew(k) to dnsnew(k+1) be preserved.
c
      do 4 k=1,knew-1
      do 4 i=ifp(j,l),ilp(j,l)
      pinteg=0.
      sigb=dnsnew(i,j,k)
      do 5 ko=1,kold
      siga=sigb
      sigb=min(dnsnew(i,j,k+1),max(dnsnew(i,j,k),oldsig(i,ko)))
      pinteg=pinteg+pold(i,j,ko)*(sigb-siga)
      if (oldsig(i,ko).ge.dnsnew(i,j,k+1)) go to 25
 5    continue
      siga=sigb
      sigb=dnsnew(i,j,k+1)
      pinteg=pinteg+pold(i,j,kold+1)*(sigb-siga)
c
 25   pnew(i,j,k+1)=pinteg/(dnsnew(i,j,k+1)-dnsnew(i,j,k))
      cloutr(i)=cloutr(i)+dnsnew(i,j,k)*(pnew(i,j,k+1)-pnew(i,j,k))
c --- remove effect of roundoff errors on monotonicity
ccc   pnew(i,j,k+1)=max(pnew(i,j,k),min(pnew(i,j,k+1),pnew(i,j,knew+1)))
 4    continue
c
      do 6 i=ifp(j,l),ilp(j,l)
      cloutr(i)=cloutr(i)+dnsnew(i,j,knew)*
     .   (pnew(i,j,knew+1)-pnew(i,j,knew))
      if (abs(cloutr(i)-colinr(i)).gt.acurcy*abs(colinr(i)))
     .  write (lp,100) i,j,' restep - dens.column intgl.error',
     .   colinr(i),cloutr(i),(cloutr(i)-colinr(i))/colinr(i)
 6    continue
c
c --- ----------------------------------------
c --- integrate -t,s- over new depth intervals
c --- ----------------------------------------
c
      do 10 i=ifp(j,l),ilp(j,l)
      at_top(i)=.true.
      clouts(i)=0.
 10   cloutt(i)=0.
c
      do 8 k=1,knew
      do 8 i=ifp(j,l),ilp(j,l)
      tinteg=0.
      sinteg=0.
      qinteg=0.
      phi=pnew(i,j,k+1)
      plo=pnew(i,j,k  )
      pb=plo
      do 7 ko=1,kold
      pa=pb
      pb=min(phi,max(plo,pold(i,j,ko+1)))
      tinteg=tinteg+temold(i,j,ko)*(pb-pa)
      sinteg=sinteg+salold(i,j,ko)*(pb-pa)
      qinteg=qinteg+trcold(i,j,ko)*(pb-pa)
      if (pa.ge.phi) go to 28
 7    continue
      pa=pb
      pb=phi
      tinteg=tinteg+temold(i,j,kold)*(pb-pa)
      sinteg=sinteg+salold(i,j,kold)*(pb-pa)
      qinteg=qinteg+trcold(i,j,kold)*(pb-pa)
 28   cloutt(i)=cloutt(i)+tinteg
      clouts(i)=clouts(i)+sinteg
c
      q=phi-plo
      if (q.gt.0.) then
        temnew(i,j,k)=tinteg/q
        salnew(i,j,k)=sinteg/q
        trcnew(i,j,k)=qinteg/q
        at_top(i)=.false.
      else if (at_top(i)) then
        temnew(i,j,k)=temold(i,j, 1)
        salnew(i,j,k)=salold(i,j, 1)
        trcnew(i,j,k)=trcold(i,j, 1)
      else
        temnew(i,j,k)=temold(i,j,kold)
        salnew(i,j,k)=salold(i,j,kold)
        trcnew(i,j,k)=trcold(i,j,kold)
      end if
 8    continue
c
      do 1 i=ifp(j,l),ilp(j,l)
      if (abs(cloutt(i)-colint(i)).gt.acurcy*5.*pold(i,j,kold+1))
     .  write (lp,100) i,j,' restep - t column intgl.error',
     .   colint(i),cloutt(i),(cloutt(i)-colint(i))/colint(i)
      if (abs(clouts(i)-colins(i)).gt.acurcy*5.*pold(i,j,kold+1))
     .  write (lp,100) i,j,' restep - s column intgl.error',
     .   colins(i),clouts(i),(clouts(i)-colins(i))/colins(i)
 100  format (2i5,a,1p,2e14.6,e9.1)
c
 1    continue
c
c     write (lp,103) itest,jtest,' old density profile:',
c    .   (dnsold(itest,jtest,k),k=1,kold)
c     write (lp,103) itest,jtest,' new density profile:',
c    .   (dnsnew(itest,jtest,k),k=1,knew)
 103  format (2i5,a/(8f9.3))
c
      do 21 j=1,jdm
      ja=mod(j-2+jdm,jdm)+1
c
      do 22 l=1,isu(j)
c
c --- --------------------------------------
c --- integrate -u- over new depth intervals
c --- --------------------------------------
c
      do 17 i=ifu(j,l),ilu(j,l)
      if (puold(i,i-1,j,kold+1).ne.punew(i,i-1,j,knew+1)) then
 104    format (3i4,2a,2i3/2(5x,3f9.2))
        write (*,104) i-1,i,j,'  error restep - ',
     .   'puold/new mismatch, kold/new =',kold,knew,
     .    pold(i-1,j,kold+1),pold(i,j,kold+1),puold(i,i-1,j,kold+1),
     .    pnew(i-1,j,knew+1),pnew(i,j,knew+1),punew(i,i-1,j,knew+1)
        abort=.true.
      end if
 17   colinu(i)=0.
c
      do 18 k=1,kold
      do 18 i=ifu(j,l),ilu(j,l)
 18   colinu(i)=colinu(i)+
     .   uold(i,j,k)*(puold(i,i-1,j,k+1)-puold(i,i-1,j,k))
c
      do 11 i=ifu(j,l),ilu(j,l)
      at_top(i)=.true.
 11   cloutu(i)=0.
c
      do 12 k=1,knew
      do 12 i=ifu(j,l),ilu(j,l)
      uinteg=0.
      phi=punew(i,i-1,j,k+1)
      plo=punew(i,i-1,j,k  )
      pb=plo
      do 13 ko=1,kold
      pa=pb
      pb=min(phi,max(plo,puold(i,i-1,j,ko+1)))
      uinteg=uinteg+uold(i,j,ko)*(pb-pa)
      if (pa.ge.phi) go to 26
 13   continue
      pa=pb
      pb=phi
      uinteg=uinteg+uold(i,j,kold)*(pb-pa)
 26   cloutu(i)=cloutu(i)+uinteg
c
      q=phi-plo
      if (q.gt.0.) then
        unew(i,j,k)=uinteg/q
        at_top(i)=.false.
      else if (at_top(i)) then
        unew(i,j,k)=uold(i,j,   1)
      else
        unew(i,j,k)=uold(i,j,kold)
      end if
 12   continue
c
      do 22 i=ifu(j,l),ilu(j,l)
      if (abs(cloutu(i)-colinu(i)).gt.acurcy*max(abs(colinu(i)),
     .  10.*(pold(i,j,kold+1)+pold(i-1,j,kold+1))))
     .  write (lp,100) i,j,' restep - u column intgl.error',
     .   colinu(i),cloutu(i),(cloutu(i)-colinu(i))/colinu(i)
 22   continue
c
c --- --------------------------------------
c --- integrate -v- over new depth intervals
c --- --------------------------------------
c
      do 23 l=1,isv(j)
c
      do 19 i=ifv(j,l),ilv(j,l)
      if (pvold(i,j,ja ,kold+1).ne.pvnew(i,j,ja ,knew+1)) then
        write (*,104) i,ja,j,'  error restep - ',
     .   'pvold/new mismatch, kold/new =',kold,knew,
     .    pold(i,ja ,kold+1),pold(i,j,kold+1),pvold(i,j,ja ,kold+1),
     .    pnew(i,ja ,knew+1),pnew(i,j,knew+1),pvnew(i,j,ja ,knew+1)
        abort=.true.
      end if
 19   colinv(i)=0.
c
      do 20 k=1,kold
      do 20 i=ifv(j,l),ilv(j,l)
 20   colinv(i)=colinv(i)+
     .   vold(i,j,k)*(pvold(i,j,ja ,k+1)-pvold(i,j,ja ,k))
c
      do 14 i=ifv(j,l),ilv(j,l)
      at_top(i)=.true.
 14   cloutv(i)=0.
c
      do 15 k=1,knew
      do 15 i=ifv(j,l),ilv(j,l)
      vinteg=0.
      phi=pvnew(i,j,ja ,k+1)
      plo=pvnew(i,j,ja ,k  )
      pb=plo
      do 16 ko=1,kold
      pa=pb
      pb=min(phi,max(plo,pvold(i,j,ja ,ko+1)))
      vinteg=vinteg+vold(i,j,ko)*(pb-pa)
      if (pa.ge.phi) go to 27
 16   continue
      pa=pb
      pb=phi
      vinteg=vinteg+vold(i,j,kold)*(pb-pa)
 27   cloutv(i)=cloutv(i)+vinteg
c
      q=phi-plo
      if (q.gt.0.) then
        vnew(i,j,k)=vinteg/q
        at_top(i)=.false.
      else if (at_top(i)) then
        vnew(i,j,k)=vold(i,j,   1)
      else
        vnew(i,j,k)=vold(i,j,kold)
      end if
 15   continue
c
      do 23 i=ifv(j,l),ilv(j,l)
      if (abs(cloutv(i)-colinv(i)).gt.acurcy*max(abs(colinv(i)),
     .  10.*(pold(i,j,kold+1)+pold(i,ja ,kold+1))))
     .  write (lp,100) i,j,' restep - v column intgl.error',
     .   colinv(i),cloutv(i),(cloutv(i)-colinv(i))/colinv(i)
 23   continue
c
      do 21 l=1,isp(j)
      do 21 i=ifp(j,l),ilp(j,l)
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,101) itest,jtest,
     .  '  restep -- new profile:   dens      dpth     temp     saln',
     .  (k,dnsnew(i,j,k),pnew(i,j,k+1),temnew(i,j,k),salnew(i,j,k),
     .  k=1,knew)
        write (lp,102) itest,jtest,
     .  '  restep -- new profile:    dpthu     u         dpthv     v',
     .  (k,punew(i,i-1,j,k+1),unew(i,j,k),
     .     pvnew(i,j,ja ,k+1),vnew(i,j,k),k=1,knew)
      end if
c
 21   continue
      if (abort) stop '(reflux error)'
c
      return
      end
