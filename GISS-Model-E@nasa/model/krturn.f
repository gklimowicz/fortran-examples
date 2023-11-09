#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
      real function krturn(pres,dens,thkold,buoyfl,u_star,corpar)
c
c --- compute Kraus-Turner mixed layer depth in a single column.
c --- sign convention: buoyancy flux is positive into the ocean, i.e.,
c --- surface density increases (column is destabilized) if buoyfl < 0
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS
      USE HYCOM_ARRAYS

      implicit none
c
      real,intent(IN) :: pres(kdm+1)	! pressure
     .                  ,dens(kdm)	! density
     .                  ,thkold		! old mixed layer depth
     .                  ,buoyfl		! buoyancy flux
     .                  ,u_star		! friction velocity, sqrt(tau/rho)
     .                  ,corpar		! coriolis parameter

      integer k
c
      real turgen,dpth,ekminv,obuinv,ex,alf1,alf2,cp1,cp3,ape,cc4,spe,
     .     em,en,ea1,ea2,em1,em2,em3,em4,em5,thknew,q,thermg,sum1,sum2,
     .     pnew,ustar3
      data ea1, ea2, em1, em2, em3, em4, em5
     .   /0.60,0.30,0.45,2.60,1.90,2.30,0.60/           ! Gaspar coefficients
      real,parameter :: athird=1./3.
c
      ustar3=u_star**3
c
c --- buoyancy flux = w_prime_buoyancy_prime_bar (m^2/sec^3)
ccc      tem=.5*(temp(i,j,k1m)+temp(i,j,k1n))
ccc      sal=.5*(saln(i,j,k1m)+saln(i,j,k1n))
ccc      buoyfl=-g*thref**2*(dsigds(tem,sal)*salflx(i,j)
ccc     .                   +dsigdt(tem,sal)*surflx(i,j)/spcifh)
c
c --- determine turb.kin.energy generation due to wind stirring (ustar)
c --- and diabatic forcing (buoyfl).
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 1 :   k r a u s  -  t u r n e r    mixed-layer t.k.e.  closure
c
      em=0.8*exp(-pres(2)/(50.*onem))   !   hadley centre choice (orig.: 1.25)
      en=0.15                           !   hadley centre choice (orig.: 0.4)
      thermg=-.5*((en+1.)*buoyfl+(en-1.)*abs(buoyfl))
      turgen=delt1*(2.*em*g*ustar3/thref+thkold*thermg)/thref**2
c
c --- find monin-obukhov length in case of receding mixed layer (turgen < 0).
c --- the monin-obukhov length is found by stipulating turgen = 0.
c
      if (turgen.lt.0.) then
        thknew=-2.*em*g*ustar3/min(-epsil,thref*thermg)
      else
        thknew=thkold
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 2 :    g a s p a r   mixed-layer t.k.e. closure
c
cc       dpth=thkold/onem
cc       ekminv=abs(corpar)/max(epsil,u_star)
cc       obuinv=buoyfl/max(epsil,ustar3)
cc       ex=exp(min(50.,dpth*obuinv))
cc       alf1=ea1+ea2*max(1.,2.5*dpth*ekminv)*ex
cc       alf2=ea1+ea2*ex
cc       cp1=((1.-em5)*(alf1/alf2)+.5*em4)*athird
cc       cp3=max(0.,(em4*(em2+em3)-(alf1/alf2)*(em2+em3-em3*em5))*athird)
cc       ape=cp3*ustar3-cp1*dpth*buoyfl
cc c
cc       if(ape.lt.0.) then					! detrainment
cc         turgen=(g*delt1/thref**3)*ape
cc         thknew=min(thkold,g*cp3/(thref*cp1*max(epsil,obuinv)))
cc c
cc       else							! entrainment
cc         cc4=2.*em4/(em1*em1) * alf1*alf1
cc         spe=(em2+em3)*ustar3-0.5*dpth*buoyfl
cc         turgen=(g*delt1/thref**3)*(sqrt((.5*ape-cp1*spe)**2
cc      .          +2.*cc4*ape*spe)-(.5*ape+cp1*spe))/(cc4-cp1)
cc         thknew=thkold
cc       end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- sum1,sum2 are used to evaluate pot.energy changes during entrainment
      sum1=dens(1)*thkold
      sum2=dens(1)*thkold**2
c
c --- find thknew in case of mixed layer deepening (turgen>0).
c --- entrain as many layers as needed to deplete -turgen-.
c
      if (turgen.ge.0.) then
        do 85 k=2,kk
        if (pres(k+1).gt.thkold) then
          q=max(pres(k),thkold)
          pnew=min(pres(k+1),(2.*turgen+dens(k)*q**2-sum2)/
     .                        max(epsil,dens(k)*q   -sum1))
c --- stop iterating for 'thknew' as soon as pnew < k-th interface pressure
          if (pnew.lt.pres(k)) exit
          thknew=pnew
          sum1=sum1+dens(k)*(pres(k+1)   -q   )
          sum2=sum2+dens(k)*(pres(k+1)**2-q**2)
        end if
 85     continue
      end if
c
c --- move mxlyr depth closer to deepest interface embedded in mixed layer.
c --- this is to artificially reduce mixing across mixed layer base.
c
ccc      do 84 k=2,kk
ccc      if (thknew.gt.pres(k) .and. thknew.lt.pres(k+1)) then
ccc        q=(thknew-pres(k))/(pres(k+1)-pres(k))
ccc        q=q*q
ccc        q=q*q
ccc        q=0.			! reduce mxlyr depth to nearest interface
ccc        thknew=pres(k)*(1.-q)+pres(k+1)*q
ccc      end if
ccc 84   continue

c --- don't allow mixed layer to get too deep or too shallow.
      krturn=min(pres(kk+1),max(pres(2),thknew))
c
      return
      end function krturn
