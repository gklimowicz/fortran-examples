! WHIZARD 3.1.0 Dec 14 2022

! TOPPIK code by M. Jezabek, T. Teubner (v1.1, 1992), T. Teubner (1998)
!
! FB: -commented out numerical recipes code for hypergeometric 2F1
!      included in hypgeo.f90;
!     -commented out unused function 'ZAPVQ1';
!     -replaced function 'cdabs' by 'abs';
!     -replaced function 'dimag' by 'aimag';
!     -replaced function 'dcmplx(,)' by 'cmplx(,,kind=kind(0d0))';
!     -replaced function 'dreal' by 'real';
!     -replaced function 'cdlog' by 'log';
!     -replaced PAUSE by PRINT statement to avoid compiler warning;
!     -initialized 'idum' explicitly as real to avoid compiler warning.
!     -modified 'adglg1', 'adglg2' and 'tttoppik' to catch unstable runs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c *********************************************************************
c
c Working version with all the different original potentials
c  like (p^2+q^2)/|p-q|^2, not transformed in terms of delta and 1/r^2;
c accuracy eps=1.d-3 possible (only), but should be save, 13.8.'98, tt.
c
c *********************************************************************
c
        subroutine tttoppik(xenergy,xtm,xtg,xalphas,xscale,xcutn,xcutv,
     u     xc0,xc1,xc2,xcdeltc,xcdeltl,xcfullc,xcfulll,xcrm2,
     u     xkincm,xkinca,jknflg,jgcflg,
     u     xkincv,jvflg,xim,xdi,np,xpp,xww,xdsdp,zvfct)
c
c *********************************************************************
c
c !! THIS IS NOT A PUBLIC VERSION !!
c
c -- Calculation of the Green function in momentum space by solving the
c     Lippmann-Schwinger equation
c     G(p) = G_0(p) + G_0(p) int_0^xcutn V(p,q) G(q) dq
c
c -- Written by Thomas Teubner, Hamburg, November 1998
c     * Based on TOPPIK Version 1.1
c        from M. Jezabek and TT, Karlsruhe, June 1992
c     * Version originally for non-constant top-width
c     * Constant width supplied here
c     * No generator included
c
c -- Use of double precision everywhere
c
c -- All masses, momenta, energies, widths in GeV
c
c -- Input parameters:
c
c    xenergy  :  E=Sqrt[s]-2*topmass
c    xtm      :  topmass (in the Pole scheme)
c    xtg      :  top-width
c    xalphas  :  alpha_s^{MSbar,n_f=5}(xscale)
c    xscale   :  soft scale  mu_{soft}
c    xcutn    :  numerical UV cutoff on all momenta
c                (UV cutoff of the Gauss-Legendre grid)
c    xcutv    :  renormalization cutoff on the
c                 delta-, the (p^2+q^2)/(p-q)^2-, and the
c                  1/r^2-[1/|p-q|]-potential:
c                 if (max(p,q).ge.xcutv) then the three potentials
c                  are set to zero in the Lippmann-Schwinger equation
c    xc0      :  0th order coefficient for the Coulomb potential,
c                 see calling example above
c    xc1      :  1st order coefficient for the Coulomb potential
c    xc2      :  2nd order coefficient for the Coulomb potential
c    xcdeltc  :  constant of the delta(r)-
c                 [= constant in momentum space-] potential
c    xcdeltl  :  constant for the additional log(q^2/mu^2)-part of the
c                 delta-potential:
c                  xcdeltc*1 + xcdeltl*log(q^2/mu^2)
c    xcfullc  :  constant of the (p^2+q^2)/(p-q)^2-potential
c    xcfulll  :  constant for the additional log(q^2/mu^2)-part of the
c                 (p^2+q^2)/(p-q)^2-potential
c    xcrm2    :  constant of the 1/r^2-[1/|p-q|]-potential
c    xkincm   :  } kinetic corrections in the 0th order Green-function:
c    xkinca   :  }  G_0(p):=1/[E+iGamma_t-p^2/m_t]*(1+xkincm)+xkinca
c     !!! WATCH THE SIGN IN G_0 !!!
c    jknflg   :  flag for these kinetic corrections:
c                 0 : no kinetic corrections applied
c                 1 : kinetic corrections applied with cutoff xcutv
c                      for  xkinca  only
c                 2 : kinetic corrections applied with cutoff xcutv
c                      for  xkinca  AND  xkincm
c    jgcflg   :  flag for G_0(p) in the LS equation:
c                 0 (standard choice) : G_0(p) as given above
c                 1 (for TIPT)        : G_0(p) = G_c^{0}(p) the 0th
c                                        order Coulomb-Green-function
c                                        in analytical form; not for
c                                        momenta  p > 1000*topmass
c    xkincv   :  additional kinematic vertexcorrection in G_0, see below:
c    jvflg    :  flag for the additional vertexcorrection  xkincv  in the
c                 ``zeroth order'' G_0(p) in the LS-equation:
c                 0 : no correction, means  G = G_0 + G_0 int V G
c                      with G_0=1/[E+iGamma_t-p^2/m_t]*(1+xkincm)+xkinca
c                 1 : apply the correction in the LS equation as
c                      G = G_0 + xkincv*p^2/m_t^2/[E+iGamma_t-p^2/m_t] +
c                          G_0 int V G
c                     and correct the integral over Im[G(p)] to get sigma_tot
c                     from the optical theorem by the same factor.
c                     The cutoff  xcutv  is applied for these corrections.
c
c -- Output:
c
c    xim      :  R_{ttbar} from the imaginary part of the green
c                 function
c    xdi      :  R_{ttbar} form the integral over the momentum
c                 distribution (no cutoff but the numerical one here!!)
c    np       :  number of points used for the grid; fixed in tttoppik
c    xpp      :  1-dim array (max. 900 elements) giving the momenta of
c                 the Gauss-Legendre grid (pp(i) in the code)
c    xww      :  1-dim array (max. 900 elements) giving the corresponding
c                 Gauss-Legendre weights for the grid
c    xdsdp    :  1-dim array (max. 900 elements) giving the
c                 momentum distribution of top: d\sigma/dp,
c                  normalized to R,
c                  at the momenta of the Gauss-Legendre grid xpp(i)
c    zvfct    :  1-dim array (max. 900 elements) of COMPLEX*16 numbers
c                 giving the vertex function K(p), G(p)=K(p)*G_0(p)
c                 at the momenta of the grid
c
c *********************************************************************
c
c
           implicit none
           real*8
     u        pi,energy,vzero,eps,
     u        pp,
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        xx,critp,consde,
     u        w1,w2,sig1,sig2,const,
     u        gtpcor,etot,
     u        xenergy,xtm,xtg,xalphas,xscale,xc0,xc1,xc2,xim,xdi,
     u        xdsdp,xpp,xww,
     u        cplas,scale,c0,c1,c2,cdeltc,cdeltl,cfullc,cfulll,crm2,
     u        xcutn,dcut,xcutv,
     u        xp,xpmax,hmass,
     u        kincom,kincoa,kincov,xkincm,xkinca,xkincv,
     u        xcdeltc,xcdeltl,xcfullc,xcfulll,xcrm2,chiggs
           complex*16 bb,gg,a1,a,g0,g0c,zvfct
           integer i,n,nmax,npot,np,gcflg,kinflg,jknflg,jgcflg,
     u             jvflg,vflag
           parameter (nmax=900)
           dimension pp(nmax), bb(nmax), xx(nmax), gg(nmax),
     u               w1(nmax), w2(nmax), a1(nmax),
     u               xdsdp(nmax),xpp(nmax),xww(nmax),zvfct(nmax)
c
           external a,gtpcor,g0,g0c
c
           common/ovalco/ pi, energy, vzero, eps, npot
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/cplcns/cplas,scale,c0,c1,c2,
     u                   cdeltc,cdeltl,cfullc,cfulll,crm2,chiggs
           common/mom/ xp,xpmax,dcut
           common/g0inf/kincom,kincoa,kincov,kinflg,gcflg,vflag
c
           pi=3.141592653589793238d0
c
c Number of points to evaluate on the integral equation
c  (<=900 and n mod 3 = 0 !!):
c          n=66
           n=600
           np=n
c
c For second order potential with free parameters:
c
           npot=5
c Internal accuracy for TOPPIK, the reachable limit may be smaller,
c  depending on the parameters. But increase in real accuracy only
c  in combination with large number of points.
           eps=1.d-3
c Some physical parameters:
           wgamma=2.07d0
           zmass=91.187d0
           wmass=80.33d0
           bmass=4.7d0
c
c Input:
           energy=xenergy
           tmass=xtm
           tgamma=xtg
           cplas=xalphas
           scale=xscale
           c0=xc0
           c1=xc1
           c2=xc2
           cdeltc=xcdeltc
           cdeltl=xcdeltl
           cfullc=xcfullc
           cfulll=xcfulll
           crm2=xcrm2
           kincom=xkincm
           kincoa=xkinca
           kincov=xkincv
           kinflg=jknflg
           gcflg=jgcflg
           vflag=jvflg
c
           alphas=xalphas
c
c Cut for divergent potential-terms for large momenta in the function vhat
c  and in the integrals a(p):
           dcut=xcutv
c
c Numerical Cutoff of all momenta (maximal momenta of the grid):
           xpmax=xcutn
           if (dcut.gt.xpmax) then
              write(*,*) ' dcut > xpmax  makes no sense! Stop.'
              stop
           endif
c
c Not needed for the fixed order potentials:
           alamb5=0.2d0
c
c      WRITE(*,*) 'INPUT TGAMMA=',TGAMMA
c Needed in subroutine GAMMAT:
           GFERMI=1.16637d-5
c           CALL GAMMAT
c           WRITE(*,*) 'CALCULATED TGAMMA=',TGAMMA
c
           etot=2.d0*tmass+energy
c
           if ((npot.eq.1).or.(npot.eq.3).or.(npot.eq.4).or.
     u         (npot.eq.5)) then
c For pure coulomb and fixed order potentials there is no delta-part:
              consde = 0.d0
           else if (npot.eq.2) then
c Initialize QCD-potential common-blocks and calculate constant multiplying
c  the delta-part of the 'qcutted' potential in momentum-space:
              call iniphc(1)
              call vqdelt(consde)
           else
              write (*,*) ' Potential not implemented! Stop.'
              stop
           endif
c Delta-part of potential is absorbed by subtracting vzero from the
c  original energy (shift from the potential to the free Hamiltonian):
           vzero = consde / (2.d0*pi)**3
c          write (*,*) 'vzero=', vzero
c
c Find x-values pp(i) and weigths w1(i) for the gaussian quadrature;
c  care about large number of points in the important intervals:
c       if (energy-vzero.le.0.d0) then
cc         call gauleg(0.d0, 1.d0, pp, w1, n/3)
cc         call gauleg(1.d0, 5.d0, pp(n/3+1), w1(n/3+1), n/3)
cc         call gauleg(0.d0, 0.2d0, pp(2*n/3+1), w1(2*n/3+1), n/3)
c          call gauleg(0.d0, 5.d0, pp, w1, n/3)
c          call gauleg(5.d0, 20.d0, pp(n/3+1), w1(n/3+1), n/3)
c          call gauleg(0.d0, 0.05d0, pp(2*n/3+1), w1(2*n/3+1), n/3)
c       else
cc Avoid numerical singular points in the inner of the intervals:
c          critp = dsqrt((energy-vzero)*tmass)
c          if (critp.le.1.d0) then
cc Gauss-Legendre is symmetric => automatically principal-value prescription:
c             call gauleg(0.d0, 2.d0*critp, pp, w1, n/3)
c             call gauleg(2.d0*critp, 20.d0, pp(n/3+1),
c     u                    w1(n/3+1), n/3)
c             call gauleg(0.d0, 0.05d0, pp(2*n/3+1), w1(2*n/3+1), n/3)
c          else
cc Better behaviour at the border of the intervals:
c             call gauleg(0.d0, critp, pp, w1, n/3)
c             call gauleg(critp, 2.d0*critp, pp(n/3+1),
c     u                    w1(n/3+1), n/3)
c             call gauleg(0.d0, 1.d0/(2.d0*critp), pp(2*n/3+1),
c     u                    w1(2*n/3+1), n/3)
c          endif
c       endif
c
c Or different (simpler) method, good for V_JKT:
           if (energy.le.0.d0) then
              critp=tmass/3.d0
           else
              critp=max(tmass/3.d0,2.d0*dsqrt(energy*tmass))
           endif
           call gauleg(0.d0, critp, pp, w1, 2*n/3)
           call gauleg(1.d0/xpmax, 1.d0/critp, pp(2*n/3+1),
     u                 w1(2*n/3+1), n/3)
c
c Do substitution p => 1/p for the last interval explicitly:
           do 10 i=2*n/3+1,n
              pp(i) = 1.d0/pp(i)
10         continue
c
c Reorder the arrays for the third interval:
           do 20 i=1,n/3
              xx(i) = pp(2*n/3+i)
              w2(i) = w1(2*n/3+i)
20         continue
           do 30 i=1,n/3
              pp(n-i+1) = xx(i)
              w1(n-i+1) = w2(i)
30         continue
c
c Calculate the integrals a(p) for the given momenta pp(i)
c  and store weights and momenta for the output arrays:
           do 40 i=1,n
              a1(i) = a(pp(i)) !!! FB: can get stuck in original Toppik!
              !!! FB: abuse 'np' as a flag to communicate unstable runs
              if ( abs(a1(i)) .gt. 1d10 ) then
                np = -1
                return
              endif
              xpp(i)=pp(i)
              xww(i)=w1(i)
40         continue
           do 41 i=n+1,nmax
              xpp(i)=0.d0
              xww(i)=0.d0
41         continue
c
c Solve the integral-equation by solving a system of algebraic equations:
           call sae(pp, w1, bb, a1, n)
c
c (The substitution for the integration to infinity  pp => 1/pp
c  is done already.)
           do 50 i=1,n
              zvfct(i)=bb(i)
              gg(i) = bb(i)*g0c(pp(i))
cc            gg(i) = (1.d0 + bb(i))*g0c(pp(i))
cc Urspruenglich anderes (Minus) VZ hier, dafuer kein Minus mehr bei der
cc  Definition des WQs ueber Im G, 2.6.1998, tt.
cc            gg(i) = - (1.d0 + bb(i))*g0c(pp(i))
50         continue
c
c Normalisation on R:
           const = 8.d0*pi/tmass**2
c
c Proove of the optical theorem for the output values of sae:
c  Simply check if sig1 = sig2.
           sig1 = 0.d0
           sig2 = 0.d0
           do 60 i=1,n*2/3
c             write(*,*) 'check! p(',i,') = ',pp(i)
cvv
              if (pp(i).lt.dcut.and.vflag.eq.1) then
                 sig1 = sig1 + w1(i)*pp(i)**2*aimag(gg(i)
cc     u                 *(1.d0+kincov*(pp(i)/tmass)**2)
     u   *(1.d0+kincov*g0(pp(i))*(pp(i)/tmass)**2/g0c(pp(i)))
     u                  )
              else
                 sig1 = sig1 + w1(i)*pp(i)**2*aimag(gg(i))
              endif
              if (pp(i).lt.dcut.and.kinflg.ne.0) then
                 sig2 = sig2 + w1(i)*pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
     u                  *(1.d0-pp(i)**2/2.d0/tmass**2)
cc     u                  *tmass/dsqrt(tmass**2+pp(i)**2)
                 xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
     u                  *(1.d0-pp(i)**2/2.d0/tmass**2)
     u                  /(2.d0*pi**2)*const
              else
                 sig2 = sig2 + w1(i)*pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
                 xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
     u                  /(2.d0*pi**2)*const
              endif
c             write(*,*) 'xdsdp = ',xdsdp(i)
c             write(*,*) 'zvfct = ',zvfct(i)
60         continue
c '*p**2' because of substitution p => 1/p in the integration of p**2*G(p)
c  to infinity
           do 70 i=n*2/3+1,n
c             write(*,*) 'check! p(',i,') = ',pp(i)
cvv
              if (pp(i).lt.dcut.and.vflag.eq.1) then
                 sig1 = sig1 + w1(i)*pp(i)**4*aimag(gg(i)
cc     u                 *(1.d0+kincov*(pp(i)/tmass)**2)
     u   *(1.d0+kincov*g0(pp(i))*(pp(i)/tmass)**2/g0c(pp(i)))
     u                  )
              else
                 sig1 = sig1 + w1(i)*pp(i)**4*aimag(gg(i))
              endif
              if (pp(i).lt.dcut.and.kinflg.ne.0) then
                 sig2 = sig2 + w1(i)*pp(i)**4*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
     u                  *(1.d0-pp(i)**2/2.d0/tmass**2)
cc     u                  *tmass/dsqrt(tmass**2+pp(i)**2)
                 xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
     u                  *(1.d0-pp(i)**2/2.d0/tmass**2)
     u                  /(2.d0*pi**2)*const
              else
                 sig2 = sig2 + w1(i)*pp(i)**4*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
                 xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
     u                  /(2.d0*pi**2)*const
              endif
c
c             write(*,*) 'xdsdp = ',xdsdp(i)
c             write(*,*) 'zvfct = ',zvfct(i)
70         continue
           do 71 i=n+1,nmax
             xdsdp(i)=0.d0
             zvfct(i)=(0.d0,0.d0)
71         continue
c
c Normalisation on R:
           sig1  = sig1 / (2.d0*pi**2) * const
           sig2  = sig2 / (2.d0*pi**2) * const
c
c The results from the momentum space approach finally are:
cc Jetzt Minus hier, 2.6.98, tt.
           xim=-sig1
           xdi=sig2
c
        end
c
c
        complex*16 function g0(p)
c
           implicit none
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi,energy,vzero,eps,
     u        p,gtpcor,hmass
           integer npot
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           external gtpcor
           save
           g0=1.d0/cmplx(energy-vzero-p**2/tmass,
     u                    tgamma*gtpcor(p,2.d0*tmass+energy),
     u                    kind=kind(0d0))
        end
c
        complex*16 function g0c(p)
c
           implicit none
           complex*16 hypgeo,green,zk,zi,amd2k,aa,bb,cc,zzp,zzm,
     u                hypp,hypm,g0
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi,energy,vzero,eps,
     u        p,gtpcor,hmass,
     u        kincom,kincoa,kincov,xp,xpmax,dcut
           integer npot,kinflg,gcflg,vflag
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           common/g0inf/kincom,kincoa,kincov,kinflg,gcflg,vflag
           common/mom/ xp,xpmax,dcut
           external hypgeo,gtpcor,g0
           save
c
           if (gcflg.eq.0) then
              if (kinflg.eq.0) then
                 g0c=g0(p)
              else if (kinflg.eq.1.and.p.lt.dcut) then
                 g0c=g0(p)*(1.d0+kincom)+kincoa
              else if (kinflg.eq.1.and.p.ge.dcut) then
                 g0c=g0(p)*(1.d0+kincom)
              else if (kinflg.eq.2.and.p.lt.dcut) then
                 g0c=g0(p)*(1.d0+kincom)+kincoa
              else if (kinflg.eq.2.and.p.ge.dcut) then
                 g0c=g0(p)
              else
                 write(*,*) ' kinflg wrong! Stop.'
                 stop
              endif
           else if (gcflg.eq.1) then
              zi=(0.d0,1.d0)
              zk=-tmass*cmplx(energy,tgamma
     u                         *gtpcor(p,2.d0*tmass+energy),
     u                         kind=kind(0d0))
              zk=sqrt(zk)
              amd2k=4.d0/3.d0*alphas*tmass/2.d0/zk
              aa=(2.d0,0.d0)
              bb=(1.d0,0.d0)
              cc=2.d0-amd2k
              zzp=(1.d0+zi*p/zk)/2.d0
              zzm=(1.d0-zi*p/zk)/2.d0
              if (abs(zzp).gt.20.d0) then
                 hypp=(1.d0-zzp)**(-aa)*
     u                hypgeo(aa,cc-bb,cc,zzp/(zzp-1.d0))
              else
                 hypp=hypgeo(aa,bb,cc,zzp)
              endif
              if (abs(zzm).gt.20.d0) then
                 hypm=(1.d0-zzm)**(-aa)*
     u                hypgeo(aa,cc-bb,cc,zzm/(zzm-1.d0))
              else
                 hypm=hypgeo(aa,bb,cc,zzm)
              endif
              green=-zi*tmass/(4.d0*p*zk)/(1.d0-amd2k)*(hypp-hypm)
c VZ anders herum als in Andres Konvention, da bei ihm G_0=1/[-E-i G+p^2/m]:
              g0c=-green
              if (p.gt.1.d3*tmass) then
                 write(*,*) ' g0cana = ',g0c,' not reliable. Stop.'
                 stop
              endif
           else
              write(*,*) ' gcflg wrong! Stop.'
              stop
           endif
c
        end
c
c
        complex*16 function a(p)
c
           implicit none
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy,ETOT,vzero, eps,
     $        QCUT,QMAT1,ALR,PCUT,
     u        p,
     u        xp,xpmax, xb1,xb2,dcut,ddcut,
     u        a1, a2, a3, a4,a5,a6,
     u        adglg1, fretil1, fretil2, fimtil1, fimtil2,
     u        ALEFVQ, gtpcor, ad8gle, buf,adglg2,
c     u        xerg,
     u        kincom,kincoa,kincov,hmass
!          complex*16 zapvq1,ZAPVGP
           complex*16 ZAPVGP !!! FB
c     u                ,acomp
           integer npot,ILFLAG,kinflg,gcflg,vflag
c
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      COMMON/PARFLG/ QCUT,QMAT1,ALR,ILFLAG
           common/ovalco/ pi, energy, vzero, eps, npot
           common/mom/ xp,xpmax,dcut
           common/g0inf/kincom,kincoa,kincov,kinflg,gcflg,vflag
c
           external adglg1, fretil1, fretil2, fimtil1, fimtil2,
!     u              zapvq1, ALEFVQ, gtpcor,ZAPVGP,ad8gle,adglg2
     u              ALEFVQ, gtpcor,ZAPVGP,ad8gle,adglg2 !!! FB
c
           if ((npot.eq.1).or.(npot.eq.3).or.(npot.eq.4).or.
     u         (npot.eq.5)) then
c
              xp=p
              buf=0.d0
c
              a1=0.d0
              a2=0.d0
              a3=0.d0
              a4=0.d0
              a5=0.d0
              a6=0.d0
              if (gcflg.eq.0) then
                 ddcut=xpmax
              else if (gcflg.eq.1) then
                 ddcut=dcut
              else
                 write(*,*) ' gcflg wrong! Stop.'
                 stop
              endif
c
              if (2.d0*xp.lt.ddcut) then
                 xb1=xp
                 xb2=2.d0*xp
c
c More stable for logarithmically divergent fixed order potentials:
c
                 a1=adglg1(fretil1, buf, xb1, eps) !!! FB: can get stuck!
                 a2=adglg1(fimtil1, buf, xb1, eps)
c Slightly unstable:
                 a3=adglg2(fretil1,xb1,xb2,eps) !!! FB: can get stuck!
c No good:
c                a3=adglg1(fretil1,xb1,xb2,eps)
c Not better:
c                call adqua(xb1,xb2,fretil1,xerg,eps)
c                a3=xerg
c Also not better:
c                a1=adglg1(fretil1, buf, xb2, eps)
c
                 a4=adglg2(fimtil1,xb1,xb2,eps)
c                a5 = adglg2(fretil1, xb2, ddcut, eps)
c                a6 = adglg2(fimtil1, xb2, ddcut, eps)
                 a5 = adglg2(fretil2, 1.d0/ddcut, 1.d0/xb2, eps)
                 a6 = adglg2(fimtil2, 1.d0/ddcut, 1.d0/xb2, eps)
              else if (xp.lt.ddcut) then
                 xb1=xp
                 xb2=ddcut
                 a1=adglg1(fretil1, buf, xb1, eps)
                 a2=adglg1(fimtil1, buf, xb1, eps)
                 a3=adglg2(fretil1,xb1,xb2,eps)
                 a4=adglg2(fimtil1,xb1,xb2,eps)
              else if (ddcut.le.xp) then
              else
                 write(*,*) ' Constellation not possible! Stop.'
                 stop
              endif
c
              a  = 1.d0/(4.d0*pi**2)*cmplx(a1+a3+a5,a2+a4+a6,
     u                    kind=kind(0d0))
c
           else if (npot.eq.2) then
      PCUT=QCUT
      ETOT=ENERGY+2*TMASS
              a  = ZAPVGP(P,ETOT,VZERO-ENERGY,PCUT,EPS)
c             acomp = zapvq1(ALEFVQ, p, vzero-energy, gtpcor, eps)
c             a = zapvq1(ALEFVQ, p, vzero-energy, gtpcor, eps)
c             acomp = acomp/a
c             if (abs(acomp-1.d0).gt.1.d-3) then
c                write (*,*) 'p=', p
c                write (*,*) 'acomp/a=', acomp
c             endif
           else
              write (*,*) ' Potential not implemented! Stop.'
              stop
           endif
c
        end
c
        real*8 function fretil1(xk)
           implicit none
           real*8 xk, freal
           external freal
           fretil1 = freal(xk)
        end
c
        real*8 function fretil2(xk)
           implicit none
           real*8 xk, freal
           external freal
           fretil2 = freal(1.d0/xk) * xk**(-2)
        end
c
        real*8 function fimtil1(xk)
           implicit none
           real*8 xk, fim
           external fim
           fimtil1 = fim(xk)
        end
c
        real*8 function fimtil2(xk)
           implicit none
           real*8 xk, fim
           external fim
           fimtil2 = fim(1.d0/xk) * xk**(-2)
        end
c
        real*8 function freal(xk)
           implicit none
           complex*16 vhat
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        p,pmax, xk, gtpcor,dcut,hmass
           complex*16 g0,g0c
           integer npot
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           common/mom/ p,pmax,dcut
           external vhat, g0, g0c, gtpcor
c
           freal = real(g0c(xk)*vhat(p, xk)) !!! FB: NaN?
        end
c
        real*8 function fim(xk)
           implicit none
           complex*16 vhat
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        p,pmax, xk, gtpcor,dcut,hmass
           complex*16 g0,g0c
           integer npot
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           common/mom/ p,pmax,dcut
           external vhat, g0, g0c, gtpcor
           fim = aimag(g0c(xk)*vhat(p, xk))
        end
c
c
        complex*16 function vhat(p, xk)
c
           implicit none
           complex*16 zi
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        p, xk,
     u        cnspot, phiint, phfqcd, AD8GLE,
     u        pm, xkm, ALPHEF,
     u        zeta3,cf,ca,tf,xnf,a1,a2,b0,b1,
     u        cplas,scale,c0,c1,c2,
     u        cdeltc,cdeltl,cfullc,cfulll,crm2,
     u        xkpln1st,xkpln2nd,xkpln3rd,
     u        pp,pmax,dcut,hmass,chiggs
           integer npot
           parameter(zi=(0.d0,1.d0))
           parameter(zeta3=1.20205690316d0,
     u               cf=4.d0/3.d0,ca=3.d0,tf=1.d0/2.d0,
     u               xnf=5.d0)
c
           external AD8GLE, phfqcd, ALPHEF
c
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           common/pmaxkm/ pm, xkm
           common/mom/ pp,pmax,dcut
           common/cplcns/cplas,scale,c0,c1,c2,
     u                   cdeltc,cdeltl,cfullc,cfulll,crm2,chiggs
c
           b0=11.d0-2.d0/3.d0*xnf
           b1=102.d0-38.d0/3.d0*xnf
c
           a1=31.d0/9.d0*ca-20.d0/9.d0*tf*xnf
           a2=(4343.d0/162.d0+4.d0*pi**2-pi**4/4.d0+
     u         22.d0/3.d0*zeta3)*ca**2-
     u        (1798.d0/81.d0+56.d0/3.d0*zeta3)*ca*tf*xnf-
     u        (55.d0/3.d0-16.d0*zeta3)*cf*tf*xnf+
     u        (20.d0/9.d0*tf*xnf)**2
c
           pm=p
           xkm=xk
           cnspot=-4.d0/3.d0*4.d0*pi
c
           if (p/xk.le.1.d-5.and.p.le.1.d-5) then
              xkpln1st=2.d0
              xkpln2nd=-4.d0*dlog(scale/xk)
              xkpln3rd=-6.d0*dlog(scale/xk)**2
           else if (xk/p.le.1.d-5.and.xk.le.1.d-5) then
              xkpln1st=2.d0*(xk/p)**2
              xkpln2nd=-4.d0*(xk/p)**2*dlog(scale/p)
              xkpln3rd=-6.d0*(xk/p)**2*dlog(scale/p)**2
           else
c             xkpln1st=xk/p*dlog(dabs((p+xk)/(p-xk)))
              xkpln1st=xk/p*(dlog(p+xk)-dlog(dabs(p-xk)))
              xkpln2nd=xk/p*(-1.d0)*(dlog(scale/(p+xk))**2-
     u                               dlog(scale/dabs(p-xk))**2)
              xkpln3rd=xk/p*(-4.d0/3.d0)*(dlog(scale/(p+xk))**3-
     u                                    dlog(scale/dabs(p-xk))**3)
           endif
c
           if (npot.eq.2) then
              if (p/xk.le.1.d-5.and.p.le.1.d-5) then
                 vhat = 2.d0 * cnspot * ALPHEF(xk)
              else if (xk/p.le.1.d-5.and.xk.le.1.d-5) then
                 vhat = 2.d0 * cnspot * xk**2 / p**2 * ALPHEF(p)
              else
                 phiint = cnspot * (AD8GLE(phfqcd, 0.d0, 0.3d0, 1.d-5)
     u                            +AD8GLE(phfqcd, 0.3d0, 1.d0, 1.d-5))
                 vhat   = xk / p * dlog(dabs((p+xk)/(p-xk))) * phiint
              endif
           else
              if (npot.eq.1) then
                 c0=1.d0
                 c1=0.d0
                 c2=0.d0
              else if (npot.eq.3) then
                 c0=1.d0+alphas/(4.d0*pi)*a1
                 c1=alphas/(4.d0*pi)*b0
                 c2=0
              else if (npot.eq.4) then
                 c0=1.d0+alphas/(4.d0*pi)*a1+(alphas/(4.d0*pi))**2*a2
                 c1=alphas/(4.d0*pi)*b0+
     u             (alphas/(4.d0*pi))**2*(b1+2.d0*b0*a1)
                 c2=(alphas/(4.d0*pi))**2*b0**2
              else if (npot.eq.5) then
              else
                 write (*,*) ' Potential not implemented! Stop.'
                 stop
              endif
              phiint=cnspot*alphas
c
c             if ((xk+p).le.dcut) then
c                vhat=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c     u               -1.d0/2.d0*(1.d0+2.d0*ca/cf)
c     u                *(pi*cf*alphas)**2/tmass
c     u                *xk/p*(p+xk-dabs(xk-p))
c             else if (dabs(xk-p).lt.dcut) then
c                vhat=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c     u               -1.d0/2.d0*(1.d0+2.d0*ca/cf)
c     u                *(pi*cf*alphas)**2/tmass
c     u                *xk/p*(dcut-dabs(xk-p))
c             else if (dcut.le.dabs(xk-p)) then
c                vhat=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c             else
c                write(*,*) ' Not possible! Stop.'
c                stop
c             endif
c
              if (max(xk,p).lt.dcut) then
c Coulomb + first + second order corrections:
                 vhat=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c All other potentials:
     u               +cdeltc*2.d0*xk**2
     u               +cdeltl*xk/p/2.d0*(
     u                (p+xk)**2*(dlog(((p+xk)/scale)**2)-1.d0)-
     u                (p-xk)**2*(dlog(((p-xk)/scale)**2)-1.d0))
     u               +cfullc*(p**2+xk**2)*xkpln1st
     u               +cfulll*(p**2+xk**2)*xk/p/4.d0*
     u                 (dlog(((p+xk)/scale)**2)**2-
     u                  dlog(((p-xk)/scale)**2)**2)
     u               +crm2*xk/p*(p+xk-dabs(xk-p))
              else
                 vhat=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
              endif
           endif
c
        end
c
c
c
c --- Routines needed for use of phenomenological potentials ---
c
      SUBROUTINE INIPHC(INIFLG)
      implicit real*8(a-h,o-z)
      save
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      common/ovalco/ pi, energy, vzero, eps, npot
      COMMON/PARFLG/ QCUT,QMAT1,ALR,ILFLAG
      CHARACTER QCTCHR,QMTCHR,ALFCHR
      DATA QCUT0/.100d0/,QMT1S/5.0d0/
c
      zmass= 91.187d0
      if(INIFLG.eq.0) then
c     standard set of parameters
        ilflag= 1
        alphas=.12d0
        qcut= qcut0
        qmat1= qmt1s
      else
c     Parameters of QCD potential specified by USER
  5     write(*,*) 'QCD coupling at M_z:   ALPHAS  or  LAMBDA  ?'
        write(*,*) 'A/L  :'
        read(*,895) ALFCHR
          if(ALFCHR.eq.'A'.or.ALFCHR.eq.'a') then
              ilflag= 1
              write(*,*) 'alpha_s(M_z)= ?'
              read(*,*) alphas
          elseif(ALFCHR.eq.'L'.or.ALFCHR.eq.'l') then
              write(*,*) 'Lambda(nf=5) =?'
              read(*,*) alamb5
              ilflag= 0
          else
              write(*,*) '!!!  PLEASE TYPE: A OR L  !!!'
              goto 5
          endif
   10   write(*,896) qcut0
        read(*,895) QCTCHR
          if(QCTCHR.eq.'Y'.or.QCTCHR.eq.'y') then
              qcut=qcut0
          elseif(QCTCHR.eq.'N'.or.QCTCHR.eq.'n') then
              write(*,*) 'QCUT (GeV) = ?'
              read(*,*) qcut
          else
              write(*,*) '!!!   PLEASE TYPE: Y OR N   !!!'
              goto 10
          endif
   15   write(*,902) qmt1s
        read(*,895) QMTCHR
          if(QMTCHR.eq.'Y'.or.QMTCHR.eq.'y') then
              qmat1=qmt1s
          elseif(QMTCHR.eq.'N'.or.QMTCHR.eq.'n') then
              write(*,*) 'QMAT1 (GeV) = ?'
              read(*,*) qmat1
          else
              write(*,*) '!!!   PLEASE TYPE: Y OR N   !!!'
              goto 15
          endif
      endif
  895 format(1A)
  896 format(1x,'Long distance cut off for QCD potential'/
     $ 1x,'QCUT = ',f5.4,' GeV.  OK ? Y/N')
  902 format(1x,
     $ 'Matching QCD for NF=5 and Richardson for NF=3 at QMAT1 =',
     $  f5.2,' GeV.'/1x,'  OK ? Y/N')
      end
c
c
        real*8 function phfqcd(x)
c     integrand over k   ?
           real*8 pm, xkm, x, ALPHEF
           external ALPHEF
           common/pmaxkm/ pm, xkm
           phfqcd = ALPHEF((pm+xkm)*(dabs(pm-xkm)/(pm+xkm))**x)
        end
c
c
      FUNCTION ALEFVQ(x)
      implicit real*8(a-h,o-z)
      external ALPHEF
      common/xtr101/ p0
      data pi/3.1415926535897930d0/
      q= p0*x
      ALEFVQ= - 4d0/3* 4*pi*ALPHEF(q)
      return
      end
C
C
C
C
      COMPLEX*16 FUNCTION ZAPVGP(P,ETOT,VME,PCUT,ACC)
C
C     A(p,E)= ZAPVGP(P,ETOT,VME,PCUT,ACC)
C     for QCD potential VQQBAR(q) and GAMTPE(P,E)  - momentum
C     dependent width of top quark in t-tbar system.
C     2-dimensional integration
C     P - intrinsic momentum of t quark, ETOT - total energy of t-tbar,
C     VME=V0-E, where V0-potential at spatial infinity, E=ETOT-2*TMASS,
C     PCUT - cut off in momentum space; e.g. for QCD potential
C     given by ALPHEF  PCUT=QCUT in COMMON/parflg/,
C     ACC - accuracy
C     external functions: VQQBAR,GAMTPE,ADQUA,AD8GLE,ADGLG1,ADGLG2
C
      IMPLICIT REAL*8(A-Z)
      EXTERNAL FIN01P,FIN02P,FIN03P,FIN04P,AD8GLE,ADGLG1,ADGLG2
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      COMMON/XTR102/ P0,E0,VMEM,TM,ACC0
      DATA PI/3.14159265/,BUF/1D-10/,SMALL/1D-2/
C For Testing only
      small = 1.d-1
C
      CONST= -TMASS/(8*PI**2*P)
      TM= TMASS
      ACC0=ACC*SMALL
      P0=P
      E0=ETOT
      VMEM=VME*TMASS
      IF(PCUT.LE.P) THEN
         XXRE=AD8GLE(FIN01P,BUF,PCUT,ACC)+ADGLG1(FIN01P,PCUT,P,ACC)+
     $        ADGLG1(FIN02P,BUF,1/P,ACC)
         XXIM=AD8GLE(FIN03P,BUF,PCUT,ACC)+ADGLG1(FIN03P,PCUT,P,ACC)+
     $        ADGLG1(FIN04P,BUF,1/P,ACC)
      ELSE
         XXRE=ADGLG1(FIN01P,BUF,P,ACC)+ADGLG2(FIN01P,P,PCUT,ACC)+
     $        AD8GLE(FIN02P,BUF,1/PCUT,ACC)
         XXIM=ADGLG1(FIN03P,BUF,P,ACC)+ADGLG2(FIN03P,P,PCUT,ACC)+
     $        AD8GLE(FIN04P,BUF,1/PCUT,ACC)
      ENDIF
      ZAPVGP=CONST*CMPLX(XXRE,XXIM,KIND=KIND(0d0))
      END
C
      REAL*8 FUNCTION FIN01P(Q)
C     this segment contains FIN01P,FIN02P,FIN03P,FIN04P
      IMPLICIT REAL*8(A-C,D-H,O-Z)
      EXTERNAL VQQBAR,FIN11P, FIN12P
      COMMON/XTR102/ P0,E0,VMEM,TM,ACC0
      DATA PI/3.14159265/,BUF/1d-10/
      Q0=Q
      XL=(P0-Q0)**2
      XU=(P0+Q0)**2
      CALL ADQUA(XL,XU,FIN11P,Y,ACC0)
      FIN01P= VQQBAR(Q0)*Q0*Y
      RETURN
      ENTRY FIN02P(Q)
      Q0=1/Q
      XL=(P0-Q0)**2
      XU=(P0+Q0)**2
      CALL ADQUA(XL,XU,FIN11P,Y,ACC0)
      FIN02P= VQQBAR(Q0)*Q0**3*Y
      RETURN
      ENTRY FIN03P(Q)
      Q0=Q
      XL=(P0-Q0)**2
      XU=(P0+Q0)**2
      CALL ADQUA(XL,XU,FIN12P,Y,ACC0)
      FIN03P= VQQBAR(Q0)*Q0*Y
      RETURN
      ENTRY FIN04P(Q)
      Q0=1/Q
      XL=(P0-Q0)**2
      XU=(P0+Q0)**2
      CALL ADQUA(XL,XU,FIN12P,Y,ACC0)
      FIN04P= VQQBAR(Q0)*Q0**3*Y
      END
      REAL*8 FUNCTION FIN11P(T)
C     this segment contains FIN11P,FIN12P
      IMPLICIT REAL*8(A-C,D-H,O-Z)
      EXTERNAL GAMTPE
      COMMON/XTR102/ P0,E0,VMEM,TM,ACC0
      T1= T+VMEM
      TSQRT= SQRT(T)
      GAMMA= TM*GAMTPE(TSQRT,E0)
      FIN11P= T1/(T1**2+GAMMA**2)
      RETURN
      ENTRY FIN12P(T)
      T1= T+VMEM
      TSQRT= SQRT(T)
      GAMMA= TM*GAMTPE(TSQRT,E0)
      FIN12P= GAMMA/(T1**2+GAMMA**2)
      END
C
c
      SUBROUTINE VQDELT(VQ)
c
c     evaluates constants multiplying Dirac delta in potentials VQCUT
c     calls: ADQUA
c
      implicit real*8(a-h,o-z)
      external alphef,fncqct
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      COMMON/PARFLG/ QCUT,QMAT1,ALR,ILFLAG
      data pi/3.141592653589793238D0/
c
      call adqua(1d-8,1d4,fncqct,y,1d-4)
      v=-4d0/3*2/pi*y
      VQ=(-.25-v)*(2*pi)**3
      end
c
      function fncqct(q)
      implicit real*8(a-h,o-z)
      fncqct=sin(q)/q*alphef(q)
      end
c
C
      REAL*8 FUNCTION VQQBAR(P)
C
C     interquark potential for q- qbar singlet state
C
      IMPLICIT REAL*8(A-C,D-H,O-Z)
      EXTERNAL ALPHEF
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      DATA PI/3.14159265/
      VQQBAR = -4D0/3*4*PI*ALPHEF(P)/P**2
      END
C
      FUNCTION ALPHEF(q)
c
c     V(q) = -4/3 * 4*pi*ALPHEF(q)/q**2
c     input: alphas or alamb5 in COMMON/PHCONS/.  If:
c     ILFLAG.EQ.0   alamb5= \Lambda_\{\bar MS}^{(5)} at M_z
c     ILFLAG.EQ.1   alphas = alpha_{strong} at M_z (91.161)
c
c     effective coupling ALPHEF is defined as follows:
c     for q > qmat1=m_b:
c       alphas*( 1 +(31/3-10*nf/9)*alphas/(4*pi) )
c       where alphas=\alpha_\bar{MS} for nf=5, i.e.
c        alpha=4*pi/( b0(nf=5)*x + b1(5)/b0(5)*ln(x) )
c        and x = ln(q**2/alamb5**2)
c     for qmat1 > q > qcut:
c       4*pi/b0(nefr=3)*(alfmt+1/log(1+q**2/alr**2))
c       where alr=.4 GeV, nefr=3, and continuity --> alfmt
c     below qcut:  alphrc*2*q**2/(q**2+qcut**2)  (cont.-->alphrc)
c
      implicit real*8(a-h,o-z)
      SAVE
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      COMMON/PARFLG/ QCUT,QMAT1,ALR,ILFLAG
      common/parpot/ a5,b5,c5,alfmt,d,alphrc
      data pi/3.141592653589793238D0/,
     $ zold/-1d0/,qctold/-1d0/,alfold/-1d0/,
     $olmbd/-1d0/
c
      if(zmass.le.0d0 .or. qcut.le.0d0) STOP 10001
      if(zold.ne.zmass .or. qcut.ne.qctold) num=0
      if(ilflag.eq.0 .and. olmbd.ne.alamb5) num=0
      if(ilflag.eq.1 .and. alfold.ne.alphas) num=0
      if(num.eq.0)then
          num=num+1
          zold=zmass
          qctold=qcut
          call potpar
          alfold= alphas
          olmbd= alamb5
      endif
      if(q.le.qcut) then
         alphef=alphrc*(2*q**2)/(qcut**2+q**2)
      elseif(q.le.qmat1) then
         alphef=alfmt+d/log(1+q**2/alr**2)
      else
         x=2*log(q/alamb5)
         alfas5=1/(a5*x+b5*log(x))
         alphef=alfas5*(1+c5*alfas5)
      endif
      end
c
c Only called by ALPHEF:
      SUBROUTINE POTPAR
      implicit real*8(a-h,o-z)
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      COMMON/PARFLG/ QCUT,QMAT1,ALR,ILFLAG
      common/parpot/ a5,b5,c5,alfmt,d,alphrc
      data pi/3.141592653589793238D0/,nefr/3/
      b0(nf)=11-2./3*nf
      b1(nf)=102-38./3*nf
      cn(nf)=31./3-10./9*nf
      alr=400d-3
      a5=b0(5)/(4*pi)
      b5=b1(5)/b0(5)/(4*pi)
      c5=cn(5)/(4*pi)
      d=4*pi/b0(nefr)
      if(ilflag.eq.0) then
         if(alamb5.le.0d0) STOP 10002
         xa=2*log(zmass/alamb5)
         alphas= 1/(a5*xa + b5*log(xa))
      else
         if(alphas.le.0d0) STOP 10003
         t0=0
         t1=max(1d0,alphas*a5)
  10     tm=(t0+t1)/2
         fm=tm/alphas+b5*tm*log(tm)-a5
         if(fm.lt.-1d-10) then
           t0=tm
           goto 10
         elseif(fm.gt.1d-10) then
           t1=tm
           goto 10
         endif
         alamb5=zmass*exp(-5d-1/tm)
      endif
      x=2*log(qmat1/alamb5)
      alfas=1/(a5*x+b5*log(x))
      alfmt=alfas*(1+c5*alfas)-d/log(1+qmat1**2/alr**2)
      alphrc=alfmt+ d/log(1+qcut**2/alr**2)
      return
      end
c
c --- End of routines for phenomenological potentials ---
c
c
c --- Routines for Gamma_top ---
C
      SUBROUTINE GAMMAT
C
C     on shell width of top quark including QCD corrections, c.f.
C     M.Jezabek and J.H. Kuhn, Nucl. Phys. B314(1989)1
C
      IMPLICIT REAL*8(A-C,D-H,O-Z)
      EXTERNAL DILOGG
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      DATA PI/3.14159265/
      F(X)= PI**2+2*DILOGG(X)-2*DILOGG(1-X)+( 4*X*(1-X-2*X**2)*LOG(X)+
     $2*(1-X)**2*(5+4*X)*LOG(1-X) - (1-X)*(5+9*X-6*X**2) ) /
     $(2*(1-X)**2*(1+2*X))
      Y= (WMASS/TMASS)**2
cc alpha_s(M_t) corresponding to alpha_s(M_Z)=0.118:
cc      alphas=0.107443d0
cc      write(*,*) 'alphas=',alphas
c Usage of alpha_s as given as input for the potential.. better use
c alpha_s at a scale close to m_t..
      TGAMMA= GFERMI*TMASS**3/(8*SQRT(2D0)*PI)*(1-Y)**2*(1+2*Y)*
     $(1- 2D0/3*ALPHAS/PI*F(Y))
      END
C
C
      REAL*8 FUNCTION GAMTPE(P,ETOT)
C
C     momentum dependent width of top quark in t-tbar system
C     GAMTPE = TGAMMA*GTPCOR(P,E), where TGAMMA includes
C     QCD corrections, see JKT, eq.(8), and
C     GTPCOR - correction factor for bound t quark
C
      IMPLICIT REAL*8(A-C,D-H,O-Z)
      EXTERNAL GTPCOR
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      GAMTPE= TGAMMA*GTPCOR(P,ETOT)
      END
C
C
C    GTPCOR and GTPCOR1 should be merged  (M.J.) !!!!
c
        real*8 function gtpcor(topp,etot)
        real*8 topp,etot,
     u         tmass,tgamma,zmass,alphas,alamb5,
     u         wmass,wgamma,bmass,GFERMI,hmass
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
c       if (topp.ge.tmass/2.d0) then
c          gtpcor1=0.001d0
c       else
           gtpcor=1.d0
c       endif
        end
c
c
c Correction function for non-constant (energy and momentum dependent) width:
      FUNCTION GTPCOR1(TOPP,ETOT)
c
c     TOPP - momentum of t quark = - momentum of tbar
c     ETOT - total energy of t-tbar system
c     calls: GENWDS, RAN2
c
c     Evaluates a correction factor to the width of t-tbar system.
c     in future has to be replaced by a function evaluating
c     width including radiative corrections and GTPCOR.
c     I include two factors reducing the width:
c     a - time dilatation: for decay in flight lifetime
c         increased accordingly to relativistic kinematics
c     b - overall energy-momentum conservation: I assume that
c         t and tbar decay in flight and in this decays energies
c         of Ws follow from 2-body kinematics. Then I calculate
c         effective mass squared of b-bar system (it may be
c         negative!) from en-momentum conservation.
c         If effective mass is < 2*Mb + 2 GeV configuration
c         is rejected. The weight is acceptance.
c
      IMPLICIT REAL*8(A-H,O-Z)
      real ran2
      external ran2
      PARAMETER(NG=20,NC=4)
      dimension gamma(0:NG),pw1(0:3),pw2(0:3),AIJ(NC,NC),BJ(NC),
     $AI(NC),SIG2IN(0:NG),XIK(0:NG,NC),INDX(NC)
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      SAVE NUM,EOLD,TOLD,AI
      data nevent/10000/, num/0/, eold/-1d5/, told/-1d0/
c
C     for test runs!!
C      nevent=1000
C
      if(etot.ne.eold) num=0
      if(tmass.ne.told) num=0
   5  if(num.eq.0) then
c            xdumm= ran2(-2)
            do 10 itp=0,NG
            tp=itp*tmass/NG*2
            gamma(itp)=0
            do 10 ix=1,nevent
            call GENWDS(tp,etot,pw1,pw2,efmsq)
            if(efmsq.gt.0d0) then
               efms=sqrt(efmsq)
               if(efms.ge. 2*bmass+2) gamma(itp)=gamma(itp)+1
            endif
  10        continue
            do 15 ix=0,NG
  15        SIG2IN(IX)= MAX(1D0,GAMMA(IX))
            DO 17 JX=1,NC
              IF(JX.EQ.1)THEN
                XIK(0,JX)= .5D0
              ELSE
                XIK(0,JX)= 0D0
              ENDIF
            DO 17 IX=1,NG
            tp= 2D0*ix/NG
  17        XIK(IX,JX)= tp**(JX-1)/(1+EXP(tp*3))
            DO 20 I=1,NC
            BJ(I)=0
            DO 20 J=1,NC
  20        AIJ(I,J)=0
            DO 30 I=1,NC
            DO 25 IX=0,NG
  25        BJ(I)= BJ(I)+GAMMA(IX)*XIK(IX,I)*SIG2IN(IX)
            DO 30 J=1,I
            DO 30 IX=0,NG
  30        AIJ(I,J)= AIJ(I,J)+XIK(IX,I)*XIK(IX,J)*SIG2IN(IX)
            DO 35 I=1,NC
            DO 35 J=I,NC
  35        AIJ(I,J)= AIJ(J,I)
            CALL LUDCMP(AIJ,NC,NC,INDX,D)
            CALL LUBKSB(AIJ,NC,NC,INDX,BJ)
            DO 40 I=1,NC
  40        AI(I)= BJ(I)/NEVENT
      do 42 i=1,nc
  42  write(*,*)'a(',i,')=',ai(i)
            do 100 ix=0,NG
  100       gamma(ix)= gamma(ix)/nevent
            eold=etot
            told=tmass
            num= 1
      endif
      SUM=AI(1)
      DO 110 I=2,NC
  110 SUM= SUM+AI(I)*(TOPP/TMASS)**(I-1)
C      CORRF2= SUM/(1+ EXP(TOPP/TMASS*3))
      CORRF2= SUM/(1+ EXP(MIN(1d1,TOPP/TMASS*3)))
C      if(topp.gt. 2d0*tmass) then
C          corrf1= 0.001d0
C      else
C           ip= NG*topp/tmass/2
C          corrf1= gamma(ip)
C      endif
C      write(*,*)'ratio=',corrf1/corrf2
C      GTPCOR1 = CORRF2
      GTPCOR1 = CORRF2*SQRT(1-TOPP**2/(TOPP**2+TMASS**2))
      END
c
c Generator: only called by GTPCOR1
      SUBROUTINE GENWDS(tp,etot,pw1,pw2,efm2)
c
c     generates 4-momenta of W's  and effective mass of b-bbar
c     from t and tbar quarks decays at flight (tp = momentum of t
c     = - momentum of tbar (in GeV) ) in Oz direction
c
      implicit real*8(a-h,o-z)
c      real ran2
      real ranf
c      external ran2
      external ranf
      dimension pw1(0:3),pw2(0:3)
      save
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
      data PI/3.141592653589793238D0/
      real idum
c  3   s1= wmass**2+wmass*wgamma*TAN((2*ran2(idum)-1)*pi/2)
  3   s1= wmass**2+wmass*wgamma*TAN((2*ranf(idum)-1)*pi/2)
      if(s1.le.0d0) goto 3
      wmass1= sqrt(s1)
      if(abs(wmass1-wmass).ge.3*wgamma) goto 3
c  4   s2= wmass**2+wmass*wgamma*TAN((2*ran2(idum)-1)*pi/2)
  4   s2= wmass**2+wmass*wgamma*TAN((2*ranf(idum)-1)*pi/2)
      if(s2.le.0d0) goto 4
      wmass2= sqrt(s2)
      if(abs(wmass2-wmass).ge.3*wgamma) goto 4
      ew1= (tmass**2+wmass1**2-bmass**2)/(2*tmass)
      pwt1= sqrt(ew1**2-wmass1**2)
      ew2= (tmass**2+wmass2**2-bmass**2)/(2*tmass)
      pwt2= sqrt(ew2**2-wmass2**2)
  5   p=tp
c      u1= 2*ran2(idum)-1
      u1= 2*ranf(idum)-1
      pw1z= pwt1*u1
c      u2= 2*ran2(idum)-1
      u2= 2*ranf(idum)-1
      pw2z= pwt2*u2
      et= sqrt(tmass**2+p**2)
      bet= p/et
      gam= et/tmass
      pw1(0)= gam*(ew1+bet*pw1z)
      pw1(3)= gam*(pw1z+bet*ew1)
      pw2(0)= gam*(ew2-bet*pw2z)
      pw2(3)= gam*(pw2z-bet*ew2)
      pw1tr= sqrt(pw1(0)**2-pw1(3)**2-wmass1**2)
      pw2tr= sqrt(pw2(0)**2-pw2(3)**2-wmass2**2)
c      phi1= 2*pi*ran2(idum)
      phi1= 2*pi*ranf(idum)
c      phi2= 2*pi*ran2(idum)
      phi2= 2*pi*ranf(idum)
      pw1(1)= pw1tr*cos(phi1)
      pw1(2)= pw1tr*sin(phi1)
      pw2(1)= pw2tr*cos(phi2)
      pw2(2)= pw2tr*sin(phi2)
      prec2= (pw1(1)+pw2(1))**2+(pw1(2)+pw2(2))**2+(pw1(3)+pw2(3))**2
      erest=etot-pw1(0)-pw2(0)
c
      efm2= erest*abs(erest)-prec2
      END
c
c --- End of routines for Gamma_top ---
c
c --- Routines for solving linear equations and matrix inversion (complex) ---
c
        subroutine sae(pp, w1, bb, a1, n)
c
           implicit none
           complex*16 vhat
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        d, pp, w1, gtpcor,hmass,
     u        xp,xpmax,dcut,kincom,kincoa,kincov
           complex*16 a, a1, bb, ff, cw, svw, g0, g0c
           integer i, j, npot, n, nmax, indx,kinflg,gcflg,vflag
           parameter (nmax=900)
           dimension bb(nmax), ff(nmax,nmax), pp(nmax), w1(nmax),
     u               indx(nmax), cw(nmax), a1(nmax)
c
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           common/mom/ xp,xpmax,dcut
           common/g0inf/kincom,kincoa,kincov,kinflg,gcflg,vflag
c
           external a, vhat, gtpcor, g0, g0c
c
           do 10 i=1,n*2/3
              cw(i) = w1(i) / (4.d0*pi**2) * g0c(pp(i))
c             cw(i) = w1(i) / (4.d0*pi**2 *
c     u                (cmplx(energy-vzero, tgamma*
c     u                 gtpcor(pp(i),2.d0*tmass+energy),
c     u                    kind=kind(0d0))-pp(i)**2/tmass))
10         continue
           do 20 i=n*2/3+1,n
              cw(i) = w1(i) / (4.d0*pi**2) * g0c(pp(i)) * pp(i)**2
c             cw(i) = w1(i) / (4.d0*pi**2 *
c     u          (cmplx(energy-vzero, tgamma*
c     u           gtpcor(pp(i),2.d0*tmass+energy),
c     u                    kind=kind(0d0)) /
c     u           pp(i)**2 - 1.d0/tmass))
20         continue
c
           do 30 i=1,n
cc            bb(i) = a1(i)
cvv
              if (pp(i).lt.dcut.and.vflag.eq.1) then
c                bb(i) = cmplx(1.d0+kincov*(pp(i)/tmass)**2,0.d0,
c     u                    kind=kind(0d0))
                 bb(i)=1.d0+kincov*
     u                       g0(pp(i))*(pp(i)/tmass)**2/g0c(pp(i))
              else
                 bb(i) = (1.d0,0.d0)
              endif
              svw = (0.d0,0.d0)
              do 40 j=1,n
                 if (i.ne.j) then
                    ff(i,j) = - vhat(pp(i),pp(j)) * cw(j)
                    svw = svw + ff(i,j)
                 endif
40            continue
              ff(i,i) = 1.d0 - a1(i) - svw
30         continue
c
           call zldcmp(ff, n, nmax, indx, d)
           call zlbksb(ff, n, nmax, indx, bb)
c
        end
c
c
      SUBROUTINE ZLBKSB(A,N,NP,INDX,B)
C complex version of lubksb
      IMPLICIT NONE
      INTEGER I, II, INDX, J, LL, N, NP
      COMPLEX*16 A, B, SUM
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.(0.D0,0.D0)) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
c
      SUBROUTINE ZLDCMP(A,N,NP,INDX,D)
C complex version of ludcmp
      IMPLICIT NONE
      INTEGER I, IMAX, INDX, J, K, N, NP, NMAX
      REAL*8 AAMAX, D, TINY, VV
      COMPLEX*16 A, DUM, SUM
      PARAMETER (NMAX=900)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)
c
        tiny=1.d-5
c
      D=1.D0
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
c        IF (AAMAX.EQ.0.D0) PAUSE 'Singular matrix.'
        IF (AAMAX.EQ.0.D0) print *, "Singular matrix."
        VV(I)=1.D0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (ABS(DUM).GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX) THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF (J.NE.N) THEN
          IF (A(J,J).EQ.(0.D0,0.D0)) A(J,J)=cmplx(TINY, 0.d0,
     u                    kind=kind(0d0))
          DUM=1.D0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.(0.D0,0.D0)) A(N,N)=cmplx(TINY, 0.d0,
     u                    kind=kind(0d0))
      RETURN
      END
C
C
C *** TOOLS ***
C
C
C     ******* ROUTINES FOR GAUSSIAN INTEGRATIONS
C
C
      SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C     Given the lower and upper limits of integration X1 and X2
C     and given N, this routine returns arrays X(N) and W(N)
C     containing the abscissas and weights of the Gauss-Legendre
C     N-point quadrature formula
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      save
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=DCOS(3.141592653589793238D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(DABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION AD8GLE(F,A,B,EPS)
      implicit double precision (a-h,o-z)
      EXTERNAL F
      DIMENSION W(12),X(12)
c      SAVE W, X
      SAVE
C
C     ******************************************************************
C
C     ADAPTIVE GAUSSIAN QUADRATURE.
C
C     AD8GLE IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
C     EPS.
C
C     ******************************************************************
C
      DATA W / 0.10122 85362 90376 25915 25313 543D0,
     $         0.22238 10344 53374 47054 43559 944D0,
     $         0.31370 66458 77887 28733 79622 020D0,
     $         0.36268 37833 78361 98296 51504 493D0,
     $         0.27152 45941 17540 94851 78057 246D-1,
     $         0.62253 52393 86478 92862 84383 699D-1,
     $         0.95158 51168 24927 84809 92510 760D-1,
     $         0.12462 89712 55533 87205 24762 822D0,
     $         0.14959 59888 16576 73208 15017 305D0,
     $         0.16915 65193 95002 53818 93120 790D0,
     $         0.18260 34150 44923 58886 67636 680D0,
     $         0.18945 06104 55068 49628 53967 232D0/
C
      DATA X / 0.96028 98564 97536 23168 35608 686D0,
     $         0.79666 64774 13626 73959 15539 365D0,
     $         0.52553 24099 16328 98581 77390 492D0,
     $         0.18343 46424 95649 80493 94761 424D0,
     $         0.98940 09349 91649 93259 61541 735D0,
     $         0.94457 50230 73232 57607 79884 155D0,
     $         0.86563 12023 87831 74388 04678 977D0,
     $         0.75540 44083 55003 03389 51011 948D0,
     $         0.61787 62444 02643 74844 66717 640D0,
     $         0.45801 67776 57227 38634 24194 430D0,
     $         0.28160 35507 79258 91323 04605 015D0,
     $         0.95012 50983 76374 40185 31933 543D-1/
C
C     ******************************************************************
C
      GAUSS=0.0D0
      AD8GLE=GAUSS
      IF(B.EQ.A) RETURN
      CONST=EPS/(B-A)
      BB=A
C
C  COMPUTATIONAL LOOP.
    1 AA=BB
      BB=B
    2    C1=0.5D0*(BB+AA)
         C2=0.5D0*(BB-AA)
         S8=0.0D0
         DO 3 I=1,4
            U=C2*X(I)
            S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3    CONTINUE
         S8=C2*S8
         S16=0.0D0
         DO 4 I=5,12
            U=C2*X(I)
            S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4    CONTINUE
         S16=C2*S16
         IF( ABS(S16-S8) .LE. EPS*(abs(s8)+ABS(S16))*0.5D0 ) GO TO 5
         BB=C1
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2
      AD8GLE=0.0D0
      write(*,*)'too high accuracy required in function ad8gle!'
      RETURN
    5 GAUSS=GAUSS+S16
      IF(BB.NE.B) GO TO 1
      AD8GLE=GAUSS
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION ADGLG1(F,A,B,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F,AD8GLE,adqua
      DIMENSION W(6),X(6),xx(6)
c      SAVE W, XX, NUM
      SAVE
C
C     ******************************************************************
C
C     ADAPTIVE GAUSSIAN QUADRATURE.
C     For x->b   f(x) = O (ln^k (b-x) )
C     A - lower limit, B - upper limit (integrable singularity)
C     AD8GLE IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
C     EPS.
C
C     ******************************************************************
      DATA W / 4.58964 673950d-1,
     $         4.17000 830772d-1,
     $         1.13373 382074d-1,
     $         1.03991 974531d-2,
     $         2.61017 202815d-4,
     $         8.98547 906430d-7/
C
      DATA X / 0.22284 66041 79d0,
     $         1.18893 21016 73d0,
     $         2.99273 63260 59d0,
     $         5.77514 35691 05d0,
     $         9.83746 74183 83d0,
     $        15.98287 39806 02d0/
      DATA NUM/0/
      IF(NUM.eq.0d0) then
      do 1 ix=1,6
  1   xx(ix)= EXP(-x(ix))
      ENDIF
      num=num+1
      sum=0d0
      c=b-a
      sum6=0d0
      do 10 in=1,6
 10   sum6= sum6+ w(in)*f(b-c*xx(in))
      sum6=sum6*c
      a1=a
 15   a2= (a1+b)/2
      c=b-a2
      sumn=0d0
      do 20 in=1,6
      !!! FB: catch NaN
      if ( c/b .lt. 1d-9 ) then
        adglg1 = 1d15
        return
      endif
 20   sumn= sumn+ w(in)*f(b-c*xx(in)) !!! FB: f(b) = NaN !
      sumn=sumn*c
ctt
c      call adqua(a1,a2,f,sum1,eps)
c      sum1=sum1+sum
      sum1=AD8GLE(F,A1,A2,eps)+sum
      IF(ABS( (sum+sum6)/(sum1+sumn)-1d0 ).lt.EPS) THEN
ctt
c      call adqua(a,a2,f,sum2,eps)
         sum2=AD8GLE(F,A,A2,eps)
         IF(ABS( (sum2+sumn)/(sum1+sumn)-1d0 ).gt.EPS) THEN
            sum=sum2
            a1=a2
            sum6=sumn
            goto 15
         ENDIF
         ADGLG1= SUM1+SUMN
         RETURN
      ELSE
         sum=sum1
         a1=a2
         sum6=sumn
         goto 15
      ENDIF
      END
C
      DOUBLE PRECISION FUNCTION ADGLG2(F,A,B,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F,AD8GLE
      DIMENSION W(6),X(6),xx(6)
c      SAVE W,XX,NUM
      SAVE
C
C     ******************************************************************
C
C     ADAPTIVE GAUSSIAN QUADRATURE.
C     For x->A   f(x) = O (ln^k (x-a) )
C     A - lower limit  (integrable singularity), B - upper limit
C     AD8GLE IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
C     EPS.
C
C     ******************************************************************
      DATA W / 4.58964 673950d-1,
     $         4.17000 830772d-1,
     $         1.13373 382074d-1,
     $         1.03991 974531d-2,
     $         2.61017 202815d-4,
     $         8.98547 906430d-7/
C
      DATA X / 0.22284 66041 79d0,
     $         1.18893 21016 73d0,
     $         2.99273 63260 59d0,
     $         5.77514 35691 05d0,
     $         9.83746 74183 83d0,
     $        15.98287 39806 02d0/
      DATA NUM/0/
      IF(NUM.eq.0d0) then
      do 1 ix=1,6
  1   xx(ix)= EXP(-x(ix))
      ENDIF
      num=num+1
      sum=0d0
      c=b-a
      sum6=0d0
      do 10 in=1,6
 10   sum6= sum6+ w(in)*f(A+c*xx(in))
      sum6=sum6*c
      b1=b
 15   b2= (a+b1)/2
      c=b2-a
      sumn=0d0
      do 20 in=1,6
      !!! FB: catch NaN
      if ( c/a .lt. 1d-9 ) then
        adglg2 = 1d15
        return
      endif
 20   sumn= sumn+ w(in)*f(a+c*xx(in)) !!! FB: f(a) = NaN !
      sumn=sumn*c
      sum1=AD8GLE(F,b2,b1,eps)+sum
      IF(ABS( (sum+sum6)/(sum1+sumn)-1d0 ).lt.EPS) THEN
         sum2=AD8GLE(F,b2,b,eps)
         IF(ABS( (sum2+sumn)/(sum1+sumn)-1d0 ).gt.EPS) THEN
            sum=sum2
            b1=b2
            sum6=sumn
            goto 15
         ENDIF
         ADGLG2= SUM1+SUMN
         RETURN
      ELSE
         sum=sum1
         b1=b2
         sum6=sumn
         goto 15
      ENDIF
      END
C
C
C------------------------------------------------------------------
C INTEGRATION ROUTINE ADQUA written by M. Jezabek            ------
C------------------------------------------------------------------
C
      SUBROUTINE ADQUA(XL,XU,F,Y,ACC)
C
C     ADAPTIVE GAUSS-LEGENDRE + SIMPSON'S RULE QUADRATURE
C     XL - LOWER LIMIT, XU - UPPER LIMIT, F - FUNCTION TO INTEGRATE
C     Y - INTEGRAL
C     ACC - ACCURACY (IF .LE. 0.  ACC=1.D-6)
c     ****** new constants,  1 error removed, Oct '92
C
C     CALLS: SIMPSA
C
C     PARAMETERS: NSUB > NO OF SUBDIVISION LEVELS IN GAUSS INTEGRATION
C          100*2**IMAX > NO OF POINTS IN SIMPSON INTEGRATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      DIMENSION VAL(25,2), BOUND(25,2,2), LEV(25),SING(25,3)
      DIMENSION W8(4),X8(4)
      DATA W8
     $/0.101228536290376D0, 0.222381034453374D0, 0.313706645877887D0,
     $ 0.362683783378362D0/
      DATA X8
     $/0.960289856497536D0, 0.796666477413627D0, 0.525532409916329D0,
     $ 0.183434642495650D0/
      save
C
      IF(ACC.LE.0.D0) ACC=1.D-6
      NSUB=24
      NSG=25
      NSC=0
      A=XL
      B=XU
      C1=0.5d0*(A+B)
      C2=C1-A
      S8=0d0
      DO 1 I=1,4
      U=X8(I)*C2
    1 S8=S8+W8(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      XM=(XL+XU)/2.d0
      BOUND(1,1,1)=XL
      BOUND(1,1,2)=XM
      BOUND(1,2,1)=XM
      BOUND(1,2,2)=XU
      NC=1
      DO 3 IX=1,2
      A=BOUND(NC,IX,1)
      B=BOUND(NC,IX,2)
      C1=0.5d0*(A+B)
      C2=C1-A
      VAL(NC,IX)=0.d0
      DO 2 I=1,4
      U=X8(I)*C2
    2 VAL(NC,IX)=VAL(NC,IX)+W8(I)*(F(C1+U)+F(C1-U))
    3 VAL(NC,IX)=VAL(NC,IX)*C2
      S16=VAL(NC,1)+VAL(NC,2)
      IF(DABS(S8-S16).GT.ACC*DABS(S16)) GOTO 4
      Y=S16
      RETURN
    4 DO 5 I=1,NSUB
    5 LEV(I)=0
      NC1= NC+1
   11 XM=(BOUND(NC,1,1)+BOUND(NC,1,2))/2.d0
      BOUND(NC1,1,1)=BOUND(NC,1,1)
      BOUND(NC1,1,2)=XM
      BOUND(NC1,2,1)=XM
      BOUND(NC1,2,2)=BOUND(NC,1,2)
      DO 13 IX=1,2
      A=BOUND(NC1,IX,1)
      B=BOUND(NC1,IX,2)
      C1=0.5d0*(A+B)
      C2=C1-A
      VAL(NC1,IX)=0.d0
      DO 12 I=1,4
      U=X8(I)*C2
   12 VAL(NC1,IX)=VAL(NC1,IX)+W8(I)*(F(C1+U)+F(C1-U))
   13 VAL(NC1,IX)=VAL(NC1,IX)*C2
      S16=VAL(NC1,1)+VAL(NC1,2)
      S8=VAL(NC,1)
      IF(DABS(S8-S16).LE.ACC*DABS(S16)) GOTO 20
      NC=NC1
      NC1= NC+1
      IF(NC1.LE.NSUB) GOTO 11
C     NC=NSUB   USE SIMPSON'S RULE
      NSC=NSC+1
      IF(NSC.LE.NSG) GOTO 15
      WRITE(*,911)
  911 FORMAT(1X,'ADQUA: TOO MANY SINGULARITIES')
      STOP
   15 SING(NSC,1)=BOUND(NC,1,1)
      SING(NSC,2)=BOUND(NC,2,2)
      SING(NSC,3)=S16
      S16=0.d0
      NC=NC-1
   20 VAL(NC,1)= S16
  121 LEV(NC)=1
   21 XM=(BOUND(NC,2,1)+BOUND(NC,2,2))/2.d0
      BOUND(NC1,1,1)=BOUND(NC,2,1)
      BOUND(NC1,1,2)=XM
      BOUND(NC1,2,1)=XM
      BOUND(NC1,2,2)=BOUND(NC,2,2)
      DO 23 IX=1,2
      A=BOUND(NC1,IX,1)
      B=BOUND(NC1,IX,2)
      C1=0.5d0*(A+B)
      C2=C1-A
      VAL(NC1,IX)=0.d0
      DO 22 I=1,4
      U=X8(I)*C2
   22 VAL(NC1,IX)=VAL(NC1,IX)+W8(I)*(F(C1+U)+F(C1-U))
   23 VAL(NC1,IX)=VAL(NC1,IX)*C2
      S16=VAL(NC1,1)+VAL(NC1,2)
      S8=VAL(NC,2)
      IF(DABS(S8-S16).LE.ACC*DABS(S16)) GOTO 40
      NC=NC+1
      NC1=NC+1
      IF(NC1.LE.NSUB) GOTO 11
C     NC=NSUB   USE SIMPSON'S RULE
      NSC=NSC+1
      IF(NSC.LE.NSG) GOTO 35
      WRITE(*,911)
      STOP
   35 SING(NSC,1)=BOUND(NC,1,1)
      SING(NSC,2)=BOUND(NC,2,2)
      SING(NSC,3)=S16
      S16=0.d0
      NC=NC-1
   40 VAL(NC,2)= S16
   45 IF(NC.GT.1) GOTO 50
      Y1=VAL(1,1)+VAL(1,2)
      GOTO 100
   50 NC0=NC-1
      IF(LEV(NC0).EQ.0) IX=1
      IF(LEV(NC0).EQ.1) IX=2
      LEV(NC)=0
      NC1=NC
      VAL(NC0,IX)=VAL(NC,1)+VAL(NC,2)
      NC=NC0
      IF(IX.EQ.1) GOTO 121
      GOTO 45
  100 CONTINUE
      IF(NSC.GT.0) GOTO 101
      Y=Y1
      RETURN
  101 FSUM=0.d0
      DO 102 IK=1,NSC
  102 FSUM=FSUM+DABS(SING(IK,3))
      ACCR=ACC*DMAX1(FSUM,DABS(Y1))/FSUM/10.d0
      DO 104 IK=1,NSC
  104 CALL SIMPSA(SING(IK,1),SING(IK,2),F,SING(IK,3),ACCR)
      DO 106 IK=1,NSC
  106 Y1=Y1+SING(IK,3)
      Y=Y1
      RETURN
      END
C
      SUBROUTINE SIMPSA(A,B,F,F0,ACC)
C     SIMPSON'S ADAPTIVE QUADRATURE
      IMPLICIT REAL*8 (A-H,O-Z)
      save
      EXTERNAL F
      IMAX=5
      N0=100
      H=(B-A)/N0
      N02=N0/2
      S2=0.d0
      IC=1
      S0=F(A)+F(B)
      DO 5 K=1,N02
    5 S2=S2+F(A+2.d0*K*H)
    7 S1=0.d0
      DO 10 K=1,N02
   10 S1=S1+F(A+(2.d0*K-1.d0)*H)
      Y=H/3.d0*(S0+4.d0*S1+2.d0*S2)
      IF(DABS(F0/Y-1.d0).GT.ACC) GOTO 20
      RETURN
   20 N02=N0
      N0=2*N0
      S2=S1+S2
      H=H/2.d0
      IF(IC.GT.IMAX) GOTO 30
      F0=Y
      IC=IC+1
      GOTO 7
   30 ACC0=DABS(Y/F0-1.d0)
      WRITE(*,900) A,B,ACC0
      STOP
  900 FORMAT(1H ,'SIMPSA: TOO HIGH ACCURACY REQUIRED'/
     /1X,   29HSINGULARITY IN THE INTERVAL  ,D20.12,1X,D20.12/
     /1X,   29HACCURACY ACHIEVED            ,D20.12)
      END
C
C
C  ******* matrix-inversion-routines
C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
!        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        IF (AAMAX.EQ.0.) print *, 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
c
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
C
C
C     *******  RANDOM NUMBER GENERATORS
C
C
      FUNCTION RANF(DUMMY)
C
C   RANDOM NUMBER FUNCTION TAKEN FROM KNUTH
C   (SEMINUMERICAL ALGORITHMS).
C   METHOD IS X(N)=MOD(X(N-55)-X(N-24),1/FMODUL)
C   NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.
C
C   RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.
C   IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.
C   IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE.
C
      PARAMETER (FMODUL=1.E-09)
      SAVE /CIRN55/
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
      CALL RANDAT
      IF( NCALL.EQ.0 ) THEN
          CALL IN55 ( IA,234612947 )
          MCALL = 55
          NCALL = 1
      ENDIF
      IF ( MCALL.EQ.0 ) THEN
          CALL IRN55(IA)
          MCALL=55
      ENDIF
      RANF=IA(MCALL)*FMODUL
      MCALL=MCALL-1
      RETURN
      END
C
      SUBROUTINE RANDAT
C
C  INITIALISES THE NUMBER NCALL TO 0 TO FLAG THE FIRST CALL
C  OF THE RANDOM NUMBER GENERATOR
C
C      SAVE /CIRN55/
C      SAVE FIRST
      SAVE
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      IF(FIRST)THEN
         FIRST=.FALSE.
         NCALL=0
      ENDIF
      RETURN
      END
C
      SUBROUTINE IN55(IA,IX)
      PARAMETER (MODULO=1000000000)
      INTEGER IA(55)
C
      IA(55)=IX
      J=IX
      K=1
      DO 10 I=1,54
         II=MOD(21*I,55)
         IA(II)=K
         K=J-K
         IF(K.LT.0)K=K+MODULO
         J=IA(II)
10    CONTINUE
      DO 20 I=1,10
         CALL IRN55(IA)
20    CONTINUE
      RETURN
      END
C
      SUBROUTINE IRN55(IA)
      PARAMETER (MODULO=1000000000)
      INTEGER IA(55)
      DO 10 I=1,24
         J=IA(I)-IA(I+31)
         IF(J.LT.0)J=J+MODULO
         IA(I)=J
10    CONTINUE
      DO 20 I=25,55
         J=IA(I)-IA(I-24)
         IF(J.LT.0)J=J+MODULO
         IA(I)=J
20    CONTINUE
      RETURN
      END
C
C
      FUNCTION RAN2(IDUM)
C     *******************
      REAL RDM(31)
      DATA  IWARM/0/
C
      IF (IDUM.LT.0.OR.IWARM.EQ.0) THEN
C INITIALIZATION OR REINITIALISATION
      IWARM=1
      IA1=         1279
      IC1=       351762
      M1=       1664557
      IA2=         2011
      IC2=       221592
      M2=       1048583
      IA3=        15091
      IC3=         6171
      M3=         29201
      IX1=MOD(-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      RM1=1./FLOAT(M1)
      RM2=1./FLOAT(M2)
      DO 10 J=1,31
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
10    RDM(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      ENDIF
C
C GENERATE NEXT NUMBER IN SEQUENCE
      IF(IWARM.EQ.0) GOTO 901
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(31*IX3)/M3
      RAN2=RDM(J)
      RDM(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
901   PRINT 9010
9010  FORMAT('   RAN2: LACK OF ITINIALISATION')
      STOP
      END
C
C
C     *******    SPECIAL FUNCTIONS
C
C
      DOUBLE PRECISION FUNCTION DILOGG(X)
C
C     SPENCE'S DILOGARITHM IN DOUBLE PRECISION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      Z=-1.644934066848226
      IF(X .LT.-1.0) GO TO 1
      IF(X .LE. 0.5) GO TO 2
      IF(X .EQ. 1.0) GO TO 3
      IF(X .LE. 2.0) GO TO 4
      Z=3.289868133696453
    1 T=1.0/X
      S=-0.5
      Z=Z-0.5*DLOG(DABS(X))**2
      GO TO 5
    2 T=X
      S=0.5
      Z=0.
      GO TO 5
    3 DILOGG=1.644934066848226
      RETURN
    4 T=1.0-X
      S=-0.5
      Z=1.644934066848226-DLOG(X)*DLOG(DABS(T))
    5 Y=2.666666666666667*T+0.666666666666667
      B=      0.00000 00000 00001
      A=Y*B  +0.00000 00000 00004
      B=Y*A-B+0.00000 00000 00011
      A=Y*B-A+0.00000 00000 00037
      B=Y*A-B+0.00000 00000 00121
      A=Y*B-A+0.00000 00000 00398
      B=Y*A-B+0.00000 00000 01312
      A=Y*B-A+0.00000 00000 04342
      B=Y*A-B+0.00000 00000 14437
      A=Y*B-A+0.00000 00000 48274
      B=Y*A-B+0.00000 00001 62421
      A=Y*B-A+0.00000 00005 50291
      B=Y*A-B+0.00000 00018 79117
      A=Y*B-A+0.00000 00064 74338
      B=Y*A-B+0.00000 00225 36705
      A=Y*B-A+0.00000 00793 87055
      B=Y*A-B+0.00000 02835 75385
      A=Y*B-A+0.00000 10299 04264
      B=Y*A-B+0.00000 38163 29463
      A=Y*B-A+0.00001 44963 00557
      B=Y*A-B+0.00005 68178 22718
      A=Y*B-A+0.00023 20021 96094
      B=Y*A-B+0.00100 16274 96164
      A=Y*B-A+0.00468 63619 59447
      B=Y*A-B+0.02487 93229 24228
      A=Y*B-A+0.16607 30329 27855
      A=Y*A-B+1.93506 43008 69969
      DILOGG=S*T*(A-B)+Z
      RETURN
      END
c
      SUBROUTINE pzext0(iest,xest,yest,yz,dy,nv)
      implicit none
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
c
c
        complex*16 function zdigamma(z)
        implicit none
        complex*16 z,psi,psipr1,psipr2
        call mkpsi(z,psi,psipr1,psipr2)
        zdigamma=psi
        end
c
      subroutine mkpsi(z,psi,psipr1,psipr2)
      implicit none
      complex*16 tmp,tmps2,tmps3,tmp0,tmp1,tmp2,ser0,ser1,ser2,ser3,
     .           zz,z,psi,psipr1,psipr2,off0,off1,off2,zcf,ser02,ser12,
     .           z1,z2
      real*8 cof(6),re1
      integer i
      data cof/76.18009173d0,-86.50532033d0,24.01409822d0,
     .    -1.231739516d0,.120858003d-2,-.536382d-5/
      save
      zz=z
      off0=cmplx(0.d0,0.d0,kind=kind(0d0))
      off1=cmplx(0.d0,0.d0,kind=kind(0d0))
      off2=cmplx(0.d0,0.d0,kind=kind(0d0))
    5 re1=real(zz)
      if (re1.le.0.d0) then
         off0=off0+1.d0/zz
         z1=zz*zz
         off1=off1-1.d0/z1
         z2=z1*zz
         off2=off2+2.d0/z2
         zz=zz+(1.d0,0.d0)
         goto 5
      endif
      tmp=zz+cmplx(4.5d0,0.d0,kind=kind(0d0))
      tmps2=tmp*tmp
      tmps3=tmp*tmps2
      tmp0=(zz-cmplx(0.5d0,0.d0,kind=kind(0d0)))/tmp+log(tmp)
     u     -cmplx(1.d0,0.d0,kind=kind(0d0))
      tmp1=(5.d0,0.d0)/tmps2+1.d0/tmp
      tmp2=(-10.0d0,0.d0)/tmps3-1.d0/tmps2
      ser0=cmplx(1.d0,0.d0,kind=kind(0d0))
      ser1=cmplx(0.d0,0.d0,kind=kind(0d0))
      ser2=cmplx(0.d0,0.d0,kind=kind(0d0))
      ser3=cmplx(0.d0,0.d0,kind=kind(0d0))
      do 10 i=1,6
         zcf=cof(i)/zz
         ser0=ser0+zcf
         zcf=zcf/zz
         ser1=ser1+zcf
         zcf=zcf/zz
         ser2=ser2+zcf
         zcf=zcf/zz
         ser3=ser3+zcf
         zz=zz+(1.d0,0.d0)
   10 continue
      ser1=-ser1
      ser2=2.d0*ser2
      ser3=-6.d0*ser3
      ser02=ser0*ser0
      ser12=ser1*ser1
      psi=tmp0+ser1/ser0-off0
      psipr1=tmp1+(ser2*ser0-ser12)/ser02-off1
      psipr2=tmp2+(ser3*ser02-3.d0*ser2*ser1*ser0+2.d0*ser12*ser1)
     .            /ser02/ser0-off2
      return
      end
