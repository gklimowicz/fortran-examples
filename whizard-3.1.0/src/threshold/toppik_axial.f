! WHIZARD 3.1.0 Dec 14 2022

! TOPPIK code by M. Jezabek, T. Teubner (v1.1, 1992), T. Teubner (1998)
!
! NOTE: axial part (p-wave) only
!
! FB: -commented out numerical recipes code for hypergeometric 2F1
!      included in hypgeo.f90;
!     -replaced function 'cdabs' by 'abs';
!     -replaced function 'dabs' by 'abs';
!     -replaced function 'dimag' by 'aimag';
!     -replaced function 'dcmplx(,)' by 'cmplx(,,kind=kind(0d0))';
!     -replaced function 'dreal' by 'real';
!     -replaced function 'dlog' by 'log';
!     -replaced function 'dsqrt' by 'sqrt';
!     -renamed function 'a' to 'aax'
!     -renamed function 'fretil1' to 'fretil1ax'
!     -renamed function 'fretil2' to 'fretil2ax'
!     -renamed function 'fimtil1' to 'fimtil1ax'
!     -renamed function 'fimtil2' to 'fimtil2ax'
!     -renamed function 'freal' to 'frealax'
!     -renamed function 'fim' to 'fimax'
!     -renamed subroutine 'vhat' to 'vhatax'
!     -renamed subroutine 'sae' to 'saeax'
!     -commented out many routines identically defined in 'toppik.f'
!     -modified 'tttoppikaxial' to catch unstable runs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c ************************************************************************
c      Version tuned to provide O(1%) relative accuracy for Coulomb axial
c      vertex function at first and second order (search for `cctt'):
c      - integrals A(p), Vhat, Vhhat provided analytically w/out cut-off
c      - grid range fixed to 0.1 ... 10**6 absolut
c      - and grid size enhanced to 600 points (900 foreseen in arrays).
c
c       This provides a compromise between stability and accuracy:
c       We need a relatively high momentum resolution and large maximal
c       momenta to achieve a ~1 percent accuracy, but the method of
c       direct inversion of the discretised integral equation for objects
c       whose integral is divergent induces instabilities at small
c       momenta. As the behaviour there is known, they can be cut off and
c       the vertex function fixed by hand; but limiting the grid
c       further would impact on the accuracy.
c      22.3.2017, tt
c ************************************************************************
c
c Working version with all the different original potentials
c  like (p^2+q^2)/|p-q|^2, not transformed in terms of delta and 1/r^2;
c accuracy eps=1.d-3 possible (only), but should be save, 13.8.'98, tt.
c cleaned up a bit, 24.2.1999, tt.
c
c *********************************************************************
c
c
        subroutine tttoppikaxial(xenergy,xtm,xtg,xalphas,xscale,xcutn,
     u     xcutv,
     u     xc0,xc1,xc2,xcdeltc,xcdeltl,xcfullc,xcfulll,xcrm2,
     u     xkincm,xkinca,jknflg,jgcflg,xkincv,jvflg,
     u     xim,xdi,np,xpp,xww,xdsdp,zftild)
c
c *********************************************************************
c
c !! THIS IS NOT A PUBLIC VERSION !!
c
c !!! Only P wave result given as output!!! 9.4.1999, tt.
c
c -- Calculation of the Green function in momentum space by solving the
c     Lippmann-Schwinger equation
c     F(p) = G_0(p) + G_0(p) int_0^xcutn V(p,q) q.p/p^2 F(q) dq
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
c    xkincm   :  } kinetic corrections in the 0th order Green function:
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
c                                        order Coulomb Green function
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
c    xim      :  R^{P wave}_{ttbar} from the imaginary part of the Green
c                 function
c    xdi      :  R^{P wave}_{ttbar} from the integral over the momentum
c                 distribution: int_0^xcutv dp p^3/m_t*|F(p,E)|^2
c    np       :  number of points used for the grid; fixed in tttoppik
c    xpp      :  1-dim array (max. 900 elements) giving the momenta of
c                 the Gauss-Legendre grid (pp(i) in the code)
c    xww      :  1-dim array (max. 900 elements) giving the corresponding
c                 Gauss-Legendre weights for the grid
c    xdsdp    :  1-dim array (max. 900 elements) giving the
c                 momentum distribution of top: d\sigma^{P wave}/dp,
c                  normalized to R,
c                  at the momenta of the Gauss-Legendre grid xpp(i)
c    zftild   :  1-dim array (max. 900 elements) of COMPLEX*16 numbers
c                 giving the vertex function K_A for the P-wave
c                 at the momenta of the grid.
c                 Then F(p)=K_A (p)*G_0(p) corresponding to G=K_V*G_0.
c
c *********************************************************************
c
c
           implicit none
           real*8
     u        pi,energy,vzero,eps,
     u        pp,
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,hmass,
     u        xx,critp,consde,
     u        w1,w2,sig1,sig2,const,
     u        gtpcor,etot,
     u        xenergy,xtm,xtg,xalphas,xscale,xc0,xc1,xc2,xim,xdi,
     u        xaai,xaad,xdsdp,xpp,xww,
     u        cplas,scale,c0,c1,c2,cdeltc,cdeltl,cfullc,cfulll,crm2,
     u        chiggs,xcutn,dcut,xcutv,
     u        xp,xpmax,
     u        kincom,kincoa,kincov,xkincm,xkinca,xkincv,
     u        xcdeltc,xcdeltl,xcfullc,xcfulll,xcrm2
           complex*16 bb,vec,gg,a1,aax,g0,g0c,zvfct,zftild
           integer i,n,nmax,npot,np,gcflg,kinflg,jknflg,jgcflg,
     u             jvflg,vflag
           parameter (nmax=900)
           dimension pp(nmax),bb(nmax),vec(nmax),xx(nmax),gg(nmax),
     u               w1(nmax),w2(nmax),a1(nmax),
     u               xdsdp(nmax),xpp(nmax),xww(nmax),
     u               zvfct(nmax),zftild(nmax)
c
           external aax,gtpcor,g0,g0c
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
           tmass=xtm
           energy=xenergy
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
c Cut for divergent potential-terms for large momenta in the function vhatax
c  and in the integrals aax(p):
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
c             call iniphc(1)
c             call vqdelt(consde)
              write(*,*) ' Not supplied with this version. Stop.'
              stop
           else
              write (*,*) ' Potential not implemented! Stop. 1'
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
c          critp = sqrt((energy-vzero)*tmass)
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
              critp=max(tmass/3.d0,2.d0*sqrt(energy*tmass))
           endif
c          call gauleg(0.d0, critp, pp, w1, 2*n/3)
c          call gauleg(1.d0/xpmax, 1.d0/critp, pp(2*n/3+1),
c     u                 w1(2*n/3+1), n/3)
cctt Tuned March 2017 for best possible numerical behaviour of P-wave
           call gauleg(0.1d0, 2.d0, pp, w1, 10)
           call gauleg(2.d0, critp, pp(11), w1(11), 2*n/3-10)
           call gauleg(1.d-6, 1.d0/critp, pp(2*n/3+1),
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
c Calculate the integrals aax(p) for the given momenta pp(i)
c  and store weights and momenta for the output arrays:
           do 40 i=1,n
              a1(i) = aax(pp(i)) !!! FB: can get stuck in original Toppik!
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
           call saeax(pp, w1, bb, vec, a1, n)
c
c (The substitution for the integration to infinity  pp => 1/pp
c  is done already.)
           do 50 i=1,n
              zvfct(i)=bb(i)
              zftild(i)=vec(i)
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
c Proove of the optical theorem for the output values of saeax:
c  Simply check if sig1 = sig2.
           sig1 = 0.d0
           sig2 = 0.d0
           xaai = 0.d0
           xaad = 0.d0
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
cc     u                  *tmass/sqrt(tmass**2+pp(i)**2)
c                xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
c     u                  tgamma*gtpcor(pp(i),etot)
c     u                  *(1.d0-pp(i)**2/2.d0/tmass**2)
c     u                  /(2.d0*pi**2)*const
              else
                 sig2 = sig2 + w1(i)*pp(i)**2*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
c                xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
c     u                  tgamma*gtpcor(pp(i),etot)
c     u                  /(2.d0*pi**2)*const
              endif
              xdsdp(i)=pp(i)**4/tmass**2*abs(zftild(i)*g0c(pp(i)))**2
     u                 *tgamma*gtpcor(pp(i),etot)
     u                 /(2.d0*pi**2)*const
              xaai=xaai+w1(i)*pp(i)**4/tmass**2*
     u                  aimag(zftild(i)*g0c(pp(i)))
              xaad=xaad+w1(i)*pp(i)**4/tmass**2*
     u                  abs(zftild(i)*g0c(pp(i)))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
c             write(*,*) 'xdsdp = ',xdsdp(i)
c             write(*,*) 'zvfct = ',zvfct(i)
c             write(*,*) 'zftild = ',zftild(i)
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
cc     u                  *tmass/sqrt(tmass**2+pp(i)**2)
c                xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
c     u                  tgamma*gtpcor(pp(i),etot)
c     u                  *(1.d0-pp(i)**2/2.d0/tmass**2)
c     u                  /(2.d0*pi**2)*const
              else
                 sig2 = sig2 + w1(i)*pp(i)**4*abs(gg(i))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
c                 xdsdp(i)=pp(i)**2*abs(gg(i))**2 *
c     u                  tgamma*gtpcor(pp(i),etot)
c     u                  /(2.d0*pi**2)*const
              endif
              xdsdp(i)=pp(i)**4/tmass**2*abs(zftild(i)*g0c(pp(i)))**2
     u                 *tgamma*gtpcor(pp(i),etot)
     u                 /(2.d0*pi**2)*const
              xaai=xaai+w1(i)*pp(i)**6/tmass**2*
     u                  aimag(zftild(i)*g0c(pp(i)))
              xaad=xaad+w1(i)*pp(i)**6/tmass**2*
     u                  abs(zftild(i)*g0c(pp(i)))**2 *
     u                  tgamma*gtpcor(pp(i),etot)
c             write(*,*) 'xdsdp = ',xdsdp(i)
c             write(*,*) 'zvfct = ',zvfct(i)
c             write(*,*) 'zftild = ',zftild(i)
70         continue
           do 71 i=n+1,nmax
             xdsdp(i)=0.d0
             zvfct(i)=(0.d0,0.d0)
             zftild(i)=(0.d0,0.d0)
71         continue
c
c Normalisation on R:
           sig1  = sig1 / (2.d0*pi**2) * const
           sig2  = sig2 / (2.d0*pi**2) * const
c
c The results from the momentum space approach finally are:
cc Jetzt Minus hier, 2.6.98, tt.
c          xim=-sig1
c          xdi=sig2
           xaai=-xaai / (2.d0*pi**2) * const
           xaad=xaad / (2.d0*pi**2) * const
c Output of P wave part only:
           xim=xaai
           xdi=xaad
c          write(*,*) 'vvi = ',-sig1,' .  vvd = ',sig2
c          write(*,*) 'aai = ',xim,' .  aad = ',xdi
c
        end
c
c
c
c
        complex*16 function aax(p)
c
c Neue Funktion fuer die Integrale aax(p), die hier im Falle Cutoff -> infinity
c  fuer reine Coulombpotentiale vollstaendig analytisch loesbar sind.
c  22.3.2001, tt.
c
           implicit none
           complex*16 zi,zb,zlp,zlm,zalo,zanlo,zannlo,zahig,za
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,hmass,
     u        pi,energy,vzero,eps,
     u        p,zeta3,cf,ca,tf,xnf,b0,b1,a1,a2,cnspot,phiint,
     u        cplas,scale,c0,c1,c2,
     u        cdeltc,cdeltl,cfullc,cfulll,crm2,chiggs
           integer npot
           parameter(zi=(0.d0,1.d0),zeta3=1.20205690316d0,
     u               cf=4.d0/3.d0,ca=3.d0,tf=1.d0/2.d0,xnf=5.d0)
c
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
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
           cnspot=-4.d0/3.d0*4.d0*pi
           phiint=cnspot*alphas
c
           zb=sqrt(tmass*cmplx(energy,tgamma,kind=kind(0d0)))
           zlp=log(zb+p)
           zlm=log(zb-p)
c LO: no log in z-integral
           zalo=zi*pi/2.d0/p*(zlp-zlm)
c from NL0: log in the z-integral
           zanlo=pi/2.d0/p*(zlp-zlm)*(pi+zi*(zlp+zlm))
c from NNLO: log**2 in the z-integral
           zannlo=pi/3.d0/p*(zlp-zlm)
     u           *(3.d0*pi*(zlp+zlm)+2.d0*zi*(zlm**2+zlm*zlp+zlp**2))
c Sum of the Coulomb contributions:
           za=c0*zalo-c1*(zanlo-2.d0*dlog(scale)*zalo)
     u       +c2*(zannlo-4.d0*dlog(scale)*zanlo
     u                  +4.d0*dlog(scale)**2*zalo)
c (Higgs) Yukawa contribution
cctt       zahig=zi*pi/2.d0/p*log((zb+p+zi*hmass)/(zb-p+zi*hmass))
c Alltogether:
cctt       aax=-tmass/(4.d0*pi**2)*(phiint*za+chiggs*zahig)  
           aax=-tmass/(4.d0*pi**2)*phiint*za
c
c          write(*,*) 'aax(',p,')= ',aax
        end
c
        real*8 function fretil1ax(xk)
           implicit none
           real*8 xk, frealax
           external frealax
           fretil1ax = frealax(xk)
        end
c
        real*8 function fretil2ax(xk)
           implicit none
           real*8 xk, frealax
           external frealax
           fretil2ax = frealax(1.d0/xk) * xk**(-2)
        end
c
        real*8 function fimtil1ax(xk)
           implicit none
           real*8 xk, fimax
           external fimax
           fimtil1ax = fimax(xk)
        end
c
        real*8 function fimtil2ax(xk)
           implicit none
           real*8 xk, fimax
           external fimax
           fimtil2ax = fimax(1.d0/xk) * xk**(-2)
        end
c
        real*8 function frealax(xk)
           implicit none
           complex*16 vhatax
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
           external vhatax, g0, g0c, gtpcor
c
           frealax = real(g0c(xk)*vhatax(p, xk))
        end
c
        real*8 function fimax(xk)
           implicit none
           complex*16 vhatax
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
           external vhatax, g0, g0c, gtpcor
           fimax = aimag(g0c(xk)*vhatax(p, xk))
        end
c
c
        complex*16 function vhatax(p, xk)
c
           implicit none
           complex*16 zi
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        p, xk,
     u        cnspot, phiint, AD8GLE,
     u        pm, xkm,
c     u        phfqcd, ALPHEF,
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
           external AD8GLE
c     u            , phfqcd, ALPHEF
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
              xkpln2nd=-4.d0*log(scale/xk)
              xkpln3rd=-6.d0*log(scale/xk)**2
           else if (xk/p.le.1.d-5.and.xk.le.1.d-5) then
              xkpln1st=2.d0*(xk/p)**2
              xkpln2nd=-4.d0*(xk/p)**2*log(scale/p)
              xkpln3rd=-6.d0*(xk/p)**2*log(scale/p)**2
           else
c             xkpln1st=xk/p*log(abs((p+xk)/(p-xk)))
              xkpln1st=xk/p*(log(p+xk)-log(abs(p-xk)))
cctt sign checked again, 2.2.2017, tt.
              xkpln2nd=xk/p*(-1.d0)*(log(scale/(p+xk))**2-
     u                               log(scale/abs(p-xk))**2)
              xkpln3rd=xk/p*(-4.d0/3.d0)*(log(scale/(p+xk))**3-
     u                                    log(scale/abs(p-xk))**3)
           endif
c
c          if (npot.eq.2) then
c             if (p/xk.le.1.d-5.and.p.le.1.d-5) then
c                vhatax = 2.d0 * cnspot * ALPHEF(xk)
c             else if (xk/p.le.1.d-5.and.xk.le.1.d-5) then
c                vhatax = 2.d0 * cnspot * xk**2 / p**2 * ALPHEF(p)
c             else
c                phiint = cnspot * (AD8GLE(phfqcd, 0.d0, 0.3d0, 1.d-5)
c     u                            +AD8GLE(phfqcd, 0.3d0, 1.d0, 1.d-5))
c                vhatax   = xk / p * log(abs((p+xk)/(p-xk))) * phiint
c             endif
c          else
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
                 write (*,*) ' Potential not implemented! Stop. 3'
                 stop
              endif
              phiint=cnspot*alphas
c
c             if ((xk+p).le.dcut) then
c                vhatax=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c     u               -1.d0/2.d0*(1.d0+2.d0*ca/cf)
c     u                *(pi*cf*alphas)**2/tmass
c     u                *xk/p*(p+xk-abs(xk-p))
c             else if (abs(xk-p).lt.dcut) then
c                vhatax=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c     u               -1.d0/2.d0*(1.d0+2.d0*ca/cf)
c     u                *(pi*cf*alphas)**2/tmass
c     u                *xk/p*(dcut-abs(xk-p))
c             else if (dcut.le.abs(xk-p)) then
c                vhatax=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c             else
c                write(*,*) ' Not possible! Stop.'
c                stop
c             endif
c
c       ctt
c Cut not applied here, should be left hard-wired in gauleg for stability of axial part. March 2017, tt.
c             if (max(xk,p).lt.dcut) then
c Coulomb + first + second order corrections:
                 vhatax=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c All other potentials:
c     u               +cdeltc*2.d0*xk**2
c     u               +cdeltl*xk/p/2.d0*(
c     u                (p+xk)**2*(log(((p+xk)/scale)**2)-1.d0)-
c     u                (p-xk)**2*(log(((p-xk)/scale)**2)-1.d0))
c     u               +cfullc*(p**2+xk**2)*xkpln1st
c     u               +cfulll*(p**2+xk**2)*xk/p/4.d0*
c     u                 (log(((p+xk)/scale)**2)**2-
c     u                  log(((p-xk)/scale)**2)**2)
c     u               +crm2*xk/p*(p+xk-abs(xk-p))
c             else
c                vhatax=phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd)
c             endif
c          endif
c
        end
c
c
        complex*16 function vhhat(p, xk)
c
           implicit none
           complex*16 zi
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        p, xk,
     u        cnspot, phiint, AD8GLE,
     u        pm, xkm,
     u        zeta3,cf,ca,tf,xnf,a1,a2,b0,b1,
     u        cplas,scale,c0,c1,c2,
     u        cdeltc,cdeltl,cfullc,cfulll,crm2,
     u        xkpln1st,xkpln2nd,
     u        pp,pmax,dcut,hmass,chiggs
           integer npot
           parameter(zi=(0.d0,1.d0))
           parameter(zeta3=1.20205690316d0,
     u               cf=4.d0/3.d0,ca=3.d0,tf=1.d0/2.d0,
     u               xnf=5.d0)
c
           external AD8GLE
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
              if (npot.eq.1) then
                 c0=1.d0
                 c1=0.d0
                 c2=0.d0
              else if (npot.eq.3) then
                 c0=1.d0+alphas/(4.d0*pi)*a1
                 c1=alphas/(4.d0*pi)*b0
                 c2=0
              else if (npot.eq.4) then
         write(*,*) '2nd order Coulomb in Vhhat not implemented yet.'
         stop
         c0=1.d0+alphas/(4.d0*pi)*a1+(alphas/(4.d0*pi))**2*a2
         c1=alphas/(4.d0*pi)*b0+
     u             (alphas/(4.d0*pi))**2*(b1+2.d0*b0*a1)
         c2=(alphas/(4.d0*pi))**2*b0**2
              else if (npot.eq.5) then
              else
                 write (*,*) ' Potential not implemented! Stop. 4'
                 stop
              endif
              phiint=cnspot*alphas
c
cctt No cut-off description used here either.
c             if (max(xk,p).lt.dcut) then
cctt Pure Coulomb in first order and second order only:
c
              xkpln1st=-(xk/p)**2*(1.d0+(xk**2+p**2)/(2.d0*xk*p)*
     u                  (dlog(dabs(p-xk))-dlog(p+xk)))
c             xkpln1st=-(xk/p)**2*(1.d0+(xk**2+p**2)/(4.d0*xk*p)*
c     u                  (dlog((p-xk)**2)-2.d0*dlog(p+xk)))
c
              xkpln2nd=((xk/p)**2/2.d0+xk*(xk**2+p**2)/8.d0/p**3*
     u                (dlog((p-xk)**2)-2.d0*dlog(p+xk)))*
     u                (-2.d0+dlog((xk-p)**2/scale**2)
     u                      +dlog((xk+p)**2/scale**2))
c
cctt 3rd order not yet.       xkpln3rd=
              if (c2.ne.0.d0) then
         write(*,*) ' Vhhat: 2nd order not implemented yet. Stop.'
         stop
              endif
c       
cctt          vhhat=dcmplx(phiint*(c0*xkpln1st+c1*xkpln2nd+c2*xkpln3rd),
cctt     u                     0.d0)
              vhhat=cmplx(phiint*(c0*xkpln1st+c1*xkpln2nd),
     u                      0.d0,kind=kind(0d0))
c             else
c                vhhat=(0.d0,0.d0)
c             endif
c
        end
c
c
c
c
c --- Routines for solving linear equations and matrix inversion (complex) ---
c
        subroutine saeax(pp, w1, bb, vec, a1, n)
c
           implicit none
           complex*16 vhatax,vhhat
           real*8
     u        tmass,tgamma,zmass,alphas,alamb5,
     u        wmass,wgamma,bmass,GFERMI,
     u        pi, energy, vzero, eps,
     u        d, d1, pp, w1, gtpcor,hmass,
     u        xp,xpmax,dcut,kincom,kincoa,kincov
           complex*16 aax, a1, bb, vec, ff, kk, cw, svw, g0, g0c
           integer i, j, npot, n, nmax, indx,kinflg,gcflg,vflag
           parameter (nmax=900)
           dimension bb(nmax),vec(nmax),ff(nmax,nmax),kk(nmax,nmax),
     u               pp(nmax),w1(nmax),indx(nmax),cw(nmax),a1(nmax)
c
      COMMON/PHCONS/TMASS,TGAMMA,ZMASS,ALPHAS,ALAMB5,
     $ WMASS,WGAMMA,BMASS,GFERMI,hmass
           common/ovalco/ pi, energy, vzero, eps, npot
           common/mom/ xp,xpmax,dcut
           common/g0inf/kincom,kincoa,kincov,kinflg,gcflg,vflag
c
           external aax, vhatax, gtpcor, g0, g0c, vhhat
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
c     u           gtpcor(pp(i),2.d0*tmass+energy),kind=kind(0d0)) /
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
c
c Without extra kinematic corrections:
              vec(i)=(1.d0,0.d0)
c
              svw = (0.d0,0.d0)
              do 40 j=1,n
                 if (i.ne.j) then
                    ff(i,j) = - vhatax(pp(i),pp(j)) * cw(j)
                    kk(i,j) = - vhhat(pp(i),pp(j)) * cw(j)
                    svw = svw + ff(i,j)
                 endif
40            continue
              ff(i,i) = 1.d0 - a1(i) - svw
              kk(i,i) = ff(i,i)
30         continue
c
           call zldcmp(ff, n, nmax, indx, d)
           call zldcmp(kk, n, nmax, indx, d1)
           call zlbksb(ff, n, nmax, indx, bb)
           call zlbksb(kk, n, nmax, indx, vec)
c
        end
c
c
