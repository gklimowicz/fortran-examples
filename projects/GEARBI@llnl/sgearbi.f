      subroutine drivbi (n, t0, h0, y0, tout, eps, fly, mf, indx, matp)
c
c this is the 2 August 2002 version of
c gearbi, a package for the solution of the initial value
c problem for systems of ordinary differential equations,
c          dy/dt = f(y,t),    y = (y(1),y(2),...,y(n)).
c This is the single precision version.
c gearbi is a variant of the gear package.  here, in the stiff case,
c the jacobian matrix df/dy is treated in block form, with blocks of
c size mp by mp, and mq blocks in each direction (n = mp*mq).
c linear systems are solved a block-iterative (block-sor) method.
c subroutine drivbi is a driver routine for the gearbi package.
c
c this version differs from the december 1977 version in that..
c (a) the history array, renamed yh, is in a common block /gear0/,
c (b) dimensions were altered, using parameters maxp, maxq, mxaaln.
c (c) error control now has ymax(i) = abs(y(i)) + fly(ip), where
c     ip = i(mod mp), i.e. i = ip + (iq - 1)*mp, ip.le.mp, iq.le.mq,
c     and fly is input,
c (d) off-diagonal blocks are accounted for in solbi by calls to
c     the user-supplied routine offdi (na and ida eliminated),
c (e) the sign of the scalar con used in psetbi/solbi was changed.
c
c                      reference
c     a. c. hindmarsh, preliminary documentation of gearbi.. solution
c       of ode systems with block-iterative treatment of the jacobian,
c       ucid-30149, lawrence livermore laboratory, p.o. box 808,
c       livermore, ca 94550, december 1976.
c-----------------------------------------------------------------------
c drivbi is to be called once for each output value of t, and
c in turn makes repeated calls to the core integrator, stifbi.
c
c the input parameters are..
c   n     =  the number of first-order differential equations.
c              n must be equal to mp*mq (see matp below)
c   t0    =  the initial value of t, the independent variable
c              (used only on first call).
c   h0    =  the next step size in t (used for input only on the
c              first call).
c   y0    =  a vector of length n containing the initial values of
c              y (used for input only on first call).
c   tout  =  the value of t at which output is desired next.
c              integration will normally go slightly beyond tout
c              and the package will interpolate to t = tout.
c   eps   =  the relative error bound  (used only on the
c              first call, unless indx = -1).  single step error
c              estimates divided by ymax(i) will be kept less than
c              eps in root-mean-square norm (i.e. euclidean norm
c              divided by sqrt(n) ).  the vector ymax of
c              weights is computed in drivbi, according to..
c                ymax(i) = abs(y(i)) + fly(ip), where i = ip(mod mp),
c              i.e. i = ip + (iq - 1)*mp, ip.le.mp, iq.le.mq.  this
c              gives a mixture of relative and absolute error control.
c   fly   =  an array of length mp containing floor values for y,
c              for purposes of error control (see above under eps).
c   mf    =  the method flag  (used only on first call, unless
c              indx = -1).  allowed values are 10, 11, 13,
c              20, 21, 23.  mf has two decimal digits, meth
c              and miter  (mf = 10*meth + miter).
c              meth is the basic method indicator..
c                meth = 1  means the adams methods.
c                meth = 2  means the backward differentiation
c                          formulas (bdf), or stiff methods of gear.
c              miter is the iteration method indicator..
c                miter = 0 means functional iteration (no partial
c                          derivatives needed).
c                miter = 1 means chord method with jacobian
c                          constructed with user-supplied routines
c                          pdbd and aset.  (see below.)
c                miter = 3 means chord method with jacobian replaced
c                          by a diagonal approximation based on a
c                          directional derivative.
c   indx  =  integer used on input to indicate type of call,
c              with the following values and meanings..
c                 1    this is the first call for this problem.
c                 0    this is not the first call for this problem,
c                      and integration is to continue.
c                -1    this is not the first call for the problem,
c                      and the user has reset eps and/or mf.
c                 2    same as 0 except that tout is to be hit
c                      exactly (no interpolation is done).
c                      assumes tout .ge. the current t.
c                 3    same as 0 except control returns to calling
c                      program after one step.  tout is ignored.
c              since the normal output value of indx is 0,
c              it need not be reset for normal continuation.
c   matp  =  an integer vector containing structure information
c              about the jacobian matrix j = df/dy.  j is assumed
c              to be in block form, with a block size of mp by mp,
c              and mq blocks in each direction.  off-diagonal blocks
c              have coefficients stored in an array aa.
c              the individual elements of matp are as follows..
c              matp(1) = mp = size of blocks in j.
c              matp(2) = mq = number of blocks in each direction in j.
c
c after the initial call, if a normal return occurred and a normal
c continuation is desired, simply reset tout and call again.
c all other parameters will be ready for the next call.
c a change of parameters with indx = -1 can be made after
c either a successful or an unsuccessful return.
c
c the output parameters are..
c   h0    =  the step size h used last, whether successfully or not.
c   y0    =  the computed values of y at t = tout.
c   tout  =  the output value of t.  if integration was successful,
c              and the input value of indx was not 3, tout is
c              unchanged from its input value.  otherwise, tout
c              is the current value of t to which integration
c              has been completed.
c   indx  =  integer used on output to indicate results,
c              with the following values and meanings..
c         0    integration was completed to tout or beyond.
c        -1    the integration was halted after failing to pass the
c              error test even after reducing h by a factor of
c              1.e10 from its initial value.
c        -2    after some initial success, the integration was
c              halted either by repeated error test failures or by
c              a test on eps.  too much accuracy has been requested.
c        -3    the integration was halted after failing to achieve
c              corrector convergence or block-iterative convergence
c              even after reducing h by a factor of 1.e10 from
c              its initial value.  see printed messages.
c        -4    immediate halt because of illegal values of input
c              parameters.  see printed message.
c        -5    indx was -1 on input, but the desired changes of
c              parameters were not implemented because tout
c              was not beyond t.  interpolation to t = tout was
c              performed as on a normal return.  to try again,
c              simply call again with indx = -1 and a new tout.
c
c  in addition to drivbi, the following routines are provided in
c  the package..
c    interp(tout,y,n0,y0)  interpolates to get the output values
c            at t = tout, from the data in the y array.
c    stifbi(yh,n0)  is the core integrator routine.  it performs a
c            single step and associated error control.
c    coset(meth,nq,el,tq,maxder)  sets coefficients for use in
c            the core integrator.
c    psetbi(y,n0,con,miter,ier)  computes and processes the jacobian
c            matrix of partial derivatives, j = df/dy.
c    solbi(b,x,con,z,ier)  performs the block-iterative (block-sor)
c            solution of a linear system.
c    omega(nu)  computes the block-sor overrelaxation parameter.
c    dec and sol  solve dense linear algebraic systems.
c  note.. the last five routines are called only if miter = 1.
c
c  the following routines are to be supplied by the user..
c    diffun(n,t,y,ydot)  computes the function ydot = f(y,t), the
c            right-hand side of the o.d.e.  here y and ydot are
c            vectors of length n.
c    pdbd(t,yi,i,di,mp)   computes the mp by mp block of the jacobian
c            j in the i-th diagonal block position, as a function of
c            t and the vector yi, which is the i-th segment of length
c            mp in the y vector.  the block is to be stored in di with
c            a column length of mp.
c    aset(t,y,aa,mq)  computes the nonzero elements in the off-diagonal
c            blocks of the jacobian, and stores them in aa.
c            the length and structure of the array aa are defined by
c            the user, with the length .le. the parameter mxaaln below.
c            if the array aa is constant, the call to aset in psetbi
c            should be moved to drivbi, following statement 15 below,
c            for maximum efficiency.
c            mq is passed for possible dynamic dimensioning use.
c    offdi(z,x,k,c,aa,mp,mq)  replaces z by z + c*j(k)*x, where
c              z is a given vector of length mp,
c              x is a given vector of length n,
c              c is a given scalar, and
c              j(k) is the k-th block-row of the jacobian matrix j,
c                with the diagonal block removed.
c            that is, the nonzero off-diagonal blocks in block
c            positions (k,l) are to be multiplied by the l-th block
c            of x, and the result, multiplied by c, is to be added to z.
c            all needed coefficients are to be taken from aa, as set by
c            aset.  mp and mq are passed for dynamic dimensioning.
c            offdi is called by solbi with k = 1,2,...,mq.
c            offdi should not alter x, k, c, aa, mp, or mq.
c  note.. if miter is not 1, pdbd, aset, and offdi are not called
c         and dummy routines can be substituted.
c
c the lengths of the common blocks below are set according to the
c parameters maxp, maxq, and mxaaln, which are upper bounds on mp, mq,
c and the size of aa, respectively.  for larger sizes, change the
c first two parameter statements below.
c the second dimension of yh can be reduced, to 1 plus the maximum
c order, if that maximum is reduced from 5 in coset and stifbi.
c the column length of the yh array as used elsewhere is n0, not maxn.
c if miter = 0, then dd, aa, and ipiv are not used.
c if miter = 3, then aa and ipiv are not used, and dd need only
c be of length n.
c
c the common block gear9 can be accessed externally by the user
c if desired.  it contains the following..
c hused  = the step size last used (successfully).
c nqused = the method order last used (successfully).
c nstep  = the number of steps taken for the problem so far.
c nfe    = the number of f evaluations for the problem so far.
c nje    = the number of jacobian evaluations for the problem so far.
c nii    = the number of inner (block-sor) iterations so far.
c
c in the following data statement, set..
c   uround =  the unit roundoff of the machine, i.e. the smallest
c             positive u such that 1. + u .ne. 1. on the machine.
c   lout   =  the logical unit number for the output of messages
c             during the integration.
c-----------------------------------------------------------------------
      parameter (maxp = 2, maxq = 400)
      parameter (mxaaln = 4*maxq)
      parameter (maxn = maxp*maxq, mxddln = maxp*maxn)
      dimension y0(n), fly(*), matp(2)
      common /gear0/ yh(maxn,6)
      common /gear1/ t,h,hmin,hmax,epsc,uround,nc,mfc,kflag,jstart
      common /gear2/ ymax(maxn)
      common /gear3/ error(maxn)
      common /gear4/ save1(maxn)
      common /gear5/ save2(maxn)
      common /gear5a/ save3(maxp)
      common /gear6/ dd(mxddln)
      common /gear6a/ aa(mxaaln)
      common /gear7/ ipiv(maxn)
      common /gear8/ epinr,srbnd,mxcor,mxnu,mp,mq,mpsq
      common /gear9/ hused,nqused,nstep,nfe,nje,nii
      data uround/7.1e-15/, lout/6/
c The following line is needed to save internal variables between calls.
      save toutp, n0
c
      if (indx .eq. 0) go to 20
      if (indx .eq. 2) go to 25
      if (indx .eq. -1) go to 30
      if (indx .eq. 3) go to 40
      if (indx .ne. 1) go to 430
      if (eps .le. 0.) go to 400
      if (n .le. 0) go to 410
      if ((t0-tout)*h0 .ge. 0.) go to 420
c-----------------------------------------------------------------------
c if initial values of ymax other than those set below are desired,
c they should be set here.  all ymax(i) must be positive.
c if values for hmin or hmax, the bounds on abs(h), other than
c those below are desired, they should be set below.
c-----------------------------------------------------------------------
      mp = matp(1)
      mq = matp(2)
      if (n .ne. mp*mq) go to 410
      if (mp .gt. maxp) go to 410
      if (mq .gt. maxq) go to 410
      do 12 iq = 1,mq
        do 10 ip = 1,mp
          i = ip + (iq - 1)*mp
          ymax(i) = abs(y0(i)) + fly(ip)
          yh(i,1) = y0(i)
 10       continue
 12     continue
      nc = n
      t = t0
      h = h0
      if ((t+h) .eq. t) write(lout,15)
 15   format(' warning..  t + h = t on next step.')
      hmin = abs(h0)
      hmax = abs(t0-tout)*10.
      epsc = eps
      mfc = mf
      jstart = 0
      n0 = n
      epinr = .01
      mxcor = 5
      mxnu = 20
      mpsq = mp**2
c aset called here if the coefficients aa are constant. ----------------
c     call aset (t, yh, aa, mq)
      nhcut = 0
      go to 50
c
c toutp is the previous value of tout for use in hmax. -----------------
 20   hmax = abs(tout-toutp)*10.
      go to 80
c
 25   hmax = abs(tout-toutp)*10.
      if ((t-tout)*h .ge. 0.) go to 500
      go to 85
c
 30   if ((t-tout)*h .ge. 0.) go to 440
      jstart = -1
      nc = n
      epsc = eps
      mfc = mf
c
 40   if ((t+h) .eq. t) write(lout,15)
c
 50   call stifbi (yh, n0)
c
      kgo = 1 - kflag
      go to (60, 100, 200, 300, 310), kgo
c kflag  =   0,  -1,  -2,  -3,  -4
c
 60   continue
c-----------------------------------------------------------------------
c normal return from integrator.
c
c the weights ymax(i) are updated.  if different values are desired,
c they should be set here.  a test is made for eps being too small
c for the machine precision.
c
c any other tests or calculations that are required after every
c step should be inserted here.
c
c if indx = 3, y0 is set to the current y values on return.
c if indx = 2, h is controlled to hit tout (within roundoff
c error), and then the current y values are put in y0 on return.
c for any other value of indx, control returns to the integrator
c unless tout has been reached.  then interpolated values of y are
c computed and stored in y0 on return.
c if interpolation is not desired, the call to interp should be
c removed and control transferred to statement 500 instead of 520.
c-----------------------------------------------------------------------
      d = 0.
      do 72 iq = 1,mq
        do 70 ip = 1,mp
          i = ip + (iq - 1)*mp
          ayi = abs(yh(i,1))
          ymax(i) = ayi + fly(ip)
          d = d + (ayi/ymax(i))**2
 70       continue
 72     continue
      d = d*(uround/eps)**2
      if (d .gt. float(n)) go to 250
      if (indx .eq. 3) go to 500
      if (indx .eq. 2) go to 85
 80   if ((t-tout)*h .lt. 0.) go to 40
      call interp (tout, yh, n0, y0)
      go to 520
 85   if (((t+h)-tout)*h .le. 0.) go to 40
      if (abs(t-tout) .le. 100.*uround*hmax) go to 500
      if ((t-tout)*h .ge. 0.) go to 500
      h = (tout - t)*(1. - 4.*uround)
      jstart = -1
      go to 40
c-----------------------------------------------------------------------
c on an error return from integrator, an immediate return occurs if
c kflag = -2, and recovery attempts are made otherwise.
c to recover, h and hmin are reduced by a factor of .1 up to 10
c times before giving up.
c-----------------------------------------------------------------------
 100  write (lout,105) t
 105  format(//' kflag = -1 from integrator at t = ',e16.8/
     1       '  error test failed with abs(h) = hmin'/)
 110  if (nhcut .eq. 10) go to 150
      nhcut = nhcut + 1
      hmin = .1*hmin
      h = .1*h
      write (lout,115) h
 115  format('  h has been reduced to ',e16.8,
     1       '  and step will be retried'//)
      jstart = -1
      go to 40
c
 150  write (lout,155)
 155  format(//' problem appears unsolvable with given input'//)
      go to 500
c
 200  write (lout,205) t,h
 205  format(//' kflag = -2 from integrator at t =',e16.8,'  h =',e16.8/
     1       '  the requested error is smaller than can be handled'//)
      go to 500
c
 250  write (lout,255) t
 255  format(//' integration halted by driver at t = ',e16.8/
     1      '  eps too small to be attained for the machine precision'/)
      kflag = -2
      go to 500
c
 300  write (lout,305) t
 305  format(//' kflag = -3 from integrator at t = ',e16.8/
     1       '  corrector convergence could not be achieved'/)
      go to 110
c
 310  write (lout,315) t
 315  format(//' kflag = -4 from integrator at t = ',e16.8/
     1       '  block iterative convergence could not be achieved'/)
      kflag = -3
      go to 110
c
 400  write (lout,405)
 405  format(//' illegal input.. eps .le. 0.'//)
      indx = -4
      return
c
 410  write (lout,415) n,mp,mq
 415  format(//' illegal input.. n =',i8,'  mp =',i8,'  mq =',i8//)
      indx = -4
      return
c
 420  write (lout,425)
 425  format(//' illegal input.. (t0-tout)*h .ge. 0.'//)
      indx = -4
      return
c
 430  write (lout,435) indx
 435  format(//' illegal input.. indx =',i5//)
      indx = -4
      return
c
 440  write (lout,445) t,tout,h
 445  format(//' indx = -1 on input with (t-tout)*h .ge. 0.'/
     1       ' t =',e16.8,'   tout =',e16.8,'   h =',e16.8/
     2       ' interpolation was done as on normal return.'/
     3       ' desired parameter changes were not made.')
      call interp (tout, yh, n0, y0)
      indx = -5
      return
c
 500  tout = t
      do 510 i = 1,n
 510    y0(i) = yh(i,1)
 520  indx = kflag
      toutp = tout
      h0 = hused
      if (kflag .ne. 0) h0 = h
      return
c----------------------- end of subroutine drivbi ----------------------
      end
      subroutine interp (tout, yh, n0, y0)
      dimension y0(n0),yh(n0,*)
      common /gear1/ t,h,dummy(4),n,idummy(2),jstart
c-----------------------------------------------------------------------
c subroutine interp computes interpolated values of the dependent
c variable y and stores them in y0.  the interpolation is to the
c point t = tout, and uses the nordsieck history array yh, as follows..
c                             nq
c                  y0(i)  =  sum  yh(i,j+1)*s**j ,
c                            j=0
c where s = -(t-tout)/h.
c-----------------------------------------------------------------------
      do 10 i = 1,n
 10     y0(i) = yh(i,1)
      l = jstart + 1
      s = (tout - t)/h
      s1 = 1.
      do 30 j = 2,l
        s1 = s1*s
        do 20 i = 1,n
 20       y0(i) = y0(i) + s1*yh(i,j)
 30     continue
      return
c----------------------- end of subroutine interp ----------------------
      end
      subroutine stifbi (yh, n0)
      dimension yh(n0,*)
      common /gear1/ t,h,hmin,hmax,eps,uround,n,mf,kflag,jstart
      common /gear2/ ymax(1)
      common /gear3/ error(1)
      common /gear4/ save1(1)
      common /gear5/ save2(1)
      common /gear5a/ save3(1)
      common /gear6/ dd(1)
      common /gear7/ ipiv(1)
      common /gear8/ epinr,srbnd,mxcor,mxnu,mp,mq,mpsq
      common /gear9/ hused,nqused,nstep,nfe,nje,nii
c-----------------------------------------------------------------------
c stifbi performs one step of the integration of an initial value
c problem for a system of ordinary differential equations, with
c block-iterative treatment of the jacobian matrix.
c communication with stifbi is done with the following variables..
c
c   yh      an n0 by lmax array containing the dependent variables
c             and their scaled derivatives.  lmax - 1 = maxder
c             is the maximum order available.  it is currently set to
c             5 by coset and stifbi.
c             yh(i,j+1) contains the j-th derivative of y(i), scaled by
c             h**j/factorial(j)  (j = 0,1,...,nq).
c   n0      a constant integer .ge. n, used for dimensioning purposes.
c   t       the independent variable. t is updated on each step taken.
c   h       the step size to be attempted on the next step.
c             h is altered by the error control algorithm during the
c             problem.  h can be either positive or negative, but its
c             sign must remain constant throughout the problem.
c   hmin,   the minimum and maximum absolute value of the step size
c    hmax     to be used for the step.  these may be changed at any
c             time, but will not take effect until the next h change.
c   eps     the relative error bound.  see description in driver.
c   n       the number of first-order differential equations.
c   mf      the method flag.  see description in driver.
c   kflag   a completion code with the following meanings..
c                     0  the step was succesful.
c                    -1  the requested error could not be achieved
c                          with abs(h) = hmin.
c                    -2  the requested error is smaller than can
c                          be handled for this problem.
c                    -3  corrector convergence could not be
c                          achieved for abs(h) = hmin.
c                    -4  block iterative convergence could not be
c                          achieved for abs(h) = hmin.
c             on a return with kflag negative, the values of t and
c             the yh array are as of the beginning of the last
c             step, and h is the last step size attempted.
c   jstart  an integer used on input and output.
c             on input, it has the following values and meanings..
c                     0  perform the first step.
c                 .gt.0  take a new step continuing from the last.
c                 .lt.0  take the next step with a new value of
c                          h, eps, and/or mf.
c             on exit, jstart is nq, the current order of the method.
c   ymax    an array of n elements with which the estimated local
c             errors in y are compared.
c   error   an array of n elements.  error(i)/tq(2) is the estimated
c             one-step error in y(i).
c   save1,  two arrays of working storage,
c    save2    each of length n.
c   save3   a working array of length mp.
c   dd      a block of locations used for partial derivatives if
c             miter is not 0.  see description in driver.
c   ipiv    an integer array of length n used for pivot
c             information if miter = 1 or 2.
c   mxcor   the maximum number of corrector iterations per step.
c-----------------------------------------------------------------------
      dimension el(13),tq(4)
      data el(2)/1./, oldl0/1./
c The following is needed to save internal variables between calls.
      save rc, miter, nstepj, nq, iweval, crate, bnd, mfold, meth, e, epsold,
     1  nold, l, idoub, lmax, eup, rmax, con, edn
c
      kflag = 0
      told = t
      if (jstart .gt. 0) go to 200
      if (jstart .ne. 0) go to 120
c-----------------------------------------------------------------------
c on the first call, the order is set to 1 and the initial ydot is
c calculated.  rmax is the maximum ratio by which h can be increased
c in a single step.  it is initially 1.e4 to compensate for the small
c initial h, but then is normally equal to 10.  if a failure
c occurs (in corrector convergence or error test), rmax is set at 2
c for the next increase.
c-----------------------------------------------------------------------
      call diffun (n, t, yh, save1)
      do 110 i = 1,n
 110    yh(i,2) = h*save1(i)
      meth = mf/10
      miter = mf - 10*meth
      nq = 1
      l = 2
      idoub = 3
      rmax = 1.e4
      rc = 0.
      crate = 1.
      hold = h
      mfold = mf
      ier = 0
      nstep = 0
      nstepj = 0
      nfe = 1
      nje = 0
      nii = 0
      iret = 3
      go to 130
c-----------------------------------------------------------------------
c if the caller has changed meth, coset is called to set
c the coefficients of the method.  if the caller has changed
c n, eps, or meth, the constants e, edn, eup, and bnd must be reset.
c e is a comparison for errors of the current order nq. eup is
c to test for increasing the order, edn for decreasing the order.
c bnd is used to test for convergence of the corrector iterates.
c if the caller has changed h, yh must be rescaled.
c if h or meth has been changed, idoub is reset to l + 1 to prevent
c further changes in h for that many steps.
c-----------------------------------------------------------------------
 120  if (mf .eq. mfold) go to 150
      meo = meth
      mio = miter
      meth = mf/10
      miter = mf - 10*meth
      mfold = mf
      if (miter .ne. mio) iweval = miter
      if (meth .eq. meo) go to 150
      idoub = l + 1
      iret = 1
 130  call coset (meth, nq, el, tq, maxder)
      if (maxder .gt. 5) maxder = 5
      lmax = maxder + 1
      rc = rc*el(1)/oldl0
      oldl0 = el(1)
 140  fn = float(n)
      edn = fn*(tq(1)*eps)**2
      e   = fn*(tq(2)*eps)**2
      eup = fn*(tq(3)*eps)**2
      bnd = fn*(tq(4)*eps)**2
      srbnd = sqrt(bnd)
      epsold = eps
      nold = n
      go to (160, 170, 200), iret
 150  if ((eps .eq. epsold) .and. (n .eq. nold)) go to 160
      if (n .ne. nold) iweval = miter
      iret = 1
      go to 140
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = amax1(rh,hmin/abs(h))
 175  rh = amin1(rh,hmax/abs(h),rmax)
      r1 = 1.
      do 180 j = 2,l
        r1 = r1*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r1
      h = h*rh
      rc = rc*rh
      idoub = l + 1
      if (iredo .eq. 0) go to 690
c-----------------------------------------------------------------------
c this section computes the predicted values by effectively
c multiplying the yh array by the pascal triangle matrix.
c rc is the ratio of new to old values of the coefficient  h*el(1).
c when rc differs from 1 by more than 30 percent, or the caller has
c changed miter, iweval is set to miter to force the partials to be
c updated, if partials are used.  in any case, the partials
c are updated at least every 20-th step.
c-----------------------------------------------------------------------
 200  if (abs(rc-1.) .gt. 0.3) iweval = miter
      if (nstep .ge. nstepj+20) iweval = miter
      t = t + h
      do 210 j1 = 1,nq
        do 210 j2 = j1,nq
          j = (nq + j1) - j2
          do 210 i = 1,n
 210        yh(i,j) = yh(i,j) + yh(i,j+1)
c-----------------------------------------------------------------------
c up to mxcor corrector iterations are taken.  a convergence test is
c made on the r.m.s. norm of each correction, using bnd, which
c is dependent on eps.  the sum of the corrections is accumulated
c in the vector error(i).  the yh array is not altered in the corrector
c loop.  the updated y vector is stored temporarily in save1.
c-----------------------------------------------------------------------
 220  m = 0
      call diffun (n, t, yh, save2)
      nfe = nfe + 1
      if (iweval .le. 0) go to 250
c-----------------------------------------------------------------------
c if indicated, the matrix p = i - h*el(1)*j is reevaluated before
c starting the corrector iteration.  iweval is set to 0 to indicate
c that this has been done.  p is computed and processed in psetbi.
c-----------------------------------------------------------------------
      iweval = 0
      rc = 1.
      nstepj = nstep
      con = h*el(1)
      call psetbi (yh, n0, con, miter, ier)
      if (ier .ne. 0) go to 420
      hl0 = h*el(1)
 250  do 260 i = 1,n
 260    error(i) = 0.
 270  if (miter .ne. 0) go to (350, 350, 310), miter
c-----------------------------------------------------------------------
c in the case of functional iteration, update y directly from
c the result of the last diffun call.
c-----------------------------------------------------------------------
      d = 0.
      do 300 i = 1,n
        r = h*save2(i) - yh(i,2)
        d = d + ( (r-error(i))/ymax(i) )**2
        save1(i) = yh(i,1) + el(1)*r
 300    error(i) = r
      go to 400
c-----------------------------------------------------------------------
c in the case of the chord method, compute the corrector error,
c f sub (m), and solve the linear system with that as right-hand
c side and p as coefficient matrix, using the lu decomposition
c if miter = 1.  if miter = 3, the coefficient h*el(1)
c in p is updated.
c-----------------------------------------------------------------------
 310  phl0 = hl0
      hl0 = h*el(1)
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        d = 1. - r*(1. - 1./dd(i))
        if (abs(d) .eq. 0.) go to 440
 320    dd(i) = 1./d
 330  do 340 i = 1,n
 340    save1(i) = dd(i)*(h*save2(i) - (yh(i,2) + error(i)))
      go to 370
 350  do 360 i = 1,n
 360    save2(i) = h*save2(i) - (yh(i,2) + error(i))
      call solbi (save2, save1, con, save3, ier)
      if (ier .ne. 0) go to 410
 370  d = 0.
      do 380 i = 1,n
        error(i) = error(i) + save1(i)
        d = d + (save1(i)/ymax(i))**2
 380    save1(i) = yh(i,1) + el(1)*error(i)
c-----------------------------------------------------------------------
c test for convergence.  if m.gt.0, the square of the convergence
c rate constant is estimated as crate, and this is used in the test.
c-----------------------------------------------------------------------
 400  if (m .ne. 0) crate = amax1(.9*crate,d/d1)
      if ((d*amin1(1.,2.*crate)) .le. bnd) go to 450
      d1 = d
      m = m + 1
      if (m .eq. mxcor) go to 410
      call diffun (n, t, save1, save2)
      nfe = nfe + 1
      go to 270
c-----------------------------------------------------------------------
c corrector iteration or block-sor failed to converge.  if partials
c are involved but are not up to date, they are reevaluated for the
c next try.  otherwise the yh array is retracted to its values
c before prediction, and h is reduced, if possible.  if not, a
c no-convergence exit is taken.
c-----------------------------------------------------------------------
 410  if (iweval .eq. -1) go to 440
 420  t = told
      rmax = 2.
      do 430 j1 = 1,nq
        do 430 j2 = j1,nq
          j = (nq + j1) - j2
          do 430 i = 1,n
 430        yh(i,j) = yh(i,j) - yh(i,j+1)
      if (abs(h) .le. hmin*1.00001) go to 680
      rh = .25
      iredo = 1
      go to 170
 440  iweval = miter
      go to 220
c-----------------------------------------------------------------------
c the corrector has converged.  iweval is set to -1 if partial
c derivatives were used, to signal that they may need updating on
c subsequent steps.  the error test is made and control passes to
c statement 500 if it fails.
c-----------------------------------------------------------------------
 450  if (miter .ne. 0) iweval = -1
      d = 0.
      do 460 i = 1,n
 460    d = d + (error(i)/ymax(i))**2
      if (d .gt. e) go to 500
c-----------------------------------------------------------------------
c after a successful step, update the yh array.
c consider changing h if idoub = 1.  otherwise decrease idoub by 1.
c if idoub is then 1 and nq .lt. maxder, then error is saved for
c use in a possible order increase on the next step.
c if a change in h is considered, an increase or decrease in order
c by one is considered also.  a change in h is made only if it is by a
c factor of at least 1.1.  if not, idoub is set to 10 to prevent
c testing for that many steps.
c-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nstep = nstep + 1
      hused = h
      nqused = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*error(i)
      if (idoub .eq. 1) go to 520
      idoub = idoub - 1
      if (idoub .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = error(i)
      go to 700
c-----------------------------------------------------------------------
c the error test failed.  kflag keeps track of multiple failures.
c restore t and the yh array to their previous values, and prepare
c to try the step again.  compute the optimum step size for this or
c one lower order.
c-----------------------------------------------------------------------
 500  kflag = kflag - 1
      t = told
      do 510 j1 = 1,nq
        do 510 j2 = j1,nq
          j = (nq + j1) - j2
          do 510 i = 1,n
 510        yh(i,j) = yh(i,j) - yh(i,j+1)
      rmax = 2.
      if (abs(h) .le. hmin*1.00001) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      pr3 = 1.e+20
      go to 540
c-----------------------------------------------------------------------
c regardless of the success or failure of the step, factors
c pr1, pr2, and pr3 are computed, by which h could be divided
c at order nq - 1, order nq, or order nq + 1, respectively.
c in the case of failure, pr3 = 1.e20 to avoid an order increase.
c the smallest of these is determined and the new order chosen
c accordingly.  if the order is to be increased, we compute one
c additional scaled derivative.
c-----------------------------------------------------------------------
 520  pr3 = 1.e+20
      if (l .eq. lmax) go to 540
      d1 = 0.
      do 530 i = 1,n
 530    d1 = d1 + ((error(i) - yh(i,lmax))/ymax(i))**2
      enq3 = .5/float(l+1)
      pr3 = ((d1/eup)**enq3)*1.4 + 1.4e-6
 540  enq2 = .5/float(l)
      pr2 = ((d/e)**enq2)*1.2 + 1.2e-6
      pr1 = 1.e+20
      if (nq .eq. 1) go to 560
      d = 0.
      do 550 i = 1,n
 550    d = d + (yh(i,l)/ymax(i))**2
      enq1 = .5/float(nq)
      pr1 = ((d/edn)**enq1)*1.3 + 1.3e-6
 560  if (pr2 .le. pr3) go to 570
      if (pr3 .lt. pr1) go to 590
      go to 580
 570  if (pr2 .gt. pr1) go to 580
      newq = nq
      rh = 1./pr2
      go to 620
 580  newq = nq - 1
      rh = 1./pr1
      go to 620
 590  newq = l
      rh = 1./pr3
      if (rh .lt. 1.1) go to 610
      do 600 i = 1,n
 600    yh(i,newq+1) = error(i)*(el(l)/float(l))
      go to 630
 610  idoub = 10
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1)) go to 610
c-----------------------------------------------------------------------
c if there is a change of order, reset nq, l, and the coefficients.
c in any case h is reset according to rh and the yh array is rescaled.
c then exit from 690 if the step was ok, or redo the step otherwise.
c-----------------------------------------------------------------------
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 130
c-----------------------------------------------------------------------
c control reaches this section if 3 or more failures have occured.
c it is assumed that the derivatives that have accumulated in the
c yh array have errors of the wrong order.  hence the first
c derivative is recomputed, and the order is set to 1.  then
c h is reduced by a factor of 10, and the step is retried.
c after a total of 7 failures, an exit is taken with kflag = -2.
c-----------------------------------------------------------------------
 640  if (kflag .eq. -7) go to 670
      rh = .1
      rh = amax1(hmin/abs(h),rh)
      h = h*rh
      call diffun (n, t, yh, save1)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*save1(i)
      iweval = miter
      idoub = 10
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 130
c-----------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c-----------------------------------------------------------------------
 660  kflag = -1
      go to 700
 670  kflag = -2
      go to 700
 680  kflag = -3
      if (ier .eq. -1) kflag = -4
      go to 700
 690  rmax = 10.
 700  hold = h
      jstart = nq
      return
c----------------------- end of subroutine stifbi ----------------------
      end
      subroutine coset (meth, nq, el, tq, maxder)
c-----------------------------------------------------------------------
c coset is called by the integrator and sets coefficients used there.
c the vector el, of length nq + 1, determines the basic method.
c the vector tq, of length 4, is involved in adjusting the step size
c in relation to truncation error.  its values are given by the
c pertst array.
c the vectors el and tq depend on meth and nq.
c coset also sets maxder, the maximum order of the method available.
c currently it is 12 for the adams methods and 5 for the bdf methods.
c lmax = maxder + 1 is the number of columns in the yh array.
c the maximum order used may be reduced simply by decreasing the
c numbers in statements 1 and/or 2 below.
c
c the coefficients in pertst need be given to only about
c one percent accuracy.  the order in which the groups appear below
c is..  coefficients for order nq - 1, coefficients for order nq,
c coefficients for order nq + 1.  within each group are the
c coefficients for the adams methods, followed by those for the
c bdf methods.
c-----------------------------------------------------------------------
      dimension pertst(12,2,3),el(13),tq(4)
      data  pertst / 1.,1.,2.,1.,.3158,.07407,.01391,.002182,
     1                 .0002945,.00003492,.000003692,.0000003524,
     2               1.,1.,.5,.1667,.04167,1.,1.,1.,1.,1.,1.,1.,
     3               2.,12.,24.,37.89,53.33,70.08,87.97,106.9,
     4                 126.7,147.4,168.8,191.0,
     5               2.0,4.5,7.333,10.42,13.7,1.,1.,1.,1.,1.,1.,1.,
     6               12.0,24.0,37.89,53.33,70.08,87.97,106.9,
     7                 126.7,147.4,168.8,191.0,1.,
     8               3.0,6.0,9.167,12.5,1.,1.,1.,1.,1.,1.,1.,1. /
c
      go to (1,2),meth
 1    maxder = 12
      go to (101,102,103,104,105,106,107,108,109,110,111,112),nq
 2    maxder = 5
      go to (201,202,203,204,205),nq
c-----------------------------------------------------------------------
c the following coefficients should be defined to machine accuracy.
c for a given order nq, they can be calculated by use of the
c generating polynomial l(t), whose coefficients are el(i)..
c      l(t) = el(1) + el(2)*t + ... + el(nq+1)*t**nq.
c for the implicit adams methods, l(t) is given by
c      dl/dt = (t+1)*(t+2)* ... *(t+nq-1)/k,    l(-1) = 0,
c where                 k = factorial(nq-1).
c for the bdf methods,
c      l(t) = (t+1)*(t+2)* ... *(t+nq)/k,
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the order in which the groups appear below is..
c implicit adams methods of orders 1 to 12,
c bdf methods of orders 1 to 5.
c-----------------------------------------------------------------------
 101  el(1) = 1.0
      go to 900
 102  el(1) = 0.5
      el(3) = 0.5
      go to 900
 103  el(1) = 4.1666666666667e-01
      el(3) = 0.75
      el(4) = 1.6666666666667e-01
      go to 900
 104  el(1) = 0.375
      el(3) = 9.1666666666667e-01
      el(4) = 3.3333333333333e-01
      el(5) = 4.1666666666667e-02
      go to 900
 105  el(1) = 3.4861111111111e-01
      el(3) = 1.0416666666667
      el(4) = 4.8611111111111e-01
      el(5) = 1.0416666666667e-01
      el(6) = 8.3333333333333e-03
      go to 900
 106  el(1) = 3.2986111111111e-01
      el(3) = 1.1416666666667
      el(4) = 0.625
      el(5) = 1.7708333333333e-01
      el(6) = 0.025
      el(7) = 1.3888888888889e-03
      go to 900
 107  el(1) = 3.1559193121693e-01
      el(3) = 1.225
      el(4) = 7.5185185185185e-01
      el(5) = 2.5520833333333e-01
      el(6) = 4.8611111111111e-02
      el(7) = 4.8611111111111e-03
      el(8) = 1.9841269841270e-04
      go to 900
 108  el(1) = 3.0422453703704e-01
      el(3) = 1.2964285714286
      el(4) = 8.6851851851852e-01
      el(5) = 3.3576388888889e-01
      el(6) = 7.7777777777778e-02
      el(7) = 1.0648148148148e-02
      el(8) = 7.9365079365079e-04
      el(9) = 2.4801587301587e-05
      go to 900
 109  el(1) = 2.9486800044092e-01
      el(3) = 1.3589285714286
      el(4) = 9.7655423280423e-01
      el(5) = 0.4171875
      el(6) = 1.1135416666667e-01
      el(7) = 0.01875
      el(8) = 1.9345238095238e-03
      el(9) = 1.1160714285714e-04
      el(10)= 2.7557319223986e-06
      go to 900
 110  el(1) = 2.8697544642857e-01
      el(3) = 1.4144841269841
      el(4) = 1.0772156084656
      el(5) = 4.9856701940035e-01
      el(6) = 0.1484375
      el(7) = 2.9060570987654e-02
      el(8) = 3.7202380952381e-03
      el(9) = 2.9968584656085e-04
      el(10)= 1.3778659611993e-05
      el(11)= 2.7557319223986e-07
      go to 900
 111  el(1) = 2.8018959644394e-01
      el(3) = 1.4644841269841
      el(4) = 1.1715145502646
      el(5) = 5.7935819003527e-01
      el(6) = 1.8832286155203e-01
      el(7) = 4.1430362654321e-02
      el(8) = 6.2111441798942e-03
      el(9) = 6.2520667989418e-04
      el(10)= 4.0417401528513e-05
      el(11)= 1.5156525573192e-06
      el(12)= 2.5052108385442e-08
      go to 900
 112  el(1) = 2.7426554003160e-01
      el(3) = 1.5099386724387
      el(4) = 1.2602711640212
      el(5) = 6.5923418209877e-01
      el(6) = 2.3045800264550e-01
      el(7) = 5.5697246105232e-02
      el(8) = 9.4394841269841e-03
      el(9) = 1.1192749669312e-03
      el(10)= 9.0939153439153e-05
      el(11)= 4.8225308641975e-06
      el(12)= 1.5031265031265e-07
      el(13)= 2.0876756987868e-09
      go to 900
 201  el(1) = 1.0
      go to 900
 202  el(1) = 6.6666666666667e-01
      el(3) = 3.3333333333333e-01
      go to 900
 203  el(1) = 5.4545454545455e-01
      el(3) = el(1)
      el(4) = 9.0909090909091e-02
      go to 900
 204  el(1) = 0.48
      el(3) = 0.7
      el(4) = 0.2
      el(5) = 0.02
      go to 900
 205  el(1) = 4.3795620437956e-01
      el(3) = 8.2116788321168e-01
      el(4) = 3.1021897810219e-01
      el(5) = 5.4744525547445e-02
      el(6) = 3.6496350364964e-03
c
 900  do 910 k = 1,3
 910    tq(k) = pertst(nq,meth,k)
      tq(4) = .5*tq(2)/float(nq+2)
      return
c----------------------- end of subroutine coset -----------------------
      end
      subroutine psetbi (yh, n0, con, miter, ier)
      dimension yh(n0,*)
      common /gear1/ t,h,dummy(3),uround,n,idummy(3)
      common /gear2/ ymax(1)
      common /gear4/ save1(1)
      common /gear5/ save2(1)
      common /gear6/ dd(1)
      common /gear6a/ aa(1)
      common /gear7/ ipiv(1)
      common /gear8/ epinr,srbnd,mxcor,mxnu,mp,mq,mpsq
      common /gear9/ hused,nqused,nstep,nfe,nje,nii
c-----------------------------------------------------------------------
c psetbi is called by stifbi to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c j is computed, either by making calls to the user-supplied routine
c pdbd if miter = 1, or as a diagonal divided difference approximation
c if miter = 3.  if miter .ne. 3, the diagonal blocks of p
c are subjected to lu decomposition in preparation for
c later solution of linear systems with p as coefficient matrix.
c
c in addition to variables described previously, communication
c with psetbi uses the following variables..
c   con   = h*el(1).
c   ier   = 0 if no singularity was found in the diagonal blocks, or
c         = the index (1 to n) at which a singularity was found.
c-----------------------------------------------------------------------
      nje = nje + 1
      go to (100, 100, 300), miter
c if miter = 1, make mq calls to pdbd and do lu on each diagonal block.
 100  iy = 1
      id = 1
      conn = -con
      do 130 iq = 1,mq
        call pdbd (t, yh(iy,1), iq, dd(id), mp)
        do 110 i = 1,mpsq
 110      dd(id+i-1) = dd(id+i-1)*conn
        i1 = 0
        do 120 i = 1,mp
          dd(id+i1) = dd(id+i1) + 1.
 120      i1 = i1 + (mp + 1)
        call dec (mp, mp, dd(id), ipiv(iy), ier)
        if (ier .ne. 0) go to 190
        iy = iy + mp
        id = id + mpsq
 130    continue
c if the off-diagonal coefficients aa are not constant, call aset here.
      call aset (t, yh, aa, mq)
      return
 190  ier = (iq-1)*mp + ier
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  con = con/h
      ier = 0
      r = con*.1
      do 310 i = 1,n
 310    dd(i) = yh(i,1) + r*(h*save2(i) - yh(i,2))
      call diffun (n, t, dd, save1)
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*save2(i) - yh(i,2)
        dd(i) = 1.
        if (abs(r0) .lt. uround*ymax(i)) go to 320
        d = .1*r0 - h*(save1(i) - save2(i))
        if (abs(d) .eq. 0.) go to 330
        dd(i) = .1*r0/d
 320    continue
      return
 330  ier = i
      return
c----------------------- end of subroutine psetbi ----------------------
      end
      subroutine solbi (b, x, con, z, ier)
      dimension b(*),x(*),z(*)
      common /gear1/ dummy(6),n,idummy(3)
      common /gear2/ ymax(1)
      common /gear6/ dd(1)
      common /gear6a/ aa(1)
      common /gear7/ ipiv(1)
      common /gear8/ epinr,srbnd,mxcor,mxnu,mp,mq,mpsq
      common /gear8a/ er,erlast,relax,is1,is2,icnt,iw
      common /gear9/ hused,nqused,nstep,nfe,nje,nii
c-----------------------------------------------------------------------
c solbi is called by stifbi and solves a blocked linear system by the
c block-sor iterative method.  the coefficient matrix has diagonal
c blocks whose lu decompositions are stored in dd, and
c off-diagonal elements which are the contents of aa multiplied by
c the scalar con.
c
c communication with solbi uses the following variables..
c b     = right-hand side vector (input).
c x     = solution vector (output).
c con   = input scalar multiplying the off-diagonal elements,
c         the value of h*el(1) as of the last jacobian update.
c z     = working vector of length mp.
c ier   = 0 on return if convergence is achieved, or -1 if not.
c epinr = constant used in the convergence test, set in driver.
c srbnd = constant used in the convergence test, set in stifbi.
c mxnu  = maximum number of block-sor iterations allowed.
c-----------------------------------------------------------------------
c initial guess for x is solution of block-diagonal system. ------------
      do 5 i = 1,n
 5      x(i) = b(i)
      ix = 1
      id = 1
      do 10 iq = 1,mq
        call sol (mp, mp, dd(id), x(ix), ipiv(ix))
        id = id + mpsq
 10     ix = ix + mp
c initialize parameters and loop over iteration index nu. --------------
      dxold = 0.
      relax = 1.
      is1 = 6
      is2 = 5
      icnt = 0
      iw = 0
      do 70 nu = 1,mxnu
c loop over block index iq = 1,...,mq. ---------------------------------
        dxnorm = 0.
        ix = 0
        id = 1
        do 50 iq = 1,mq
c load segment of b into z. --------------------------------------------
          ix1 = ix + 1
          do 20 ip = 1,mp
 20         z(ip) = b(ix+ip)
c call offdi to subtract terms from off-diagonal elements. -------------
          kq = iq
          call offdi (z, x, kq, con, aa, mp, mq)
c solve mp by mp system for block of new values. -----------------------
          call sol (mp, mp, dd(id), z, ipiv(ix1))
c apply correction to x and form needed norms. -------------------------
          do 40 ip = 1,mp
            z(ip) = relax*(z(ip) - x(ix+ip))
            dxnorm = dxnorm + ( z(ip)/ymax(ix+ip) )**2
 40         x(ix+ip) = x(ix+ip) + z(ip)
          ix = ix + mp
          id = id + mpsq
 50       continue
        dxnorm = sqrt(dxnorm)
c update relaxation parameters and test for convergence. ---------------
        er = .99
        if (dxold .ne. 0.) er = dxnorm/dxold
        dxold = dxnorm
        relbnd = srbnd*epinr*(1. - er)
        if (dxnorm .lt. relbnd) go to 80
        if (nu .ge. is1  .or.  nu .ge. is2) call omega(nu)
        erlast = er
 70     continue
c failed to converge.  error return. -----------------------------------
      nii = nii + mxnu
      ier = -1
      return
c normal return. -------------------------------------------------------
 80   nii = nii + nu
      ier = 0
      return
c----------------------- end of subroutine solbi -----------------------
      end
      subroutine omega (nu)
      common /gear8a/ er,erlast,relax,is1,is2,icnt,iw
c-----------------------------------------------------------------------
c omega is called by solbi to compute the overrelaxation parameter
c used there, when miter = 1.
c programmed by n. k. madsen, based on the work of l. a. hageman.
c
c communication with omega uses the following variables..
c   nu      = iteration index = 1,2,...
c   er      = ratio of successive norms of iterate differences.
c   erlast  = previous value of er.
c   relax   = the overrelaxation parameter.
c   is1,is2 = lower limits on nu for which omega is called.
c   icnt,iw = internal counters, initialized by solbi.
c-----------------------------------------------------------------------
      dimension rlxmax(8)
      data rlxmax/1.75,1.83,1.9,1.94,1.98,1.99,1.995,1.999/
      data xc1/0.5/, xc2/0.75/
c The following is needed to save an internal variable between calls.
      save thrate
c-----------------------------------------------------------------------
c first check error reduction factor er for possible oscillatory
c behavior which indicates that relax is already too large.
c-----------------------------------------------------------------------
      if ( nu .lt. is1 ) go to 30
        temp = relax - 1.
        if ( er .ge. temp ) go to 10
          icnt = icnt + 1
          er = temp
          return
 10   if ( (er+.00005) .lt. erlast ) return
        if ( er .lt. 1. ) go to 20
          er = .9999999
          return
 20   if ( icnt .gt. 0 ) return
c-----------------------------------------------------------------------
c compute new values of relax and/or is1 and is2 if
c conditions are appropriate.  eff measures the efficiency of the
c iterations with the current value of relax.  if eff is at least xc2
c (currently .75), the current value of relax is considered adequate.
c-----------------------------------------------------------------------
 30     if ( nu .lt. is2 ) return
          if ( relax .ne. 1. ) go to 40
            rhosq = er
            if ( er .lt. 1. ) go to 50
              is1 = is1 + 4
              is2 = is2 + 4
              er = .9999999
              return
 40   eff = -alog(er) / thrate
      if ( eff .ge. xc2 ) return
        rhosq = (er+relax-1.)**2 / (relax*relax*er)
 50     relax = 2. / ( 1. + sqrt(1.-rhosq) )
        iw = min0(8,iw+1)
        if ( relax .gt. rlxmax(iw) ) relax = rlxmax(iw)
        inc = 6
 60   if ( float(inc)*(relax-1.)**inc .le. xc1 ) go to 70
        inc = inc + 1
        go to 60
 70   is1 = nu + 5
      is2 = is1 + inc
      icnt = 0
      thrate = -alog(relax-1.)
      return
c----------------------- end of subroutine omega -----------------------
      end
      subroutine dec (n, ndim, a, ip, ier)
      dimension a(ndim,n), ip(n)
c-----------------------------------------------------------------------
c matrix triangularization by gauss elimination with partial pivoting.
c input..
c    n = order of matrix.
c    ndim = declared first dimension of array  a.
c    a = matrix to be triangularized.
c output..
c    a(i,j), i.le.j = upper triangular factor, u .
c    a(i,j), i.gt.j = multipliers = lower triangular factor, i - l.
c    ip(k), k.lt.n = index of k-th pivot row.
c    ier = 0 if matrix a is nonsingular, or k if found to be
c          singular at stage k.
c row interchanges are finished in u, only partly in l.
c use  sol  to obtain solution of linear system.
c if ier .ne. 0, a is singular, sol will divide by zero.
c-----------------------------------------------------------------------
      ier = 0
      if (n .eq. 1) go to 70
      nm1 = n - 1
      do 60 k = 1,nm1
        kp1 = k + 1
c find the pivot in column k.  search rows k to n. ---------------------
        m = k
        do 10 i = kp1,n
 10       if (abs(a(i,k)) .gt. abs(a(m,k))) m = i
        ip(k) = m
c interchange elements in rows k and m. --------------------------------
        t = a(m,k)
        if (m .eq. k) go to 20
        a(m,k) = a(k,k)
        a(k,k) = t
 20     if (t .eq. 0.) go to 80
c store multipliers in a(i,k), i = k+1,...,n. --------------------------
        t = 1./t
        do 30 i = kp1,n
 30       a(i,k) = -a(i,k)*t
c apply multipliers to other columns of a. -----------------------------
        do 50 j = kp1,n
          t = a(m,j)
          a(m,j) = a(k,j)
          a(k,j) = t
          if (t .eq. 0.) go to 50
          do 40 i = kp1,n
 40         a(i,j) = a(i,j) + a(i,k)*t
 50       continue
 60     continue
 70   k = n
      if (a(n,n) .eq. 0.) go to 80
      return
 80   ier = k
      return
c----------------------- end of subroutine dec -------------------------
      end
      subroutine sol (n, ndim, a, b, ip)
      dimension a(ndim,n), b(n), ip(n)
c-----------------------------------------------------------------------
c solution of linear system a*x = b using output of dec.
c input..
c    n = order of matrix.
c    ndim = declared first dimension of array  a.
c    a = triangularized matrix obtained from dec.
c    b = right hand side vector.
c    ip = pivot information vector obtained from dec.
c do not use if dec has set ier .ne. 0.
c output..
c    b = solution vector, x .
c-----------------------------------------------------------------------
      if (n .eq. 1) go to 50
      nm1 = n - 1
c apply row permutations and multipliers to b. -------------------------
      do 20 k = 1,nm1
        kp1 = k + 1
        m = ip(k)
        t = b(m)
        b(m) = b(k)
        b(k) = t
        do 10 i = kp1,n
 10       b(i) = b(i) + a(i,k)*t
 20     continue
c back solve. ----------------------------------------------------------
      do 40 kb = 1,nm1
        km1 = n - kb
        k = km1 + 1
        b(k) = b(k)/a(k,k)
        t = -b(k)
        do 30 i = 1,km1
 30       b(i) = b(i) + a(i,k)*t
 40     continue
 50   b(1) = b(1)/a(1,1)
      return
c----------------------- end of subroutine sol -------------------------
      end
