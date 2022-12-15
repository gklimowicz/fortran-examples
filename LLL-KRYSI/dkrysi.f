      subroutine krysi(n,f,jac,psol,y,t,tend,aeps,reps,apre,ipre,
     &                  work,wk,iwk,iopt,jpre,iflag,output,ipt)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     krysi solver for integrating y'=f(t,y) , y(t)=given as y
c     over the interval (t,tend) using
c     a semi-implicit rk-method of norsett and thomsen,
c     in combination with a preconditioned krylov subspace iteration
c     method (scaled incomplete generalized minimum residual method).
c     see comments in subroutines news and spigmr.
c
c     Reference:
c       Alan C. Hindmarsh and Syvert Paul Norsett,
c       KRYSI, An ODE Solver Combining a Semi-Implicit 
c       Runge-Kutta Method and a Preconditioned Krylov Method,
c       NTH, Trondheim, 1987 (Mathematics of Computation Report 
c       No. 8/87); also available as LLNL Report UCID-21422, 
c       May 1988.
c
c     Version of 6 August 2002.
c
c     this version is in double precision.
c
c     input parameters are:
c              n    : integer giving the number of equations
c
c              f    : subroutine for righthand side.
c                     subroutine f(t,y,r)
c                     double precision t,y(n),r(n)
c                     r=f(t,y)
c
c              jac  : subroutine for the computing jacobian elements
c                     needed for preconditioning.
c                     subroutine jac(f,n,t,y,ewt,fty,ftem,gah,
c                                    apre,ipre,jflag)
c                     double precision y(n),ewt(n),fty(n),ftem(n),
c                                      apre(1),ipre(1)
c                     this routine must evaluate and preprocess any
c                     parts of the jacobian matrix df/dy used in the
c                     left and right preconditioner matrices, p1 and p2.
c                     fty is space for the current value of f(t,y),
c                     and should be loaded by a call to f if needed.
c                     ftem is work space, e.g. for values of f(t,y+e)
c                     for use in difference quotients approximations.
c                     on computing jacobian elements, jac must multiply
c                     all computed elements by -gah and add the
c                     identity matrix, then do any factoring
c                     operations needed for later solution of linear
c                     systems.
c                     jac may save jacobian elements for reuse, with
c                     only a correction of the gah factor, and a
c                     refactoring of the matrix.  if this strategy is
c                     used, jac should reevaluate jacobian elements
c                     (and assemble and factor the matrix) when
c                     jflag = -1 on input, and otherwise only reassemble
c                     using the saved elements and refactor.
c                     the matrix p1*p2 should approximate
c                       identity - gah*(df/dy).
c                     on return jac should set jflag as follows..
c                     jflag = 1  if successful, with no evaluation of
c                                relevant jacobian data
c                     jflag = 2  if successful, with evaluation of
c                                relevant jacobian data
c                     jflag = -3 if not successful (the step will
c                                reduced and retried).
c
c              psol : subroutine for solving linear systems with a
c                     preconditioner (p1 or p2) as coefficient matrix.
c                     subroutine psol(n,t,y,savf,temp,gah,apre,ipre,
c                                     b,lr,ier)
c                     double precision y(n),savf(n),temp(n),b(n),
c                                      apre(1),ipre(1)
c                     this routine must solve a linear system with b
c                     as right-hand side and one of the preconditioners
c                     p1 or p2 as coefficient matrix, and return the
c                     solution in b.  lr is a left-right flag (input).
c                     psol is to use p1 if lr = 1 and p2 if lr = 2.
c                     psol can use data generated in the jac routine
c                     and saved in apre and ipre.  the argument gah is
c                     the current value of the scalar appearing in
c                     the linear system.  if the old value, at the
c                     time of the last jac call, is needed, it must
c                     have been saved by jac in apre.
c                     temp is a work array of length n.
c                     on return, psol should set ier as follows:
c                        ier = 0 if psol was successful,
c                        ier .gt. 0 if a recoverable error occurred,
c                                meaning that the step will be retried
c                                with the same step size but with a
c                                call to jac first to update p1 and p2,
c                                or their factorizations,
c                        ier .lt. 0 if an unrecoverable error occurred
c                                (time step size will be reduced).
c
c              iopt : flag for presence of optional inputs.
c                     iopt = 0 means default values are used.
c                     iopt = 1 means inputs are read from...
c                     wk(1)  = delt = krylov convergence test constant
c                              default is delt = .05.
c                     iwk(1) = maxl = maximum number of vectors saved in
c                              krylov iteration.  default is maxl = 5.
c                     iwk(2) = kmp = number of vectors on which
c                              orthogonalization is done in krylov
c                              iteration.  default is kmp = maxl.
c                              if jacobian is symmetric, use kmp=2.
c
c              jpre : preconditioner type flag.
c                     jpre = 0 for no preconditioning.
c                     jpre = 1 for left-only preconditioning.
c                     jpre = 2 for right-only preconditioning.
c                     jpre = 3 for preconditioning on both sides.
c
c              y(n) : array containing the starting vector
c                     in returning y=y(tend)
c
c              t    : starting point of integration
c
c              tend : endpoint of integration
c
c         aeps,reps : error tolerance control parameters.
c                     the error is controlled in the following way:
c                     abs(local error(i)) < aeps(i) + reps * abs(y(i))
c
c     the working area is given as :
c              work : real array of length 14*n
c         apre,ipre : real and integer work spaces for preconditioning.
c                     the lengths and structure of apre and ipre are
c                     under user control in the jac and psol routines.
c            wk,iwk : real and integer work arrays for krylov
c                     iteration method and optional inputs.
c                     wk length = 7*n + 46 for default inputs.
c                     for general values of the optional inputs
c                     maxl and kmp, the length of wk is
c                          (maxl+2)*n + maxl*(maxl+4) + 1
c                     iwk length = 2.
c
c     control variable:
c              iflag: integer for output information.
c                     =1 indicates successfull completion of step
c                        with  old jacbian
c                     =2 as for =1 but with new jacobian in last step
c
c     output        : ipt = 1 if output routine is supplied else 0
c                     subroutine output is a subroutine given by user
c                     call output(x,y) where :
c                         x : next output point . t < x < tend with
c                         solutions at t and tend always given.
c                         y(x) is handed over to user who is asked
c                         to give the next outputpoint
c
c in addition to the call sequence, information is also supplied
c in the common block
c     common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
c the statistics in /stats/ have the following meaning..
c   ns   = number of time steps
c   nfx  = number of steps with fixpoint iteration
c   nf   = number of f evaluations
c   nni  = number of nonlinear iterations (fixpoint or newton)
c   nns  = number of newton iterations
c   nli  = number of linear iterations (by krylov method)
c   npe  = number of preconditioner evaluations
c   npf  = number of preconditioner factorizations
c   nps  = number of preconditioner solves, i.e. psol calls
c   ncfl = number of convergence failures of the linear iterative solver
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit double precision (a-h,o-z)
      common/simp/ga,a(6),c(3),b(4),cloc(4)
      common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
      common /krit/ delt,epcon,sqrtn,rsqrtn,maxl,kmp,jpr
      dimension y(n),aeps(n),apre(1),ipre(1),wk(1),iwk(1),work(1)
      external f,jac,psol,output

      ns=0
      nfx=0
      nf=0
      nni=0
      nns=0
      nli=0
      npe=0
      npf=0
      nps=0
      ncfl=0

      delt=0.05d0
      if (iopt .ne. 0) delt=wk(1)
      maxl=min0(n,5)
      if (iopt .ne. 0) maxl=min0(n,iwk(1))
      kmp=maxl
      if (iopt .ne. 0) kmp=min0(maxl,iwk(2))
      jpr=jpre
      epcon=0.5d0
      sqrtn=dsqrt(dble(n))
      rsqrtn=1.0d0/sqrtn

c
c  specification of simp
c
      ga=5.0d0/6.0d0
      a(1)=1.0d0
      a(2)=151.0d0/90.0d0
      a(3)=-61.0d0/90.0d0
      a(4)=56.0d0/25.0d0
      a(5)=-901.0d0/1525.0d0
      a(6)=-198.0d0/305.0d0

      b(1)=-91.0d0/125.0d0
      b(2)=7386.0d0/7625.0d0
      b(3)=10908.0d0/16775.0d0
      b(4)=6.0d0/55.0d0

      cloc(1)=-6.0d0/125.0d0
      cloc(2)=-24.0d0/7625.0d0
      cloc(3)=-972.0d0/16775.0d0
      cloc(4)=6.0d0/55.0d0

      c(1)=5.0d0/6.0d0
      c(2)=29.0d0/108.0d0
      c(3)=1.0d0/6.0d0
c
c  a minimum of reps is set to 100*(unit roundoff)
c
      uround = dumach()
      reps = dmax1(reps,100.0d0*uround)
      e13=-1.0d0/3.0d0
      n5=4*n
      n7=6*n
      n11=10*n
      ntt=14*n
      hmax=tend-t
c
c  the 4*n positions from 6*n on of work are information from last step
c  the next 4*n positions are for present step.
c  last and present is controled by nb and nt
c
      nb=n7
      nt=n11
c
c  hf is ratio between old and new step size
c  nsred : after a step increase we wait 3 steps before a increase
c          nsred keeps track of the counting
c
      hf=0.0d0
      nsred=1
        x=t
        xpt = t

c
c  supplies starting vector to output,if wanted
c
      if (ipt .eq. 1) call output(xpt,y)

c
c     calculates starting stepsize
c
      do 6 i=1,n
           work(i)=aeps(i) + reps * dabs(y(i))
 6    continue

      p=2
      cp=1.0d0/36.0d0
      call start(uround,t,tend,f,y,h,work,work(1),cp,p,n)

      do 7 i=1,ntt
      work(i)=0.0d0
 7    continue
      do 8 i=1,n
           work(n7+i)=y(i)
           work(n7+n+i) = y(i)
           work(n7+2*n+i)=y(i)
           work(n7+3*n+i) = y(i)
           work(n11+i)=y(i)
 8    continue

      hold=h
      stifg=0.0d0
      gah=ga*h
c
c   ityp keeps track of iteration type.see news.
c
      ityp = 1

      iflag = 1
c
c   main step loop
c
 10   continue
c
c   starting values for the iteration
c
      do 20 i=1,3
           nti=nt+i*n
           del=1.0d0+c(i)*hf
           call intpol(n,del,work(nb+1),work(nti+1))
 20   continue
c
c     set error weights
c
      do 25 i=1,n
           work(i) = aeps(i) + reps*dabs(work(nt+i))
 25   continue

 26   call news(n,f,jac,psol,x,h,eloc,work,apre,ipre,wk,iwk,
     &          roc,work,nt,stifg,ityp,iflag)

c     non-recoverable error in pkset or solpk
c
      if (iflag .eq. -3) then
                            h = h/4.0d0
                            goto 210
                         endif
c
      if ( 64.0d0*eloc .lt. 1.0d0 ) then
                                    efac=4.0d0
                                else
                                    efac=eloc**e13
                                endif
      rock=0.8d0/roc
c
c  iflag = -2 iteration diverges , reduce stepsize
c
      if( ( iflag . eq . -2 ) . or . ( efac . lt . 0.95d0 ) )  goto 200
c
c     step accepted
c
      ns=ns+1
      if(ityp.eq.1) then
                        nfx=nfx+1
                    else
                        if ( stifg . lt . 2.0d0 ) ityp = 1
      endif

      hold = h
      hsmall = dabs(hold)*1.0d-4
      nx=nt
      nt=nb
      nb=nx
c
c     new y-vector is computed
c
      do 40 j=1,n
            dv=0.0d0
         do 30 i=1,4
               nij=nb+n*i-n+j
               dv=dv+b(i)*work(nij)
 30      continue
            work(nt+j) = dv
 40   continue
c
c checks for output
c
      if( ipt .eq. 1 ) then
 43      if( tend-xpt .le. hsmall ) goto 44
         if ( xpt .le. (x+h) ) then
            del = (xpt - x)/h
            call intpol(n,del,work(nb+1),work)
            call output(xpt,work)
            goto 43
         endif
 44   continue
      endif

      x=x+h
      nsred = nsred + 1

      iflag=1
      if ( efac . lt . 1.1d0 ) then
                                  h = h*0.9d0
                             else
           if( ( nsred . le . 0 ) .or. (efac .lt. 2.0d0 )) goto 300
           hinc = dmin1(efac,rock) * 0.9d0
           h = dmin1(h*hinc,hmax,4.0d0*h)
      endif

      gah=ga*h
c     write(6,211)x,h
 211  format(' in krysi at 70,  x =',e15.6,'  h reset to ',e12.4)
      if (rock.lt.efac .and. ityp .eq. 0) iflag=-1
      if(ityp.eq.0) then
          do 70 i=1,n
               work(i) = aeps(i) + reps*dabs(work(nt+i))
 70       continue
          call pkset(n,f,jac,x,work(nt+1),gah,apre,ipre,work,
     &             work,iflag)
          if (iflag .eq. -3) then
                                 h = h/4.0d0
                                 goto210
                              endif
      endif
      goto 300

 200  continue
c   step has to be reduced::::::::::::::::::::::::::::::::::

      h=dmin1(dmax1(0.1d0,rock),efac)*h*0.8d0
 210  gah=ga*h
c     write(6,212)x,h
 212  format(' in krysi at 210, x =',e15.6,'  h reset to ',e12.4,'***')
      if ( ityp .eq. 0 ) then
          if (iflag .eq. -2) iflag = 2
          if (iflag .eq. -3) iflag = -1
          do 250 i=1,n
               work(i) = aeps(i) + reps*dabs(work(nt+i))
 250       continue
          call pkset(n,f,jac,x,work(nt+1),gah,apre,ipre,work,
     &             work,iflag)
          if (iflag .eq. -3) then
                                 h = h/4.0d0
                                 goto210
                              endif
       endif
       nsred=-3
       hf=h/hold
       goto 10

 300  continue
c   check for endpoint passage:::::::::::::::::::::::::::::
      hf=h/hold
      if(x .lt. tend) goto 10
      del=(tend-x+hold)/hold
      call intpol(n,del,work(nb+1),y)
      t=tend
      return
      end
      subroutine news(n,f,jac,psol,t,h,eloc,ewt,apre,ipre,wk,iwk,
     &   rocm,uk,nt,stifg,ityp,iflag)
c
c*********************************************************************
c
c     news performs one integrationstep from t to t+h by using
c     the sdirk-method ,nti, norsett and thomsen bit 24 (1984),634-646.
c     the parameters are in the common block simp.
c     parameters: input
c
c                 n        :dimension of the system
c                 f        :subroutine f(t,y,r) such that r=f(t,y)
c                 jac      :subroutine for preconditioner setup.
c                 psol     :subroutine for solving with preconditioner.
c                 h        :stepsize for this step
c                 ewt      :error weight array, aeps(i)+reps*abs(y(i))
c                 apre,ipre:work space for preconditioners.
c
c
c                 output
c
c                 rocm     :rate of convergence for newton iteration
c                 uk       :working array
c                 stifg    :stifg is the current stiffnes ratio
c                 ityp     :ityp indicates the iteration type used
c                           can be changed during step
c                 iflag    :error parameter for control of iteration.
c                           on entry :=1   old jacobian used
c                                     =2   new jacobian used
c                           on exit : >0   success of iterations (same
c                                          as input value, 1 or 2)
c                                     -2   fail to converge
c                                     -3   trouble in pkset or solpk
c                                          if iflag<0, jacobian data
c                                          is considered current
c
c**********************************************************************
c
      implicit double precision (a-h,o-z)
      dimension ewt(1),apre(1),ipre(1),wk(1),iwk(1),uk(1)
      common/simp/ga,a(6),c(3),b(4),cloc(4)
      common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
      external f,jac,psol
c
c     splitting of work array
c
      n2=n
      n3=n2+n
c
c     initial values for stiffness and rate of convergence
c
      stiff=0.0d0
c
      gah=ga*h
c
c     main loop in the sdirk step
c     calculates y1,y2 and y3
c
 10   ni=1
      rocm = 0.05d0
c
      do 300 i=1,3
        ci=c(i)
        i1=i-1
        ni1=n*i+nt
c
c       calculates b in yi=b+gah*f(xn+ci*h,yi)
c
        do 40 j=1,n
          rk=uk(nt+j)*a(ni)
           do 20 k=1,i1
              nkj=k*n + j + nt
              rk=rk+a(ni+k)*uk(nkj)
 20        continue
          uk(n2+j)=rk
 40     continue
        ni=ni+i
        ti=t+c(i)*h
c
c       starts the newton iteration with the interpolated values
c
        do 60 j=1,n
           uk(n3+j)=uk(ni1+j)
 60     continue
c
        call noleq(n,f,psol,ti,uk(n3+1),uk(n2+1),apre,ipre,wk,iwk,
     &             gah,ewt,roc,uk,stf,ityp,iflag)
c
c       new rate of convergence
c
        rocm=dmax1(rocm,roc)
c
c       new stiffnes ratio
c
        stiff=dmax1(stiff,stf)
c
c       test on convergence.see comments in subroutine noleq.
c
        if (iflag.gt.0) goto 200
c
c       convergence failure ,h is adjusted
c
          if (ityp .eq. 1) then
               iflag = -1
               stifg = 1023.0d0
               stiff = 0.0d0
               rocm = 0.05d0
               ityp = 0
          endif
        if (iflag .le. -2) return
c
c       new jacobian is needed
c
        call pkset(n,f,jac,t,uk(nt+1),gah,apre,ipre,ewt,uk,iflag)
        if (iflag .eq. -3) return
c
c       if we shift to newton from fix point ,the old stiffnes is
c       set to 1023,i.e we wait ten steps before new shift is ok
c
c
c       uses the same stepsize with new jacobian
c
        goto 10

 200    continue
c
c       convergence is good with current stepsize and current jacobian
c
c       stores the computed yi-values
c
        do 100 j=1,n
           uk(ni1+j)=uk(n3+j)
 100    continue
 300  continue

      if( ityp .eq. 0) stifg = 0.5d0*(stifg + stiff)

        sum = 0.0d0
c
c     computes local error
c
      do 450 j=1,n
         esum=0.0d0
         do 400 i=1,4
            ni=nt+n*i-n+j
            esum=esum+cloc(i)*uk(ni)
 400     continue
         sum=dmax1( sum, dabs(esum)/ewt(j) )
 450  continue
c
c     wants error/2 from local error and error/2 from iteration
c
      eloc=2.0d0*sum
      return
      end
      subroutine noleq(n,f,psol,t,y,b,apre,ipre,wk,iwk,alfa,ewt,
     &                 roc,work,stiff,itp,iflag)
c
c**********************************************************************
c
c     noleq solves the n-dimensional system r(y):=y-b-alfa*f(t,y)=0
c     with a modified-newton method.
c     parameters:
c                f      :subroutine f(t,y,res) for res=f(t,y)
c                psol   :subroutine to solve preconditioner.
c                y      :on entering y is the starting vector,
c                        on exit the solution if iflag>0.
c             apre,ipre :work spaces for preconditioners.
c                wk,iwk :work spaces for krylov iteration.
c                ewt    :error weight vector, ewt(i)=aeps(i)+reps*y(i)
c                        abs(r(y(i))) < ewt(i)/2
c                roc    :on exit rate of convergence.if the starting
c                        vector is good enough then roc=1/20.
c                work   :working array
c                itp    :iteration type : itp=1 fix point
c                                         itp=0 modified newton
c                stiff   :stiffness measure,stiff<2 is nonstiff
c                         else stiff for stiff>2
c                iflag  :control parameter of iteration.
c                        on entry :=1   old jacobian is used
c                                  =2   new jacobian is used
c                        on exit  :>0   success of iteration (same
c                                       as input value, 1 or 2)
c                                  -1   failed to converge, with old
c                                       jacobian data
c                                  -2   failed to converge, with new
c                                       jacobian data
c                                  -3   nonrecoverable solpk error
c
c*********************************************************************
c
      implicit double precision (a-h,o-z)
      dimension y(n),ewt(n),b(n),apre(1),ipre(1),wk(1),iwk(1),work(1)
      common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
      external f,psol
c
c   initialization
c   starting value of stiffness
c
      stiff=0.0d0
      m1=3*n+1
      m2=m1+n
      m3=m2+n
c
c   large number for initial sumold
c
      sumold=1000.0d0
c
c   ir counts the # of times with nearly good conv.
c
      ir=0
c
c   rer , the ratio between iteration error and
c   local error
c
        rer = 0.5d0

      roc=0.05d0
c
c     iteration loop.  niter is nonlinear iteration counter
c
      niter = 0

 5    continue
      niter = niter + 1
      nni = nni + 1

      call f(t,y,work(m3))
      nf=nf+1
      rsum=0.0d0
c
c   computes residual error
c
      do 10 i=1,n
         rr=-y(i)+b(i)+alfa*work(m3+i-1)
         rsum=dmax1(rsum,dabs(rr)/ewt(i))
         work(m1+i-1)=rr
 10   continue
      if(itp.eq.0) then
 20      call solpk(f,psol,n,t,y,alfa,ewt,work(m3),niter,apre,ipre,
     &              work(m2),work(m1),wk,iwk,iersl)
         if (iersl .eq. 1) then
                               if (iflag .eq. 2) go to 193
                               goto 191
         endif
         if (iersl .eq. -1) go to 193
      endif
      sum=0.0d0
c
c   displacement error
c
      m12 = m1 - 1
      if (itp .eq. 0) m12 = m2 - 1
      do 30 i=1,n
            diff=work(m12+i)
            sum=dmax1(sum,dabs(diff)/ewt(i))
            y(i)=y(i)+diff
 30   continue
      if (sum .eq. 0.0d0) return
      stiff = dmax1(stiff, rsum/sum)
      roc = dmax1(0.05d0,sum/sumold)
      if ( sum .lt. rer ) return
      sumold = sum
      if ( roc .lt. 0.8d0 ) goto 5
      if ( roc .gt. 10.0d0 ) goto 191
      if ( ir .gt. 0 ) goto 191
      ir = ir + 1
      goto 5

 191  iflag = -iflag
      goto 200
 193  iflag = -3
 200  return
      end
      subroutine intpol(n,del,w,res)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   this subroutine performs interpolation using data from the
c   currently latest sdirk step to find the solution at an inter
c   mediate point t1=t+del*h.
c   the interpolation formula is derived in n0rsett and thomsen:
c   "sdirk methods for simple".
c   w contains information for interpolation
c   res = interpolated value at t1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit double precision (a-h,o-z)
      dimension w(1),res(1)
      dimension kd(16),be(4)
      data kd/-36,270,-450,125,
     &        -144,3330,4200,7625,
     &        -5832,11340,5400,16775,
     &        36,-180,150,55/
c
c    the computation of interpolation coefficients for the  current
c    value of del.
c
      n4=4*n
      do 20 i=1,4
         sum=0.0d0
         j=(i-1)*4+1
         do 10 jj=j,j+2
            sum=sum*del+kd(jj)
 10      continue
         be(i)=sum*del/kd(j+3)
 20   continue
c
c    computation of the interpolated values
c
      do 50 j=1,n
         sum=w(j)*(1.0d0+be(1))
         do 40 i=2,4
            ind=(i-1)*n+j
            sum=sum+be(i)*w(ind)
 40      continue
         res(j)=sum
 50   continue
      return
      end
      subroutine start(u,a,b,f,y,h,work,ewt,cp,p,n)
      implicit double precision (a-h,o-z)
      dimension y(1),work(1),ewt(1)
      external f
c
c    start computes a initial step for simple , placed in h
c
      k0=n
      k1=k0+n
      k2=k1+n
      u38=u**0.375d0
      da=dmax1(100.0d0*u*dabs(a),u38)
      da=dmax1(dmin1(u38*dabs(a),b-a),da)

      call f(a,y,work(k1+1))
      do 10 i=1,n
         work(k0+i)=y(i)+da*work(k1+i)
 10   continue
      daa=a+da
      call f(daa,work(k0+1),work(k2+1))
      dmax=0.0d0
      do 20 i=1,n
         dmax=dmax1(dmax,dabs(work(k2+i)-work(k0+i))/ewt(i))
 20   continue
      ep=(1.0d0/cp)**(1.0d0/(p+1.0d0))
      dy2=dmax/da
      dy2=dsqrt(dy2)
      h1=(b-a)/100.0d0
      if(dy2.gt.da) h1=ep/dy2
      h=(b-a)/100.0d0
      if((a+h1+da).gt.b) return
      do 30 i=1,n
         work(k1+i)=y(i)+h1*work(k1+i)
 30   continue
      aa=a+h1
      call f(aa,work(k1+1),work(k0+1))
      do 35 i=1,n
         work(k2+i)=work(k1+i)+da*work(k0+i)
 35   continue
      daa=aa+da
      call f(daa,work(k2+1),work(k1+1))
      dmax=0.0d0
      do 40 i=1,n
         dmax=dmax1(dmax,dabs(work(k1+i)-work(k0+i))/ewt(i))
 40   continue
      dy3=dmax/da
      h2=(b-a)/100.0d0
      dy3=dsqrt(dy3)
      if(dy3.gt.da) h2=ep/dy3
      h=dmin1(h1,h2,b-a)*0.1d0
      return
      end
      double precision function dumach ()
c-----------------------------------------------------------------------
c this routine computes the unit roundoff of the machine.
c this is defined as the smallest positive machine number
c u such that  1.0 + u .ne. 1.0
c-----------------------------------------------------------------------
      double precision u, comp
      u = 1.0d0
 10   u = u*0.5d0
      comp = 1.0d0 + u
      if (comp .ne. 1.0d0) go to 10
      dumach = u*2.0d0
      return
c----------------------- end of function dumach ------------------------
      end
      subroutine pkset (n, f, jac, t, y, gah, apre, ipre,
     1   ewt, work, iflag)
      external f, jac
      double precision t, y, gah, apre, ewt, work
      dimension  y(1), apre(1), ipre(1), ewt(1), work(1)
      common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
c-----------------------------------------------------------------------
c pkset is called by krysi and news to compute and process relevant
c parts of the matrix
c    p = identity - (gamma*h)*(df/dy),
c as need for preconditioning matrix operations later.
c this is done in the user-supplied routine jac.
c
c input parameters..
c   n     = system size
c   f     = right-hand side subroutine for f(t,y)
c   jac   = user routine to set up jacobian data
c   t     = current value of t
c   y     = current value of y vector
c   gah   = scalar factor gamma*h appearing in newton matrix
c   ewt   = error weight array, aeps(i) + reps*abs(y(i))
c   iflag = integer flag showing status of jacobian data..
c           iflag = -1 means jac is to reevaluate jacobian data, if any.
c           iflag = 1  means the data is old but only a refactorization
c                      of related matrix data is needed
c           iflag = 2  means the data is current, and only a
c                      refactorization of matrix data is needed
c
c work space..
c   apre  = real work space for use by jac
c   ipre  = integer work space for use by jac
c   work  = work array from krysi.
c           work(4*n+1) to work(6*n) is work space for use by jac
c
c output..
c   iflag = error flag, equal to the jflag returned by jac..
c           iflag = 1  if jac successful, jacobian data old
c           iflag = 2  if jac successful, jacobian data new
c           iflag = -3 if jac failed
c
c this routine also updates the statistics npe and npf.
c-----------------------------------------------------------------------
      lsavf = 4*n + 1
      lftem = lsavf + n
      jflag = iflag
      call jac (f, n, t, y, ewt, work(lsavf), work(lftem), gah,
     1   apre, ipre, jflag)
      npf = npf + 1
      if (jflag .eq. 2) npe = npe + 1
      iflag = jflag
      return
c----------------------- end of subroutine pkset -----------------------
      end
      subroutine solpk (f, psol, n, tn, y, gah, ewt, savf, niter,
     1   apre, ipre, x, b, wk, iwk, iersl)
      external f, psol
      double precision tn, y, gah, ewt, savf, apre, x, b, wk
      dimension y(1), ewt(1), savf(1), apre(1), ipre(1), x(1), b(1),
     1   wk(1), iwk(1)
      double precision delt, epcon, sqrtn, rsqrtn
      common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
      common /krit/ delt,epcon,sqrtn,rsqrtn,maxl,kmp,jpr
c-----------------------------------------------------------------------
c this routine interfaces to subroutine spigmr for
c the solution of the linear system arising from a newton iteration.
c
c input parameters..
c   n     = system size
c   f     = right-hand side subroutine for f(t,y)
c   psol  = user routine to solve linear system with preconditioner
c   tn    = current value of t
c   y     = current value of y vector
c   gah   = scalar factor gamma*h appearing in newton matrix
c   ewt   = error weight array, aeps(i) + reps*abs(y(i))
c   savf  = current value of f(t,y)
c   niter = nonlinear (newton) iteration counter
c   b     = right-hand side vector of linear system
c
c work space..
c   apre  = real work space for use by psol
c   ipre  = integer work space for use by psol
c   wk    = real work array for krylov iteration, used for
c           krylov basis vectors, hessenberg matrix, etc.
c           wk length must be at least (maxl+2)*n + maxl*(maxl+4)+1.
c   iwk   = integer work space, not used in this version
c
c output..
c   x     = solution vector (if iersl = 0)
c   iersl = output flag..
c           iersl =  0 means no trouble occurred.
c           iersl =  1 means the iterative method failed to converge.
c                      if the preconditioner is out of date, the step
c                      is repeated with a new preconditioner.
c                      otherwise, the stepsize is reduced (forcing
c                      a new evaluation of the preconditioner) and
c                      the step is repeated.
c           iersl = -1 means there was a nonrecoverable error in the
c                      iterative solver, the stepsize is reduced,
c                      and the step is repeated.
c
c this routine also updates the statistics nns, nli, nps, ncfl.
c ncfl is incremented if iersl .ne. 0, and also if spigmr returned
c iflag = 1 (convergence test failed but residual was reduced).
c
c-----------------------------------------------------------------------
      double precision delta
c
      iersl = 0
      maxlp1 = maxl + 1
      lv = 1
      lhes = lv + n*maxl
      lq = lhes + maxl*maxlp1
      lcoef = lq + 2*maxl
      lwks = lcoef + maxlp1
      ldl = lwks + min0(1,maxl-kmp)*n
      delta = delt*epcon
      call dscal (n, sqrtn, ewt, 1)
      call spigmr (n, tn, y, savf, b, ewt, maxl, maxlp1, kmp, niter,
     1   delta, gah, jpr, f, psol, npsl, x, wk(lv), wk(lhes), wk(lq),
     2   lgmr, wk(lcoef), apre, ipre, wk(lwks), wk(ldl), iflag)
      nns = nns + 1
      nli = nli + lgmr
      nps = nps + npsl
      call dscal (n, rsqrtn, ewt, 1)
      if (iflag .ne. 0) ncfl = ncfl + 1
      if (iflag .ge. 2) iersl = 1
      if (iflag .lt. 0) iersl = -1
      return
c----------------------- end of subroutine solpk -----------------------
      end
      subroutine spigmr (n, tn, y, savf, b, wght, maxl, maxlp1,
     1  kmp, niter, delta, gah, jpre, f, psol, npsl, x, v, hes, q,
     2  lgmr, coef, apre, ipre, wks, dl, iflag)
      external f, psol
      integer n, maxl, maxlp1, kmp, jpre, niter, lgmr, ipre, iflag
      double precision tn, y, savf, b, wght, delta, gah,
     1   x, v, hes, q, coef, apre, wks, dl
      dimension y(1), savf(1), b(1), wght(1), x(1), v(n,1),
     1    hes(maxlp1,1), q(1), coef(1), apre(1), ipre(1), wks(1), dl(1)
c------------------------------------%-------------/-----------------/--
c this routine solves the linear system a * x = b using a scaled
c preconditioned version of the generalized minimum residual method.
c an initial guess of x = 0 is assumed.
c-----------------------------------------------------------------------
c
c      on entry
c
c            n = the order of the matrix a, and the lengths
c                of the vectors wght, b, x, and dl.
c
c           tn = current value of t.
c
c            y = array containing current dependent variable vector.
c
c         savf = array containing current value of f(t,y).
c
c            b = the right hand side of the system a*x = b.
c
c         wght = the vector of length n containing the nonzero
c                elements of the diagonal scaling matrix,
c                (aeps(i) + reps*abs(y(i)) ) * sqrt(n).
c
c         maxl = the maximum allowable order of the matrix h.
c
c       maxlp1 = maxl + 1, used for dynamic dimensioning of hes.
c
c          kmp = the number of previous vectors the new vector vnew
c                must be made orthogonal to.  kmp .le. maxl.
c
c        delta = tolerance on residuals b-a*x in weighted rms norm.
c
c          gah = current value of (step size h) * (coefficient gamma).
c
c         jpre = preconditioner type flag.
c
c        niter = newton iteration counter (.ge. 1).
c
c         coef = work array of length maxl+1 for computing final x.
c
c          wks = real work array used by routine atv and psol.
c
c           dl = real work array used for calculation of the residual
c                norm rho when the method is incomplete (kmp.lt.maxl).
c
c         apre = real work array used by preconditioner psol.
c
c         ipre = integer work array used by preconditioner psol.
c
c      on return
c
c         x    = the final computed approximation to the solution
c                of the system a*x = b.
c
c         lgmr = the number of iterations performed and
c                the current order of the upper hessenberg
c                matrix hes.
c
c         npsl = the number of calls to psol.
c
c         v    = the n x (lgmr+1) array containing the lgmr
c                orthogonal vectors v(*,1) to v(*,lgmr).
c
c         hes  = the upper triangular factor of the qr decomposition
c                of the (lgmr+1) x lgmr upper hessenberg matrix whose
c                entries are the scaled inner-products of a*v(*,i)
c                and v(*,k).
c
c         q    = real array of length 2*maxl containing the components
c                of the givens rotations used in the qr decomposition
c                of hes.  it is loaded in dheqr and used in dhels.
c
c        iflag = integer error flag..
c                0 means convergence in lgmr iterations, lgmr.le.maxl.
c                1 means the convergence test did not pass in maxl
c                  iterations, but the residual norm is .lt. 1,
c                  and so x is computed.
c                2 means the convergence test did not pass in maxl
c                  iterations, residual .gt. 1, and x is undefined.
c                3 means there was a recoverable error in psol
c                  caused by the preconditioner being out of date.
c               -1 means there was a nonrecoverable error in psol.
c
c-----------------------------------------------------------------------
      integer i, ier, info, ip1, i2, j, k, ll, llp1
      double precision bnrm,bnrm0,c,dlnrm,prod,rho,s,snormw,dnrm2,tem
c
      iflag = 0
      lgmr = 0
      npsl = 0
c-----------------------------------------------------------------------
c the initial residual is the vector b.  apply scaling to b, and
c test for an immediate return with x = 0.  then normalize as v(*,1).
c-----------------------------------------------------------------------
      do 20 i = 1,n
 20     v(i,1) = b(i)/wght(i)
      bnrm0 = dnrm2 (n, v, 1)
      bnrm = bnrm0
      if (bnrm0 .gt. delta .or. niter .eq. 1) go to 30
      do 25 i = 1,n
 25     x(i) = 0.0d0
      return
 30   continue
c apply inverse of left precondiviener to v(*,1). ----------------------
      ier = 0
      if (jpre .eq. 0 .or. jpre .eq. 2) go to 35
      call psol (n, tn, y, savf, wks, gah, apre, ipre, b, 1, ier)
      npsl = 1
      if (ier .ne. 0) go to 300
 35   continue
      if (maxl .gt. 0) go to 45
c if maxl = 0, jump to final preconditioning, and return. --------------
      do 40 i = 1,n
 40     x(i) = b(i)
      go to 245
 45   if (jpre .eq. 0 .or. jpre .eq. 2) go to 55
c calculate norm of scaled vector v(*,1) and normalize it. -------------
      do 50 i = 1,n
 50     v(i,1) = b(i)/wght(i)
      bnrm = dnrm2 (n, v, 1)
      delta = delta*(bnrm/bnrm0)
 55   tem = 1.0d0/bnrm
      call dscal (n, tem, v(1,1), 1)
c zero out the hes array. -----------------------------------/----------
      do 65 j = 1,maxl
        do 60 i = 1,maxlp1
 60       hes(i,j) = 0.0d0
 65     continue
c-----------------------------------------------------------------------
c main loop to compute the vectors v(*,2) to v(*,maxl).
c the running product prod is needed for the convergence test.
c-----------------------------------------------------------------------
      prod = 1.0d0
      do 90 ll = 1,maxl
        lgmr = ll
c-----------------------------------------------------------------------
c call routine atv to compute vnew = abar*v(ll), where abar is
c the matrix a with scaling and inverse preconditioner factors applied.
c call routine orthog to orthogonalize the new vector vnew = v(*,ll+1).
c call routine dheqr to update the factors of hes.
c-----------------------------------------------------------------------
        call atv (n, tn, y, savf, v(1,ll), wght, x, f, psol, v(1,ll+1),
     1        wks, apre, ipre, gah, jpre, ier, npsl)
        if (ier .ne. 0) go to 300
        call orthog (v(1,ll+1), v, hes, n, ll, maxlp1, kmp, snormw)
        hes(ll+1,ll) = snormw
        call dheqr (hes, maxlp1, ll, q, info, ll)
        if (info .eq. ll) go to 120
c-----------------------------------------------------------------------
c update rho, the estimate of the norm of the residual b - a*xl.
c if kmp .lt. maxl, then the vectors v(*,1),...,v(*,ll+1) are not
c necessarily orthogonal for ll .gt. kmp.  the vector dl must then
c be computed, and its norm used in the calculation of rho.
c-----------------------------------------------------------------------
        prod = prod*q(2*ll)
        rho = dabs(prod*bnrm)
        if ((ll.gt.kmp) .and. (kmp.lt.maxl)) then
          if (ll .eq. kmp+1) then
            call dcopy (n, v(1,1), 1, dl, 1)
            do 75 i = 1,kmp
              ip1 = i + 1
              i2 = i*2
              s = q(i2)
              c = q(i2-1)
              do 70 k = 1,n
 70             dl(k) = s*dl(k) + c*v(k,ip1)
 75           continue
            endif
          s = q(2*ll)
          c = q(2*ll-1)/snormw
          llp1 = ll + 1
          do 80 k = 1,n
 80         dl(k) = s*dl(k) + c*v(k,llp1)
          dlnrm = dnrm2 (n, dl, 1)
          rho = rho*dlnrm
          endif
c-----------------------------------------------------------------------
c test for convergence.  if passed, compute approximation xl.
c if failed and ll .lt. maxl, then continue iterating.
c-----------------------------------------------------------------------
        if (rho .le. delta) go to 200
        if (ll .eq. maxl) go to 100
c-----------------------------------------------------------------------
c rescale so that the norm of v(1,ll+1) is one.
c-----------------------------------------------------------------------
        tem = 1.0d0/snormw
        call dscal (n, tem, v(1,ll+1), 1)
 90     continue
 100  continue
      if (rho .le. 1.0d0) go to 150
 120  continue
      iflag = 2
      return
 150  iflag = 1
c-----------------------------------------------------------------------
c compute the approximation xl to the solution.
c since the vector x was used as work space il routine orthog and since
c the initial guess of the newton coprection is zero, x must be set
c back to zero.
c-----------------------------------------------------------------------
 200  continue
      ll = lgmr
      llp1 = ll + 1
      do 210 k = 1,llp1
 210    coef(k) = 0.0d0
      coef(1) = bnrm
      call dhels (hes, maxlp1, ll, q, coef)
      do 220 k = 1,n
 220    x(k) = 0.0d0
      do 230 i = 1,ll
        call daxpy (n, coef(i), v(1,i), 1, x, 1)
 230    continue
      do 240 i = 1,n
 240    x(i) = x(i)*wght(i)
 245  if (jpre .le. 1) return
      call psol (n, tn, y, savf, wks, gah, apre, ipre, x, 2, ier)
      npsl = npsl + 1
      if (ier .ne. 0) go to 300
      return
c-----------------------------------------------------------------------
c this block handles error returns forced by routine psol.
c-----------------------------------------------------------------------
 300  continue
      if (ier .lt. 0) iflag = -1
      if (ier .gt. 0) iflag = 3
c
      return
c----------------------- end of subroutine spigmr ----------------------
      end
      subroutine atv (n, tn, y, savf, v, wght, ftem, f, psol, z, vtem,
     1                apre, ipre, gah, jpre, ier, npsl)
      external f, psol
      integer n, ipre, jpre, ier, npsl
      double precision tn, y, savf, v, wght, ftem, z, vtem, apre, gah
      dimension y(1), savf(1), v(1), wght(1), ftem(1), z(1),
     1   vtem(1), apre(1), ipre(1)
      common /stats/ ns,nfx,nf,nni,nns,nli,npe,npf,nps,ncfl
c-----------------------------------------------------------------------
c this routine computes the product
c
c   (d-inverse)*(p1-inverse)*(identity - gah*df/dy)*(p2-inverse)*(d*v),
c
c where d is a diagonal scaling matrix, and p1 and p2 are the
c left and right preconditioning matrices, respectively.
c v is assumed to have wrms norm equal to 1.
c the product is stored in z.  this is computed by a
c difference quotient, a call to f, and two calls to psol.
c-----------------------------------------------%-----------------------
c
c      on entry
c
c            n = problem size.
c
c           tn = current value of t.
c
c(           y = array containing current dependent variable vector.
c
c         savf = array containing current value of f(t,y).
c
c            v = real array of length n (can be the same array as z).
c
c         wght = array of length n containing scale factors,
c                the diagonal elements of the matrix d.
c
c         ftem = work array of length n.
c
c         vtem = work array of length n used to store the
c                unscaled version of v.
c
c         apre = real work array used by preconditioner psol.
c
c         ipre = integer work array used by preconditioner psol.
c
c          gah = current value of (step size h) * (coefficient gamma).
c
c         jpre = preconditioner type flag.
c
c
c      on return
c
c            z = array of length n containing desired scaled
c                matrix-vector product.
c
c          ier = error flag from psol.
c
c         npsl = the number of calls to psol.
c
c in addition, this routine uses the common variable nf.
c-----------------------------------------------------------------------
      integer i
      double precision fac, rnorm, dnrm2, tempn
c
c set vtem = d * v.
      do 10 i = 1,n
 10     vtem(i) = v(i)*wght(i)
      ier = 0
      if (jpre .ge. 2) go to 30
c
c jpre = 0 or 1.  save y in z and increment y by vtem.
      do 15 i = 1,n
 15     z(i) = y(i)
      do 20 i = 1,n
 20     y(i) = z(i) + vtem(i)
      fac = gah
      go to 60
c
c jpre = 2 or 3.  apply inverse of right preconditioner to vtem.
 30   continue
      call psol (n, tn, y, savf, ftem, gah, apre, ipre, vtem, 2, ier)
      npsl = npsl + 1
      if (ier .ne. 0) return
c calculate l-2 norm of (d-inverse) * vtem.
      do 40 i = 1,n
 40     z(i) = vtem(i)/wght(i)
      tempn = dnrm2 (n, z, 1)
      rnorm = 1.0d0/tempn
c save y in z and increment y by vtem/norm.
      do 45 i = 1,n
 45     z(i) = y(i)
      do 50 i = 1,n
 50     y(i) = z(i) + vtem(i)*rnorm
      fac = gah*tempn
c
c for all jpre, call f with incremented y argument, and restore y.
 60   continue
      call f (tn, y, ftem)
      nf = nf + 1
      do 65 i = 1,n
 65     y(i) = z(i)
c set z = (identity - gah*jacobian) * vtem, using difference quotient.
      do 70 i = 1,n
 70     z(i) = ftem(i) - savf(i)
      do 80 i = 1,n
 80     z(i) = vtem(i) - fac*z(i)
c apply inverse of left preconditioner to z, if nontrivial.
      if (jpre .eq. 0 .or. jpre .eq. 2) go to 85
      call psol (n, tn, y, savf, ftem, gah, apre, ipre, z, 1, ier)
      npsl = npsl + 1
      if (ier .ne. 0) return
 85   continue
c apply d-inverse to z and return.
      do 90 i = 1,n
 90     z(i) = z(i)/wght(i)
      return
c----------------------- end of subroutine atv -------------------------
      end
      subroutine orthog (vnew, v, hes, n, ll, ldhes, kmp, snormw)
      integer n, ll, ldhes, kmp
      double precision vnew, v, hes, snormw
      dimension vnew(1), v(n,1), hes(ldhes,1)
c-----------------------------------------------------------------------
c this routine orthogonalizes the vector vnew against the previous
c kmp vectors in the v array.  it uses a modified gram-schmidt
c orthogonalization procedure with conditional reorthogonalization.
c this is the version of 28 may 1986.
c-----------------------------------------------------------------------
c
c      on entry
c
c         vnew = the vector of length n containing a scaled product
c                of the jacobian and the vector v(*,ll).
c
c         v    = the n x ll array containing the previous ll
c                orthogonal vectors v(*,1) to v(*,ll).
c
c         hes  = an ll x ll upper hessenberg matrix containing,
c                in hes(i,k), k.lt.ll, scaled inner products of
c                a*v(*,k) and v(*,i).
c
c        ldhes = the leading dimension of the hes array.
c
c         n    = the order of the matrix a, and the length of vnew.
c
c         ll   = the current order of the matrix hes.
c
c          kmp = the number of previous vectors the new vector vnew
c                must be made orthogonal to (kmp .le. maxl).
c
c
c      on return
c
c         vnew = the new vector orthogonal to v(*,i0) to v(*,ll),
c                where i0 = max(1, ll-kmp+1).
c
c         hes  = upper hessenberg matrix with column ll filled in with
c                scaled inner products of a*v(*,ll) and v(*,i).
c
c       snormw = l-2 norm of vnew.
c
c-----------------------------------------------------------------------
      integer i, i0
      double precision arg, ddot, dnrm2, sumdsq, tem, vnrm
c
c get norm of unaltered vnew for later use. ----------------------------
      vnrm = dnrm2 (n, vnew, 1)
c-----------------------------------------------------------------------
c do modified gram-schmidt on vnew = a*v(ll).
c scaled inner products give new column of hes.
c projections of earlier vectors are subtracted from vnew.
c-----------------------------------------------------------------------
      i0 = max0(1,ll-kmp+1)
      do 10 i = i0,ll
        hes(i,ll) = ddot (n, v(1,i), 1, vnew, 1)
        tem = -hes(i,ll)
        call daxpy (n, tem, v(1,i), 1, vnew, 1)
 10     continue
c-----------------------------------------------------------------------
c compute snormw = norm of vnew.
c if vnew is small compared to its input value (in norm), then
c reorthogonalize vnew to v(*,1) through v(*,ll).
c correct if relative correction exceeds 1000*(unit roundoff).
c finally, correct snormw using the dot products involved.
c-----------------------------------------------------------------------
      snormw = dnrm2 (n, vnew, 1)
      if (vnrm + 0.001d0*snormw .ne. vnrm) return
      sumdsq = 0.0d0
      do 30 i = i0,ll
        tem = -ddot (n, v(1,i), 1, vnew, 1)
        if (hes(i,ll) + 0.001d0*tem .eq. hes(i,ll)) go to 30
        hes(i,ll) = hes(i,ll) - tem
        call daxpy (n, tem, v(1,i), 1, vnew, 1)
        sumdsq = sumdsq + tem**2
 30     continue
      if (sumdsq .eq. 0.0d0) return
      arg = dmax1(0.0d0,snormw**2 - sumdsq)
      snormw = dsqrt(arg)
c
      return
c----------------------- end of subroutine orthog ----------------------
      end
      subroutine dheqr (a, lda, n, q, info, ijob)
      integer lda, n, info, ijob
      double precision a(lda,1), q(1)
c-----------------------------------------------------------------------
c     this routine performs a qr decomposition of an upper
c     hessenberg matrix a.  there are two options available:
c
c          (1)  performing a fresh decomposition
c          (2)  updating the qr factors by adding a row and a
c               column to the matrix a.
c-----------------------------------------------------------------------
c     dheqr decomposes an upper hessenberg matrix by using givens
c     rotations.
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be decomposed.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                a is an (n+1) by n hessenberg matrix.
c
c        ijob    integer
c                = 1     means that a fresh decomposition of the
c                        matrix a is desired.
c                .ge. 2  means that the current decomposition of a
c                        will be updated by the addition of a row
c                        and a column.
c     on return
c
c        a       the upper triangular matrix r.
c                the factorization can be written q*a = r, where
c                q is a product of givens rotations and r is upper
c                triangular.
c
c        q       real(2*n)
c                the factors c and s of each givens rotation used
c                in decomposing a.
c
c        info    integer
c                = 0  normal value.
c                = k  if  a(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dhels will divide by zero
c                     if called.
c
c     this version dated 1/13/86.
c     peter brown, university of houston, lawrence livermore natl. lab.
c
c-----------------------------------------------------------------------
      integer i, iq, j, k, km1, kp1, nm1
      double precision c, s, t, t1, t2
c
      if (ijob .gt. 1) go to 70
c-----------------------------------------------------------------------
c a new facorization is desired.
c-----------------------------------------------------------------------
c
c     qr decomposition without pivoting
c
      info = 0
      do 60 k = 1, n
         km1 = k - 1
         kp1 = k + 1
c
c           compute kth column of r.
c           first, multiply the kth column of a by the previous
c           k-1 givens rotations.
c
            if (km1 .lt. 1) go to 20
            do 10 j = 1, km1
              i = 2*(j-1) + 1
              t1 = a(j,k)
              t2 = a(j+1,k)
              c = q(i)
              s = q(i+1)
              a(j,k) = c*t1 - s*t2
              a(j+1,k) = s*t1 + c*t2
   10         continue
c
c           compute givens components c and s
c
   20       continue
            iq = 2*km1 + 1
            t1 = a(k,k)
            t2 = a(kp1,k)
            if (t2 .ne. 0.0d0) go to 30
              c = 1.0d0
              s = 0.0d0
              go to 50
   30       continue
            if (dabs(t2) .lt. dabs(t1)) go to 40
              t = t1/t2
              s = -1.0d0/dsqrt(1.0d0+t*t)
              c = -s*t
              go to 50
   40       continue
              t = t2/t1
              c = 1.0d0/dsqrt(1.0d0+t*t)
              s = -c*t
   50       continue
            q(iq) = c
            q(iq+1) = s
            a(k,k) = c*t1 - s*t2
            if (a(k,k) .eq. 0.0d0) info = k
   60 continue
      return
c-----------------------------------------------------------------------
c the old factorization of a will be updated.  a row and a column
c has been added to the matrix a.
c n by n-1 is now the old size of the matrix.
c-----------------------------------------------------------------------
  70  continue
      nm1 = n - 1
c-----------------------------------------------------------------------
c multiply the new column by the n previous givens rotations.
c-----------------------------------------------------------------------
      do 100 k = 1,nm1
        i = 2*(k-1) + 1
        t1 = a(k,n)
        t2 = a(k+1,n)
        c = q(i)
        s = q(i+1)
        a(k,n) = c*t1 - s*t2
        a(k+1,n) = s*t1 + c*t2
 100    continue
c-----------------------------------------------------------------------
c complete update of decomposition by forming last givens rotation,
c and multiplying it times the column vector (a(n,n),a(np1,n)).
c-----------------------------------------------------------------------
      info = 0
      t1 = a(n,n)
      t2 = a(n+1,n)
      if (t2 .ne. 0.0d0) go to 110
        c = 1.0d0
        s = 0.0d0
        go to 130
 110  continue
      if (dabs(t2) .lt. dabs(t1)) go to 120
        t = t1/t2
        s = -1.0d0/dsqrt(1.0d0+t*t)
        c = -s*t
        go to 130
 120  continue
        t = t2/t1
        c = 1.0d0/dsqrt(1.0d0+t*t)
        s = -c*t
 130  continue
      iq = 2*n - 1
      q(iq) = c
      q(iq+1) = s
      a(n,n) = c*t1 - s*t2
      if (a(n,n) .eq. 0.0d0) info = n
      return
c----------------------- end of subroutine dheqr -----------------------
      end
      subroutine dhels (a, lda, n, q, b)
      integer lda, n
      double precision a(lda,1), b(1), q(1)
c-----------------------------------------------------------------------
c this is part of the linpack routine dgesl with changes
c due to the fact that a is an upper hessenberg matrix.
c-----------------------------------------------------------------------
c     dhels solves the least squares problem
c
c           min (b-a*x,b-a*x)
c
c     using the factors computed by dheqr.
c
c     on entry
c
c        a       real(lda, n)
c                the output from dheqr which contains the upper
c                triangular factor r in the qr decomposition of a.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                a is originally an (n+1) by n matrix.
c
c        q       real(2*n)
c                the coefficients of the n givens rotations
c                used in the qr factorization of a.
c
c        b       real(n+1)
c                the right hand side vector.
c
c
c     on return
c
c        b       the solution vector  x .
c
c
c     modification of linpack. this version dated 1/13/86.
c     peter brown, university of houston, lawrence livermore natl. lab.
c
c     blas daxpy
c-----------------------------------------------------------------------
      integer iq, k, kb, kp1
      double precision c, s, t, t1, t2
c
c        minimize (b-a*x,b-a*x)
c        first form q*b.
c
         do 20 k = 1, n
            kp1 = k + 1
            iq = 2*(k-1) + 1
            c = q(iq)
            s = q(iq+1)
            t1 = b(k)
            t2 = b(kp1)
            b(k) = c*t1 - s*t2
            b(kp1) = s*t1 + c*t2
   20    continue
c
c        now solve  r*x = q*b.
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy (k-1, t, a(1,k), 1, b(1), 1)
   40    continue
      return
c----------------------- end of subroutine dhels -----------------------
      end
*DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C***BEGIN PROLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGESL
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGESL
C***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
C            factors computed by DGECO or DGEFA.
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DCOPY
C***PURPOSE  Copy a vector.
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCOPY
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2 (N, DX, INCX)
C***BEGIN PROLOGUE  DNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNRM2
      INTEGER NEXT
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,
     +                 ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C***FIRST EXECUTABLE STATEMENT  DNRM2
      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J = I,NN,INCX
      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = SQRT(SUM)
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
*DECK DSCAL
      SUBROUTINE DSCAL (N, DA, DX, INCX)
C***BEGIN PROLOGUE  DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSCAL
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
*DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***BEGIN PROLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***CATEGORY  D1A2
C***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  IDAMAX
      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
