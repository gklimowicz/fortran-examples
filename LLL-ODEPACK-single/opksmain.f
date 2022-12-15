*DECK SLSODE
      SUBROUTINE SLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C***BEGIN PROLOGUE  SLSODE
C***PURPOSE  Livermore Solver for Ordinary Differential Equations.
C            SLSODE solves the initial-value problem for stiff or
C            nonstiff systems of first-order ODE's,
C               dy/dt = f(t,y),   or, in component form,
C               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
C***CATEGORY  I1A
C***TYPE      SINGLE PRECISION (SLSODE-S, DLSODE-D)
C***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM,
C             STIFF, NONSTIFF
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C             Center for Applied Scientific Computing, L-561
C             Lawrence Livermore National Laboratory
C             Livermore, CA 94551.
C***DESCRIPTION
C
C     NOTE: The "Usage" and "Arguments" sections treat only a subset of
C           available options, in condensed fashion.  The options
C           covered and the information supplied will support most
C           standard uses of SLSODE.
C
C           For more sophisticated uses, full details on all options are
C           given in the concluding section, headed "Long Description."
C           A synopsis of the SLSODE Long Description is provided at the
C           beginning of that section; general topics covered are:
C           - Elements of the call sequence; optional input and output
C           - Optional supplemental routines in the SLSODE package
C           - internal COMMON block
C
C *Usage:
C     Communication between the user and the SLSODE package, for normal
C     situations, is summarized here.  This summary describes a subset
C     of the available options.  See "Long Description" for complete
C     details, including optional communication, nonstandard options,
C     and instructions for special situations.
C
C     A sample program is given in the "Examples" section.
C
C     Refer to the argument descriptions for the definitions of the
C     quantities that appear in the following sample declarations.
C
C     For MF = 10,
C        PARAMETER  (LRW = 20 + 16*NEQ,           LIW = 20)
C     For MF = 21 or 22,
C        PARAMETER  (LRW = 22 +  9*NEQ + NEQ**2,  LIW = 20 + NEQ)
C     For MF = 24 or 25,
C        PARAMETER  (LRW = 22 + 10*NEQ + (2*ML+MU)*NEQ,
C       *                                         LIW = 20 + NEQ)
C
C        EXTERNAL F, JAC
C        INTEGER  NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK(LIW),
C       *         LIW, MF
C        REAL Y(NEQ), T, TOUT, RTOL, ATOL(ntol), RWORK(LRW)
C
C        CALL SLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
C       *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
C
C *Arguments:
C     F     :EXT    Name of subroutine for right-hand-side vector f.
C                   This name must be declared EXTERNAL in calling
C                   program.  The form of F must be:
C
C                   SUBROUTINE  F (NEQ, T, Y, YDOT)
C                   INTEGER  NEQ
C                   REAL T, Y(*), YDOT(*)
C
C                   The inputs are NEQ, T, Y.  F is to set
C
C                   YDOT(i) = f(i,T,Y(1),Y(2),...,Y(NEQ)),
C                                                     i = 1, ..., NEQ .
C
C     NEQ   :IN     Number of first-order ODE's.
C
C     Y     :INOUT  Array of values of the y(t) vector, of length NEQ.
C                   Input:  For the first call, Y should contain the
C                           values of y(t) at t = T. (Y is an input
C                           variable only if ISTATE = 1.)
C                   Output: On return, Y will contain the values at the
C                           new t-value.
C
C     T     :INOUT  Value of the independent variable.  On return it
C                   will be the current value of t (normally TOUT).
C
C     TOUT  :IN     Next point where output is desired (.NE. T).
C
C     ITOL  :IN     1 or 2 according as ATOL (below) is a scalar or
C                   an array.
C
C     RTOL  :IN     Relative tolerance parameter (scalar).
C
C     ATOL  :IN     Absolute tolerance parameter (scalar or array).
C                   If ITOL = 1, ATOL need not be dimensioned.
C                   If ITOL = 2, ATOL must be dimensioned at least NEQ.
C
C                   The estimated local error in Y(i) will be controlled
C                   so as to be roughly less (in magnitude) than
C
C                   EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C                   EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C
C                   Thus the local error test passes if, in each
C                   component, either the absolute error is less than
C                   ATOL (or ATOL(i)), or the relative error is less
C                   than RTOL.
C
C                   Use RTOL = 0.0 for pure absolute error control, and
C                   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative
C                   error control.  Caution:  Actual (global) errors may
C                   exceed these local tolerances, so choose them
C                   conservatively.
C
C     ITASK :IN     Flag indicating the task SLSODE is to perform.
C                   Use ITASK = 1 for normal computation of output
C                   values of y at t = TOUT.
C
C     ISTATE:INOUT  Index used for input and output to specify the state
C                   of the calculation.
C                   Input:
C                    1   This is the first call for a problem.
C                    2   This is a subsequent call.
C                   Output:
C                    1   Nothing was done, as TOUT was equal to T.
C                    2   SLSODE was successful (otherwise, negative).
C                        Note that ISTATE need not be modified after a
C                        successful return.
C                   -1   Excess work done on this call (perhaps wrong
C                        MF).
C                   -2   Excess accuracy requested (tolerances too
C                        small).
C                   -3   Illegal input detected (see printed message).
C                   -4   Repeated error test failures (check all
C                        inputs).
C                   -5   Repeated convergence failures (perhaps bad
C                        Jacobian supplied or wrong choice of MF or
C                        tolerances).
C                   -6   Error weight became zero during problem
C                        (solution component i vanished, and ATOL or
C                        ATOL(i) = 0.).
C
C     IOPT  :IN     Flag indicating whether optional inputs are used:
C                   0   No.
C                   1   Yes.  (See "Optional inputs" under "Long
C                       Description," Part 1.)
C
C     RWORK :WORK   Real work array of length at least:
C                   20 + 16*NEQ                    for MF = 10,
C                   22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
C                   22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
C
C     LRW   :IN     Declared length of RWORK (in user's DIMENSION
C                   statement).
C
C     IWORK :WORK   Integer work array of length at least:
C                   20        for MF = 10,
C                   20 + NEQ  for MF = 21, 22, 24, or 25.
C
C                   If MF = 24 or 25, input in IWORK(1),IWORK(2) the
C                   lower and upper Jacobian half-bandwidths ML,MU.
C
C                   On return, IWORK contains information that may be
C                   of interest to the user:
C
C            Name   Location   Meaning
C            -----  ---------  -----------------------------------------
C            NST    IWORK(11)  Number of steps taken for the problem so
C                              far.
C            NFE    IWORK(12)  Number of f evaluations for the problem
C                              so far.
C            NJE    IWORK(13)  Number of Jacobian evaluations (and of
C                              matrix LU decompositions) for the problem
C                              so far.
C            NQU    IWORK(14)  Method order last used (successfully).
C            LENRW  IWORK(17)  Length of RWORK actually required.  This
C                              is defined on normal returns and on an
C                              illegal input return for insufficient
C                              storage.
C            LENIW  IWORK(18)  Length of IWORK actually required.  This
C                              is defined on normal returns and on an
C                              illegal input return for insufficient
C                              storage.
C
C     LIW   :IN     Declared length of IWORK (in user's DIMENSION
C                   statement).
C
C     JAC   :EXT    Name of subroutine for Jacobian matrix (MF =
C                   21 or 24).  If used, this name must be declared
C                   EXTERNAL in calling program.  If not used, pass a
C                   dummy name.  The form of JAC must be:
C
C                   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C                   INTEGER NEQ, ML, MU, NROWPD
C                   REAL T, Y(*), PD(NROWPD,*)
C
C                   See item c, under "Description" below for more
C                   information about JAC.
C
C     MF    :IN     Method flag.  Standard values are:
C                   10  Nonstiff (Adams) method, no Jacobian used.
C                   21  Stiff (BDF) method, user-supplied full Jacobian.
C                   22  Stiff method, internally generated full
C                       Jacobian.
C                   24  Stiff method, user-supplied banded Jacobian.
C                   25  Stiff method, internally generated banded
C                       Jacobian.
C
C *Description:
C     SLSODE solves the initial value problem for stiff or nonstiff
C     systems of first-order ODE's,
C
C        dy/dt = f(t,y) ,
C
C     or, in component form,
C
C        dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ))
C                                                  (i = 1, ..., NEQ) .
C
C     SLSODE is a package based on the GEAR and GEARB packages, and on
C     the October 23, 1978, version of the tentative ODEPACK user
C     interface standard, with minor modifications.
C
C     The steps in solving such a problem are as follows.
C
C     a. First write a subroutine of the form
C
C           SUBROUTINE  F (NEQ, T, Y, YDOT)
C           INTEGER  NEQ
C           REAL T, Y(*), YDOT(*)
C
C        which supplies the vector function f by loading YDOT(i) with
C        f(i).
C
C     b. Next determine (or guess) whether or not the problem is stiff.
C        Stiffness occurs when the Jacobian matrix df/dy has an
C        eigenvalue whose real part is negative and large in magnitude
C        compared to the reciprocal of the t span of interest.  If the
C        problem is nonstiff, use method flag MF = 10.  If it is stiff,
C        there are four standard choices for MF, and SLSODE requires the
C        Jacobian matrix in some form.  This matrix is regarded either
C        as full (MF = 21 or 22), or banded (MF = 24 or 25).  In the
C        banded case, SLSODE requires two half-bandwidth parameters ML
C        and MU. These are, respectively, the widths of the lower and
C        upper parts of the band, excluding the main diagonal.  Thus the
C        band consists of the locations (i,j) with
C
C           i - ML <= j <= i + MU ,
C
C        and the full bandwidth is ML + MU + 1 .
C
C     c. If the problem is stiff, you are encouraged to supply the
C        Jacobian directly (MF = 21 or 24), but if this is not feasible,
C        SLSODE will compute it internally by difference quotients (MF =
C        22 or 25).  If you are supplying the Jacobian, write a
C        subroutine of the form
C
C           SUBROUTINE  JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C           INTEGER  NEQ, ML, MU, NRWOPD
C           REAL T, Y(*), PD(NROWPD,*)
C
C        which provides df/dy by loading PD as follows:
C        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
C          the partial derivative of f(i) with respect to y(j).  (Ignore
C          the ML and MU arguments in this case.)
C        - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
C          df(i)/dy(j); i.e., load the diagonal lines of df/dy into the
C          rows of PD from the top down.
C        - In either case, only nonzero elements need be loaded.
C
C     d. Write a main program that calls subroutine SLSODE once for each
C        point at which answers are desired.  This should also provide
C        for possible use of logical unit 6 for output of error messages
C        by SLSODE.
C
C        Before the first call to SLSODE, set ISTATE = 1, set Y and T to
C        the initial values, and set TOUT to the first output point.  To
C        continue the integration after a successful return, simply
C        reset TOUT and call SLSODE again.  No other parameters need be
C        reset.
C
C *Examples:
C     The following is a simple example problem, with the coding needed
C     for its solution by SLSODE. The problem is from chemical kinetics,
C     and consists of the following three rate equations:
C
C        dy1/dt = -.04*y1 + 1.E4*y2*y3
C        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
C        dy3/dt = 3.E7*y2**2
C
C     on the interval from t = 0.0 to t = 4.E10, with initial conditions
C     y1 = 1.0, y2 = y3 = 0. The problem is stiff.
C
C     The following coding solves this problem with SLSODE, using 
C     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses 
C     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2 
C     has much smaller values.  At the end of the run, statistical 
C     quantities of interest are printed.
C
C        EXTERNAL  FEX, JEX
C        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
C       *         MF, NEQ
C        REAL  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
C        NEQ = 3
C        Y(1) = 1.
C        Y(2) = 0.
C        Y(3) = 0.
C        T = 0.
C        TOUT = .4
C        ITOL = 2
C        RTOL = 1.E-4
C        ATOL(1) = 1.E-6
C        ATOL(2) = 1.E-10
C        ATOL(3) = 1.E-6
C        ITASK = 1
C        ISTATE = 1
C        IOPT = 0
C        LRW = 58
C        LIW = 23
C        MF = 21
C        DO 40 IOUT = 1,12
C          CALL SLSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
C       *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
C          WRITE(6,20)  T, Y(1), Y(2), Y(3)
C    20    FORMAT(' At t =',E12.4,'   y =',3E14.6)
C          IF (ISTATE .LT. 0)  GO TO 80
C    40    TOUT = TOUT*10.
C        WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
C    60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
C        STOP
C    80  WRITE(6,90)  ISTATE
C    90  FORMAT(///' Error halt.. ISTATE =',I3)
C        STOP
C        END
C
C        SUBROUTINE  FEX (NEQ, T, Y, YDOT)
C        INTEGER  NEQ
C        REAL  T, Y(3), YDOT(3)
C        YDOT(1) = -.04*Y(1) + 1.E4*Y(2)*Y(3)
C        YDOT(3) = 3.E7*Y(2)*Y(2)
C        YDOT(2) = -YDOT(1) - YDOT(3)
C        RETURN
C        END
C
C        SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
C        INTEGER  NEQ, ML, MU, NRPD
C        REAL  T, Y(3), PD(NRPD,3)
C        PD(1,1) = -.04
C        PD(1,2) = 1.E4*Y(3)
C        PD(1,3) = 1.E4*Y(2)
C        PD(2,1) = .04
C        PD(2,3) = -PD(1,3)
C        PD(3,2) = 6.E7*Y(2)
C        PD(2,2) = -PD(1,2) - PD(3,2)
C        RETURN
C        END
C
C     The output from this program (on a Cray-1 in single precision)
C     is as follows.
C
C     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
C     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
C     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
C     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
C     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
C     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
C     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
C     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
C     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
C     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01
C     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
C     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00
C
C     No. steps = 330,  No. f-s = 405,  No. J-s = 69
C
C *Accuracy:
C     The accuracy of the solution depends on the choice of tolerances
C     RTOL and ATOL.  Actual (global) errors may exceed these local
C     tolerances, so choose them conservatively.
C
C *Cautions:
C     The work arrays should not be altered between calls to SLSODE for
C     the same problem, except possibly for the conditional and optional
C     inputs.
C
C *Portability:
C     Since NEQ is dimensioned inside SLSODE, some compilers may object
C     to a call to SLSODE with NEQ a scalar variable.  In this event, 
C     use DIMENSION NEQ(1).  Similar remarks apply to RTOL and ATOL.
C
C     Note to Cray users:
C     For maximum efficiency, use the CFT77 compiler.  Appropriate
C     compiler optimization directives have been inserted for CFT77.
C
C *Reference:
C     Alan C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE
C     Solvers," in Scientific Computing, R. S. Stepleman, et al., Eds.
C     (North-Holland, Amsterdam, 1983), pp. 55-64.
C
C *Long Description:
C     The following complete description of the user interface to
C     SLSODE consists of four parts:
C
C     1.  The call sequence to subroutine SLSODE, which is a driver
C         routine for the solver.  This includes descriptions of both
C         the call sequence arguments and user-supplied routines.
C         Following these descriptions is a description of optional
C         inputs available through the call sequence, and then a
C         description of optional outputs in the work arrays.
C
C     2.  Descriptions of other routines in the SLSODE package that may
C         be (optionally) called by the user.  These provide the ability
C         to alter error message handling, save and restore the internal
C         COMMON, and obtain specified derivatives of the solution y(t).
C
C     3.  Descriptions of COMMON block to be declared in overlay or
C         similar environments, or to be saved when doing an interrupt
C         of the problem and continued solution later.
C
C     4.  Description of two routines in the SLSODE package, either of
C         which the user may replace with his own version, if desired.
C         These relate to the measurement of errors.
C
C
C                         Part 1.  Call Sequence
C                         ----------------------
C
C     Arguments
C     ---------
C     The call sequence parameters used for input only are
C
C        F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
C
C     and those used for both input and output are
C
C        Y, T, ISTATE.
C
C     The work arrays RWORK and IWORK are also used for conditional and
C     optional inputs and optional outputs.  (The term output here
C     refers to the return from subroutine SLSODE to the user's calling
C     program.)
C
C     The legality of input parameters will be thoroughly checked on the
C     initial call for the problem, but not checked thereafter unless a
C     change in input parameters is flagged by ISTATE = 3 on input.
C
C     The descriptions of the call arguments are as follows.
C
C     F        The name of the user-supplied subroutine defining the ODE
C              system.  The system must be put in the first-order form
C              dy/dt = f(t,y), where f is a vector-valued function of
C              the scalar t and the vector y. Subroutine F is to compute
C              the function f. It is to have the form
C
C                 SUBROUTINE F (NEQ, T, Y, YDOT)
C                 REAL T, Y(*), YDOT(*)
C
C              where NEQ, T, and Y are input, and the array YDOT =
C              f(T,Y) is output.  Y and YDOT are arrays of length NEQ.
C              Subroutine F should not alter Y(1),...,Y(NEQ).  F must be
C              declared EXTERNAL in the calling program.
C
C              Subroutine F may access user-defined quantities in
C              NEQ(2),... and/or in Y(NEQ(1)+1),..., if NEQ is an array
C              (dimensioned in F) and/or Y has length exceeding NEQ(1).
C              See the descriptions of NEQ and Y below.
C
C              If quantities computed in the F routine are needed
C              externally to SLSODE, an extra call to F should be made
C              for this purpose, for consistent and accurate results.
C              If only the derivative dy/dt is needed, use SINTDY
C              instead.
C
C     NEQ      The size of the ODE system (number of first-order
C              ordinary differential equations).  Used only for input.
C              NEQ may be decreased, but not increased, during the
C              problem.  If NEQ is decreased (with ISTATE = 3 on input),
C              the remaining components of Y should be left undisturbed,
C              if these are to be accessed in F and/or JAC.
C
C              Normally, NEQ is a scalar, and it is generally referred
C              to as a scalar in this user interface description.
C              However, NEQ may be an array, with NEQ(1) set to the
C              system size.  (The SLSODE package accesses only NEQ(1).)
C              In either case, this parameter is passed as the NEQ
C              argument in all calls to F and JAC.  Hence, if it is an
C              array, locations NEQ(2),... may be used to store other
C              integer data and pass it to F and/or JAC.  Subroutines
C              F and/or JAC must include NEQ in a DIMENSION statement
C              in that case.
C
C     Y        A real array for the vector of dependent variables, of
C              length NEQ or more.  Used for both input and output on
C              the first call (ISTATE = 1), and only for output on
C              other calls.  On the first call, Y must contain the
C              vector of initial values.  On output, Y contains the
C              computed solution vector, evaluated at T. If desired,
C              the Y array may be used for other purposes between
C              calls to the solver.
C
C              This array is passed as the Y argument in all calls to F
C              and JAC.  Hence its length may exceed NEQ, and locations
C              Y(NEQ+1),... may be used to store other real data and
C              pass it to F and/or JAC.  (The SLSODE package accesses
C              only Y(1),...,Y(NEQ).)
C
C     T        The independent variable.  On input, T is used only on
C              the first call, as the initial point of the integration.
C              On output, after each call, T is the value at which a
C              computed solution Y is evaluated (usually the same as
C              TOUT).  On an error return, T is the farthest point
C              reached.
C
C     TOUT     The next value of T at which a computed solution is
C              desired.  Used only for input.
C
C              When starting the problem (ISTATE = 1), TOUT may be equal
C              to T for one call, then should not equal T for the next
C              call.  For the initial T, an input value of TOUT .NE. T
C              is used in order to determine the direction of the
C              integration (i.e., the algebraic sign of the step sizes)
C              and the rough scale of the problem.  Integration in
C              either direction (forward or backward in T) is permitted.
C
C              If ITASK = 2 or 5 (one-step modes), TOUT is ignored
C              after the first call (i.e., the first call with
C              TOUT .NE. T).  Otherwise, TOUT is required on every call.
C
C              If ITASK = 1, 3, or 4, the values of TOUT need not be
C              monotone, but a value of TOUT which backs up is limited
C              to the current internal T interval, whose endpoints are
C              TCUR - HU and TCUR.  (See "Optional Outputs" below for
C              TCUR and HU.)
C
C
C     ITOL     An indicator for the type of error control.  See
C              description below under ATOL.  Used only for input.
C
C     RTOL     A relative error tolerance parameter, either a scalar or
C              an array of length NEQ.  See description below under
C              ATOL.  Input only.
C
C     ATOL     An absolute error tolerance parameter, either a scalar or
C              an array of length NEQ.  Input only.
C
C              The input parameters ITOL, RTOL, and ATOL determine the
C              error control performed by the solver.  The solver will
C              control the vector e = (e(i)) of estimated local errors
C              in Y, according to an inequality of the form
C
C                 rms-norm of ( e(i)/EWT(i) ) <= 1,
C
C              where
C
C                 EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C
C              and the rms-norm (root-mean-square norm) here is
C
C                 rms-norm(v) = SQRT(sum v(i)**2 / NEQ).
C
C              Here EWT = (EWT(i)) is a vector of weights which must
C              always be positive, and the values of RTOL and ATOL
C              should all be nonnegative.  The following table gives the
C              types (scalar/array) of RTOL and ATOL, and the
C              corresponding form of EWT(i).
C
C              ITOL    RTOL      ATOL      EWT(i)
C              ----    ------    ------    -----------------------------
C              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
C              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
C              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
C              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C              When either of these parameters is a scalar, it need not
C              be dimensioned in the user's calling program.
C
C              If none of the above choices (with ITOL, RTOL, and ATOL
C              fixed throughout the problem) is suitable, more general
C              error controls can be obtained by substituting
C              user-supplied routines for the setting of EWT and/or for
C              the norm calculation.  See Part 4 below.
C
C              If global errors are to be estimated by making a repeated
C              run on the same problem with smaller tolerances, then all
C              components of RTOL and ATOL (i.e., of EWT) should be
C              scaled down uniformly.
C
C     ITASK    An index specifying the task to be performed.  Input
C              only.  ITASK has the following values and meanings:
C              1   Normal computation of output values of y(t) at
C                  t = TOUT (by overshooting and interpolating).
C              2   Take one step only and return.
C              3   Stop at the first internal mesh point at or beyond
C                  t = TOUT and return.
C              4   Normal computation of output values of y(t) at
C                  t = TOUT but without overshooting t = TCRIT.  TCRIT
C                  must be input as RWORK(1).  TCRIT may be equal to or
C                  beyond TOUT, but not behind it in the direction of
C                  integration.  This option is useful if the problem
C                  has a singularity at or beyond t = TCRIT.
C              5   Take one step, without passing TCRIT, and return.
C                  TCRIT must be input as RWORK(1).
C
C              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C              (within roundoff), it will return T = TCRIT (exactly) to
C              indicate this (unless ITASK = 4 and TOUT comes before
C              TCRIT, in which case answers at T = TOUT are returned
C              first).
C
C     ISTATE   An index used for input and output to specify the state
C              of the calculation.
C
C              On input, the values of ISTATE are as follows:
C              1   This is the first call for the problem
C                  (initializations will be done).  See "Note" below.
C              2   This is not the first call, and the calculation is to
C                  continue normally, with no change in any input
C                  parameters except possibly TOUT and ITASK.  (If ITOL,
C                  RTOL, and/or ATOL are changed between calls with
C                  ISTATE = 2, the new values will be used but not
C                  tested for legality.)
C              3   This is not the first call, and the calculation is to
C                  continue normally, but with a change in input
C                  parameters other than TOUT and ITASK.  Changes are
C                  allowed in NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C                  ML, MU, and any of the optional inputs except H0.
C                  (See IWORK description for ML and MU.)
C
C              Note:  A preliminary call with TOUT = T is not counted as
C              a first call here, as no initialization or checking of
C              input is done.  (Such a call is sometimes useful for the
C              purpose of outputting the initial conditions.)  Thus the
C              first call for which TOUT .NE. T requires ISTATE = 1 on
C              input.
C
C              On output, ISTATE has the following values and meanings:
C               1  Nothing was done, as TOUT was equal to T with
C                  ISTATE = 1 on input.
C               2  The integration was performed successfully.
C              -1  An excessive amount of work (more than MXSTEP steps)
C                  was done on this call, before completing the
C                  requested task, but the integration was otherwise
C                  successful as far as T. (MXSTEP is an optional input
C                  and is normally 500.)  To continue, the user may
C                  simply reset ISTATE to a value >1 and call again (the
C                  excess work step counter will be reset to 0).  In
C                  addition, the user may increase MXSTEP to avoid this
C                  error return; see "Optional Inputs" below.
C              -2  Too much accuracy was requested for the precision of
C                  the machine being used.  This was detected before
C                  completing the requested task, but the integration
C                  was successful as far as T. To continue, the
C                  tolerance parameters must be reset, and ISTATE must
C                  be set to 3. The optional output TOLSF may be used
C                  for this purpose.  (Note:  If this condition is
C                  detected before taking any steps, then an illegal
C                  input return (ISTATE = -3) occurs instead.)
C              -3  Illegal input was detected, before taking any
C                  integration steps.  See written message for details.
C                  (Note:  If the solver detects an infinite loop of
C                  calls to the solver with illegal input, it will cause
C                  the run to stop.)
C              -4  There were repeated error-test failures on one
C                  attempted step, before completing the requested task,
C                  but the integration was successful as far as T.  The
C                  problem may have a singularity, or the input may be
C                  inappropriate.
C              -5  There were repeated convergence-test failures on one
C                  attempted step, before completing the requested task,
C                  but the integration was successful as far as T. This
C                  may be caused by an inaccurate Jacobian matrix, if
C                  one is being used.
C              -6  EWT(i) became zero for some i during the integration.
C                  Pure relative error control (ATOL(i)=0.0) was
C                  requested on a variable which has now vanished.  The
C                  integration was successful as far as T.
C
C              Note:  Since the normal output value of ISTATE is 2, it
C              does not need to be reset for normal continuation.  Also,
C              since a negative input value of ISTATE will be regarded
C              as illegal, a negative output value requires the user to
C              change it, and possibly other inputs, before calling the
C              solver again.
C
C     IOPT     An integer flag to specify whether any optional inputs
C              are being used on this call.  Input only.  The optional
C              inputs are listed under a separate heading below.
C              0   No optional inputs are being used.  Default values
C                  will be used in all cases.
C              1   One or more optional inputs are being used.
C
C     RWORK    A real working array (single precision).  The length of
C              RWORK must be at least
C
C                 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM
C
C              where
C                 NYH = the initial value of NEQ,
C              MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                       smaller value is given as an optional input),
C                 LWM = 0           if MITER = 0,
C                 LWM = NEQ**2 + 2  if MITER = 1 or 2,
C                 LWM = NEQ + 2     if MITER = 3, and
C                 LWM = (2*ML + MU + 1)*NEQ + 2
C                                   if MITER = 4 or 5.
C              (See the MF description below for METH and MITER.)
C
C              Thus if MAXORD has its default value and NEQ is constant,
C              this length is:
C              20 + 16*NEQ                    for MF = 10,
C              22 + 16*NEQ + NEQ**2           for MF = 11 or 12,
C              22 + 17*NEQ                    for MF = 13,
C              22 + 17*NEQ + (2*ML + MU)*NEQ  for MF = 14 or 15,
C              20 +  9*NEQ                    for MF = 20,
C              22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
C              22 + 10*NEQ                    for MF = 23,
C              22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
C
C              The first 20 words of RWORK are reserved for conditional
C              and optional inputs and optional outputs.
C
C              The following word in RWORK is a conditional input:
C              RWORK(1) = TCRIT, the critical value of t which the
C                         solver is not to overshoot.  Required if ITASK
C                         is 4 or 5, and ignored otherwise.  See ITASK.
C
C     LRW      The length of the array RWORK, as declared by the user.
C              (This will be checked by the solver.)
C
C     IWORK    An integer work array.  Its length must be at least
C              20       if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
C              20 + NEQ otherwise (MF = 11, 12, 14, 15, 21, 22, 24, 25).
C              (See the MF description below for MITER.)  The first few
C              words of IWORK are used for conditional and optional
C              inputs and optional outputs.
C
C              The following two words in IWORK are conditional inputs:
C              IWORK(1) = ML   These are the lower and upper half-
C              IWORK(2) = MU   bandwidths, respectively, of the banded
C                              Jacobian, excluding the main diagonal.
C                         The band is defined by the matrix locations
C                         (i,j) with i - ML <= j <= i + MU. ML and MU
C                         must satisfy 0 <= ML,MU <= NEQ - 1. These are
C                         required if MITER is 4 or 5, and ignored
C                         otherwise.  ML and MU may in fact be the band
C                         parameters for a matrix to which df/dy is only
C                         approximately equal.
C
C     LIW      The length of the array IWORK, as declared by the user.
C              (This will be checked by the solver.)
C
C     Note:  The work arrays must not be altered between calls to SLSODE
C     for the same problem, except possibly for the conditional and
C     optional inputs, and except for the last 3*NEQ words of RWORK.
C     The latter space is used for internal scratch space, and so is
C     available for use by the user outside SLSODE between calls, if
C     desired (but not for use by F or JAC).
C
C     JAC      The name of the user-supplied routine (MITER = 1 or 4) to
C              compute the Jacobian matrix, df/dy, as a function of the
C              scalar t and the vector y.  (See the MF description below
C              for MITER.)  It is to have the form
C
C                 SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C                 REAL T, Y(*), PD(NROWPD,*)
C
C              where NEQ, T, Y, ML, MU, and NROWPD are input and the
C              array PD is to be loaded with partial derivatives
C              (elements of the Jacobian matrix) on output.  PD must be
C              given a first dimension of NROWPD.  T and Y have the same
C              meaning as in subroutine F.
C
C              In the full matrix case (MITER = 1), ML and MU are
C              ignored, and the Jacobian is to be loaded into PD in
C              columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
C
C              In the band matrix case (MITER = 4), the elements within
C              the band are to be loaded into PD in columnwise manner,
C              with diagonal lines of df/dy loaded into the rows of PD.
C              Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).  ML
C              and MU are the half-bandwidth parameters (see IWORK).
C              The locations in PD in the two triangular areas which
C              correspond to nonexistent matrix elements can be ignored
C              or loaded arbitrarily, as they are overwritten by SLSODE.
C
C              JAC need not provide df/dy exactly. A crude approximation
C              (possibly with a smaller bandwidth) will do.
C
C              In either case, PD is preset to zero by the solver, so
C              that only the nonzero elements need be loaded by JAC.
C              Each call to JAC is preceded by a call to F with the same
C              arguments NEQ, T, and Y. Thus to gain some efficiency,
C              intermediate quantities shared by both calculations may
C              be saved in a user COMMON block by F and not recomputed
C              by JAC, if desired.  Also, JAC may alter the Y array, if
C              desired.  JAC must be declared EXTERNAL in the calling
C              program.
C
C              Subroutine JAC may access user-defined quantities in
C              NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C              (dimensioned in JAC) and/or Y has length exceeding
C              NEQ(1).  See the descriptions of NEQ and Y above.
C
C     MF       The method flag.  Used only for input.  The legal values
C              of MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24,
C              and 25.  MF has decimal digits METH and MITER:
C                 MF = 10*METH + MITER .
C
C              METH indicates the basic linear multistep method:
C              1   Implicit Adams method.
C              2   Method based on backward differentiation formulas
C                  (BDF's).
C
C              MITER indicates the corrector iteration method:
C              0   Functional iteration (no Jacobian matrix is
C                  involved).
C              1   Chord iteration with a user-supplied full (NEQ by
C                  NEQ) Jacobian.
C              2   Chord iteration with an internally generated
C                  (difference quotient) full Jacobian (using NEQ
C                  extra calls to F per df/dy value).
C              3   Chord iteration with an internally generated
C                  diagonal Jacobian approximation (using one extra call
C                  to F per df/dy evaluation).
C              4   Chord iteration with a user-supplied banded Jacobian.
C              5   Chord iteration with an internally generated banded
C                  Jacobian (using ML + MU + 1 extra calls to F per
C                  df/dy evaluation).
C
C              If MITER = 1 or 4, the user must supply a subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For other values of MITER, a dummy argument can be used.
C
C     Optional Inputs
C     ---------------
C     The following is a list of the optional inputs provided for in the
C     call sequence.  (See also Part 2.)  For each such input variable,
C     this table lists its name as used in this documentation, its
C     location in the call sequence, its meaning, and the default value.
C     The use of any of these inputs requires IOPT = 1, and in that case
C     all of these inputs are examined.  A value of zero for any of
C     these optional inputs will cause the default value to be used.
C     Thus to use a subset of the optional inputs, simply preload
C     locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively,
C     and then set those of interest to nonzero values.
C
C     Name    Location   Meaning and default value
C     ------  ---------  -----------------------------------------------
C     H0      RWORK(5)   Step size to be attempted on the first step.
C                        The default value is determined by the solver.
C     HMAX    RWORK(6)   Maximum absolute step size allowed.  The
C                        default value is infinite.
C     HMIN    RWORK(7)   Minimum absolute step size allowed.  The
C                        default value is 0.  (This lower bound is not
C                        enforced on the final step before reaching
C                        TCRIT when ITASK = 4 or 5.)
C     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
C                        is 12 if METH = 1, and 5 if METH = 2. (See the
C                        MF description above for METH.)  If MAXORD
C                        exceeds the default value, it will be reduced
C                        to the default value.  If MAXORD is changed
C                        during the problem, it may cause the current
C                        order to be reduced.
C     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
C                        allowed during one call to the solver.  The
C                        default value is 500.
C     MXHNIL  IWORK(7)   Maximum number of messages printed (per
C                        problem) warning that T + H = T on a step
C                        (H = step size).  This must be positive to
C                        result in a nondefault value.  The default
C                        value is 10.
C
C     Optional Outputs
C     ----------------
C     As optional additional output from SLSODE, the variables listed
C     below are quantities related to the performance of SLSODE which 
C     are available to the user.  These are communicated by way of the
C     work arrays, but also have internal mnemonic names as shown. 
C     Except where stated otherwise, all of these outputs are defined on
C     any successful return from SLSODE, and on any return with ISTATE =
C     -1, -2, -4, -5, or -6.  On an illegal input return (ISTATE = -3),
C     they will be unchanged from their existing values (if any), except
C     possibly for TOLSF, LENRW, and LENIW.  On any error return,
C     outputs relevant to the error will be defined, as noted below.
C
C     Name   Location   Meaning
C     -----  ---------  ------------------------------------------------
C     HU     RWORK(11)  Step size in t last used (successfully).
C     HCUR   RWORK(12)  Step size to be attempted on the next step.
C     TCUR   RWORK(13)  Current value of the independent variable which
C                       the solver has actually reached, i.e., the
C                       current internal mesh point in t. On output,
C                       TCUR will always be at least as far as the
C                       argument T, but may be farther (if interpolation
C                       was done).
C     TOLSF  RWORK(14)  Tolerance scale factor, greater than 1.0,
C                       computed when a request for too much accuracy
C                       was detected (ISTATE = -3 if detected at the
C                       start of the problem, ISTATE = -2 otherwise).
C                       If ITOL is left unaltered but RTOL and ATOL are
C                       uniformly scaled up by a factor of TOLSF for the
C                       next call, then the solver is deemed likely to
C                       succeed.  (The user may also ignore TOLSF and
C                       alter the tolerance parameters in any other way
C                       appropriate.)
C     NST    IWORK(11)  Number of steps taken for the problem so far.
C     NFE    IWORK(12)  Number of F evaluations for the problem so far.
C     NJE    IWORK(13)  Number of Jacobian evaluations (and of matrix LU
C                       decompositions) for the problem so far.
C     NQU    IWORK(14)  Method order last used (successfully).
C     NQCUR  IWORK(15)  Order to be attempted on the next step.
C     IMXER  IWORK(16)  Index of the component of largest magnitude in
C                       the weighted local error vector ( e(i)/EWT(i) ),
C                       on an error return with ISTATE = -4 or -5.
C     LENRW  IWORK(17)  Length of RWORK actually required.  This is
C                       defined on normal returns and on an illegal
C                       input return for insufficient storage.
C     LENIW  IWORK(18)  Length of IWORK actually required.  This is
C                       defined on normal returns and on an illegal
C                       input return for insufficient storage.
C
C     The following two arrays are segments of the RWORK array which may
C     also be of interest to the user as optional outputs.  For each
C     array, the table below gives its internal name, its base address
C     in RWORK, and its description.
C
C     Name  Base address  Description
C     ----  ------------  ----------------------------------------------
C     YH    21            The Nordsieck history array, of size NYH by
C                         (NQCUR + 1), where NYH is the initial value of
C                         NEQ.  For j = 0,1,...,NQCUR, column j + 1 of
C                         YH contains HCUR**j/factorial(j) times the jth
C                         derivative of the interpolating polynomial
C                         currently representing the solution, evaluated
C                         at t = TCUR.
C     ACOR  LENRW-NEQ+1   Array of size NEQ used for the accumulated
C                         corrections on each step, scaled on output to
C                         represent the estimated local error in Y on
C                         the last step.  This is the vector e in the
C                         description of the error control.  It is
C                         defined only on successful return from SLSODE.
C
C
C                    Part 2.  Other Callable Routines
C                    --------------------------------
C
C     The following are optional calls which the user may make to gain
C     additional capabilities in conjunction with SLSODE.
C
C     Form of call              Function
C     ------------------------  ----------------------------------------
C     CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                               output of messages from SLSODE, if the
C                               default is not desired.  The default
C                               value of LUN is 6. This call may be made
C                               at any time and will take effect
C                               immediately.
C     CALL XSETF(MFLAG)         Set a flag to control the printing of
C                               messages by SLSODE.  MFLAG = 0 means do
C                               not print.  (Danger:  this risks losing
C                               valuable information.)  MFLAG = 1 means
C                               print (the default).  This call may be
C                               made at any time and will take effect
C                               immediately.
C     CALL SSRCOM(RSAV,ISAV,JOB)  Saves and restores the contents of the
C                               internal COMMON blocks used by SLSODE
C                               (see Part 3 below).  RSAV must be a
C                               real array of length 218 or more, and
C                               ISAV must be an integer array of length
C                               37 or more.  JOB = 1 means save COMMON
C                               into RSAV/ISAV.  JOB = 2 means restore
C                               COMMON from same.  SSRCOM is useful if
C                               one is interrupting a run and restarting
C                               later, or alternating between two or
C                               more problems solved with SLSODE.
C     CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C     (see below)               orders, at a specified point t, if
C                               desired.  It may be called only after a
C                               successful return from SLSODE.  Detailed
C                               instructions follow.
C
C     Detailed instructions for using SINTDY
C     --------------------------------------
C     The form of the CALL is:
C
C           CALL SINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C     The input parameters are:
C
C     T          Value of independent variable where answers are
C                desired (normally the same as the T last returned by
C                SLSODE).  For valid results, T must lie between
C                TCUR - HU and TCUR.  (See "Optional Outputs" above
C                for TCUR and HU.)
C     K          Integer order of the derivative desired.  K must
C                satisfy 0 <= K <= NQCUR, where NQCUR is the current
C                order (see "Optional Outputs").  The capability
C                corresponding to K = 0, i.e., computing y(t), is
C                already provided by SLSODE directly.  Since
C                NQCUR >= 1, the first derivative dy/dt is always
C                available with SINTDY.
C     RWORK(21)  The base address of the history array YH.
C     NYH        Column length of YH, equal to the initial value of NEQ.
C
C     The output parameters are:
C
C     DKY        Real array of length NEQ containing the computed value
C                of the Kth derivative of y(t).
C     IFLAG      Integer flag, returned as 0 if K and T were legal,
C                -1 if K was illegal, and -2 if T was illegal.
C                On an error return, a message is also written.
C
C
C                          Part 3.  Common Blocks
C                          ----------------------
C
C     If SLSODE is to be used in an overlay situation, the user must
C     declare, in the primary overlay, the variables in:
C     (1) the call sequence to SLSODE,
C     (2) the internal COMMON block /SLS001/, of length 255 
C         (218 single precision words followed by 37 integer words).
C
C     If SLSODE is used on a system in which the contents of internal
C     COMMON blocks are not preserved between calls, the user should
C     declare the above COMMON block in his main program to insure that
C     its contents are preserved.
C
C     If the solution of a given problem by SLSODE is to be interrupted
C     and then later continued, as when restarting an interrupted run or
C     alternating between two or more problems, the user should save,
C     following the return from the last SLSODE call prior to the
C     interruption, the contents of the call sequence variables and the
C     internal COMMON block, and later restore these values before the
C     next SLSODE call for that problem.   In addition, if XSETUN and/or
C     XSETF was called for non-default handling of error messages, then
C     these calls must be repeated.  To save and restore the COMMON
C     block, use subroutine SSRCOM (see Part 2 above).
C
C
C              Part 4.  Optionally Replaceable Solver Routines
C              -----------------------------------------------
C
C     Below are descriptions of two routines in the SLSODE package which
C     relate to the measurement of errors.  Either routine can be
C     replaced by a user-supplied version, if desired.  However, since
C     such a replacement may have a major impact on performance, it
C     should be done only when absolutely necessary, and only with great
C     caution.  (Note:  The means by which the package version of a
C     routine is superseded by the user's version may be system-
C     dependent.)
C
C     SEWSET
C     ------
C     The following subroutine is called just before each internal
C     integration step, and sets the array of error weights, EWT, as
C     described under ITOL/RTOL/ATOL above:
C
C           SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C
C     where NEQ, ITOL, RTOL, and ATOL are as in the SLSODE call
C     sequence, YCUR contains the current dependent variable vector,
C     and EWT is the array of weights set by SEWSET.
C
C     If the user supplies this subroutine, it must return in EWT(i)
C     (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C     in Y(i) to.  The EWT array returned by SEWSET is passed to the
C     SVNORM routine (see below), and also used by SLSODE in the
C     computation of the optional output IMXER, the diagonal Jacobian
C     approximation, and the increments for difference quotient
C     Jacobians.
C
C     In the user-supplied version of SEWSET, it may be desirable to use
C     the current values of derivatives of y. Derivatives up to order NQ
C     are available from the history array YH, described above under
C     optional outputs.  In SEWSET, YH is identical to the YCUR array,
C     extended to NQ + 1 columns with a column length of NYH and scale
C     factors of H**j/factorial(j).  On the first call for the problem,
C     given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C     NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C     can be obtained by including in SEWSET the statements:
C           REAL RLS
C           COMMON /SLS001/ RLS(218),ILS(37)
C           NQ = ILS(33)
C           NST = ILS(34)
C           H = RLS(212)
C     Thus, for example, the current value of dy/dt can be obtained as
C     YCUR(NYH+i)/H (i=1,...,NEQ) (and the division by H is unnecessary
C     when NST = 0).
C
C     SVNORM
C     ------
C     SVNORM is a real function routine which computes the weighted
C     root-mean-square norm of a vector v:
C
C        d = SVNORM (n, v, w)
C
C     where:
C     n = the length of the vector,
C     v = real array of length n containing the vector,
C     w = real array of length n containing weights,
C     d = SQRT( (1/n) * sum(v(i)*w(i))**2 ).
C
C     SVNORM is called with n = NEQ and with w(i) = 1.0/EWT(i), where
C     EWT is as set by subroutine SEWSET.
C
C     If the user supplies this function, it should return a nonnegative
C     value of SVNORM suitable for use in the error control in SLSODE.
C     None of the arguments should be altered by SVNORM.  For example, a
C     user-supplied SVNORM routine might:
C     - Substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
C     - Ignore some components of v in the norm, with the effect of
C       suppressing the error control on those components of Y.
C  ---------------------------------------------------------------------
C***ROUTINES CALLED  SEWSET, SINTDY, RUMACH, SSTODE, SVNORM, XERRWV
C***COMMON BLOCKS    SLS001
C***REVISION HISTORY  (YYYYMMDD)
C 19791129  DATE WRITTEN
C 19791213  Minor changes to declarations; DELP init. in STODE.
C 19800118  Treat NEQ as array; integer declarations added throughout;
C           minor changes to prologue.
C 19800306  Corrected TESCO(1,NQP1) setting in CFODE.
C 19800519  Corrected access of YH on forced order reduction;
C           numerous corrections to prologues and other comments.
C 19800617  In main driver, added loading of SQRT(UROUND) in RWORK;
C           minor corrections to main prologue.
C 19800923  Added zero initialization of HU and NQU.
C 19801218  Revised XERRWV routine; minor corrections to main prologue.
C 19810401  Minor changes to comments and an error message.
C 19810814  Numerous revisions: replaced EWT by 1/EWT; used flags
C           JCUR, ICF, IERPJ, IERSL between STODE and subordinates;
C           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF;
C           reorganized returns from STODE; reorganized type decls.;
C           fixed message length in XERRWV; changed default LUNIT to 6;
C           changed Common lengths; changed comments throughout.
C 19870330  Major update by ACH: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODE;
C           in STODE, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 19890426  Modified prologue to SLATEC/LDOC format.  (FNF)
C 19890501  Many improvements to prologue.  (FNF)
C 19890503  A few final corrections to prologue.  (FNF)
C 19890504  Minor cosmetic changes.  (FNF)
C 19890510  Corrected description of Y in Arguments section.  (FNF)
C 19890517  Minor corrections to prologue.  (FNF)
C 19920514  Updated with prologue edited 891025 by G. Shaw for manual.
C 19920515  Converted source lines to upper case.  (FNF)
C 19920603  Revised XERRWV calls using mixed upper-lower case.  (ACH)
C 19920616  Revised prologue comment regarding CFT.  (ACH)
C 19921116  Revised prologue comments regarding Common.  (ACH).
C 19930326  Added comment about non-reentrancy.  (FNF)
C 19930723  Changed R1MACH to RUMACH. (FNF)
C 19930801  Removed ILLIN and NTREP from Common (affects driver logic);
C           minor changes to prologue and internal comments;
C           changed Hollerith strings to quoted strings; 
C           changed internal comments to mixed case;
C           replaced XERRWV with new version using character type;
C           changed dummy dimensions from 1 to *. (ACH)
C 19930809  Changed to generic intrinsic names; changed names of
C           subprograms and Common blocks to SLSODE etc. (ACH)
C 19930929  Eliminated use of REAL intrinsic; other minor changes. (ACH)
C 20010412  Removed all 'own' variables from Common block /SLS001/
C           (affects declarations in 6 routines). (ACH)
C 20010509  Minor corrections to prologue. (ACH)
C 20031105  Restored 'own' variables to Common block /SLS001/, to
C           enable interrupt/restart feature. (ACH)
C 20031112  Added SAVE statements for data-loaded constants.
C
C***  END PROLOGUE  SLSODE
C
C*Internal Notes:
C
C Other Routines in the SLSODE Package.
C
C In addition to Subroutine SLSODE, the SLSODE package includes the
C following subroutines and function routines:
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODE   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPREPJ   computes and preprocesses the Jacobian matrix J = df/dy
C           and the Newton iteration matrix P = I - h*l0*J.
C  SSOLSY   manages solution of linear system in chord iteration.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted R.M.S. norm of a vector.
C  SSRCOM   is a user-callable routine to save and restore
C           the contents of the internal Common block.
C  SGEFA and SGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  SGBFA and SGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, IUMACH   handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note: SVNORM, RUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C**End
C
C  Declare externals.
      EXTERNAL SPREPJ, SSOLSY
      REAL RUMACH, SVNORM
C
C  Declare all other variables.
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*80 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following internal Common block contains
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODE, SINTDY, SSTODE,
C SPREPJ, and SSOLSY.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      DATA  MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .GT. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  SLSODE
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENWM = N*N + 2
      IF (MITER .EQ. 3) LENWM = N + 2
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 20
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to SSTODE. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO 80 I = 1,N
 80     RWORK(I+LSAVF-1) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
C if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
C                                       1
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         f(i)   = i-th component of initial value of f,
C         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODE.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODE-  Warning..internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      (H = step size). Solver will continue anyway'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODE-  Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      It will not be issued again for this problem'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C  CALL SSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,SPREPJ,SSOLSY)
C-----------------------------------------------------------------------
      CALL SSTODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, SPREPJ, SSOLSY)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODE.
C If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.  The optional outputs
C are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODE-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODE-  At T (=R1), EWT(I1) has become R2 .LE. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODE-  At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  see TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODE-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODE-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0E0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input 
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODE-  ISTATE (=I1) illegal '
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODE-  ITASK (=I1) illegal  '
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODE-  ISTATE .GT. 1 but SLSODE not initialized '
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODE-  NEQ (=I1) .LT. 1     '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODE-  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODE-  ITOL (=I1) illegal   '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODE-  IOPT (=I1) illegal   '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODE-  MF (=I1) illegal     '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 609  MSG = 'SLSODE-  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 610  MSG = 'SLSODE-  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODE-  MAXORD (=I1) .LT. 0  '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODE-  MXSTEP (=I1) .LT. 0  '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODE-  MXHNIL (=I1) .LT. 0  '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODE-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODE-  HMAX (=R1) .LT. 0.0  '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODE-  HMIN (=R1) .LT. 0.0  '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  CONTINUE
      MSG='SLSODE-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  CONTINUE
      MSG='SLSODE-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODE-  RTOL(I1) is R1 .LT. 0.0        '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODE-  ATOL(I1) is R1 .LT. 0.0        '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODE-  EWT(I1) is R1 .LE. 0.0         '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  CONTINUE
      MSG='SLSODE-  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='SLSODE-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='SLSODE-  ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='SLSODE-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODE-  At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODE-  Trouble in SINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODE-  Run aborted.. apparent infinite loop     '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- END OF SUBROUTINE SLSODE ----------------------
      END
*DECK SLSODES
      SUBROUTINE SLSODES (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 12 November 2003 version of
C SLSODES: Livermore Solver for Ordinary Differential Equations
C          with general Sparse Jacobian matrix.
C
C This version is in single precision.
C
C SLSODES solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C SLSODES is a variant of the SLSODE package, and is intended for
C problems in which the Jacobian matrix df/dy has an arbitrary
C sparse structure (when the problem is stiff).
C
C Authors:       Alan C. Hindmarsh
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Andrew H. Sherman
C                J. S. Nolen and Associates
C                Houston, TX 77084
C-----------------------------------------------------------------------
C References:
C 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C
C 2.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
C     Yale Sparse Matrix Package: I. The Symmetric Codes,
C     Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.
C
C 3.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
C     Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
C     Research Report No. 114, Dept. of Computer Sciences, Yale
C     University, 1977.
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODES package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First provide a subroutine of the form:
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Next determine (or guess) whether or not the problem is stiff.
C Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
C whose real part is negative and large in magnitude, compared to the
C reciprocal of the t span of interest.  If the problem is nonstiff,
C use a method flag MF = 10.  If it is stiff, there are two standard
C choices for the method flag, MF = 121 and MF = 222.  In both cases,
C SLSODES requires the Jacobian matrix in some form, and it treats this
C matrix in general sparse form, with sparsity structure determined
C internally.  (For options where the user supplies the sparsity
C structure, see the full description of MF below.)
C
C C. If the problem is stiff, you are encouraged to supply the Jacobian
C directly (MF = 121), but if this is not feasible, SLSODES will
C compute it internally by difference quotients (MF = 222).
C If you are supplying the Jacobian, provide a subroutine of the form:
C               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
C               DIMENSION Y(*), IAN(*), JAN(*), PDJ(*)
C Here NEQ, T, Y, and J are input arguments, and the JAC routine is to
C load the array PDJ (of length NEQ) with the J-th column of df/dy.
C I.e., load PDJ(i) with df(i)/dy(J) for all relevant values of i.
C The arguments IAN and JAN should be ignored for normal situations.
C SLSODES will call the JAC routine with J = 1,2,...,NEQ.
C Only nonzero elements need be loaded.  Usually, a crude approximation
C to df/dy, possibly with fewer nonzero elements, will suffice.
C
C D. Write a main program which calls Subroutine SLSODES once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages by
C SLSODES.  On the first call to SLSODES, supply arguments as follows:
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared External in calling program.
C NEQ    = number of first order ODEs.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable t.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          The estimated local error in Y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of Y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             20 + 16*NEQ            for MF = 10,
C             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
C                                    for MF = 121 or 222,
C          where:
C          NNZ    = the number of nonzero elements in the sparse
C                   Jacobian (if this is unknown, use an estimate), and
C          LENRAT = the real to integer wordlength ratio (usually 1 in
C                   single precision and 2 in double precision).
C          In any case, the required size of RWORK cannot generally
C          be predicted in advance if MF = 121 or 222, and the value
C          above is a rough estimate of a crude lower bound.  Some
C          experimentation with this size may be necessary.
C          (When known, the correct required length is an optional
C          output, available in IWORK(17).)
C LRW    = declared length of RWORK (in user dimension).
C IWORK  = integer work array of length at least 30.
C LIW    = declared length of IWORK (in user dimension).
C JAC    = name of subroutine for Jacobian matrix (MF = 121).
C          If used, this name must be declared External in calling
C          program.  If not used, pass a dummy name.
C MF     = method flag.  Standard values are:
C          10  for nonstiff (Adams) method, no Jacobian used
C          121 for stiff (BDF) method, user-supplied sparse Jacobian
C          222 for stiff method, internally generated sparse Jacobian
C Note that the main program must declare arrays Y, RWORK, IWORK,
C and possibly ATOL.
C
C E. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if SLSODES was successful, negative otherwise.
C          -1 means excess work done on this call (perhaps wrong MF).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of MF or tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means a fatal error return flag came from sparse solver
C             SCDRV by way of SPRJS or SSOLSS.  Should never happen.
C          A return with ISTATE = -1, -4, or -5 may result from using
C          an inappropriate sparsity structure, one that is quite
C          different from the initial structure.  Consider calling
C          SLSODES again with ISTATE = 3 to force the structure to be
C          reevaluated.  See the full description of ISTATE below.
C
C F. To continue the integration after a successful return, simply
C reset TOUT and call SLSODES again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is a simple example problem, with the coding
C needed for its solution by SLSODES.  The problem is from chemical
C kinetics, and consists of the following 12 rate equations:
C    dy1/dt  = -rk1*y1
C    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
C                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
C    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
C                + rk11*rk14*y4 + rk12*rk14*y6
C    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
C    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
C    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
C    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
C    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
C    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
C    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
C                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
C                - rk6*y10 - rk9*y10
C    dy11/dt = rk10*y8
C    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
C                - rk15*y2*y12 - rk17*y10*y12
C
C with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
C      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
C      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
C      rk15 = rk17 = 100.0.
C
C The t interval is from 0 to 1000, and the initial conditions
C are y1 = 1, y2 = y3 = ... = y12 = 0.  The problem is stiff.
C
C The following coding solves this problem with SLSODES, using MF = 121
C and printing results at t = .1, 1., 10., 100., 1000.  It uses
C ITOL = 1 and mixed relative/absolute tolerance controls.
C During the run and at the end, statistical quantities of interest
C are printed (see optional outputs in the full description below).
C
C     EXTERNAL FEX, JEX
C     REAL ATOL, RTOL, RWORK, T, TOUT, Y
C     DIMENSION Y(12), RWORK(500), IWORK(30)
C     DATA LRW/500/, LIW/30/
C     NEQ = 12
C     DO 10 I = 1,NEQ
C 10    Y(I) = 0.0E0
C     Y(1) = 1.0E0
C     T = 0.0E0
C     TOUT = 0.1E0
C     ITOL = 1
C     RTOL = 1.0E-4
C     ATOL = 1.0E-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     MF = 121
C     DO 40 IOUT = 1,5
C       CALL SLSODES (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL,
C    1     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
C       WRITE(6,30)T,IWORK(11),RWORK(11),(Y(I),I=1,NEQ)
C 30    FORMAT(//' At t =',E11.3,4X,
C    1    ' No. steps =',I5,4X,' Last step =',E11.3/
C    2    '  Y array =  ',4E14.5/13X,4E14.5/13X,4E14.5)
C       IF (ISTATE .LT. 0) GO TO 80
C       TOUT = TOUT*10.0E0
C 40    CONTINUE
C     LENRW = IWORK(17)
C     LENIW = IWORK(18)
C     NST = IWORK(11)
C     NFE = IWORK(12)
C     NJE = IWORK(13)
C     NLU = IWORK(21)
C     NNZ = IWORK(19)
C     NNZLU = IWORK(25) + IWORK(26) + NEQ
C     WRITE (6,70) LENRW,LENIW,NST,NFE,NJE,NLU,NNZ,NNZLU
C 70  FORMAT(//' Required RWORK size =',I4,'   IWORK size =',I4/
C    1   ' No. steps =',I4,'   No. f-s =',I4,'   No. J-s =',I4,
C    2   '   No. LU-s =',I4/' No. of nonzeros in J =',I5,
C    3   '   No. of nonzeros in LU =',I5)
C     STOP
C 80  WRITE(6,90)ISTATE
C 90  FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE FEX (NEQ, T, Y, YDOT)
C     REAL T, Y, YDOT
C     REAL RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
C    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
C     DIMENSION Y(12), YDOT(12)
C     DATA RK1/0.1E0/, RK2/10.0E0/, RK3/50.0E0/, RK4/2.5E0/, RK5/0.1E0/,
C    1   RK6/10.0E0/, RK7/50.0E0/, RK8/2.5E0/, RK9/50.0E0/, RK10/5.0E0/,
C    2   RK11/50.0E0/, RK12/50.0E0/, RK13/50.0E0/, RK14/30.0E0/,
C    3   RK15/100.0E0/, RK16/2.5E0/, RK17/100.0E0/, RK18/2.5E0/,
C    4   RK19/50.0E0/, RK20/50.0E0/
C     YDOT(1)  = -RK1*Y(1)
C     YDOT(2)  = RK1*Y(1) + RK11*RK14*Y(4) + RK19*RK14*Y(5)
C    1           - RK3*Y(2)*Y(3) - RK15*Y(2)*Y(12) - RK2*Y(2)
C     YDOT(3)  = RK2*Y(2) - RK5*Y(3) - RK3*Y(2)*Y(3) - RK7*Y(10)*Y(3)
C    1           + RK11*RK14*Y(4) + RK12*RK14*Y(6)
C     YDOT(4)  = RK3*Y(2)*Y(3) - RK11*RK14*Y(4) - RK4*Y(4)
C     YDOT(5)  = RK15*Y(2)*Y(12) - RK19*RK14*Y(5) - RK16*Y(5)
C     YDOT(6)  = RK7*Y(10)*Y(3) - RK12*RK14*Y(6) - RK8*Y(6)
C     YDOT(7)  = RK17*Y(10)*Y(12) - RK20*RK14*Y(7) - RK18*Y(7)
C     YDOT(8)  = RK9*Y(10) - RK13*RK14*Y(8) - RK10*Y(8)
C     YDOT(9)  = RK4*Y(4) + RK16*Y(5) + RK8*Y(6) + RK18*Y(7)
C     YDOT(10) = RK5*Y(3) + RK12*RK14*Y(6) + RK20*RK14*Y(7)
C    1           + RK13*RK14*Y(8) - RK7*Y(10)*Y(3) - RK17*Y(10)*Y(12)
C    2           - RK6*Y(10) - RK9*Y(10)
C     YDOT(11) = RK10*Y(8)
C     YDOT(12) = RK6*Y(10) + RK19*RK14*Y(5) + RK20*RK14*Y(7)
C    1           - RK15*Y(2)*Y(12) - RK17*Y(10)*Y(12)
C     RETURN
C     END
C
C     SUBROUTINE JEX (NEQ, T, Y, J, IA, JA, PDJ)
C     REAL T, Y, PDJ
C     REAL RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
C    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
C     DIMENSION Y(12), IA(*), JA(*), PDJ(12)
C     DATA RK1/0.1E0/, RK2/10.0E0/, RK3/50.0E0/, RK4/2.5E0/, RK5/0.1E0/,
C    1   RK6/10.0E0/, RK7/50.0E0/, RK8/2.5E0/, RK9/50.0E0/, RK10/5.0E0/,
C    2   RK11/50.0E0/, RK12/50.0E0/, RK13/50.0E0/, RK14/30.0E0/,
C    3   RK15/100.0E0/, RK16/2.5E0/, RK17/100.0E0/, RK18/2.5E0/,
C    4   RK19/50.0E0/, RK20/50.0E0/
C     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), J
C 1   PDJ(1) = -RK1
C     PDJ(2) = RK1
C     RETURN
C 2   PDJ(2) = -RK3*Y(3) - RK15*Y(12) - RK2
C     PDJ(3) = RK2 - RK3*Y(3)
C     PDJ(4) = RK3*Y(3)
C     PDJ(5) = RK15*Y(12)
C     PDJ(12) = -RK15*Y(12)
C     RETURN
C 3   PDJ(2) = -RK3*Y(2)
C     PDJ(3) = -RK5 - RK3*Y(2) - RK7*Y(10)
C     PDJ(4) = RK3*Y(2)
C     PDJ(6) = RK7*Y(10)
C     PDJ(10) = RK5 - RK7*Y(10)
C     RETURN
C 4   PDJ(2) = RK11*RK14
C     PDJ(3) = RK11*RK14
C     PDJ(4) = -RK11*RK14 - RK4
C     PDJ(9) = RK4
C     RETURN
C 5   PDJ(2) = RK19*RK14
C     PDJ(5) = -RK19*RK14 - RK16
C     PDJ(9) = RK16
C     PDJ(12) = RK19*RK14
C     RETURN
C 6   PDJ(3) = RK12*RK14
C     PDJ(6) = -RK12*RK14 - RK8
C     PDJ(9) = RK8
C     PDJ(10) = RK12*RK14
C     RETURN
C 7   PDJ(7) = -RK20*RK14 - RK18
C     PDJ(9) = RK18
C     PDJ(10) = RK20*RK14
C     PDJ(12) = RK20*RK14
C     RETURN
C 8   PDJ(8) = -RK13*RK14 - RK10
C     PDJ(10) = RK13*RK14
C     PDJ(11) = RK10
C 9   RETURN
C 10  PDJ(3) = -RK7*Y(3)
C     PDJ(6) = RK7*Y(3)
C     PDJ(7) = RK17*Y(12)
C     PDJ(8) = RK9
C     PDJ(10) = -RK7*Y(3) - RK17*Y(12) - RK6 - RK9
C     PDJ(12) = RK6 - RK17*Y(12)
C 11  RETURN
C 12  PDJ(2) = -RK15*Y(2)
C     PDJ(5) = RK15*Y(2)
C     PDJ(7) = RK17*Y(10)
C     PDJ(10) = -RK17*Y(10)
C     PDJ(12) = -RK15*Y(2) - RK17*Y(10)
C     RETURN
C     END
C
C The output of this program (on a Cray-1 in single precision)
C is as follows:
C
C
C At t =  1.000e-01     No. steps =   12     Last step =  1.515e-02
C  Y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
C                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
C                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
C
C
C At t =  1.000e+00     No. steps =   33     Last step =  7.880e-02
C  Y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
C                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
C                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
C
C
C At t =  1.000e+01     No. steps =   48     Last step =  1.239e+00
C  Y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
C                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
C                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
C
C
C At t =  1.000e+02     No. steps =   91     Last step =  3.764e+00
C  Y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
C                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
C                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
C
C
C At t =  1.000e+03     No. steps =  111     Last step =  4.156e+02
C  Y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
C               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
C                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
C
C
C Required RWORK size = 442   IWORK size =  30
C No. steps = 111   No. f-s = 142   No. J-s =   2   No. LU-s =  20
C No. of nonzeros in J =   44   No. of nonzeros in LU =   50
C
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSODES.
C
C The user interface to SLSODES consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODES, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODES package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the SLSODES package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODES to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C F      = the name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter y(1),...,y(NEQ).
C          F must be declared External in the calling program.
C
C          Subroutine F may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in F) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y below.
C
C          If quantities computed in the F routine are needed
C          externally to SLSODES, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use SINTDY instead.
C
C NEQ    = the size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in F and/or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODES package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to F and JAC.  Hence, if it is an array, locations
C          NEQ(2),... may be used to store other integer data and pass
C          it to F and/or JAC.  Subroutines F and/or JAC must include
C          NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          on the first call, Y must contain the vector of initial
C          values.  On output, Y contains the computed solution vector,
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to
C          F and JAC.  Hence its length may exceed NEQ, and locations
C          Y(NEQ+1),... may be used to store other real data and
C          pass it to F and/or JAC.  (The SLSODES package accesses only
C          Y(1),...,Y(NEQ).)
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          on output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C             the conditional inputs IA and JA,
C             and any of the optional inputs except H0.
C             In particular, if MITER = 1 or 2, a call with ISTATE = 3
C             will cause the sparsity structure of the problem to be
C             recomputed (or reread from IA and JA if MOSS = 0).
C          Note:  a preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix,
C              if one is being used.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means a fatal error return flag came from the sparse
C              solver SCDRV by way of SPRJS or SSOLSS (numerical
C              factorization or backsolve).  This should never happen.
C              The integration was successful as far as T.
C
C          Note: an error return with ISTATE = -1, -4, or -5 and with
C          MITER = 1 or 2 may mean that the sparsity structure of the
C          problem has changed significantly since it was last
C          determined (or input).  In that case, one can attempt to
C          complete the integration by setting ISTATE = 3 on the next
C          call, so that a new structure determination is done.
C
C          Note:  since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a work array used for a mixture of real (single precision)
C          and integer work space.
C          The length of RWORK (in real words) must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LWM = 0                                    if MITER = 0,
C          LWM = 2*NNZ + 2*NEQ + (NNZ+9*NEQ)/LENRAT   if MITER = 1,
C          LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT  if MITER = 2,
C          LWM = NEQ + 2                              if MITER = 3.
C          In the above formulas,
C          NNZ    = number of nonzero elements in the Jacobian matrix.
C          LENRAT = the real to integer wordlength ratio (usually 1 in
C                   single precision and 2 in double precision).
C          (See the MF description for METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          the minimum length of RWORK is:
C             20 + 16*NEQ        for MF = 10,
C             20 + 16*NEQ + LWM  for MF = 11, 111, 211, 12, 112, 212,
C             22 + 17*NEQ        for MF = 13,
C             20 +  9*NEQ        for MF = 20,
C             20 +  9*NEQ + LWM  for MF = 21, 121, 221, 22, 122, 222,
C             22 + 10*NEQ        for MF = 23.
C          If MITER = 1 or 2, the above formula for LWM is only a
C          crude lower bound.  The required length of RWORK cannot
C          be readily predicted in general, as it depends on the
C          sparsity structure of the problem.  Some experimentation
C          may be necessary.
C
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C             31 + NEQ + NNZ   if MOSS = 0 and MITER = 1 or 2, or
C             30               otherwise.
C          (NNZ is the number of nonzero elements in df/dy.)
C
C          In SLSODES, IWORK is used only for conditional and
C          optional inputs and optional outputs.
C
C          The following two blocks of words in IWORK are conditional
C          inputs, required if MOSS = 0 and MITER = 1 or 2, but not
C          otherwise (see the description of MF for MOSS).
C            IWORK(30+j) = IA(j)     (j=1,...,NEQ+1)
C            IWORK(31+NEQ+k) = JA(k) (k=1,...,NNZ)
C          The two arrays IA and JA describe the sparsity structure
C          to be assumed for the Jacobian matrix.  JA contains the row
C          indices where nonzero elements occur, reading in columnwise
C          order, and IA contains the starting locations in JA of the
C          descriptions of columns 1,...,NEQ, in that order, with
C          IA(1) = 1.  Thus, for each column index j = 1,...,NEQ, the
C          values of the row index i in column j where a nonzero
C          element may occur are given by
C            i = JA(k),  where   IA(j) .le. k .lt. IA(j+1).
C          If NNZ is the total number of nonzero locations assumed,
C          then the length of the JA array is NNZ, and IA(NEQ+1) must
C          be NNZ + 1.  Duplicate entries are not allowed.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to SLSODES
C for the same problem, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODES between calls, if
C desired (but not for use by F or JAC).
C
C JAC    = name of user-supplied routine (MITER = 1 or MOSS = 1) to
C          compute the Jacobian matrix, df/dy, as a function of
C          the scalar t and the vector y.  It is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
C               DIMENSION Y(*), IAN(*), JAN(*), PDJ(*)
C          where NEQ, T, Y, J, IAN, and JAN are input, and the array
C          PDJ, of length NEQ, is to be loaded with column J
C          of the Jacobian on output.  Thus df(i)/dy(J) is to be
C          loaded into PDJ(i) for all relevant values of i.
C          Here T and Y have the same meaning as in Subroutine F,
C          and J is a column index (1 to NEQ).  IAN and JAN are
C          undefined in calls to JAC for structure determination
C          (MOSS = 1).  otherwise, IAN and JAN are structure
C          descriptors, as defined under optional outputs below, and
C          so can be used to determine the relevant row indices i, if
C          desired.
C               JAC need not provide df/dy exactly.  A crude
C          approximation (possibly with greater sparsity) will do.
C               In any case, PDJ is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Calls to JAC are made with J = 1,...,NEQ, in that order, and
C          each such set of calls is preceded by a call to F with the
C          same arguments NEQ, T, and Y.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by F and not recomputed by JAC,
C          if desired.  JAC must not alter its input arguments.
C          JAC must be declared External in the calling program.
C               Subroutine JAC may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C MF     = the method flag.  Used only for input.
C          MF has three decimal digits-- MOSS, METH, MITER--
C             MF = 100*MOSS + 10*METH + MITER.
C          MOSS indicates the method to be used to obtain the sparsity
C          structure of the Jacobian matrix if MITER = 1 or 2:
C            MOSS = 0 means the user has supplied IA and JA
C                     (see descriptions under IWORK above).
C            MOSS = 1 means the user has supplied JAC (see below)
C                     and the structure will be obtained from NEQ
C                     initial calls to JAC.
C            MOSS = 2 means the structure will be obtained from NEQ+1
C                     initial calls to F.
C          METH indicates the basic linear multistep method:
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on Backward
C                     Differentiation Formulas (BDFs).
C          MITER indicates the corrector iteration method:
C            MITER = 0 means functional iteration (no Jacobian matrix
C                      is involved).
C            MITER = 1 means chord iteration with a user-supplied
C                      sparse Jacobian, given by Subroutine JAC.
C            MITER = 2 means chord iteration with an internally
C                      generated (difference quotient) sparse Jacobian
C                      (using NGP extra calls to F per df/dy value,
C                      where NGP is an optional output described below.)
C            MITER = 3 means chord iteration with an internally
C                      generated diagonal Jacobian approximation
C                      (using 1 extra call to F per df/dy evaluation).
C          If MITER = 1 or MOSS = 1, the user must supply a Subroutine
C          JAC (the name is arbitrary) as described above under JAC.
C          Otherwise, a dummy argument can be used.
C
C          The standard choices for MF are:
C            MF = 10  for a nonstiff problem,
C            MF = 21 or 22 for a stiff problem with IA/JA supplied
C                     (21 if JAC is supplied, 22 if not),
C            MF = 121 for a stiff problem with JAC supplied,
C                     but not IA/JA,
C            MF = 222 for a stiff problem with neither IA/JA nor
C                     JAC supplied.
C          The sparseness structure can be changed during the
C          problem by making a call to SLSODES with ISTATE = 3.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C SETH    RWORK(8)  the element threshhold for sparsity determination
C                   when MOSS = 1 or 2.  If the absolute value of
C                   an estimated Jacobian element is .le. SETH, it
C                   will be assumed to be absent in the structure.
C                   The default value of SETH is 0.
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODES, the variables listed
C below are quantities related to the performance of SLSODES
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODES, and on any return with
C ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NFE     IWORK(12) the number of f evaluations for the problem so far,
C                   excluding those for structure determination
C                   (MOSS = 2).
C
C NJE     IWORK(13) the number of Jacobian evaluations for the problem
C                   so far, excluding those for structure determination
C                   (MOSS = 1).
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NNZ     IWORK(19) the number of nonzero elements in the Jacobian
C                   matrix, including the diagonal (MITER = 1 or 2).
C                   (This may differ from that given by IA(NEQ+1)-1
C                   if MOSS = 0, because of added diagonal entries.)
C
C NGP     IWORK(20) the number of groups of column indices, used in
C                   difference quotient Jacobian aproximations if
C                   MITER = 2.  This is also the number of extra f
C                   evaluations needed for each Jacobian evaluation.
C
C NLU     IWORK(21) the number of sparse LU decompositions for the
C                   problem so far.
C
C LYH     IWORK(22) the base address in RWORK of the history array YH,
C                   described below in this list.
C
C IPIAN   IWORK(23) the base address of the structure descriptor array
C                   IAN, described below in this list.
C
C IPJAN   IWORK(24) the base address of the structure descriptor array
C                   JAN, described below in this list.
C
C NZL     IWORK(25) the number of nonzero elements in the strict lower
C                   triangle of the LU factorization used in the chord
C                   iteration (MITER = 1 or 2).
C
C NZU     IWORK(26) the number of nonzero elements in the strict upper
C                   triangle of the LU factorization used in the chord
C                   iteration (MITER = 1 or 2).
C                   The total number of nonzeros in the factorization
C                   is therefore NZL + NZU + NEQ.
C
C The following four arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address, and its description.
C For YH and ACOR, the base addresses are in RWORK (a real array).
C The integer arrays IAN and JAN are to be obtained by declaring an
C integer array IWK and identifying IWK(1) with RWORK(21), using either
C an equivalence statement or a subroutine call.  Then the base
C addresses IPIAN (of IAN) and IPJAN (of JAN) in IWK are to be obtained
C as optional outputs IWORK(23) and IWORK(24), respectively.
C Thus IAN(1) is IWK(IPIAN), etc.
C
C Name    Base Address      Description
C
C IAN    IPIAN (in IWK)  structure descriptor array of size NEQ + 1.
C JAN    IPJAN (in IWK)  structure descriptor array of size NNZ.
C         (see above)    IAN and JAN together describe the sparsity
C                        structure of the Jacobian matrix, as used by
C                        SLSODES when MITER = 1 or 2.
C                        JAN contains the row indices of the nonzero
C                        locations, reading in columnwise order, and
C                        IAN contains the starting locations in JAN of
C                        the descriptions of columns 1,...,NEQ, in
C                        that order, with IAN(1) = 1.  Thus for each
C                        j = 1,...,NEQ, the row indices i of the
C                        nonzero locations in column j are
C                        i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).
C                        Note that IAN(NEQ+1) = NNZ + 1.
C                        (If MOSS = 0, IAN/JAN may differ from the
C                        input IA/JA because of a different ordering
C                        in each column, and added diagonal entries.)
C
C YH      LYH            the Nordsieck history array, of size NYH by
C          (optional     (NQCUR + 1), where NYH is the initial value
C           output)      of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.  The base address LYH
C                        is another optional output, listed above.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output
C                        to represent the estimated local error in y
C                        on the last step.  This is the vector E  in
C                        the description of the error control.  It is
C                        defined only on a successful return from
C                        SLSODES.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODES.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSODES, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSODES.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCMS(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODES (see Part 3 below).
C                             RSAV must be a real array of length 224
C                             or more, and ISAV must be an integer
C                             array of length 71 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCMS is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODES.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODES.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   LYH = IWORK(22)
C   CALL SINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODES).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (See optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by SLSODES directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C LYH       = the base address of the history array YH, obtained
C             as an optional output as shown above.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODES is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODES, and
C   (2) the two internal Common blocks
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C         /SLSS01/  of length  40  (6 single precision words
C                      followed by 34 integer words),
C
C If SLSODES is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODES is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODES call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODES call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCMS (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the SLSODES package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     Subroutine SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODES call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the SVNORM
C routine (see below), and also used by SLSODES in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) SVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = SVNORM (N, V, W)
C where
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine SEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of SVNORM suitable for use in the error control in SLSODES.
C None of the arguments should be altered by SVNORM.
C For example, a user-supplied SVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19810120  DATE WRITTEN
C 19820315  Upgraded MDI in ODRV package: operates on M + M-transpose.
C 19820426  Numerous revisions in use of work arrays;
C           use wordlength ratio LENRAT; added IPISP & LRAT to Common;
C           added optional outputs IPIAN/IPJAN;
C           numerous corrections to comments.
C 19830503  Added routine CNTNZU; added NZL and NZU to /LSS001/;
C           changed ADJLR call logic; added optional outputs NZL & NZU;
C           revised counter initializations; revised PREP stmt. numbers;
C           corrections to comments throughout.
C 19870320  Corrected jump on test of umax in CDRV routine;
C           added ISTATE = -7 return.
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODE;
C           in STODE, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           converted arithmetic IF statements to logical IF statements;
C           numerous corrections to prologues and internal comments.
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODES package.
C
C In addition to Subroutine SLSODES, the SLSODES package includes the
C following subroutines and function routines:
C  SIPREP   acts as an iterface between SLSODES and SPREP, and also does
C           adjusting of work space pointers and work arrays.
C  SPREP    is called by SIPREP to compute sparsity and do sparse matrix
C           preprocessing if MITER = 1 or 2.
C  JGROUP   is called by SPREP to compute groups of Jacobian column
C           indices for use when MITER = 2.
C  ADJLR    adjusts the length of required sparse matrix work space.
C           It is called by SPREP.
C  CNTNZU   is called by SPREP and counts the nonzero elements in the
C           strict upper triangle of J + J-transpose, where J = df/dy.
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODE   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPRJS    computes and preprocesses the Jacobian matrix J = df/dy
C           and the Newton iteration matrix P = I - h*l0*J.
C  SSOLSS   manages solution of linear system in chord iteration.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted RMS-norm of a vector.
C  SSRCMS   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SODRV    constructs a reordering of the rows and columns of
C           a matrix by the minimum degree algorithm.  SODRV is a
C           driver routine which calls Subroutines MD, MDI, MDM,
C           MDP, MDU, and SSRO.  See Ref. 2 for details.  (The SODRV
C           module has been modified since Ref. 2, however.)
C  SCDRV    performs reordering, symbolic factorization, numerical
C           factorization, or linear system solution operations,
C           depending on a path argument ipath.  SCDRV is a
C           driver routine which calls Subroutines SNROC, SNSFC,
C           SNNFC, SNNSC, and SNNTC.  See Ref. 3 for details.
C           SLSODES uses SCDRV to solve linear systems in which the
C           coefficient matrix is  P = I - con*J, where I is the
C           identity, con is a scalar, and J is an approximation to
C           the Jacobian df/dy.  Because SCDRV deals with rowwise
C           sparsity descriptions, SCDRV works with P-transpose, not P.
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SVNORM, RUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL SPRJS, SSOLSS
      REAL RUMACH, SVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, I1, I2, IFLAG, IMAX, IMUL, IMXER, IPFLAG, IPGO, IREM,
     1   J, KGO, LENRAT, LENYHT, LENIW, LENRW, LF0, LIA, LJA,
     2   LRTEM, LWTEM, LYHD, LYHN, MF1, MORD, MXHNL0, MXSTP0, NCOLM
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE LENRAT, MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following two internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODES, SIPREP, SPREP,
C SINTDY, SSTODE, SPRJS, and SSOLSS.
C The block SLSS01 is declared in subroutines SLSODES, SIPREP, SPREP,
C SPRJS, and SSOLSS.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /SLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH,
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C In the Data statement below, set LENRAT equal to the ratio of
C the wordlength for a real number to that for an integer.  Usually,
C LENRAT = 1 for single precision and 2 for double precision.  If the
C true ratio is not an integer, use the next smaller integer (.ge. 1).
C-----------------------------------------------------------------------
      DATA LENRAT/1/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C If ISTATE = 1, the final setting of work space pointers, the matrix
C preprocessing, and other initializations are done in Block C.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      MOSS = MF/100
      MF1 = MF - 100*MOSS
      METH = MF1/10
      MITER = MF1 - 10*METH
      IF (MOSS .LT. 0 .OR. MOSS .GT. 2) GO TO 608
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 3) GO TO 608
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) MOSS = 0
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      SETH = 0.0E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
      SETH = RWORK(8)
      IF (SETH .LT. 0.0E0) GO TO 609
C Check RTOL and ATOL for legality. ------------------------------------
 60   RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 65 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 65     CONTINUE
C-----------------------------------------------------------------------
C Compute required work array lengths, as far as possible, and test
C these against LRW and LIW.  Then set tentative pointers for work
C arrays.  Pointers to RWORK/IWORK segments are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  WM, YH, SAVF, EWT, ACOR.
C If MITER = 1 or 2, the required length of the matrix work space WM
C is not yet known, and so a crude minimum value is used for the
C initial tests of LRW and LIW, and YH is temporarily stored as far
C to the right in RWORK as possible, to leave the maximum amount
C of space for WM for matrix preprocessing.  Thus if MITER = 1 or 2
C and MOSS .ne. 2, some of the segments of RWORK are temporarily
C omitted, as they are not needed in the preprocessing.  These
C omitted segments are: ACOR if ISTATE = 1, EWT and ACOR if ISTATE = 3
C and MOSS = 1, and SAVF, EWT, and ACOR if ISTATE = 3 and MOSS = 0.
C-----------------------------------------------------------------------
      LRAT = LENRAT
      IF (ISTATE .EQ. 1) NYH = N
      LWMIN = 0
      IF (MITER .EQ. 1) LWMIN = 4*N + 10*N/LRAT
      IF (MITER .EQ. 2) LWMIN = 4*N + 11*N/LRAT
      IF (MITER .EQ. 3) LWMIN = N + 2
      LENYH = (MAXORD+1)*NYH
      LREST = LENYH + 3*N
      LENRW = 20 + LWMIN + LREST
      IWORK(17) = LENRW
      LENIW = 30
      IF (MOSS .EQ. 0 .AND. MITER .NE. 0 .AND. MITER .NE. 3)
     1   LENIW = LENIW + N + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
      LIA = 31
      IF (MOSS .EQ. 0 .AND. MITER .NE. 0 .AND. MITER .NE. 3)
     1   LENIW = LENIW + IWORK(LIA+N) - 1
      IWORK(18) = LENIW
      IF (LENIW .GT. LIW) GO TO 618
      LJA = LIA + N + 1
      LIA = MIN(LIA,LIW)
      LJA = MIN(LJA,LIW)
      LWM = 21
      IF (ISTATE .EQ. 1) NQ = 1
      NCOLM = MIN(NQ+1,MAXORD+2)
      LENYHM = NCOLM*NYH
      LENYHT = LENYH
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENYHT = LENYHM
      IMUL = 2
      IF (ISTATE .EQ. 3) IMUL = MOSS
      IF (MOSS .EQ. 2) IMUL = 3
      LRTEM = LENYHT + IMUL*N
      LWTEM = LWMIN
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LWTEM = LRW - 20 - LRTEM
      LENWK = LWTEM
      LYHN = LWM + LWTEM
      LSAVF = LYHN + LENYHT
      LEWT = LSAVF + N
      LACOR = LEWT + N
      ISTATC = ISTATE
      IF (ISTATE .EQ. 1) GO TO 100
C-----------------------------------------------------------------------
C ISTATE = 3.  Move YH to its new location.
C Note that only the part of YH needed for the next step, namely
C MIN(NQ+1,MAXORD+2) columns, is actually moved.
C A temporary error weight array EWT is loaded if MOSS = 2.
C Sparse matrix processing is done in SIPREP/SPREP if MITER = 1 or 2.
C If MAXORD was reduced below NQ, then the pointers are finally set
C so that SAVF is identical to YH(*,MAXORD+2).
C-----------------------------------------------------------------------
      LYHD = LYH - LYHN
      IMAX = LYHN - 1 + LENYHM
C Move YH.  Move right if LYHD < 0; move left if LYHD > 0. -------------
      IF (LYHD .LT. 0) THEN
        DO 72 I = LYHN,IMAX
          J = IMAX + LYHN - I
 72       RWORK(J) = RWORK(J+LYHD)
      ENDIF
      IF (LYHD .GT. 0) THEN
        DO 76 I = LYHN,IMAX
 76       RWORK(I) = RWORK(I+LYHD)
      ENDIF
 80   LYH = LYHN
      IWORK(22) = LYH
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) GO TO 92
      IF (MOSS .NE. 2) GO TO 85
C Temporarily load EWT if MITER = 1 or 2 and MOSS = 2. -----------------
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 82 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 82     RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 85   CONTINUE
C SIPREP and SPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
      LSAVF = MIN(LSAVF,LRW)
      LEWT = MIN(LEWT,LRW)
      LACOR = MIN(LACOR,LRW)
      CALL SIPREP (NEQ, Y, RWORK, IWORK(LIA),IWORK(LJA), IPFLAG, F, JAC)
      LENRW = LWM - 1 + LENWK + LREST
      IWORK(17) = LENRW
      IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
      IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
      IPGO = -IPFLAG + 1
      GO TO (90, 628, 629, 630, 631, 632, 633), IPGO
 90   IWORK(22) = LYH
      IF (LENRW .GT. LRW) GO TO 617
C Set flag to signal parameter changes to SSTODE. ----------------------
 92   JSTART = -1
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C the sparse matrix preprocessing (MITER = 1 or 2), and the
C calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  CONTINUE
      LYH = LYHN
      IWORK(22) = LYH
      TN = T
      NST = 0
      H = 1.0E0
      NNZ = 0
      NGP = 0
      NZL = 0
      NZU = 0
C Load the initial value vector in YH. ---------------------------------
      DO 105 I = 1,N
 105    RWORK(I+LYH-1) = Y(I)
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 110 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 110    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) GO TO 120
C SIPREP and SPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
      LACOR = MIN(LACOR,LRW)
      CALL SIPREP (NEQ, Y, RWORK, IWORK(LIA),IWORK(LJA), IPFLAG, F, JAC)
      LENRW = LWM - 1 + LENWK + LREST
      IWORK(17) = LENRW
      IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
      IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
      IPGO = -IPFLAG + 1
      GO TO (115, 628, 629, 630, 631, 632, 633), IPGO
 115  IWORK(22) = LYH
      IF (LENRW .GT. LRW) GO TO 617
C Check TCRIT for legality (ITASK = 4 or 5). ---------------------------
 120  CONTINUE
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 125
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
C Initialize all remaining parameters. ---------------------------------
 125  UROUND = RUMACH()
      JSTART = 0
      IF (MITER .NE. 0) RWORK(LWM) = SQRT(UROUND)
      MSBJ = 50
      NSLJ = 0
      CCMXJ = 0.2E0
      PSMALL = 1000.0E0*UROUND
      RBIG = 0.01E0/PSMALL
      NHNIL = 0
      NJE = 0
      NLU = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( f(i)/ywt(i) )**2  )
C                                       1
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         f(i)   = i-th component of initial value of f,
C         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C ABS(H0) is made .le. ABS(TOUT-T) in any case.
C-----------------------------------------------------------------------
      LF0 = LYH + NYH
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODE.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODES- Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODES- Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C    CALL SSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,WM,F,JAC,SPRJS,SSOLSS)
C-----------------------------------------------------------------------
      CALL SSTODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), RWORK(LWM),
     2   F, JAC, SPRJS, SSOLSS)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 550), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  if TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODES.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNZ
      IWORK(20) = NGP
      IWORK(21) = NLU
      IWORK(25) = NZL
      IWORK(26) = NZU
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODES- At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODES- At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODES- At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODES- At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODES- At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C KFLAG = -3.  Fatal error flag returned by SPRJS or SSOLSS (SCDRV). ---
 550  MSG = 'SLSODES- At T (=R1) and step size H (=R2), a fatal'
      CALL XERRWV (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      error flag was returned by SCDRV (by way of '
      CALL XERRWV (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      Subroutine SPRJS or SSOLSS)       '
      CALL XERRWV (MSG, 40, 207, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0E0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNZ
      IWORK(20) = NGP
      IWORK(21) = NLU
      IWORK(25) = NZL
      IWORK(26) = NZU
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODES- ISTATE (=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODES- ITASK (=I1) illegal. '
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODES- ISTATE.gt.1 but SLSODES not initialized. '
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODES- NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODES- ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODES- ITOL (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODES- IOPT (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODES- MF (=I1) illegal.    '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 609  MSG = 'SLSODES- SETH (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 9, 0, 0, 0, 0, 1, SETH, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODES- MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODES- MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODES- MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODES- TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODES- HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODES- HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG = 'SLSODES- RWORK length is insufficient to proceed. '
      CALL XERRWV (MSG, 50, 17, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG = 'SLSODES- IWORK length is insufficient to proceed. '
      CALL XERRWV (MSG, 50, 18, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODES- RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODES- ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODES- EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODES- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODES- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODES- At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODES- Trouble in SINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
      GO TO 700
 628  MSG='SLSODES- RWORK length insufficient (for Subroutine SPREP).  '
      CALL XERRWV (MSG, 60, 28, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 28, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 629  MSG='SLSODES- RWORK length insufficient (for Subroutine JGROUP). '
      CALL XERRWV (MSG, 60, 29, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 29, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 630  MSG='SLSODES- RWORK length insufficient (for Subroutine SODRV).  '
      CALL XERRWV (MSG, 60, 30, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 30, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 631  MSG='SLSODES- Error from SODRV in Yale Sparse Matrix Package.    '
      CALL XERRWV (MSG, 60, 31, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      IMUL = (IYS - 1)/N
      IREM = IYS - IMUL*N
      MSG='      At T (=R1), SODRV returned error flag = I1*NEQ + I2.  '
      CALL XERRWV (MSG, 60, 31, 0, 2, IMUL, IREM, 1, TN, 0.0E0)
      GO TO 700
 632  MSG='SLSODES- RWORK length insufficient (for Subroutine SCDRV).  '
      CALL XERRWV (MSG, 60, 32, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 32, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 633  MSG='SLSODES- Error from SCDRV in Yale Sparse Matrix Package.    '
      CALL XERRWV (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      IMUL = (IYS - 1)/N
      IREM = IYS - IMUL*N
      MSG='      At T (=R1), SCDRV returned error flag = I1*NEQ + I2.  '
      CALL XERRWV (MSG, 60, 33, 0, 2, IMUL, IREM, 1, TN, 0.0E0)
      IF (IMUL .EQ. 2) THEN
      MSG='        Duplicate entry in sparsity structure descriptors.  '
      CALL XERRWV (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (IMUL .EQ. 3 .OR. IMUL .EQ. 6) THEN
      MSG='        Insufficient storage for SNSFC (called by SCDRV).   '
      CALL XERRWV (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODES- Run aborted.. apparent infinite loop.    '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODES ---------------------
      END
*DECK SLSODA
      SUBROUTINE SLSODA (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 12 November 2003 version of
C SLSODA: Livermore Solver for Ordinary Differential Equations, with
C         Automatic method switching for stiff and nonstiff problems.
C
C This version is in single precision.
C
C SLSODA solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C
C This a variant version of the SLSODE package.
C It switches automatically between stiff and nonstiff methods.
C This means that the user does not have to determine whether the
C problem is stiff or not, and the solver will automatically choose the
C appropriate method.  It always starts with the nonstiff method.
C
C Authors:       Alan C. Hindmarsh
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Linda R. Petzold
C                Univ. of California at Santa Barbara
C                Dept. of Computer Science
C                Santa Barbara, CA 93106
C
C References:
C 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C 2.  Linda R. Petzold, Automatic Selection of Methods for Solving
C     Stiff and Nonstiff Systems of Ordinary Differential Equations,
C     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODA package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including alternative treatment of the Jacobian matrix,
C optional inputs and outputs, nonstandard options, and
C instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First provide a subroutine of the form:
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Write a main program which calls Subroutine SLSODA once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by SLSODA.  On the first call to SLSODA, supply arguments as follows:
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared External in calling program.
C NEQ    = number of first order ODEs.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be less than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 + NEQ * MAX(16, NEQ + 9).
C          See also Paragraph E below.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least  20 + NEQ.
C LIW    = declared length of IWORK (in user's dimension).
C JAC    = name of subroutine for Jacobian matrix.
C          Use a dummy name.  See also Paragraph E below.
C JT     = Jacobian type indicator.  Set JT = 2.
C          See also Paragraph E below.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C and possibly ATOL.
C
C C. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if SLSODA was successful, negative otherwise.
C          -1 means excess work done on this call (perhaps wrong JT).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of JT or tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means work space insufficient to finish (see messages).
C
C D. To continue the integration after a successful return, simply
C reset TOUT and call SLSODA again.  No other parameters need be reset.
C
C E. Note: If and when SLSODA regards the problem as stiff, and
C switches methods accordingly, it must make use of the NEQ by NEQ
C Jacobian matrix, J = df/dy.  For the sake of simplicity, the
C inputs to SLSODA recommended in Paragraph B above cause SLSODA to
C treat J as a full matrix, and to approximate it internally by
C difference quotients.  Alternatively, J can be treated as a band
C matrix (with great potential reduction in the size of the RWORK
C array).  Also, in either the full or banded case, the user can supply
C J in closed form, with a routine whose name is passed as the JAC
C argument.  These alternatives are described in the paragraphs on
C RWORK, JAC, and JT in the full description of the call sequence below.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is a simple example problem, with the coding
C needed for its solution by SLSODA.  The problem is from chemical
C kinetics, and consists of the following three rate equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
C
C The following coding solves this problem with SLSODA,
C printing results at t = .4, 4., ..., 4.e10.  It uses
C ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
C y2 has much smaller values.
C At the end of the run, statistical quantities of interest are
C printed (see optional outputs in the full description below).
C
C     EXTERNAL FEX
C     REAL ATOL, RTOL, RWORK, T, TOUT, Y
C     DIMENSION Y(3), ATOL(3), RWORK(70), IWORK(23)
C     NEQ = 3
C     Y(1) = 1.
C     Y(2) = 0.
C     Y(3) = 0.
C     T = 0.
C     TOUT = .4
C     ITOL = 2
C     RTOL = 1.E-4
C     ATOL(1) = 1.E-6
C     ATOL(2) = 1.E-10
C     ATOL(3) = 1.E-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 70
C     LIW = 23
C     JT = 2
C     DO 40 IOUT = 1,12
C       CALL SLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
C    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
C       WRITE(6,20)T,Y(1),Y(2),Y(3)
C 20    FORMAT(' At t =',E12.4,'   Y =',3E14.6)
C       IF (ISTATE .LT. 0) GO TO 80
C 40    TOUT = TOUT*10.
C     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
C 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/
C    1   ' Method last used =',I2,'   Last switch was at t =',E12.4)
C     STOP
C 80  WRITE(6,90)ISTATE
C 90  FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE FEX (NEQ, T, Y, YDOT)
C     REAL T, Y, YDOT
C     DIMENSION Y(3), YDOT(3)
C     YDOT(1) = -.04*Y(1) + 1.E4*Y(2)*Y(3)
C     YDOT(3) = 3.E7*Y(2)*Y(2)
C     YDOT(2) = -YDOT(1) - YDOT(3)
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C   At t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
C   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02
C   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01
C   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01
C   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01
C   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01
C   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01
C   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01
C   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01
C   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01
C   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01
C   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00
C
C   No. steps = 361  No. f-s = 693  No. J-s =  64
C   Method last used = 2   Last switch was at t =  6.0092e-03
C-----------------------------------------------------------------------
C Full description of user interface to SLSODA.
C
C The user interface to SLSODA consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODA, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODA package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of a subroutine in the SLSODA package,
C      which the user may replace with his/her own version, if desired.
C      this relates to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODA to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C F      = the name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared External in the calling program.
C
C          Subroutine F may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in F) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y below.
C
C          If quantities computed in the F routine are needed
C          externally to SLSODA, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use SINTDY instead.
C
C NEQ    = the size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in F and/or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODA package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to F and JAC.  Hence, if it is an array, locations
C          NEQ(2),... may be used to store other integer data and pass
C          it to F and/or JAC.  Subroutines F and/or JAC must include
C          NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  On output, Y contains the computed solution vector,
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to
C          F and JAC.  Hence its length may exceed NEQ, and locations
C          Y(NEQ+1),... may be used to store other real data and
C          pass it to F and/or JAC.  (The SLSODA package accesses only
C          Y(1),...,Y(NEQ).)
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          on output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          on an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial t, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      max-norm of ( E(i)/EWT(i) )   .le.   1,
C          where EWT = (EWT(i)) is a vector of positive error weights.
C          The values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting a
C          user-supplied routine for the setting of EWT.
C          See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
C             and any optional inputs except H0, MXORDN, and MXORDS.
C             (See IWORK description for ML and MU.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix,
C              if one is being used.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means the length of RWORK and/or IWORK was too small to
C              proceed, but the integration was successful as far as T.
C              This happens when SLSODA chooses to switch methods
C              but LRW and/or LIW is too small for the new method.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real array (single precision) for work space, and (in the
C          first 20 words) for conditional and optional inputs and
C          optional outputs.
C          As SLSODA switches automatically between stiff and nonstiff
C          methods, the required length of RWORK can change during the
C          problem.  Thus the RWORK array passed to SLSODA can either
C          have a static (fixed) length large enough for both methods,
C          or have a dynamic (changing) length altered by the calling
C          program in response to output from SLSODA.
C
C                       --- Fixed Length Case ---
C          If the RWORK length is to be fixed, it should be at least
C               MAX (LRN, LRS),
C          where LRN and LRS are the RWORK lengths required when the
C          current method is nonstiff or stiff, respectively.
C
C          The separate RWORK length requirements LRN and LRS are
C          as follows:
C          IF NEQ is constant and the maximum method orders have
C          their default values, then
C             LRN = 20 + 16*NEQ,
C             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
C             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
C          Under any other conditions, LRN and LRS are given by:
C             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ,
C             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT,
C          where
C             NYH    = the initial value of NEQ,
C             MXORDN = 12, unless a smaller value is given as an
C                      optional input,
C             MXORDS = 5, unless a smaller value is given as an
C                      optional input,
C             LMAT   = length of matrix work space:
C             LMAT   = NEQ**2 + 2              if JT = 1 or 2,
C             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.
C
C                       --- Dynamic Length Case ---
C          If the length of RWORK is to be dynamic, then it should
C          be at least LRN or LRS, as defined above, depending on the
C          current method.  Initially, it must be at least LRN (since
C          SLSODA starts with the nonstiff method).  On any return
C          from SLSODA, the optional output MCUR indicates the current
C          method.  If MCUR differs from the value it had on the
C          previous return, or if there has only been one call to
C          SLSODA and MCUR is now 2, then SLSODA has switched
C          methods during the last call, and the length of RWORK
C          should be reset (to LRN if MCUR = 1, or to LRS if
C          MCUR = 2).  (An increase in the RWORK length is required
C          if SLSODA returned ISTATE = -7, but not otherwise.)
C          After resetting the length, call SLSODA with ISTATE = 3
C          to signal that change.
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer array for work space.
C          As SLSODA switches automatically between stiff and nonstiff
C          methods, the required length of IWORK can change during
C          problem, between
C             LIS = 20 + NEQ   and   LIN = 20,
C          respectively.  Thus the IWORK array passed to SLSODA can
C          either have a fixed length of at least 20 + NEQ, or have a
C          dynamic length of at least LIN or LIS, depending on the
C          current method.  The comments on dynamic length under
C          RWORK above apply here.  Initially, this length need
C          only be at least LIN = 20.
C
C          The first few words of IWORK are used for conditional and
C          optional inputs and optional outputs.
C
C          The following 2 words in IWORK are conditional inputs:
C            IWORK(1) = ML     these are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       banded Jacobian, excluding the main diagonal.
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if JT is 4 or 5, and
C                       ignored otherwise.  ML and MU may in fact be
C                       the band parameters for a matrix to which
C                       df/dy is only approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note: The base addresses of the work arrays must not be
C altered between calls to SLSODA for the same problem.
C The contents of the work arrays must not be altered
C between calls, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODA between calls, if
C desired (but not for use by F or JAC).
C
C JAC    = the name of the user-supplied routine to compute the
C          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine
C          is optional, but if the problem is expected to be stiff much
C          of the time, you are encouraged to supply JAC, for the sake
C          of efficiency.  (Alternatively, set JT = 2 or 5 to have
C          SLSODA compute df/dy internally by difference quotients.)
C          If and when SLSODA uses df/dy, it treats this NEQ by NEQ
C          matrix either as full (JT = 1 or 2), or as banded (JT =
C          4 or 5) with half-bandwidths ML and MU (discussed under
C          IWORK above).  In either case, if JT = 1 or 4, the JAC
C          routine must compute df/dy as a function of the scalar t
C          and the vector y.  It is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C               DIMENSION Y(*), PD(NROWPD,*)
C          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
C          PD is to be loaded with partial derivatives (elements of
C          the Jacobian matrix) on output.  PD must be given a first
C          dimension of NROWPD.  T and Y have the same meaning as in
C          Subroutine F.
C               In the full matrix case (JT = 1), ML and MU are
C          ignored, and the Jacobian is to be loaded into PD in
C          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
C               In the band matrix case (JT = 4), the elements
C          within the band are to be loaded into PD in columnwise
C          manner, with diagonal lines of df/dy loaded into the rows
C          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
C          ML and MU are the half-bandwidth parameters (see IWORK).
C          The locations in PD in the two triangular areas which
C          correspond to nonexistent matrix elements can be ignored
C          or loaded arbitrarily, as they are overwritten by SLSODA.
C               JAC need not provide df/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               In either case, PD is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to F with the same
C          arguments NEQ, T, and Y.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by F and not recomputed by JAC,
C          if desired.  Also, JAC may alter the Y array, if desired.
C          JAC must be declared External in the calling program.
C               Subroutine JAC may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C JT     = Jacobian type indicator.  Used only for input.
C          JT specifies how the Jacobian matrix df/dy will be
C          treated, if and when SLSODA requires this matrix.
C          JT has the following values and meanings:
C           1 means a user-supplied full (NEQ by NEQ) Jacobian.
C           2 means an internally generated (difference quotient) full
C             Jacobian (using NEQ extra calls to F per df/dy value).
C           4 means a user-supplied banded Jacobian.
C           5 means an internally generated banded Jacobian (using
C             ML+MU+1 extra calls to F per df/dy evaluation).
C          If JT = 1 or 4, the user must supply a Subroutine JAC
C          (the name is arbitrary) as described above under JAC.
C          If JT = 2 or 5, a dummy argument can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C IXPR    IWORK(5)  flag to generate extra printing at method switches.
C                   IXPR = 0 means no extra printing (the default).
C                   IXPR = 1 means print data on each switch.
C                   T, H, and NST will be printed on the same logical
C                   unit as used for error messages.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
C                   (Adams) method.  the default value is 12.
C                   if MXORDN exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDN is held constant during the problem.
C
C MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
C                   (BDF) method.  The default value is 5.
C                   If MXORDS exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDS is held constant during the problem.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODA, the variables listed
C below are quantities related to the performance of SLSODA
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODA, and on any return with
C ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C TSW     RWORK(15) the value of t at the time of the last method
C                   switch, if any.
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NFE     IWORK(12) the number of f evaluations for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (and of matrix
C                   LU decompositions) for the problem so far.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required, assuming
C                   that the length of RWORK is to be fixed for the
C                   rest of the problem, and that switching may occur.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required, assuming
C                   that the length of IWORK is to be fixed for the
C                   rest of the problem, and that switching may occur.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C MUSED   IWORK(19) the method indicator for the last successful step:
C                   1 means Adams (nonstiff), 2 means BDF (stiff).
C
C MCUR    IWORK(20) the current method indicator:
C                   1 means Adams (nonstiff), 2 means BDF (stiff).
C                   This is the method to be attempted
C                   on the next step.  Thus it differs from MUSED
C                   only if a method switch has just been made.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at T = TCUR.
C
C ACOR     LACOR         array of size NEQ used for the accumulated
C         (from Common   corrections on each step, scaled on output
C           as noted)    to represent the estimated local error in y
C                        on the last step.  This is the vector E in
C                        the description of the error control.  It is
C                        defined only on a successful return from
C                        SLSODA.  The base address LACOR is obtained by
C                        including in the user's program the
C                        following 2 lines:
C                           COMMON /SLS001/ RLS(218), ILS(37)
C                           LACOR = ILS(22)
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODA.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          set the logical unit number, LUN, for
C                             output of messages from SLSODA, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         set a flag to control the printing of
C                             messages by SLSODA.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCMA(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODA (see Part 3 below).
C                             RSAV must be a real array of length 240
C                             or more, and ISAV must be an integer
C                             array of length 46 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCMA is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODA.
C
C   CALL SINTDY(,,,,,)        provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODA.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   CALL SINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODA).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by SLSODA directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODA is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODA, and
C   (2) the two internal Common blocks
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C         /SLSA01/  of length  31    (22 single precision words
C                      followed by  9 integer words).
C
C If SLSODA is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODA is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODA call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODA call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCMA (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below is a description of a routine in the SLSODA package which
C relates to the measurement of errors, and can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     Subroutine SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODA call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the
C SMNORM routine, and also used by SLSODA in the computation
C of the optional output IMXER, and the increments for difference
C quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19811102  DATE WRITTEN
C 19820126  Fixed bug in tests of work space lengths;
C           minor corrections in main prologue and comments.
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODA;
C           in STODA, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODA.
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           numerous corrections to prologues and internal comments.
C 20010613  Revised excess accuracy test (to match rest of ODEPACK).
C 20010808  Fixed bug in SPRJA (matrix in SBNORM call).
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODA package.
C
C In addition to Subroutine SLSODA, the SLSODA package includes the
C following subroutines and function routines:
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODA   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPRJA    computes and preprocesses the Jacobian matrix J = df/dy
C           and the Newton iteration matrix P = I - h*l0*J.
C  SSOLSY   manages solution of linear system in chord iteration.
C  SEWSET   sets the error weight vector EWT before each step.
C  SMNORM   computes the weighted max-norm of a vector.
C  SFNORM   computes the norm of a full matrix consistent with the
C           weighted max-norm on vectors.
C  SBNORM   computes the norm of a band matrix consistent with the
C           weighted max-norm on vectors.
C  SSRCMA   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SGEFA and SGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  SGBFA and SGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SMNORM, SFNORM, SBNORM, RUMACH, IXSAV, and IUMACH are
C function routines.  All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL SPRJA, SSOLSY
      REAL RUMACH, SMNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER INSUFR, INSUFI, IXPR, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0
      INTEGER LEN1, LEN1C, LEN1N, LEN1S, LEN2, LENIWC, LENRWC
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL TSW, ROWNS2, PDNORM
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following two internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODA, SINTDY, SSTODA,
C SPRJA, and SSOLSY.
C The block SLSA01 is declared in subroutines SLSODA, SSTODA, and SPRJA.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /SLSA01/ TSW, ROWNS2(20), PDNORM,
     1   INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C JT, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      IF (JT .EQ. 3 .OR. JT .LT. 1 .OR. JT .GT. 5) GO TO 608
      JTYP = JT
      IF (JT .LE. 2) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      IXPR = 0
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      HMXI = 0.0E0
      HMIN = 0.0E0
      IF (ISTATE .NE. 1) GO TO 60
      H0 = 0.0E0
      MXORDN = MORD(1)
      MXORDS = MORD(2)
      GO TO 60
 40   IXPR = IWORK(5)
      IF (IXPR .LT. 0 .OR. IXPR .GT. 1) GO TO 611
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      MXORDN = IWORK(8)
      IF (MXORDN .LT. 0) GO TO 628
      IF (MXORDN .EQ. 0) MXORDN = 100
      MXORDN = MIN(MXORDN,MORD(1))
      MXORDS = IWORK(9)
      IF (MXORDS .LT. 0) GO TO 629
      IF (MXORDS .EQ. 0) MXORDS = 100
      MXORDS = MIN(MXORDS,MORD(2))
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C If ISTATE = 1, METH is initialized to 1 here to facilitate the
C checking of work space lengths.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
C If the lengths provided are insufficient for the current method,
C an error return occurs.  This is treated as illegal input on the
C first call, but as a problem interruption with ISTATE = -7 on a
C continuation call.  If the lengths are sufficient for the current
C method but not for both methods, a warning message is sent.
C-----------------------------------------------------------------------
 60   IF (ISTATE .EQ. 1) METH = 1
      IF (ISTATE .EQ. 1) NYH = N
      LYH = 21
      LEN1N = 20 + (MXORDN + 1)*NYH
      LEN1S = 20 + (MXORDS + 1)*NYH
      LWM = LEN1S + 1
      IF (JT .LE. 2) LENWM = N*N + 2
      IF (JT .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEN1S = LEN1S + LENWM
      LEN1C = LEN1N
      IF (METH .EQ. 2) LEN1C = LEN1S
      LEN1 = MAX(LEN1N,LEN1S)
      LEN2 = 3*N
      LENRW = LEN1 + LEN2
      LENRWC = LEN1C + LEN2
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      LENIWC = 20
      IF (METH .EQ. 2) LENIWC = LENIW
      IWORK(18) = LENIW
      IF (ISTATE .EQ. 1 .AND. LRW .LT. LENRWC) GO TO 617
      IF (ISTATE .EQ. 1 .AND. LIW .LT. LENIWC) GO TO 618
      IF (ISTATE .EQ. 3 .AND. LRW .LT. LENRWC) GO TO 550
      IF (ISTATE .EQ. 3 .AND. LIW .LT. LENIWC) GO TO 555
      LEWT = LEN1 + 1
      INSUFR = 0
      IF (LRW .GE. LENRW) GO TO 65
      INSUFR = 2
      LEWT = LEN1C + 1
      MSG='SLSODA-  Warning.. RWORK length is sufficient for now, but  '
      CALL XERRWV (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWV (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      Length needed is LENRW = I1, while LRW = I2.'
      CALL XERRWV (MSG, 50, 103, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
 65   LSAVF = LEWT + N
      LACOR = LSAVF + N
      INSUFI = 0
      IF (LIW .GE. LENIW) GO TO 70
      INSUFI = 2
      MSG='SLSODA-  Warning.. IWORK length is sufficient for now, but  '
      CALL XERRWV (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWV (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      Length needed is LENIW = I1, while LIW = I2.'
      CALL XERRWV (MSG, 50, 104, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
 70   CONTINUE
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 75 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 75     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to SSTODA. -------
      JSTART = -1
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      TSW = T
      MAXORD = MXORDN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      MUSED = 0
      MITER = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by:
C
C   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
C
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         F      = the initial value of the vector f(t,y), and
C         norm() = the weighted vector norm used throughout, given by
C                  the SMNORM function routine, and weighted by the
C                  tolerances initially loaded into the EWT array.
C The sign of H0 is inferred from the initial values of TOUT and T.
C ABS(H0) is made .le. ABS(TOUT-T) in any case.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SMNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      T = TN
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) T = TCRIT
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2 .AND. JSTART .GE. 0) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODA.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF (METH .EQ. MUSED) GO TO 255
      IF (INSUFR .EQ. 1) GO TO 550
      IF (INSUFI .EQ. 1) GO TO 555
 255  IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SMNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODA-  Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODA-  Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C   CALL SSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,SPRJA,SSOLSY)
C-----------------------------------------------------------------------
      CALL SSTODA (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, SPRJA, SSOLSY)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).
C If a method switch was just made, record TSW, reset MAXORD,
C set JSTART to -1 to signal SSTODA to complete the switch,
C and do extra printing of data if IXPR = 1.
C Then, in any case, check for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      IF (METH .EQ. MUSED) GO TO 310
      TSW = TN
      MAXORD = MXORDN
      IF (METH .EQ. 2) MAXORD = MXORDS
      IF (METH .EQ. 2) RWORK(LWM) = SQRT(UROUND)
      INSUFR = MIN(INSUFR,1)
      INSUFI = MIN(INSUFI,1)
      JSTART = -1
      IF (IXPR .EQ. 0) GO TO 310
      IF (METH .EQ. 2) THEN
      MSG='SLSODA- A switch to the BDF (stiff) method has occurred     '
      CALL XERRWV (MSG, 60, 105, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      ENDIF
      IF (METH .EQ. 1) THEN
      MSG='SLSODA- A switch to the Adams (nonstiff) method has occurred'
      CALL XERRWV (MSG, 60, 106, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      MSG='     at T = R1,  tentative step size H = R2,  step NST = I1 '
      CALL XERRWV (MSG, 60, 107, 0, 1, NST, 0, 2, TN, H)
 310  GO TO (320, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 320  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (JSTART .GE. 0) JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODA.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODA-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODA-  At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODA-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODA-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C RWORK length too small to proceed. -----------------------------------
 550  MSG = 'SLSODA-  At current T(=R1), RWORK length too small'
      CALL XERRWV (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWV (MSG, 60, 206, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 580
C IWORK length too small to proceed. -----------------------------------
 555  MSG = 'SLSODA-  At current T(=R1), IWORK length too small'
      CALL XERRWV (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWV (MSG, 60, 207, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0E0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODA-  ISTATE (=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODA-  ITASK (=I1) illegal. '
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODA-  ISTATE .gt. 1 but SLSODA not initialized.'
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODA-  NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODA-  ITOL (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODA-  IOPT (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODA-  JT (=I1) illegal.    '
      CALL XERRWV (MSG, 30, 8, 0, 1, JT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 609  MSG = 'SLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
      CALL XERRWV (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 610  MSG = 'SLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
      CALL XERRWV (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODA-  IXPR (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 11, 0, 1, IXPR, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODA-  MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODA-  MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODA-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODA-  HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODA-  HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG='SLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG='SLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODA-  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODA-  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODA-  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODA-  TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODA-  At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODA-  Trouble in SINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
      GO TO 700
 628  MSG = 'SLSODA-  MXORDN (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 28, 0, 1, MXORDN, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 629  MSG = 'SLSODA-  MXORDS (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 29, 0, 1, MXORDS, 0, 0, 0.0E0, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODA-  Run aborted.. apparent infinite loop.    '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODA ----------------------
      END
*DECK SLSODAR
      SUBROUTINE SLSODAR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT,
     2            G, NG, JROOT)
      EXTERNAL F, JAC, G
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT,
     1   NG, JROOT
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1   JROOT(NG)
C-----------------------------------------------------------------------
C This is the 12 November 2003 version of
C SLSODAR: Livermore Solver for Ordinary Differential Equations, with
C          Automatic method switching for stiff and nonstiff problems,
C          and with Root-finding.
C
C This version is in single precision.
C
C SLSODAR solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C At the same time, it locates the roots of any of a set of functions
C     g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).
C
C This a variant version of the SLSODE package.  It differs from it
C in two ways:
C (a) It switches automatically between stiff and nonstiff methods.
C This means that the user does not have to determine whether the
C problem is stiff or not, and the solver will automatically choose the
C appropriate method.  It always starts with the nonstiff method.
C (b) It finds the root of at least one of a set of constraint
C functions g(i) of the independent and dependent variables.
C It finds only those roots for which some g(i), as a function
C of t, changes sign in the interval of integration.
C It then returns the solution at the root, if that occurs
C sooner than the specified stop condition, and otherwise returns
C the solution according the specified stop condition.
C
C Authors:       Alan C. Hindmarsh,
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Linda R. Petzold
C                Univ. of California at Santa Barbara
C                Dept. of Computer Science
C                Santa Barbara, CA 93106
C
C References:
C 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C 2.  Linda R. Petzold, Automatic Selection of Methods for Solving
C     Stiff and Nonstiff Systems of Ordinary Differential Equations,
C     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
C 3.  Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
C     Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
C     February 1980.
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODAR package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including alternative treatment of the Jacobian matrix,
C optional inputs and outputs, nonstandard options, and
C instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First provide a subroutine of the form:
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Provide a subroutine of the form:
C               SUBROUTINE G (NEQ, T, Y, NG, GOUT)
C               DIMENSION Y(*), GOUT(NG)
C which supplies the vector function g by loading GOUT(i) with
C g(i), the i-th constraint function whose root is sought.
C
C C. Write a main program which calls Subroutine SLSODAR once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages by
C SLSODAR.  On the first call to SLSODAR, supply arguments as follows:
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared External in calling program.
C NEQ    = number of first order ODEs.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be less than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 + NEQ * MAX(16, NEQ + 9) + 3*NG.
C          See also Paragraph F below.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least  20 + NEQ.
C LIW    = declared length of IWORK (in user's dimension).
C JAC    = name of subroutine for Jacobian matrix.
C          Use a dummy name.  See also Paragraph F below.
C JT     = Jacobian type indicator.  Set JT = 2.
C          See also Paragraph F below.
C G      = name of subroutine for constraint functions, whose
C          roots are desired during the integration.
C          This name must be declared External in calling program.
C NG     = number of constraint functions g(i).  If there are none,
C          set NG = 0, and pass a dummy name for G.
C JROOT  = integer array of length NG for output of root information.
C          See next paragraph.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C JROOT, and possibly ATOL.
C
C D. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable.  This is
C          TOUT if ISTATE = 2, or the root location if ISTATE = 3,
C          or the farthest point reached if SLSODAR was unsuccessful.
C ISTATE = 2 or 3  if SLSODAR was successful, negative otherwise.
C           2 means no root was found, and TOUT was reached as desired.
C           3 means a root was found prior to reaching TOUT.
C          -1 means excess work done on this call (perhaps wrong JT).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of JT or tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means work space insufficient to finish (see messages).
C JROOT  = array showing roots found if ISTATE = 3 on return.
C          JROOT(i) = 1 if g(i) has a root at t, or 0 otherwise.
C
C E. To continue the integration after a successful return, proceed
C as follows:
C  (a) If ISTATE = 2 on return, reset TOUT and call SLSODAR again.
C  (b) If ISTATE = 3 on return, reset ISTATE to 2, call SLSODAR again.
C In either case, no other parameters need be reset.
C
C F. Note: If and when SLSODAR regards the problem as stiff, and
C switches methods accordingly, it must make use of the NEQ by NEQ
C Jacobian matrix, J = df/dy.  For the sake of simplicity, the
C inputs to SLSODAR recommended in Paragraph C above cause SLSODAR to
C treat J as a full matrix, and to approximate it internally by
C difference quotients.  Alternatively, J can be treated as a band
C matrix (with great potential reduction in the size of the RWORK
C array).  Also, in either the full or banded case, the user can supply
C J in closed form, with a routine whose name is passed as the JAC
C argument.  These alternatives are described in the paragraphs on
C RWORK, JAC, and JT in the full description of the call sequence below.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is a simple example problem, with the coding
C needed for its solution by SLSODAR.  The problem is from chemical
C kinetics, and consists of the following three rate equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
C In addition, we want to find the values of t, y1, y2, and y3 at which
C   (1) y1 reaches the value 1.e-4, and
C   (2) y3 reaches the value 1.e-2.
C
C The following coding solves this problem with SLSODAR,
C printing results at t = .4, 4., ..., 4.e10, and at the computed
C roots.  It uses ITOL = 2 and ATOL much smaller for y2 than y1 or y3
C because y2 has much smaller values.
C At the end of the run, statistical quantities of interest are
C printed (see optional outputs in the full description below).
C
C     EXTERNAL FEX, GEX
C     REAL ATOL, RTOL, RWORK, T, TOUT, Y
C     DIMENSION Y(3), ATOL(3), RWORK(76), IWORK(23), JROOT(2)
C     NEQ = 3
C     Y(1) = 1.
C     Y(2) = 0.
C     Y(3) = 0.
C     T = 0.
C     TOUT = .4
C     ITOL = 2
C     RTOL = 1.E-4
C     ATOL(1) = 1.E-6
C     ATOL(2) = 1.E-10
C     ATOL(3) = 1.E-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 76
C     LIW = 23
C     JT = 2
C     NG = 2
C     DO 40 IOUT = 1,12
C 10    CALL SLSODAR(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
C    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT,GEX,NG,JROOT)
C       WRITE(6,20)T,Y(1),Y(2),Y(3)
C 20    FORMAT(' At t =',E12.4,'   Y =',3E14.6)
C       IF (ISTATE .LT. 0) GO TO 80
C       IF (ISTATE .EQ. 2) GO TO 40
C       WRITE(6,30)JROOT(1),JROOT(2)
C 30    FORMAT(5X,' The above line is a root,  JROOT =',2I5)
C       ISTATE = 2
C       GO TO 10
C 40    TOUT = TOUT*10.
C     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(10),
C    1   IWORK(19),RWORK(15)
C 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4,
C    1   '  No. g-s =',I4/
C    2   ' Method last used =',I2,'   Last switch was at t =',E12.4)
C     STOP
C 80  WRITE(6,90)ISTATE
C 90  FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE FEX (NEQ, T, Y, YDOT)
C     REAL T, Y, YDOT
C     DIMENSION Y(3), YDOT(3)
C     YDOT(1) = -.04*Y(1) + 1.E4*Y(2)*Y(3)
C     YDOT(3) = 3.E7*Y(2)*Y(2)
C     YDOT(2) = -YDOT(1) - YDOT(3)
C     RETURN
C     END
C
C     SUBROUTINE GEX (NEQ, T, Y, NG, GOUT)
C     REAL T, Y, GOUT
C     DIMENSION Y(3), GOUT(2)
C     GOUT(1) = Y(1) - 1.E-4
C     GOUT(2) = Y(3) - 1.E-2
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C   At t =  2.6400e-01   y =  9.899653e-01  3.470563e-05  1.000000e-02
C        The above line is a root,  JROOT =    0    1
C   At t =  4.0000e-01   Y =  9.851712e-01  3.386380e-05  1.479493e-02
C   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02
C   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01
C   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01
C   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01
C   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01
C   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01
C   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01
C   At t =  2.0745e+07   Y =  1.000000e-04  4.000395e-10  9.999000e-01
C        The above line is a root,  JROOT =    1    0
C   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01
C   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01
C   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01
C   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00
C
C   No. steps = 361  No. f-s = 693  No. J-s =  64  No. g-s = 390
C   Method last used = 2   Last switch was at t =  6.0092e-03
C
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSODAR.
C
C The user interface to SLSODAR consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODAR, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODAR package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of a subroutine in the SLSODAR package,
C      which the user may replace with his/her own version, if desired.
C      this relates to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC,
C     JT, G, and NG,
C that used only for output is  JROOT,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODAR to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C F      = the name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared External in the calling program.
C
C          Subroutine F may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in F) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y below.
C
C          If quantities computed in the F routine are needed
C          externally to SLSODAR, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use SINTDY instead.
C
C NEQ    = the size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in F and/or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODAR package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to F, JAC, and G.  Hence, if it is an array, locations
C          NEQ(2),... may be used to store other integer data and pass
C          it to F, JAC, and G.  Each such subroutine must include
C          NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  On output, Y contains the computed solution vector,
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to F,
C          JAC, and G.  Hence its length may exceed NEQ, and locations
C          Y(NEQ+1),... may be used to store other real data and
C          pass it to F, JAC, and G.  (The SLSODAR package accesses only
C          Y(1),...,Y(NEQ).)
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution y is evaluated (usually the same as TOUT).
C          If a root was found, T is the computed location of the
C          root reached first, on output.
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      max-norm of ( E(i)/EWT(i) )   .le.   1,
C          where EWT = (EWT(i)) is a vector of positive error weights.
C          The values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting a
C          user-supplied routine for the setting of EWT.
C          See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
C             and any optional inputs except H0, MXORDN, and MXORDS.
C             (See IWORK description for ML and MU.)
C             In addition, immediately following a return with
C             ISTATE = 3 (root found), NG and G may be changed.
C             (But changing NG from 0 to .gt. 0 is not allowed.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           1  means nothing was done; TOUT = t and ISTATE = 1 on input.
C           2  means the integration was performed successfully, and
C              no roots were found.
C           3  means the integration was successful, and one or more
C              roots were found before satisfying the stop condition
C              specified by ITASK.  See JROOT.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix,
C              if one is being used.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means the length of RWORK and/or IWORK was too small to
C              proceed, but the integration was successful as far as T.
C              This happens when SLSODAR chooses to switch methods
C              but LRW and/or LIW is too small for the new method.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real array (single precision) for work space, and (in the
C          first 20 words) for conditional and optional inputs and
C          optional outputs.
C          As SLSODAR switches automatically between stiff and nonstiff
C          methods, the required length of RWORK can change during the
C          problem.  Thus the RWORK array passed to SLSODAR can either
C          have a static (fixed) length large enough for both methods,
C          or have a dynamic (changing) length altered by the calling
C          program in response to output from SLSODAR.
C
C                       --- Fixed Length Case ---
C          If the RWORK length is to be fixed, it should be at least
C               max (LRN, LRS),
C          where LRN and LRS are the RWORK lengths required when the
C          current method is nonstiff or stiff, respectively.
C
C          The separate RWORK length requirements LRN and LRS are
C          as follows:
C          If NEQ is constant and the maximum method orders have
C          their default values, then
C             LRN = 20 + 16*NEQ + 3*NG,
C             LRS = 22 + 9*NEQ + NEQ**2 + 3*NG           (JT = 1 or 2),
C             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ + 3*NG   (JT = 4 or 5).
C          Under any other conditions, LRN and LRS are given by:
C             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ + 3*NG,
C             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT + 3*NG,
C          where
C             NYH    = the initial value of NEQ,
C             MXORDN = 12, unless a smaller value is given as an
C                      optional input,
C             MXORDS = 5, unless a smaller value is given as an
C                      optional input,
C             LMAT   = length of matrix work space:
C             LMAT   = NEQ**2 + 2              if JT = 1 or 2,
C             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.
C
C                       --- Dynamic Length Case ---
C          If the length of RWORK is to be dynamic, then it should
C          be at least LRN or LRS, as defined above, depending on the
C          current method.  Initially, it must be at least LRN (since
C          SLSODAR starts with the nonstiff method).  On any return
C          from SLSODAR, the optional output MCUR indicates the current
C          method.  If MCUR differs from the value it had on the
C          previous return, or if there has only been one call to
C          SLSODAR and MCUR is now 2, then SLSODAR has switched
C          methods during the last call, and the length of RWORK
C          should be reset (to LRN if MCUR = 1, or to LRS if
C          MCUR = 2).  (An increase in the RWORK length is required
C          if SLSODAR returned ISTATE = -7, but not otherwise.)
C          After resetting the length, call SLSODAR with ISTATE = 3
C          to signal that change.
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer array for work space.
C          As SLSODAR switches automatically between stiff and nonstiff
C          methods, the required length of IWORK can change during
C          problem, between
C             LIS = 20 + NEQ   and   LIN = 20,
C          respectively.  Thus the IWORK array passed to SLSODAR can
C          either have a fixed length of at least 20 + NEQ, or have a
C          dynamic length of at least LIN or LIS, depending on the
C          current method.  The comments on dynamic length under
C          RWORK above apply here.  Initially, this length need
C          only be at least LIN = 20.
C
C          The first few words of IWORK are used for conditional and
C          optional inputs and optional outputs.
C
C          The following 2 words in IWORK are conditional inputs:
C            IWORK(1) = ML     These are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       banded Jacobian, excluding the main diagonal.
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if JT is 4 or 5, and
C                       ignored otherwise.  ML and MU may in fact be
C                       the band parameters for a matrix to which
C                       df/dy is only approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note: The base addresses of the work arrays must not be
C altered between calls to SLSODAR for the same problem.
C The contents of the work arrays must not be altered
C between calls, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODAR between calls, if
C desired (but not for use by F, JAC, or G).
C
C JAC    = the name of the user-supplied routine to compute the
C          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine
C          is optional, but if the problem is expected to be stiff much
C          of the time, you are encouraged to supply JAC, for the sake
C          of efficiency.  (Alternatively, set JT = 2 or 5 to have
C          SLSODAR compute df/dy internally by difference quotients.)
C          If and when SLSODAR uses df/dy, it treats this NEQ by NEQ
C          matrix either as full (JT = 1 or 2), or as banded (JT =
C          4 or 5) with half-bandwidths ML and MU (discussed under
C          IWORK above).  In either case, if JT = 1 or 4, the JAC
C          routine must compute df/dy as a function of the scalar t
C          and the vector y.  It is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C               DIMENSION Y(*), PD(NROWPD,*)
C          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
C          PD is to be loaded with partial derivatives (elements of
C          the Jacobian matrix) on output.  PD must be given a first
C          dimension of NROWPD.  T and Y have the same meaning as in
C          Subroutine F.
C               In the full matrix case (JT = 1), ML and MU are
C          ignored, and the Jacobian is to be loaded into PD in
C          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
C               In the band matrix case (JT = 4), the elements
C          within the band are to be loaded into PD in columnwise
C          manner, with diagonal lines of df/dy loaded into the rows
C          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
C          ML and MU are the half-bandwidth parameters (see IWORK).
C          The locations in PD in the two triangular areas which
C          correspond to nonexistent matrix elements can be ignored
C          or loaded arbitrarily, as they are overwritten by SLSODAR.
C               JAC need not provide df/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               In either case, PD is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to F with the same
C          arguments NEQ, T, and Y.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by F and not recomputed by JAC,
C          if desired.  Also, JAC may alter the Y array, if desired.
C          JAC must be declared External in the calling program.
C               Subroutine JAC may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C JT     = Jacobian type indicator.  Used only for input.
C          JT specifies how the Jacobian matrix df/dy will be
C          treated, if and when SLSODAR requires this matrix.
C          JT has the following values and meanings:
C           1 means a user-supplied full (NEQ by NEQ) Jacobian.
C           2 means an internally generated (difference quotient) full
C             Jacobian (using NEQ extra calls to F per df/dy value).
C           4 means a user-supplied banded Jacobian.
C           5 means an internally generated banded Jacobian (using
C             ML+MU+1 extra calls to F per df/dy evaluation).
C          If JT = 1 or 4, the user must supply a Subroutine JAC
C          (the name is arbitrary) as described above under JAC.
C          If JT = 2 or 5, a dummy argument can be used.
C
C G      = the name of subroutine for constraint functions, whose
C          roots are desired during the integration.  It is to have
C          the form
C               SUBROUTINE G (NEQ, T, Y, NG, GOUT)
C               DIMENSION Y(*), GOUT(NG)
C          where NEQ, T, Y, and NG are input, and the array GOUT
C          is output.  NEQ, T, and Y have the same meaning as in
C          the F routine, and GOUT is an array of length NG.
C          For i = 1,...,NG, this routine is to load into GOUT(i)
C          the value at (T,Y) of the i-th constraint function g(i).
C          SLSODAR will find roots of the g(i) of odd multiplicity
C          (i.e. sign changes) as they occur during the integration.
C          G must be declared External in the calling program.
C
C          Caution:  Because of numerical errors in the functions
C          g(i) due to roundoff and integration error, SLSODAR may
C          return false roots, or return the same root at two or more
C          nearly equal values of t.  If such false roots are
C          suspected, the user should consider smaller error tolerances
C          and/or higher precision in the evaluation of the g(i).
C
C          If a root of some g(i) defines the end of the problem,
C          the input to SLSODAR should nevertheless allow integration
C          to a point slightly past that root, so that SLSODAR can
C          locate the root by interpolation.
C
C          Subroutine G may access user-defined quantities in
C          NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in G) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C NG     = number of constraint functions g(i).  If there are none,
C          set NG = 0, and pass a dummy name for G.
C
C JROOT  = integer array of length NG.  Used only for output.
C          On a return with ISTATE = 3 (one or more roots found),
C          JROOT(i) = 1 if g(i) has a root at T, or JROOT(i) = 0 if not.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C IXPR    IWORK(5)  flag to generate extra printing at method switches.
C                   IXPR = 0 means no extra printing (the default).
C                   IXPR = 1 means print data on each switch.
C                   T, H, and NST will be printed on the same logical
C                   unit as used for error messages.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
C                   (Adams) method.  The default value is 12.
C                   If MXORDN exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDN is held constant during the problem.
C
C MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
C                   (BDF) method.  The default value is 5.
C                   If MXORDS exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDS is held constant during the problem.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODAR, the variables listed
C below are quantities related to the performance of SLSODAR
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODAR, and on any return with
C ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C TSW     RWORK(15) the value of t at the time of the last method
C                   switch, if any.
C
C NGE     IWORK(10) the number of g evaluations for the problem so far.
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NFE     IWORK(12) the number of f evaluations for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (and of matrix
C                   LU decompositions) for the problem so far.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required, assuming
C                   that the length of RWORK is to be fixed for the
C                   rest of the problem, and that switching may occur.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required, assuming
C                   that the length of IWORK is to be fixed for the
C                   rest of the problem, and that switching may occur.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C MUSED   IWORK(19) the method indicator for the last successful step:
C                   1 means Adams (nonstiff), 2 means BDF (stiff).
C
C MCUR    IWORK(20) the current method indicator:
C                   1 means Adams (nonstiff), 2 means BDF (stiff).
C                   This is the method to be attempted
C                   on the next step.  Thus it differs from MUSED
C                   only if a method switch has just been made.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21 + 3*NG      the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LACOR         array of size NEQ used for the accumulated
C         (from Common   corrections on each step, scaled on output
C           as noted)    to represent the estimated local error in y
C                        on the last step.  This is the vector E in
C                        the description of the error control.  It is
C                        defined only on a successful return from
C                        SLSODAR.  The base address LACOR is obtained by
C                        including in the user's program the
C                        following 2 lines:
C                           COMMON /SLS001/ RLS(218), ILS(37)
C                           LACOR = ILS(22)
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODAR.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSODAR, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSODAR.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCAR(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODAR (see Part 3 below).
C                             RSAV must be a real array of length 245
C                             or more, and ISAV must be an integer
C                             array of length 55 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCAR is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODAR.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODAR.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   LYH = 21 + 3*NG
C   CALL SINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODAR).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(t), is already provided
C             by SLSODAR directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C LYH       = 21 + 3*NG = base address in RWORK of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODAR is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODAR, and
C   (2) the three internal Common blocks
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C         /SLSA01/  of length  31    (22 single precision words
C                      followed by  9 integer words).
C         /SLSR01/  of length   7  (3 single precision words
C                      followed by  4 integer words).
C
C If SLSODAR is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODAR is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODAR call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODAR call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCAR (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below is a description of a routine in the SLSODAR package which
C relates to the measurement of errors, and can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     Subroutine SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODAR call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the
C SMNORM routine, and also used by SLSODAR in the computation
C of the optional output IMXER, and the increments for difference
C quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19811102  DATE WRITTEN
C 19820126  Fixed bug in tests of work space lengths;
C           minor corrections in main prologue and comments.
C 19820507  Fixed bug in RCHEK in setting HMING.
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODA;
C           in STODA, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODAR.
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           numerous corrections to prologues and internal comments.
C 20010613  Revised excess accuracy test (to match rest of ODEPACK).
C 20010808  Fixed bug in SPRJA (matrix in SBNORM call).
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODAR package.
C
C In addition to Subroutine SLSODAR, the SLSODAR package includes the
C following subroutines and function routines:
C  SRCHEK   does preliminary checking for roots, and serves as an
C           interface between Subroutine SLSODAR and Subroutine SROOTS.
C  SROOTS   finds the leftmost root of a set of functions.
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODA   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPRJA    computes and preprocesses the Jacobian matrix J = df/dy
C           and the Newton iteration matrix P = I - h*l0*J.
C  SSOLSY   manages solution of linear system in chord iteration.
C  SEWSET   sets the error weight vector EWT before each step.
C  SMNORM   computes the weighted max-norm of a vector.
C  SFNORM   computes the norm of a full matrix consistent with the
C           weighted max-norm on vectors.
C  SBNORM   computes the norm of a band matrix consistent with the
C           weighted max-norm on vectors.
C  SSRCAR   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SGEFA and SGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  SGBFA and SGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  SCOPY    is one of the basic linear algebra modules (BLAS).
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SMNORM, SFNORM, SBNORM, RUMACH, IXSAV, and IUMACH are
C function routines.  All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL SPRJA, SSOLSY
      REAL RUMACH, SMNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER INSUFR, INSUFI, IXPR, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      INTEGER LG0, LG1, LGX, IOWNR3, IRFND, ITASKC, NGC, NGE
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LENIW,
     1   LENRW, LENWM, LF0, ML, MORD, MU, MXHNL0, MXSTP0
      INTEGER LEN1, LEN1C, LEN1N, LEN1S, LEN2, LENIWC, LENRWC
      INTEGER IRFP, IRT, LENYH, LYHNEW
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL TSW, ROWNS2, PDNORM
      REAL ROWNR3, T0, TLAST, TOUTC
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following three internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODAR, SINTDY, SSTODA,
C SPRJA, and SSOLSY.
C The block SLSA01 is declared in subroutines SLSODAR, SSTODA, SPRJA.
C The block SLSR01 is declared in subroutines SLSODAR, SRCHEK, SROOTS.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /SLSA01/ TSW, ROWNS2(20), PDNORM,
     1   INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
C
      COMMON /SLSR01/ ROWNR3(2), T0, TLAST, TOUTC,
     1   LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      ITASKC = ITASK
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C JT, ML, MU, and NG.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      IF (JT .EQ. 3 .OR. JT .LT. 1 .OR. JT .GT. 5) GO TO 608
      JTYP = JT
      IF (JT .LE. 2) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
      IF (NG .LT. 0) GO TO 630
      IF (ISTATE .EQ. 1) GO TO 35
      IF (IRFND .EQ. 0 .AND. NG .NE. NGC) GO TO 631
 35   NGC = NG
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      IXPR = 0
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      HMXI = 0.0E0
      HMIN = 0.0E0
      IF (ISTATE .NE. 1) GO TO 60
      H0 = 0.0E0
      MXORDN = MORD(1)
      MXORDS = MORD(2)
      GO TO 60
 40   IXPR = IWORK(5)
      IF (IXPR .LT. 0 .OR. IXPR .GT. 1) GO TO 611
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      MXORDN = IWORK(8)
      IF (MXORDN .LT. 0) GO TO 628
      IF (MXORDN .EQ. 0) MXORDN = 100
      MXORDN = MIN(MXORDN,MORD(1))
      MXORDS = IWORK(9)
      IF (MXORDS .LT. 0) GO TO 629
      IF (MXORDS .EQ. 0) MXORDS = 100
      MXORDS = MIN(MXORDS,MORD(2))
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C If ISTATE = 1, METH is initialized to 1 here to facilitate the
C checking of work space lengths.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  G0, G1, GX, YH, WM,
C EWT, SAVF, ACOR.
C If the lengths provided are insufficient for the current method,
C an error return occurs.  This is treated as illegal input on the
C first call, but as a problem interruption with ISTATE = -7 on a
C continuation call.  If the lengths are sufficient for the current
C method but not for both methods, a warning message is sent.
C-----------------------------------------------------------------------
 60   IF (ISTATE .EQ. 1) METH = 1
      IF (ISTATE .EQ. 1) NYH = N
      LG0 = 21
      LG1 = LG0 + NG
      LGX = LG1 + NG
      LYHNEW = LGX + NG
      IF (ISTATE .EQ. 1) LYH = LYHNEW
      IF (LYHNEW .EQ. LYH) GO TO 62
C If ISTATE = 3 and NG was changed, shift YH to its new location. ------
      LENYH = L*NYH
      IF (LRW .LT. LYHNEW-1+LENYH) GO TO 62
      I1 = 1
      IF (LYHNEW .GT. LYH) I1 = -1
      CALL SCOPY (LENYH, RWORK(LYH), I1, RWORK(LYHNEW), I1)
      LYH = LYHNEW
 62   CONTINUE
      LEN1N = LYHNEW - 1 + (MXORDN + 1)*NYH
      LEN1S = LYHNEW - 1 + (MXORDS + 1)*NYH
      LWM = LEN1S + 1
      IF (JT .LE. 2) LENWM = N*N + 2
      IF (JT .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEN1S = LEN1S + LENWM
      LEN1C = LEN1N
      IF (METH .EQ. 2) LEN1C = LEN1S
      LEN1 = MAX(LEN1N,LEN1S)
      LEN2 = 3*N
      LENRW = LEN1 + LEN2
      LENRWC = LEN1C + LEN2
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      LENIWC = 20
      IF (METH .EQ. 2) LENIWC = LENIW
      IWORK(18) = LENIW
      IF (ISTATE .EQ. 1 .AND. LRW .LT. LENRWC) GO TO 617
      IF (ISTATE .EQ. 1 .AND. LIW .LT. LENIWC) GO TO 618
      IF (ISTATE .EQ. 3 .AND. LRW .LT. LENRWC) GO TO 550
      IF (ISTATE .EQ. 3 .AND. LIW .LT. LENIWC) GO TO 555
      LEWT = LEN1 + 1
      INSUFR = 0
      IF (LRW .GE. LENRW) GO TO 65
      INSUFR = 2
      LEWT = LEN1C + 1
      MSG='SLSODAR-  Warning.. RWORK length is sufficient for now, but '
      CALL XERRWV (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWV (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      Length needed is LENRW = I1, while LRW = I2.'
      CALL XERRWV (MSG, 50, 103, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
 65   LSAVF = LEWT + N
      LACOR = LSAVF + N
      INSUFI = 0
      IF (LIW .GE. LENIW) GO TO 70
      INSUFI = 2
      MSG='SLSODAR-  Warning.. IWORK length is sufficient for now, but '
      CALL XERRWV (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWV (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      Length needed is LENIW = I1, while LIW = I2.'
      CALL XERRWV (MSG, 50, 104, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
 70   CONTINUE
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 75 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 75     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C if ISTATE = 3, set flag to signal parameter changes to SSTODA. -------
      JSTART = -1
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  zero part of yh to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      TSW = T
      MAXORD = MXORDN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      MUSED = 0
      MITER = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by:
C
C   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
C
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         F      = the initial value of the vector f(t,y), and
C         norm() = the weighted vector norm used throughout, given by
C                  the SMNORM function routine, and weighted by the
C                  tolerances initially loaded into the EWT array.
C The sign of H0 is inferred from the initial values of TOUT and T.
C ABS(H0) is made .le. ABS(TOUT-T) in any case.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SMNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
C
C Check for a zero of g at T. ------------------------------------------
      IRFND = 0
      TOUTC = TOUT
      IF (NGC .EQ. 0) GO TO 270
      CALL SRCHEK (1, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT)
      IF (IRT .EQ. 0) GO TO 270
      GO TO 632
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C First, SRCHEK is called to check for a root within the last step
C taken, other than the last root found there, if any.
C If ITASK = 2 or 5, and y(TN) has not yet been returned to the user
C because of an intervening root, return through Block G.
C-----------------------------------------------------------------------
 200  NSLAST = NST
C
      IRFP = IRFND
      IF (NGC .EQ. 0) GO TO 205
      IF (ITASK .EQ. 1 .OR. ITASK .EQ. 4) TOUTC = TOUT
      CALL SRCHEK (2, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT)
      IF (IRT .NE. 1) GO TO 205
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 205  CONTINUE
      IRFND = 0
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 2) GO TO 400
C
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      T = TN
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) T = TCRIT
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 5) GO TO 400
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2 .AND. JSTART .GE. 0) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODA.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF (METH .EQ. MUSED) GO TO 255
      IF (INSUFR .EQ. 1) GO TO 550
      IF (INSUFI .EQ. 1) GO TO 555
 255  IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SMNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODAR-  Warning..Internal T(=R1) and H(=R2) are '
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODAR-  Above warning has been issued I1 times. '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C   CALL SSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,SPRJA,SSOLSY)
C-----------------------------------------------------------------------
      CALL SSTODA (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, SPRJA, SSOLSY)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).
C If a method switch was just made, record TSW, reset MAXORD,
C set JSTART to -1 to signal SSTODA to complete the switch,
C and do extra printing of data if IXPR = 1.
C Then call SRCHEK to check for a root within the last step.
C Then, if no root was found, check for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      IF (METH .EQ. MUSED) GO TO 310
      TSW = TN
      MAXORD = MXORDN
      IF (METH .EQ. 2) MAXORD = MXORDS
      IF (METH .EQ. 2) RWORK(LWM) = SQRT(UROUND)
      INSUFR = MIN(INSUFR,1)
      INSUFI = MIN(INSUFI,1)
      JSTART = -1
      IF (IXPR .EQ. 0) GO TO 310
      IF (METH .EQ. 2) THEN
      MSG='SLSODAR- A switch to the BDF (stiff) method has occurred    '
      CALL XERRWV (MSG, 60, 105, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (METH .EQ. 1) THEN
      MSG='SLSODAR- A switch to the Adams (nonstiff) method occurred   '
      CALL XERRWV (MSG, 60, 106, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      MSG='     at T = R1,  tentative step size H = R2,  step NST = I1 '
      CALL XERRWV (MSG, 60, 107, 0, 1, NST, 0, 2, TN, H)
 310  CONTINUE
C
      IF (NGC .EQ. 0) GO TO 315
      CALL SRCHEK (3, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT)
      IF (IRT .NE. 1) GO TO 315
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 315  CONTINUE
C
      GO TO (320, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 320  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (JSTART .GE. 0) JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODAR.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
 425  CONTINUE
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODAR-  At current T (=R1), MXSTEP (=I1) steps  '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODAR-  At T(=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODAR-  At T (=R1), too much accuracy requested '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODAR-  At T(=R1), step size H(=R2), the error  '
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODAR-  At T (=R1) and step size H (=R2), the   '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C RWORK length too small to proceed. -----------------------------------
 550  MSG = 'SLSODAR- At current T(=R1), RWORK length too small'
      CALL XERRWV (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWV (MSG, 60, 206, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 580
C IWORK length too small to proceed. -----------------------------------
 555  MSG = 'SLSODAR- At current T(=R1), IWORK length too small'
      CALL XERRWV (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWV (MSG, 60, 207, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0E0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODAR-  ISTATE(=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODAR-  ITASK (=I1) illegal.'
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODAR-  ISTATE.gt.1 but SLSODAR not initialized.'
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODAR-  NEQ (=I1) .lt. 1    '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODAR-  ISTATE = 3 and NEQ increased (I1 to I2).'
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODAR-  ITOL (=I1) illegal. '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODAR-  IOPT (=I1) illegal. '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODAR-  JT (=I1) illegal.   '
      CALL XERRWV (MSG, 30, 8, 0, 1, JT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 609  MSG = 'SLSODAR-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 610  MSG = 'SLSODAR-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODAR-  IXPR (=I1) illegal. '
      CALL XERRWV (MSG, 30, 11, 0, 1, IXPR, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODAR-  MXSTEP (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODAR-  MXHNIL (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODAR-  TOUT (=R1) behind T (=R2)     '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODAR-  HMAX (=R1) .lt. 0.0 '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODAR-  HMIN (=R1) .lt. 0.0 '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG='SLSODAR-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) '
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG='SLSODAR-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) '
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODAR-  RTOL(I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODAR-  ATOL(I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODAR-  EWT(I1) is R1 .le. 0.0        '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODAR- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODAR-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2) '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODAR-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODAR-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODAR-  At start of problem, too much accuracy  '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODAR-  Trouble in SINTDY. ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
      GO TO 700
 628  MSG = 'SLSODAR-  MXORDN (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 28, 0, 1, MXORDN, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 629  MSG = 'SLSODAR-  MXORDS (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 29, 0, 1, MXORDS, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 630  MSG = 'SLSODAR-  NG (=I1) .lt. 0     '
      CALL XERRWV (MSG, 30, 30, 0, 1, NG, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 631  MSG = 'SLSODAR-  NG changed (from I1 to I2) illegally,   '
      CALL XERRWV (MSG, 50, 31, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      i.e. not immediately after a root was found.'
      CALL XERRWV (MSG, 50, 31, 0, 2, NGC, NG, 0, 0.0E0, 0.0E0)
      GO TO 700
 632  MSG = 'SLSODAR-  One or more components of g has a root  '
      CALL XERRWV (MSG, 50, 32, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      too near to the initial point.    '
      CALL XERRWV (MSG, 40, 32, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODAR-  Run aborted.. apparent infinite loop.   '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODAR ---------------------
      END
*DECK SLSODPK
      SUBROUTINE SLSODPK (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL, MF)
      EXTERNAL F, JAC, PSOL
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 18 November 2003 version of
C SLSODPK: Livermore Solver for Ordinary Differential equations,
C          with Preconditioned Krylov iteration methods for the
C          Newton correction linear systems.
C
C This version is in single precision.
C
C SLSODPK solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C-----------------------------------------------------------------------
C Introduction.
C
C This is a modification of the SLSODE package which incorporates
C various preconditioned Krylov subspace iteration methods for the
C linear algebraic systems that arise in the case of stiff systems.
C
C The linear systems that must be solved have the form
C   A * x  = b ,  where  A = identity - hl0 * (df/dy) .
C Here hl0 is a scalar, and df/dy is the Jacobian matrix of partial
C derivatives of f (NEQ by NEQ).
C
C The particular Krylov method is chosen by setting the second digit,
C MITER, in the method flag MF.
C Currently, the values of MITER have the following meanings:
C
C  MITER = 1 means the preconditioned Scaled Incomplete
C            Orthogonalization Method (SPIOM).
C
C          2 means an incomplete version of the Preconditioned Scaled
C            Generalized Minimal Residual method (SPIGMR).
C            This is the best choice in general.
C
C          3 means the Preconditioned Conjugate Gradient method (PCG).
C            Recommended only when df/dy is symmetric or nearly so.
C
C          4 means the scaled Preconditioned Conjugate Gradient method
C            (PCGS).  Recommended only when D-inverse * df/dy * D is
C            symmetric or nearly so, where D is the diagonal scaling
C            matrix with elements 1/EWT(i) (see RTOL/ATOL description).
C
C          9 means that only a user-supplied matrix P (approximating A)
C            will be used, with no Krylov iteration done.  This option
C            allows the user to provide the complete linear system
C            solution algorithm, if desired.
C
C The user can apply preconditioning to the linear system A*x = b,
C by means of arbitrary matrices (the preconditioners).
C     In the case of SPIOM and SPIGMR, one can apply left and right
C preconditioners P1 and P2, and the basic iterative method is then
C applied to the matrix (P1-inverse)*A*(P2-inverse) instead of to the
C matrix A.  The product P1*P2 should be an approximation to matrix A
C such that linear systems with P1 or P2 are easier to solve than with
C A.  Preconditioning from the left only or right only means using
C P2 = identity or P1 = identity, respectively.
C     In the case of the PCG and PCGS methods, there is only one
C preconditioner matrix P (but it can be the product of more than one).
C It should approximate the matrix A but allow for relatively
C easy solution of linear systems with coefficient matrix P.
C For PCG, P should be positive definite symmetric, or nearly so,
C and for PCGS, the scaled preconditioner D-inverse * P * D
C should be symmetric or nearly so.
C     If the Jacobian J = df/dy splits in a natural way into a sum
C J = J1 + J2, then one possible choice of preconditioners is
C     P1 = identity - hl0 * J1  and  P2 = identity - hl0 * J2
C provided each of these is easy to solve (or approximately solve).
C
C-----------------------------------------------------------------------
C References:
C 1.  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage Matrix
C     Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31 (1989),
C     pp. 40-91; also  L.L.N.L. Report UCRL-95088, Rev. 1, June 1987.
C 2.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C-----------------------------------------------------------------------
C Authors:       Alan C. Hindmarsh and Peter N. Brown
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODPK package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the demonstration
C program distributed with this solver.
C
C A. First provide a subroutine of the form:
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Next determine (or guess) whether or not the problem is stiff.
C Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
C whose real part is negative and large in magnitude, compared to the
C reciprocal of the t span of interest.  If the problem is nonstiff,
C use a method flag MF = 10.  If it is stiff, MF should be between 21
C and 24, or possibly 29.  MF = 22 is generally the best choice.
C Use 23 or 24 only if symmetry is present.  Use MF = 29 if the
C complete linear system solution is to be provided by the user.
C The following four parameters must also be set.
C  IWORK(1) = LWP  = length of real array WP for preconditioning.
C  IWORK(2) = LIWP = length of integer array IWP for preconditioning.
C  IWORK(3) = JPRE = preconditioner type flag:
C                  = 0 for no preconditioning (P1 = P2 = P = identity)
C                  = 1 for left-only preconditioning (P2 = identity)
C                  = 2 for right-only preconditioning (P1 = identity)
C                  = 3 for two-sided preconditioning (and PCG or PCGS)
C  IWORK(4) = JACFLG = flag for whether JAC is called.
C                    = 0 if JAC is not to be called,
C                    = 1 if JAC is to be called.
C  Use JACFLG = 1 if JAC computes any nonconstant data for use in
C  preconditioning, such as Jacobian elements.
C  The arrays WP and IWP are work arrays under the user's control,
C  for use in the routines that perform preconditioning operations.
C
C C. If the problem is stiff, you must supply two routines that deal
C with the preconditioning of the linear systems to be solved.
C These are as follows:
C
C     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HL0, WP,IWP, IER)
C     DIMENSION Y(*), YSV(*), REWT(*), FTY(*), V(*), WP(*), IWP(*)
C        This routine must evaluate and preprocess any parts of the
C     Jacobian matrix df/dy involved in the preconditioners P1, P2, P.
C     The Y and FTY arrays contain the current values of y and f(t,y),
C     respectively, and YSV also contains the current value of y.
C     The array V is work space of length NEQ.
C     JAC must multiply all computed Jacobian elements by the scalar
C     -HL0, add the identity matrix, and do any factorization
C     operations called for, in preparation for solving linear systems
C     with a coefficient matrix of P1, P2, or P.  The matrix P1*P2 or P
C     should be an approximation to  identity - HL0 * (df/dy).
C     JAC should return IER = 0 if successful, and IER .ne. 0 if not.
C     (If IER .ne. 0, a smaller time step will be tried.)
C
C     SUBROUTINE PSOL (NEQ, T, Y, SAVF, WK, HL0, WP, IWP, B, LR, IER)
C     DIMENSION Y(*), SAVF(*), WK(*), WP(*), IWP(*), B(*)
C        This routine must solve a linear system with B as right-hand
C     side and one of the preconditioning matrices, P1, P2, or P, as
C     coefficient matrix, and return the solution vector in B.
C     LR is a flag concerning left vs right preconditioning, input
C     to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
C     In the case of the PCG or PCGS method, LR will be 3, and PSOL
C     should solve the system P*x = B with the preconditioner matrix P.
C     In the case MF = 29 (no Krylov iteration), LR will be 0,
C     and PSOL is to return in B the desired approximate solution
C     to A * x = B, where A = identity - HL0 * (df/dy).
C     PSOL can use data generated in the JAC routine and stored in
C     WP and IWP.  WK is a work array of length NEQ.
C     The argument HL0 is the current value of the scalar appearing
C     in the linear system.  If the old value, at the time of the last
C     JAC call, is needed, it must have been saved by JAC in WP.
C     On return, PSOL should set the error flag IER as follows:
C       IER = 0 if PSOL was successful,
C       IER .gt. 0 if a recoverable error occurred, meaning that the
C              time step will be retried,
C       IER .lt. 0 if an unrecoverable error occurred, meaning that the
C              solver is to stop immediately.
C
C D. Write a main program which calls Subroutine SLSODPK once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages by
C SLSODPK.  On the first call to SLSODPK, supply arguments as follows:
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared External in calling program.
C NEQ    = number of first order ODEs.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             20 + 16*NEQ           for MF = 10,
C             45 + 17*NEQ + LWP     for MF = 21,
C             61 + 17*NEQ + LWP     for MF = 22,
C             20 + 15*NEQ + LWP     for MF = 23 or 24,
C             20 + 12*NEQ + LWP     for MF = 29.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least:
C             30            for MF = 10,
C             35 + LIWP     for MF = 21,
C             30 + LIWP     for MF = 22, 23, 24, or 29.
C LIW    = declared length of IWORK (in user's dimension).
C JAC,PSOL = names of subroutines for preconditioning.
C          These names must be declared External in the calling program.
C MF     = method flag.  Standard values are:
C          10 for nonstiff (Adams) method.
C          21 for stiff (BDF) method, with preconditioned SIOM.
C          22 for stiff method, with preconditioned GMRES method.
C          23 for stiff method, with preconditioned CG method.
C          24 for stiff method, with scaled preconditioned CG method.
C          29 for stiff method, with user's PSOL routine only.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C and possibly ATOL.
C
C E. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if SLSODPK was successful, negative otherwise.
C          -1 means excess work done on this call (perhaps wrong MF).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad JAC
C             or PSOL routine supplied or wrong choice of MF or
C             tolerances, or this solver is inappropriate).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means an unrecoverable error occurred in PSOL.
C
C F. To continue the integration after a successful return, simply
C reset TOUT and call SLSODPK again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSODPK.
C
C The user interface to SLSODPK consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODPK, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODPK package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the SLSODPK package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C  F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, PSOL, MF,
C and those used for both input and output are
C  Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODPK to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C F      = the name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared External in the calling program.
C
C          Subroutine F may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in F) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y below.
C
C          If quantities computed in the F routine are needed
C          externally to SLSODPK, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use SINTDY instead.
C
C NEQ    = the size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in the user-supplied subroutines.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODPK package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to F, JAC, and PSOL.  Hence, if it is an array, locations
C          NEQ(2),... may be used to store other integer data and pass
C          it to the user-supplied subroutines.  Each such routine must
C          include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  On output, Y contains the computed solution vector,
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to F,
C          JAC, and PSOL. Hence its length may exceed NEQ, and locations
C          Y(NEQ+1),... may be used to store other real data and
C          pass it to the user-supplied subroutines.  (The SLSODPK
C          package accesses only Y(1),...,Y(NEQ).)
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          the following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C             and any of the optional inputs except H0.
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means the PSOL routine returned an unrecoverable error
C              flag (IER .lt. 0).  The integration was successful as
C              far as T.
C
C          Note:  since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real working array (single precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LENLS + LWP    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENLS = length of work space for linear system (Krylov)
C                  method, excluding preconditioning:
C            LENLS = 0                               if MITER = 0,
C            LENLS = NEQ*(MAXL+3) + MAXL**2          if MITER = 1,
C            LENLS = NEQ*(MAXL+3+MIN(1,MAXL-KMP))
C                 + (MAXL+3)*MAXL + 1                if MITER = 2,
C            LENLS = 6*NEQ                           if MITER = 3 or 4,
C            LENLS = 3*NEQ                           if MITER = 9.
C          (See the MF description for METH and MITER, and the
C          list of optional inputs for MAXL and KMP.)
C          LWP = length of real user work space for preconditioning
C          (see JAC/PSOL).
C          Thus if default values are used and NEQ is constant,
C          this length is:
C             20 + 16*NEQ           for MF = 10,
C             45 + 24*NEQ + LWP     FOR MF = 11,
C             61 + 24*NEQ + LWP     FOR MF = 12,
C             20 + 22*NEQ + LWP     FOR MF = 13 OR 14,
C             20 + 19*NEQ + LWP     FOR MF = 19,
C             20 + 9*NEQ            FOR MF = 20,
C             45 + 17*NEQ + LWP     FOR MF = 21,
C             61 + 17*NEQ + LWP     FOR MF = 22,
C             20 + 15*NEQ + LWP     FOR MF = 23 OR 24,
C             20 + 12*NEQ + LWP     for MF = 29.
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C             30                 if MITER = 0 (MF = 10 or 20),
C             30 + MAXL + LIWP   if MITER = 1 (MF = 11, 21),
C             30 + LIWP          if MITER = 2, 3, 4, or 9.
C          MAXL = 5 unless a different optional input value is given.
C          LIWP = length of integer user work space for preconditioning
C          (see conditional input list following).
C          The first few words of IWORK are used for conditional and
C          optional inputs and optional outputs.
C
C          The following 4 words in IWORK are conditional inputs,
C          required if MITER .ge. 1:
C          IWORK(1) = LWP  = length of real array WP for use in
C                     preconditioning (part of RWORK array).
C          IWORK(2) = LIWP = length of integer array IWP for use in
C                     preconditioning (part of IWORK array).
C                     The arrays WP and IWP are work arrays under the
C                     user's control, for use in the routines that
C                     perform preconditioning operations (JAC and PSOL).
C          IWORK(3) = JPRE = preconditioner type flag:
C                   = 0 for no preconditioning (P1 = P2 = P = identity)
C                   = 1 for left-only preconditioning (P2 = identity)
C                   = 2 for right-only preconditioning (P1 = identity)
C                   = 3 for two-sided preconditioning (and PCG or PCGS)
C          IWORK(4) = JACFLG = flag for whether JAC is called.
C                   = 0 if JAC is not to be called,
C                   = 1 if JAC is to be called.
C                     Use JACFLG = 1 if JAC computes any nonconstant
C                     data needed in preconditioning operations,
C                     such as some of the Jacobian elements.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to SLSODPK
C for the same problem, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODPK between calls, if
C desired (but not for use by any of the user-supplied subroutines).
C
C JAC    = the name of the user-supplied routine to compute any
C          Jacobian elements (or approximations) involved in the
C          matrix preconditioning operations (MITER .ge. 1).
C          It is to have the form
C            SUBROUTINE JAC (F,NEQ,T,Y,YSV,REWT,FTY,V,HL0,WP,IWP,IER)
C            DIMENSION Y(*),YSV(*),REWT(*),FTY(*),V(*),WP(*),IWP(*)
C          This routine must evaluate and preprocess any parts of the
C          Jacobian matrix df/dy used in the preconditioners P1, P2, P.
C          the Y and FTY arrays contain the current values of y and
C          f(t,y), respectively, and YSV also contains the current
C          value of y.  The array V is work space of length
C          NEQ for use by JAC.  REWT is the array of reciprocal error
C          weights (1/EWT).  JAC must multiply all computed Jacobian
C          elements by the scalar -HL0, add the identity matrix, and do
C          any factorization operations called for, in preparation
C          for solving linear systems with a coefficient matrix of
C          P1, P2, or P.  The matrix P1*P2 or P should be an
C          approximation to  identity - HL0 * (df/dy).  JAC should
C          return IER = 0 if successful, and IER .ne. 0 if not.
C          (If IER .ne. 0, a smaller time step will be tried.)
C          The arrays WP (of length LWP) and IWP (of length LIWP)
C          are for use by JAC and PSOL for work space and for storage
C          of data needed for the solution of the preconditioner
C          linear systems.  Their lengths and contents are under the
C          user's control.
C          The JAC routine may save relevant Jacobian elements (or
C          approximations) used in the preconditioners, along with the
C          value of HL0, and use these to reconstruct preconditioner
C          matrices later without reevaluationg those elements.
C          This may be cost-effective if JAC is called with HL0
C          considerably different from its earlier value, indicating
C          that a corrector convergence failure has occurred because
C          of the change in HL0, not because of changes in the
C          value of the Jacobian.  In doing this, use the saved and
C          current values of HL0 to decide whether to use saved
C          or reevaluated elements.
C          JAC may alter V, but may not alter Y, YSV, REWT, FTY, or HL0.
C          JAC must be declared External in the calling program.
C               Subroutine JAC may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C PSOL   = the name of the user-supplied routine for the
C          solution of preconditioner linear systems.
C          It is to have the form
C            SUBROUTINE PSOL (NEQ,T,Y,FTY,WK,HL0,WP,IWP,B,LR,IER)
C            DIMENSION Y(*),FTY(*),WK(*),WP(*),IWP(*),B(*)
C          This routine must solve a linear system with B as right-hand
C          side and one of the preconditioning matrices, P1, P2, or P,
C          as coefficient matrix, and return the solution vector in B.
C          LR is a flag concerning left vs right preconditioning, input
C          to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
C          In the case of the PCG or PCGS method, LR will be 3, and PSOL
C          should solve the system P*x = B with the preconditioner P.
C          In the case MITER = 9 (no Krylov iteration), LR will be 0,
C          and PSOL is to return in B the desired approximate solution
C          to A * x = B, where A = identity - HL0 * (df/dy).
C          PSOL can use data generated in the JAC routine and stored in
C          WP and IWP.
C          The Y and FTY arrays contain the current values of y and
C          f(t,y), respectively.  The array WK is work space of length
C          NEQ for use by PSOL.
C          The argument HL0 is the current value of the scalar appearing
C          in the linear system.  If the old value, as of the last
C          JAC call, is needed, it must have been saved by JAC in WP.
C          On return, PSOL should set the error flag IER as follows:
C            IER = 0 if PSOL was successful,
C            IER .gt. 0 on a recoverable error, meaning that the
C                   time step will be retried,
C            IER .lt. 0 on an unrecoverable error, meaning that the
C                   solver is to stop immediately.
C          PSOL may not alter Y, FTY, or HL0.
C          PSOL must be declared External in the calling program.
C               Subroutine PSOL may access user-defined quantities in
C          NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in PSOL) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C MF     = the method flag.  Used only for input.  The legal values of
C          MF are 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, and 29.
C          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
C          METH indicates the basic linear multistep method:
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on Backward
C                     Differentiation Formulas (BDFs).
C          MITER indicates the corrector iteration method:
C            MITER = 0 means functional iteration (no linear system
C                      is involved).
C            MITER = 1 means Newton iteration with Scaled Preconditioned
C                      Incomplete Orthogonalization Method (SPIOM)
C                      for the linear systems.
C            MITER = 2 means Newton iteration with Scaled Preconditioned
C                      Generalized Minimal Residual method (SPIGMR)
C                      for the linear systems.
C            MITER = 3 means Newton iteration with Preconditioned
C                      Conjugate Gradient method (PCG)
C                      for the linear systems.
C            MITER = 4 means Newton iteration with scaled Preconditioned
C                      Conjugate Gradient method (PCGS)
C                      for the linear systems.
C            MITER = 9 means Newton iteration with only the
C                      user-supplied PSOL routine called (no Krylov
C                      iteration) for the linear systems.
C                      JPRE is ignored, and PSOL is called with LR = 0.
C          See comments in the introduction about the choice of MITER.
C          If MITER .ge. 1, the user must supply routines JAC and PSOL
C          (the names are arbitrary) as described above.
C          For MITER = 0, dummy arguments can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C DELT    RWORK(8)  convergence test constant in Krylov iteration
C                   algorithm.  The default is .05.
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MAXL    IWORK(8)  maximum number of iterations in the SPIOM, SPIGMR,
C                   PCG, or PCGS algorithm (.le. NEQ).
C                   The default is MAXL = MIN(5,NEQ).
C
C KMP     IWORK(9)  number of vectors on which orthogonalization
C                   is done in SPIOM or SPIGMR algorithm (.le. MAXL).
C                   The default is KMP = MAXL.
C                   Note:  When KMP .lt. MAXL and MF = 22, the length
C                          of RWORK must be defined accordingly.  See
C                          the definition of RWORK above.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODPK, the variables listed
C below are quantities related to the performance of SLSODPK
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODPK, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NFE     IWORK(12) the number of f evaluations for the problem so far.
C
C NPE     IWORK(13) the number of calls to JAC so far (for Jacobian
C                   evaluation associated with preconditioning).
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NNI     IWORK(19) number of nonlinear iterations so far (each of
C                   which calls an iterative linear solver).
C
C NLI     IWORK(20) number of linear iterations so far.
C                   Note: A measure of the success of algorithm is
C                   the average number of linear iterations per
C                   nonlinear iteration, given by NLI/NNI.
C                   If this is close to MAXL, MAXL may be too small.
C
C NPS     IWORK(21) number of preconditioning solve operations
C                   (PSOL calls) so far.
C
C NCFN    IWORK(22) number of convergence failures of the nonlinear
C                   (Newton) iteration so far.
C                   Note: A measure of success is the overall
C                   rate of nonlinear convergence failures, NCFN/NST.
C
C NCFL    IWORK(23) number of convergence failures of the linear
C                   iteration so far.
C                   Note: A measure of success is the overall
C                   rate of linear convergence failures, NCFL/NNI.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output
C                        to represent the estimated local error in y
C                        on the last step.  This is the vector E in
C                        the description of the error control.  It is
C                        defined only on a successful return from
C                        SLSODPK.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODPK.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSODPK, if
C                             the default is not desired.
C                             The default value of lun is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSODPK.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCPK(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODPK (see Part 3 below).
C                             RSAV must be a real array of length 222
C                             or more, and ISAV must be an integer
C                             array of length 50 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCPK is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODPK.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (See below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODPK.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   CALL SINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODPK).
C             for valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by SLSODPK directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODPK is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODPK, and
C   (2) the two internal Common blocks
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C         /SLPK01/  of length  17  (4 single precision words
C                      followed by 13 integer words).
C
C If SLSODPK is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODPK is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODPK call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODPK call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCPK (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C below are descriptions of two routines in the SLSODPK package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODPK call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the SVNORM
C routine (see below), and also used by SLSODPK in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) SVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = SVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine SEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of SVNORM suitable for use in the error control in SLSODPK.
C None of the arguments should be altered by SVNORM.
C For example, a user-supplied SVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19860901  DATE WRITTEN
C 19861010  Numerous minor revisions to SPIOM and SPGMR routines;
C           minor corrections to prologues and comments.
C 19870114  Changed name SPGMR to SPIGMR; revised residual norm
C           calculation in SPIGMR (for incomplete case);
C           revised error return logic in SPIGMR;
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODPK;
C           in STODPK, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 19871130  Added option MITER = 9; shortened WM array by 2;
C           revised early return from SPIOM and SPIGMR;
C           replaced copy loops with SCOPY/DCOPY calls;
C           minor corrections/revisions to SOLPK, SPIGMR, ATV, ATP;
C           corrections to main prologue and internal comments.
C 19880304  Corrections to type declarations in SOLPK, SPIOM, USOL.
C 19891025  Added ISTATE = -7 return; minor revisions to USOL;
C           added initialization of JACFLG in main driver;
C           removed YH and NYH from PKSET call list;
C           minor revisions to SPIOM and SPIGMR;
C           corrections to main prologue and internal comments.
C 19900803  Added YSV to JAC call list; minor comment corrections.
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           numerous corrections to prologues and internal comments.
C 20030603  Corrected duplicate type declaration for RUMACH.
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C 20031117  Changed internal name NPE to NJE.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODPK package.
C
C In addition to Subroutine SLSODPK, the SLSODPK package includes the
C following subroutines and function routines:
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted RMS-norm of a vector.
C  SSTODPK  is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPKSET   interfaces between SSTODPK and the JAC routine.
C  SSOLPK   manages solution of linear system in Newton iteration.
C  SSPIOM   performs the SPIOM algorithm.
C  SATV     computes a scaled, preconditioned product (I-hl0*J)*v.
C  SORTHOG  orthogonalizes a vector against previous basis vectors.
C  SHEFA    generates an LU factorization of a Hessenberg matrix.
C  SHESL    solves a Hessenberg square linear system.
C  SSPIGMR  performs the SPIGMR algorithm.
C  SHEQR    generates a QR factorization of a Hessenberg matrix.
C  SHELS    finds the least squares solution of a Hessenberg system.
C  SPCG     performs Preconditioned Conjugate Gradient algorithm (PCG).
C  SPCGS    performs the PCGS algorithm.
C  SATP     computes the product A*p, where A = I - hl0*df/dy.
C  SUSOL    interfaces to the user's PSOL routine (MITER = 9).
C  SSRCPK   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SAXPY, SCOPY, SDOT, SNRM2, and SSCAL   are basic linear
C           algebra modules (from the BLAS collection).
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SVNORM, SDOT, SNRM2, RUMACH, IXSAV, and IUMACH are function
C routines.  All the others are subroutines.
C
C-----------------------------------------------------------------------
      REAL RUMACH, SVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0, LENIW,
     1   LENIWK, LENRW, LENWM, LENWK, LIWP, LWP, MORD, MXHNL0, MXSTP0,
     2   NCFN0, NCFL0, NLI0, NNI0, NNID, NSTD, NWARN
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL DELT, EPCON, SQRTN, RSQRTN
      REAL ATOLI, AVDIM, AYI, BIG, EWTI, H0, HMAX, HMX,
     1   RCFL, RCFN, RH, RTOLI, TCRIT,
     2   TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT, LAVD, LCFN, LCFL, LWARN
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following two internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODPK, SINTDY, SSTODPK,
C SSOLPK, and SATV.
C The block SLPK01 is declared in subroutines SLSODPK, SSTODPK, SPKSET,
C and SSOLPK.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /SLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0) GO TO 608
      IF (MITER .GT. 4 .AND. MITER .LT. 9) GO TO 608
      IF (MITER .GE. 1) JPRE = IWORK(3)
      JACFLG = 0
      IF (MITER .GE. 1) JACFLG = IWORK(4)
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      MAXL = MIN(5,N)
      KMP = MAXL
      DELT = 0.05E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
      MAXL = IWORK(8)
      IF (MAXL .EQ. 0) MAXL = 5
      MAXL = MIN(MAXL,N)
      KMP = IWORK(9)
      IF (KMP .EQ. 0 .OR. KMP .GT. MAXL) KMP = MAXL
      DELT = RWORK(8)
      IF (DELT .EQ. 0.0E0) DELT = 0.05E0
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C RWORK segments (in order) are denoted  YH, WM, EWT, SAVF, SAVX, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWK = 0
      IF (MITER .EQ. 1) LENWK = N*(MAXL+2) + MAXL*MAXL
      IF (MITER .EQ. 2)
     1   LENWK = N*(MAXL+2+MIN(1,MAXL-KMP)) + (MAXL+3)*MAXL + 1
      IF (MITER .EQ. 3 .OR. MITER .EQ. 4) LENWK = 5*N
      IF (MITER .EQ. 9) LENWK = 2*N
      LWP = 0
      IF (MITER .GE. 1) LWP = IWORK(1)
      LENWM = LENWK + LWP
      LOCWP = LENWK + 1
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LSAVX = LSAVF + N
      LACOR = LSAVX + N
      IF (MITER .EQ. 0) LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 31
      LENIWK = 0
      IF (MITER .EQ. 1) LENIWK = MAXL
      LIWP = 0
      IF (MITER .GE. 1) LIWP = IWORK(2)
      LENIW = 30 + LENIWK + LIWP
      LOCIWP = LENIWK + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 70     CONTINUE
C Load SQRT(N) and its reciprocal in Common. ---------------------------
      SQRTN = SQRT(REAL(N))
      RSQRTN = 1.0E0/SQRTN
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to SSTODPK. ------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO 80 I = 1,N
 80     RWORK(I+LSAVF-1) = RWORK(I+LWM-1)
 90   CONTINUE
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      NLI0 = 0
      NNI0 = 0
      NCFN0 = 0
      NCFL0 = 0
      NWARN = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
      NNI = 0
      NLI = 0
      NPS = 0
      NCFN = 0
      NCFL = 0
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( f(i)/ywt(i) )**2  )
C                                       1
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         f(i)   = i-th component of initial value of f,
C         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      NLI0 = NLI
      NNI0 = NNI
      NCFN0 = NCFN
      NCFL0 = NCFL
      NWARN = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODPK.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken,
C Check for poor Newton/Krylov method performance, update EWT (if not
C at start of problem), check for too much accuracy being requested,
C and check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      NSTD = NST - NSLAST
      NNID = NNI - NNI0
      IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 255
      AVDIM = REAL(NLI - NLI0)/REAL(NNID)
      RCFN = REAL(NCFN - NCFN0)/REAL(NSTD)
      RCFL = REAL(NCFL - NCFL0)/REAL(NNID)
      LAVD = AVDIM .GT. (MAXL - 0.05E0)
      LCFN = RCFN .GT. 0.9E0
      LCFL = RCFL .GT. 0.9E0
      LWARN = LAVD .OR. LCFN .OR. LCFL
      IF (.NOT.LWARN) GO TO 255
      NWARN = NWARN + 1
      IF (NWARN .GT. 10) GO TO 255
      IF (LAVD) THEN
      MSG='SLSODPK- Warning. Poor iterative algorithm performance seen '
      CALL XERRWV (MSG, 60, 111, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (LAVD) THEN
      MSG='      at T = R1 by average no. of linear iterations = R2    '
      CALL XERRWV (MSG, 60, 111, 0, 0, 0, 0, 2, TN, AVDIM)
      ENDIF
      IF (LCFN) THEN
      MSG='SLSODPK- Warning. Poor iterative algorithm performance seen '
      CALL XERRWV (MSG, 60, 112, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (LCFN) THEN
      MSG='      at T = R1 by nonlinear convergence failure rate = R2  '
      CALL XERRWV (MSG, 60, 112, 0, 0, 0, 0, 2, TN, RCFN)
      ENDIF
      IF (LCFL) THEN
      MSG='SLSODPK- Warning. Poor iterative algorithm performance seen '
      CALL XERRWV (MSG, 60, 113, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (LCFL) THEN
      MSG='      at T = R1 by linear convergence failure rate = R2     '
      CALL XERRWV (MSG, 60, 113, 0, 0, 0, 0, 2, TN, RCFL)
      ENDIF
 255  CONTINUE
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODPK-  Warning..Internal T(=R1) and H(=R2) are '
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODPK-  Above warning has been issued I1 times. '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL SSTODPK(NEQ,Y,YH,NYH,YH,EWT,SAVF,SAVX,ACOR,WM,IWM,F,JAC,PSOL)
C-----------------------------------------------------------------------
      CALL SSTODPK (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LSAVX), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), F, JAC, PSOL)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 550), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  see if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODPK.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNI
      IWORK(20) = NLI
      IWORK(21) = NPS
      IWORK(22) = NCFN
      IWORK(23) = NCFL
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODPK-  At current T (=R1), MXSTEP (=I1) steps  '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODPK-  At T (=R1), EWT(I1) has become R2.le.0. '
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODPK-  At T (=R1), too much accuracy requested '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODPK-  At T(=R1), step size H(=R2), the error  '
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODPK-  At T (=R1) and step size H (=R2), the   '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C KFLAG = -3.  Unrecoverable error from PSOL. --------------------------
 550  MSG = 'SLSODPK-  At T (=R1) an unrecoverable error return'
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      was made from Subroutine PSOL     '
      CALL XERRWV (MSG, 40, 205, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0E0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNI
      IWORK(20) = NLI
      IWORK(21) = NPS
      IWORK(22) = NCFN
      IWORK(23) = NCFL
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODPK-  ISTATE(=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODPK-  ITASK (=I1) illegal.'
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODPK-  ISTATE.gt.1 but SLSODPK not initialized.'
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODPK-  NEQ (=I1) .lt. 1    '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODPK-  ISTATE = 3 and NEQ increased (I1 to I2).'
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODPK-  ITOL (=I1) illegal. '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODPK-  IOPT (=I1) illegal. '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODPK-  MF (=I1) illegal.   '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODPK-  MAXORD (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODPK-  MXSTEP (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODPK-  MXHNIL (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODPK-  TOUT (=R1) behind T (=R2)     '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODPK-  HMAX (=R1) .lt. 0.0 '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODPK-  HMIN (=R1) .lt. 0.0 '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG='SLSODPK-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) '
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG='SLSODPK-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) '
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODPK-  RTOL(I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODPK-  ATOL(I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODPK-  EWT(I1) is R1 .le. 0.0        '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODPK- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODPK-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2) '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODPK-  At start of problem, too much accuracy  '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODPK-  Trouble in SINTDY. ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODPK-  Run aborted.. apparent infinite loop.   '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODPK ---------------------
      END
*DECK SLSODKR
      SUBROUTINE SLSODKR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL,
     2            MF, G, NG, JROOT)
      EXTERNAL F, JAC, PSOL, G
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF,
     1        NG, JROOT
      REAL Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          JROOT(*)
C-----------------------------------------------------------------------
C This is the 18 November 2003 version of
C SLSODKR: Livermore Solver for Ordinary Differential equations,
C          with preconditioned Krylov iteration methods for the
C          Newton correction linear systems, and with Rootfinding.
C
C This version is in single precision.
C
C SLSODKR solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C At the same time, it locates the roots of any of a set of functions
C     g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).
C
C-----------------------------------------------------------------------
C Introduction.
C
C This is a modification of the SLSODE package, and differs from it
C in five ways:
C (a) It uses various preconditioned Krylov subspace iteration methods
C for the linear algebraic systems that arise in the case of stiff
C systems.  See the introductory notes below.
C (b) It does automatic switching between functional (fixpoint)
C iteration and Newton iteration in the corrector iteration.
C (c) It finds the root of at least one of a set of constraint
C functions g(i) of the independent and dependent variables.
C It finds only those roots for which some g(i), as a function
C of t, changes sign in the interval of integration.
C It then returns the solution at the root, if that occurs
C sooner than the specified stop condition, and otherwise returns
C the solution according the specified stop condition.
C (d) It supplies to JAC an input flag, JOK, which indicates whether
C JAC may (optionally) bypass the evaluation of Jacobian matrix data
C and instead process saved data (with the current value of scalar hl0).
C (e) It contains a new subroutine that calculates the initial step
C size to be attempted.
C
C
C Introduction to the Krylov methods in SLSODKR:
C
C The linear systems that must be solved have the form
C   A * x  = b ,  where  A = identity - hl0 * (df/dy) .
C Here hl0 is a scalar, and df/dy is the Jacobian matrix of partial
C derivatives of f (NEQ by NEQ).
C
C The particular Krylov method is chosen by setting the second digit,
C MITER, in the method flag MF.
C Currently, the values of MITER have the following meanings:
C
C  MITER = 1 means the Scaled Preconditioned Incomplete
C            Orthogonalization Method (SPIOM).
C
C          2 means an incomplete version of the preconditioned scaled
C            Generalized Minimal Residual method (SPIGMR).
C            This is the best choice in general.
C
C          3 means the Preconditioned Conjugate Gradient method (PCG).
C            Recommended only when df/dy is symmetric or nearly so.
C
C          4 means the scaled Preconditioned Conjugate Gradient method
C            (PCGS).  Recommended only when D-inverse * df/dy * D is
C            symmetric or nearly so, where D is the diagonal scaling
C            matrix with elements 1/EWT(i) (see RTOL/ATOL description).
C
C          9 means that only a user-supplied matrix P (approximating A)
C            will be used, with no Krylov iteration done.  This option
C            allows the user to provide the complete linear system
C            solution algorithm, if desired.
C
C The user can apply preconditioning to the linear system A*x = b,
C by means of arbitrary matrices (the preconditioners).
C     In the case of SPIOM and SPIGMR, one can apply left and right
C preconditioners P1 and P2, and the basic iterative method is then
C applied to the matrix (P1-inverse)*A*(P2-inverse) instead of to the
C matrix A.  The product P1*P2 should be an approximation to matrix A
C such that linear systems with P1 or P2 are easier to solve than with
C A.  Preconditioning from the left only or right only means using
C P2 = identity or P1 = identity, respectively.
C     In the case of the PCG and PCGS methods, there is only one
C preconditioner matrix P (but it can be the product of more than one).
C It should approximate the matrix A but allow for relatively
C easy solution of linear systems with coefficient matrix P.
C For PCG, P should be positive definite symmetric, or nearly so,
C and for PCGS, the scaled preconditioner D-inverse * P * D
C should be symmetric or nearly so.
C     If the Jacobian J = df/dy splits in a natural way into a sum
C J = J1 + J2, then one possible choice of preconditioners is
C     P1 = identity - hl0 * J1  and  P2 = identity - hl0 * J2
C provided each of these is easy to solve (or approximately solve).
C
C-----------------------------------------------------------------------
C References:
C 1.  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage Matrix
C     Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31 (1989),
C     pp. 40-91; also  L.L.N.L. Report UCRL-95088, Rev. 1, June 1987.
C 2.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C-----------------------------------------------------------------------
C Authors:       Alan C. Hindmarsh and Peter N. Brown
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODKR package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the demonstration
C program distributed with this solver.
C
C A. First provide a subroutine of the form:
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Provide a subroutine of the form:
C               SUBROUTINE G (NEQ, T, Y, NG, GOUT)
C               DIMENSION Y(*), GOUT(*)
C which supplies the vector function g by loading GOUT(i) with
C g(i), the i-th constraint function whose root is sought.
C
C C. Next determine (or guess) whether or not the problem is stiff.
C Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
C whose real part is negative and large in magnitude, compared to the
C reciprocal of the t span of interest.  If the problem is nonstiff,
C use a method flag MF = 10.  If it is stiff, MF should be between 21
C and 24, or possibly 29.  MF = 22 is generally the best choice.
C Use 23 or 24 only if symmetry is present.  Use MF = 29 if the
C complete linear system solution is to be provided by the user.
C The following four parameters must also be set.
C  IWORK(1) = LWP  = length of real array WP for preconditioning.
C  IWORK(2) = LIWP = length of integer array IWP for preconditioning.
C  IWORK(3) = JPRE = preconditioner type flag:
C                  = 0 for no preconditioning (P1 = P2 = P = identity)
C                  = 1 for left-only preconditioning (P2 = identity)
C                  = 2 for right-only preconditioning (P1 = identity)
C                  = 3 for two-sided preconditioning (and PCG or PCGS)
C  IWORK(4) = JACFLG = flag for whether JAC is called.
C                    = 0 if JAC is not to be called,
C                    = 1 if JAC is to be called.
C  Use JACFLG = 1 if JAC computes any nonconstant data for use in
C  preconditioning, such as Jacobian elements.
C  The arrays WP and IWP are work arrays under the user's control,
C  for use in the routines that perform preconditioning operations.
C
C D. If the problem is stiff, you must supply two routines that deal
C with the preconditioning of the linear systems to be solved.
C These are as follows:
C
C     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY,V,HL0,JOK,WP,IWP,IER)
C     DIMENSION Y(*), YSV(*), REWT(*), FTY(*), V(*), WP(*), IWP(*)
C        This routine must evaluate and preprocess any parts of the
C     Jacobian matrix df/dy involved in the preconditioners P1, P2, P.
C     The Y and FTY arrays contain the current values of y and f(t,y),
C     respectively, and YSV also contains the current value of y.
C     The array V is work space of length NEQ.
C     JAC must multiply all computed Jacobian elements by the scalar
C     -HL0, add the identity matrix, and do any factorization
C     operations called for, in preparation for solving linear systems
C     with a coefficient matrix of P1, P2, or P.  The matrix P1*P2 or P
C     should be an approximation to  identity - hl0 * (df/dy).
C     JAC should return IER = 0 if successful, and IER .ne. 0 if not.
C     (If IER .ne. 0, a smaller time step will be tried.)
C     JAC may alter Y and V, but not YSV, REWT, FTY, or HL0.
C     The JOK argument can be ignored (or see full description below).
C
C     SUBROUTINE PSOL (NEQ, T, Y, SAVF, WK, HL0, WP, IWP, B, LR, IER)
C     DIMENSION Y(*), SAVF(*), WK(*), WP(*), IWP(*), B(*)
C        This routine must solve a linear system with B as right-hand
C     side and one of the preconditioning matrices, P1, P2, or P, as
C     coefficient matrix, and return the solution vector in B.
C     LR is a flag concerning left vs right preconditioning, input
C     to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
C     In the case of the PCG or PCGS method, LR will be 3, and PSOL
C     should solve the system P*x = B with the preconditioner matrix P.
C     In the case MF = 29 (no Krylov iteration), LR will be 0,
C     and PSOL is to return in B the desired approximate solution
C     to A * x = B, where A = identity - hl0 * (df/dy).
C     PSOL can use data generated in the JAC routine and stored in
C     WP and IWP.  WK is a work array of length NEQ.
C     The argument HL0 is the current value of the scalar appearing
C     in the linear system.  If the old value, at the time of the last
C     JAC call, is needed, it must have been saved by JAC in WP.
C     on return, PSOL should set the error flag IER as follows:
C       IER = 0 if PSOL was successful,
C       IER .gt. 0 if a recoverable error occurred, meaning that the
C              time step will be retried,
C       IER .lt. 0 if an unrecoverable error occurred, meaning that the
C              solver is to stop immediately.
C
C E. Write a main program which calls Subroutine SLSODKR once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by SLSODKR.  On the first call to SLSODKR, supply arguments as
C follows:
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared External in calling program.
C NEQ    = number of first order ODEs.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          The estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             20 + 16*NEQ + 3*NG           for MF = 10,
C             45 + 17*NEQ + 3*NG + LWP     for MF = 21,
C             61 + 17*NEQ + 3*NG + LWP     for MF = 22,
C             20 + 15*NEQ + 3*NG + LWP     for MF = 23 or 24,
C             20 + 12*NEQ + 3*NG + LWP     for MF = 29.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least:
C             30            for MF = 10,
C             35 + LIWP     for MF = 21,
C             30 + LIWP     for MF = 22, 23, 24, or 29.
C LIW    = declared length of IWORK (in user's dimension).
C JAC,PSOL = names of subroutines for preconditioning.
C          These names must be declared External in the calling program.
C MF     = method flag.  Standard values are:
C          10 for nonstiff (Adams) method.
C          21 for stiff (BDF) method, with preconditioned SIOM.
C          22 for stiff method, with preconditioned GMRES method.
C          23 for stiff method, with preconditioned CG method.
C          24 for stiff method, with scaled preconditioned CG method.
C          29 for stiff method, with user's PSOL routine only.
C G      = name of subroutine for constraint functions, whose
C          roots are desired during the integration.
C          This name must be declared External in calling program.
C NG     = number of constraint functions g(i).  If there are none,
C          set NG = 0, and pass a dummy name for G.
C JROOT  = integer array of length NG for output of root information.
C          See next paragraph.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C JROOT, and possibly ATOL.
C
C F. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2 or 3  if SLSODKR was successful, negative otherwise.
C           2 means no root was found, and TOUT was reached as desired.
C           3 means a root was found prior to reaching TOUT.
C          -1 means excess work done on this call (perhaps wrong MF).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad JAC
C             or PSOL routine supplied or wrong choice of MF or
C             tolerances, or this solver is inappropriate).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means an unrecoverable error occurred in PSOL.
C JROOT  = array showing roots found if ISTATE = 3 on return.
C          JROOT(i) = 1 if g(i) has a root at T, or 0 otherwise.
C
C G. To continue the integration after a successful return, proceed
C as follows:
C (a) If ISTATE = 2 on return, reset TOUT and call SLSODKR again.
C (b) If ISTATE = 3 on return, reset ISTATE to 2 and call SLSODKR again.
C In either case, no other parameters need be reset.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSODKR.
C
C The user interface to SLSODKR consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODKR, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODKR package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the SLSODKR package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C  F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, PSOL, MF,
C  G, and NG,
C that used only for output is  JROOT,
C and those used for both input and output are
C  Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODKR to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C F      = the name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DIMENSION Y(*), YDOT(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared External in the calling program.
C
C          Subroutine F may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in F) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y below.
C
C          If quantities computed in the F routine are needed
C          externally to SLSODKR, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use SINTDY instead.
C
C NEQ    = the size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in the user-supplied routines.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODKR package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to the user-supplied routines.  Hence, if it is an array,
C          locations NEQ(2),... may be used to store other integer data
C          and pass it to the user-supplied routines. Each such routine
C          must include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  On output, Y contains the computed solution vector,
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to F, G,
C          JAC, and PSOL.  Hence its length may exceed NEQ, and
C          locations Y(NEQ+1),... may be used to store other real data
C          and pass it to the user-supplied routines.
C          (The SLSODKR package accesses only Y(1),...,Y(NEQ).)
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution y is evaluated (usually the same as TOUT).
C          If a root was found, T is the computed location of the
C          root reached first, on output.
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at T = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C             and any of the optional inputs except H0.
C             In addition, immediately following a return with
C             ISTATE = 3 (root found), NG and G may be changed.
C             (But changing NG from 0 to .gt. 0 is not allowed.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
C           2  means the integration was performed successfully.
C           3  means the integration was successful, and one or more
C              roots were found before satisfying the stop condition
C              specified by ITASK.  See JROOT.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means the PSOL routine returned an unrecoverable error
C              flag (IER .lt. 0).  The integration was successful as
C              far as T.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real working array (single precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD+1) + 3*NEQ + 3*NG + LENLS + LWP    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENLS = length of work space for linear system (Krylov)
C                  method, excluding preconditioning:
C            LENLS = 0                               if MITER = 0,
C            LENLS = NEQ*(MAXL+3) + MAXL**2          if MITER = 1,
C            LENLS = NEQ*(MAXL+3+MIN(1,MAXL-KMP))
C                 + (MAXL+3)*MAXL + 1                if MITER = 2,
C            LENLS = 6*NEQ                           if MITER = 3 or 4,
C            LENLS = 3*NEQ                           if MITER = 9.
C          (See the MF description for METH and MITER, and the
C          list of optional inputs for MAXL and KMP.)
C          LWP = length of real user work space for preconditioning
C          (see JAC/PSOL).
C          Thus if default values are used and NEQ is constant,
C          this length is:
C             20 + 16*NEQ + 3*NG           for MF = 10,
C             45 + 24*NEQ + 3*NG + LWP     for MF = 11,
C             61 + 24*NEQ + 3*NG + LWP     for MF = 12,
C             20 + 22*NEQ + 3*NG + LWP     for MF = 13 or 14,
C             20 + 19*NEQ + 3*NG + LWP     for MF = 19,
C             20 + 9*NEQ + 3*NG            for MF = 20,
C             45 + 17*NEQ + 3*NG + LWP     for MF = 21,
C             61 + 17*NEQ + 3*NG + LWP     for MF = 22,
C             20 + 15*NEQ + 3*NG + LWP     for MF = 23 or 24,
C             20 + 12*NEQ + 3*NG + LWP     for MF = 29.
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C             30                 if MITER = 0 (MF = 10 or 20),
C             30 + MAXL + LIWP   if MITER = 1 (MF = 11, 21),
C             30 + LIWP          if MITER = 2, 3, 4, or 9.
C          MAXL = 5 unless a different optional input value is given.
C          LIWP = length of integer user work space for preconditioning
C          (see conditional input list following).
C          The first few words of IWORK are used for conditional and
C          optional inputs and optional outputs.
C
C          The following 4 words in IWORK are conditional inputs,
C          required if MITER .ge. 1:
C          IWORK(1) = LWP  = length of real array WP for use in
C                     preconditioning (part of RWORK array).
C          IWORK(2) = LIWP = length of integer array IWP for use in
C                     preconditioning (part of IWORK array).
C                     The arrays WP and IWP are work arrays under the
C                     user's control, for use in the routines that
C                     perform preconditioning operations (JAC and PSOL).
C          IWORK(3) = JPRE = preconditioner type flag:
C                   = 0 for no preconditioning (P1 = P2 = P = identity)
C                   = 1 for left-only preconditioning (P2 = identity)
C                   = 2 for right-only preconditioning (P1 = identity)
C                   = 3 for two-sided preconditioning (and PCG or PCGS)
C          IWORK(4) = JACFLG = flag for whether JAC is called.
C                   = 0 if JAC is not to be called,
C                   = 1 if JAC is to be called.
C                     Use JACFLG = 1 if JAC computes any nonconstant
C                     data needed in preconditioning operations,
C                     such as some of the Jacobian elements.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to SLSODKR
C for the same problem, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODKR between calls, if
C desired (but not for use by any of the user-supplied routines).
C
C JAC    = the name of the user-supplied routine to compute any
C          Jacobian elements (or approximations) involved in the
C          matrix preconditioning operations (MITER .ge. 1).
C          It is to have the form
C            SUBROUTINE JAC(F,NEQ,T,Y,YSV,REWT,FTY,V,HL0,JOK,WP,IWP,IER)
C            DIMENSION Y(*),YSV(*),REWT(*),FTY(*),V(*),WP(*),IWP(*)
C          This routine must evaluate and preprocess any parts of the
C          Jacobian matrix df/dy used in the preconditioners P1, P2, P.
C          The Y and FTY arrays contain the current values of y and
C          f(t,y), respectively, and the YSV array also contains
C          the current y vector.  The array V is work space of length
C          NEQ for use by JAC.  REWT is the array of reciprocal error
C          weights (1/EWT).  JAC must multiply all computed Jacobian
C          elements by the scalar -HL0, add the identity matrix, and do
C          any factorization operations called for, in preparation
C          for solving linear systems with a coefficient matrix of
C          P1, P2, or P.  The matrix P1*P2 or P should be an
C          approximation to  identity - hl0 * (df/dy).  JAC should
C          return IER = 0 if successful, and IER .ne. 0 if not.
C          (If IER .ne. 0, a smaller time step will be tried.)
C          The arrays WP (of length LWP) and IWP (of length LIWP)
C          are for use by JAC and PSOL for work space and for storage
C          of data needed for the solution of the preconditioner
C          linear systems.  Their lengths and contents are under the
C          user's control.
C               The argument JOK is an input flag for optional use
C          by JAC in deciding whether to recompute Jacobian elements
C          or use saved values.  If JOK = -1, then JAC must compute
C          any relevant Jacobian elements (or approximations) used in
C          the preconditioners.  Optionally, JAC may also save these
C          elements for later reuse.  If JOK = 1, the integrator has
C          made a judgement (based on the convergence history and the
C          value of HL0) that JAC need not recompute Jacobian elements,
C          but instead use saved values, and the current value of HL0,
C          to reconstruct the preconditioner matrices, followed by
C          any required factorizations.  This may be cost-effective if
C          Jacobian elements are costly and storage is available.
C               JAC may alter Y and V, but not YSV, REWT, FTY, or HL0.
C          JAC must be declared External in the calling program.
C               Subroutine JAC may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C PSOL   = the name of the user-supplied routine for the
C          solution of preconditioner linear systems.
C          It is to have the form
C            SUBROUTINE PSOL (NEQ,T,Y,FTY,WK,HL0,WP,IWP,B,LR,IER)
C            DIMENSION Y(*),FTY(*),WK(*),WP(*),IWP(*),B(*)
C          This routine must solve a linear system with B as right-hand
C          side and one of the preconditioning matrices, P1, P2, or P,
C          as coefficient matrix, and return the solution vector in B.
C          LR is a flag concerning left vs right preconditioning, input
C          to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
C          In the case of the PCG or PCGS method, LR will be 3, and PSOL
C          should solve the system P*x = B with the preconditioner P.
C          In the case MITER = 9 (no Krylov iteration), LR will be 0,
C          and PSOL is to return in B the desired approximate solution
C          to A * x = B, where A = identity - hl0 * (df/dy).
C          PSOL can use data generated in the JAC routine and stored in
C          WP and IWP.
C          The Y and FTY arrays contain the current values of y and
C          f(t,y), respectively.  The array WK is work space of length
C          NEQ for use by PSOL.
C          The argument HL0 is the current value of the scalar appearing
C          in the linear system.  If the old value, as of the last
C          JAC call, is needed, it must have been saved by JAC in WP.
C          On return, PSOL should set the error flag IER as follows:
C	     IER = 0 if PSOL was successful,
C            IER .gt. 0 on a recoverable error, meaning that the
C                   time step will be retried,
C            IER .lt. 0 on an unrecoverable error, meaning that the
C                   solver is to stop immediately.
C          PSOL may not alter Y, FTY, or HL0.
C          PSOL must be declared External in the calling program.
C               Subroutine PSOL may access user-defined quantities in
C          NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in PSOL) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C MF     = the method flag.  Used only for input.  The legal values of
C          MF are 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, and 29.
C          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
C          METH indicates the basic linear multistep method:
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on Backward
C                     Differentiation Formulas (BDFs).
C          MITER indicates the corrector iteration method:
C            MITER = 0 means functional iteration (no linear system
C                      is involved).
C            MITER = 1 means Newton iteration with Scaled Preconditioned
C                      Incomplete Orthogonalization Method (SPIOM)
C                      for the linear systems.
C            MITER = 2 means Newton iteration with Scaled Preconditioned
C                      Incomplete Generalized Minimal Residual method
C                      (SPIGMR) for the linear systems.
C            MITER = 3 means Newton iteration with Preconditioned
C                      Conjugate Gradient method (PCG)
C                      for the linear systems.
C            MITER = 4 means Newton iteration with scaled preconditioned
C                      Conjugate Gradient method (PCGS)
C                      for the linear systems.
C            MITER = 9 means Newton iteration with only the
C                      user-supplied PSOL routine called (no Krylov
C                      iteration) for the linear systems.
C                      JPRE is ignored, and PSOL is called with LR = 0.
C          See comments in the introduction about the choice of MITER.
C          If MITER .ge. 1, the user must supply routines JAC and PSOL
C          (the names are arbitrary) as described above.
C          For MITER = 0, a dummy argument can be used.
C
C G      = the name of subroutine for constraint functions, whose
C          roots are desired during the integration.  It is to have
C          the form
C               SUBROUTINE G (NEQ, T, Y, NG, GOUT)
C               DIMENSION Y(*), GOUT(NG)
C          where NEQ, T, Y, and NG are input, and the array GOUT
C          is output.  NEQ, T, and Y have the same meaning as in
C          the F routine, and GOUT is an array of length NG.
C          For i = 1,...,NG, this routine is to load into GOUT(i)
C          the value at (t,y) of the i-th constraint function g(i).
C          SLSODKR will find roots of the g(i) of odd multiplicity
C          (i.e. sign changes) as they occur during the integration.
C          G must be declared External in the calling program.
C
C          Caution: Because of numerical errors in the functions
C          g(i) due to roundoff and integration error, SLSODKR may
C          return false roots, or return the same root at two or more
C          nearly equal values of t.  If such false roots are
C          suspected, the user should consider smaller error tolerances
C          and/or higher precision in the evaluation of the g(i).
C
C          If a root of some g(i) defines the end of the problem,
C          the input to SLSODKR should nevertheless allow integration
C          to a point slightly past that root, so that SLSODKR can
C          locate the root by interpolation.
C
C          Subroutine G may access user-defined quantities in
C          NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in G) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C NG     = number of constraint functions g(i).  If there are none,
C          set NG = 0, and pass a dummy name for G.
C
C JROOT  = integer array of length NG.  Used only for output.
C          On a return with ISTATE = 3 (one or more roots found),
C          JROOT(i) = 1 if g(i) has a root at t, or JROOT(i) = 0 if not.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C DELT    RWORK(8)  convergence test constant in Krylov iteration
C                   algorithm.  The default is .05.
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MAXL    IWORK(8)  maximum number of iterations in the SPIOM, SPIGMR,
C                   PCG, or PCGS algorithm (.le. NEQ).
C                   The default is MAXL = MIN(5,NEQ).
C
C KMP     IWORK(9)  number of vectors on which orthogonalization
C                   is done in SPIOM or SPIGMR algorithm (.le. MAXL).
C                   The default is KMP = MAXL.
C                   Note:  When KMP .lt. MAXL and MF = 22, the length
C                          of RWORK must be defined accordingly.  See
C                          the definition of RWORK above.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODKR, the variables listed
C below are quantities related to the performance of SLSODKR
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODKR, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NGE     IWORK(10) the number of g evaluations for the problem so far.
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NFE     IWORK(12) the number of f evaluations for the problem so far.
C
C NPE     IWORK(13) the number of calls to JAC so far (for evaluation
C                   of preconditioners).
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NNI     IWORK(19) number of nonlinear iterations so far (each of
C                   which calls an iterative linear solver).
C
C NLI     IWORK(20) number of linear iterations so far.
C                   Note: A measure of the success of algorithm is
C                   the average number of linear iterations per
C                   nonlinear iteration, given by NLI/NNI.
C                   If this is close to MAXL, MAXL may be too small.
C
C NPS     IWORK(21) number of preconditioning solve operations
C                   (PSOL calls) so far.
C
C NCFN    IWORK(22) number of convergence failures of the nonlinear
C                   (Newton) iteration so far.
C                   Note: A measure of success is the overall
C                   rate of nonlinear convergence failures, NCFN/NST.
C
C NCFL    IWORK(23) number of convergence failures of the linear
C                   iteration so far.
C                   Note: A measure of success is the overall
C                   rate of linear convergence failures, NCFL/NNI.
C
C NSFI    IWORK(24) number of functional iteration steps so far.
C                   Note: A measure of the extent to which the
C                   problem is nonstiff is the ratio NSFI/NST.
C
C NJEV    IWORK(25) number of JAC calls with JOK = -1 so far
C                   (number of evaluations of Jacobian data).
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21 + 3*NG      the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output
C                        to represent the estimated local error in y
C                        on the last step.  This is the vector E in
C                        the description of the error control.  It is
C                        defined only on a successful return from
C                        SLSODKR.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODKR.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSODKR, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSODKR.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCKR(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODKR (see Part 3 below).
C                             RSAV must be a real array of length 228
C                             or more, and ISAV must be an integer
C                             array of length 63 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCKR is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODKR.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODKR.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   LYH = 21 + 3*NG
C   CALL SINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODKR).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by SLSODKR directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C LYH       = 21 + 3*NG = base address in RWORK of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODKR is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODKR, and
C   (2) the four internal Common blocks
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C         /SLS002/  of length   5  (1 single precision word
C                      followed by  4 integer words),
C         /SLPK01/  of length  17  (4 single precision words
C                      followed by 13 integer words),
C         /SLSR01/  of length  14     (5 single precision words
C                      followed by  9 integer words).
C
C If SLSODKR is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODKR is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODKR call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODKR call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCKR (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the SLSODKR package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODKR call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the SVNORM
C routine (see below), and also used by SLSODKR in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) SVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = SVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine SEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of SVNORM suitable for use in the error control in SLSODKR.
C None of the arguments should be altered by SVNORM.
C For example, a user-supplied SVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19900117  DATE WRITTEN
C 19900503  Added iteration switching (functional/Newton).
C 19900802  Added flag for Jacobian-saving in user preconditioner.
C 19900910  Added new initial stepsize routine LHIN.
C 19901019  Corrected LHIN - y array restored.
C 19910909  Changed names STOPK to STOKA, PKSET to SETPK;
C           removed unused variables in driver declarations;
C           minor corrections to main prologue.
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           numerous corrections to prologues and internal comments.
C 20030603  Corrected duplicate type declaration for RUMACH.
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C 20031117  Changed internal name NPE to NJE.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODKR package.
C
C In addition to Subroutine SLSODKR, the SLSODKR package includes the
C following subroutines and function routines:
C  SLHIN    calculates a step size to be attempted initially.
C  SRCHEK   does preliminary checking for roots, and serves as an
C           interface between Subroutine SLSODKR and Subroutine SROOTS.
C  SROOTS   finds the leftmost root of a set of functions.
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted RMS-norm of a vector.
C  SSTOKA   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SSETPK   interfaces between SSTOKA and the JAC routine.
C  SSOLPK   manages solution of linear system in Newton iteration.
C  SSPIOM   performs the SPIOM algorithm.
C  SATV     computes a scaled, preconditioned product (I-hl0*J)*v.
C  SORTHOG  orthogonalizes a vector against previous basis vectors.
C  SHEFA    generates an LU factorization of a Hessenberg matrix.
C  SHESL    solves a Hessenberg square linear system.
C  SSPIGMR  performs the SPIGMR algorithm.
C  SHEQR    generates a QR factorization of a Hessenberg matrix.
C  SHELS    finds the least squares solution of a Hessenberg system.
C  SPCG     performs preconditioned conjugate gradient algorithm (PCG).
C  SPCGS    performs the PCGS algorithm.
C  SATP     computes the product A*p, where A = I - hl0*df/dy.
C  SUSOL    interfaces to the user's PSOL routine (MITER = 9).
C  SSRCKR   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SAXPY, SCOPY, SDOT, SNRM2, and SSCAL   are basic linear
C           algebra modules (from the BLAS collection).
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SVNORM, SDOT, SNRM2, RUMACH, IXSAV, and IUMACH are function
C routines.  All the others are subroutines.
C
C-----------------------------------------------------------------------
      REAL RUMACH, SVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER NEWT, NSFI, NSLJ, NJEV
      INTEGER LG0, LG1, LGX, IOWNR3, IRFND, ITASKC, NGC, NGE
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      INTEGER I, I1, I2, IER, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENIWK, LENRW, LENWM, LENWK, LIWP, LWP, MORD, MXHNL0,
     2   MXSTP0, NCFN0, NCFL0, NITER, NLI0, NNI0, NNID, NSTD, NWARN
      INTEGER IRFP, IRT, LENYH, LYHNEW
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL STIFR
      REAL ROWNR3, T0, TLAST, TOUTC
      REAL DELT, EPCON, SQRTN, RSQRTN
      REAL ATOLI, AVDIM, BIG, EWTI, H0, HMAX, HMX, RCFL,
     1   RCFN, RH, RTOLI, TCRIT, TNEXT, TOLSF, TP, SIZE
      DIMENSION MORD(2)
      LOGICAL IHIT, LAVD, LCFN, LCFL, LWARN
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following four internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODKR, SINTDY,
C SSTOKA, SSOLPK, and SATV.
C The block SLS002 is declared in subroutines SLSODKR and SSTOKA.
C The block SLSR01 is declared in subroutines SLSODKR, SRCHEK, SROOTS.
C The block SLPK01 is declared in subroutines SLSODKR, SSTOKA, SSETPK,
C and SSOLPK.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /SLS002/ STIFR, NEWT, NSFI, NSLJ, NJEV
C
      COMMON /SLSR01/ ROWNR3(2), T0, TLAST, TOUTC,
     1   LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE
C
      COMMON /SLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      ITASKC = ITASK
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF,
C and NG.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0) GO TO 608
      IF (MITER .GT. 4 .AND. MITER .LT. 9) GO TO 608
      IF (MITER .GE. 1) JPRE = IWORK(3)
      JACFLG = 0
      IF (MITER .GE. 1) JACFLG = IWORK(4)
      IF (NG .LT. 0) GO TO 630
      IF (ISTATE .EQ. 1) GO TO 35
      IF (IRFND .EQ. 0 .AND. NG .NE. NGC) GO TO 631
 35   NGC = NG
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      MAXL = MIN(5,N)
      KMP = MAXL
      DELT = 0.05E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
      MAXL = IWORK(8)
      IF (MAXL .EQ. 0) MAXL = 5
      MAXL = MIN(MAXL,N)
      KMP = IWORK(9)
      IF (KMP .EQ. 0 .OR. KMP .GT. MAXL) KMP = MAXL
      DELT = RWORK(8)
      IF (DELT .EQ. 0.0E0) DELT = 0.05E0
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C RWORK segments (in order) are denoted  G0, G1, GX, YH, WM,
C EWT, SAVF, SAVX, ACOR.
C-----------------------------------------------------------------------
 60   IF (ISTATE .EQ. 1) NYH = N
      LG0 = 21
      LG1 = LG0 + NG
      LGX = LG1 + NG
      LYHNEW = LGX + NG
      IF (ISTATE .EQ. 1) LYH = LYHNEW
      IF (LYHNEW .EQ. LYH) GO TO 62
C If ISTATE = 3 and NG was changed, shift YH to its new location. ------
      LENYH = L*NYH
      IF (LRW .LT. LYHNEW-1+LENYH) GO TO 62
      I1 = 1
      IF (LYHNEW .GT. LYH) I1 = -1
      CALL SCOPY (LENYH, RWORK(LYH), I1, RWORK(LYHNEW), I1)
      LYH = LYHNEW
 62   CONTINUE
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWK = 0
      IF (MITER .EQ. 1) LENWK = N*(MAXL+2) + MAXL*MAXL
      IF (MITER .EQ. 2)
     1   LENWK = N*(MAXL+2+MIN(1,MAXL-KMP)) + (MAXL+3)*MAXL + 1
      IF (MITER .EQ. 3 .OR. MITER .EQ. 4) LENWK = 5*N
      IF (MITER .EQ. 9) LENWK = 2*N
      LWP = 0
      IF (MITER .GE. 1) LWP = IWORK(1)
      LENWM = LENWK + LWP
      LOCWP = LENWK + 1
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LSAVX = LSAVF + N
      LACOR = LSAVX + N
      IF (MITER .EQ. 0) LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 31
      LENIWK = 0
      IF (MITER .EQ. 1) LENIWK = MAXL
      LIWP = 0
      IF (MITER .GE. 1) LIWP = IWORK(2)
      LENIW = 30 + LENIWK + LIWP
      LOCIWP = LENIWK + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 70     CONTINUE
C Load SQRT(N) and its reciprocal in Common. ---------------------------
      SQRTN = SQRT(REAL(N))
      RSQRTN = 1.0E0/SQRTN
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to SSTOKA.--------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO 80 I = 1,N
 80     RWORK(I+LSAVF-1) = RWORK(I+LWM-1)
 90   CONTINUE
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      NLI0 = 0
      NNI0 = 0
      NCFN0 = 0
      NCFL0 = 0
      NWARN = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
      NNI = 0
      NLI = 0
      NPS = 0
      NCFN = 0
      NCFL = 0
      NSFI = 0
      NJEV = 0
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
      IF (H0 .NE. 0.0E0) GO TO 180
C Call SLHIN to set initial step size H0 to be attempted. --------------
      CALL SLHIN (NEQ, N, T, RWORK(LYH), RWORK(LF0), F, TOUT, UROUND,
     1   RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0, NITER, IER)
      NFE = NFE + NITER
      IF (IER .NE. 0) GO TO 622
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
C Check for a zero of g at T. ------------------------------------------
      IRFND = 0
      TOUTC = TOUT
      IF (NGC .EQ. 0) GO TO 270
      CALL SRCHEK (1, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT)
      IF (IRT .EQ. 0) GO TO 270
      GO TO 632
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C First, SRCHEK is called to check for a root within the last step
C taken, other than the last root found there, if any.
C If ITASK = 2 or 5, and y(TN) has not yet been returned to the user
C because of an intervening root, return through Block G.
C-----------------------------------------------------------------------
 200  NSLAST = NST
C
      IRFP = IRFND
      IF (NGC .EQ. 0) GO TO 205
      IF (ITASK .EQ. 1 .OR. ITASK .EQ. 4) TOUTC = TOUT
      CALL SRCHEK (2, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT)
      IF (IRT .NE. 1) GO TO 205
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 205  CONTINUE
      IRFND = 0
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 2) GO TO 400
C
      NLI0 = NLI
      NNI0 = NNI
      NCFN0 = NCFN
      NCFL0 = NCFL
      NWARN = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) T = TCRIT
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 5) GO TO 400
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTOKA.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken,
C check for poor Newton/Krylov method performance, update EWT (if not
C at start of problem), check for too much accuracy being requested,
C and check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      NSTD = NST - NSLAST
      NNID = NNI - NNI0
      IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 255
      AVDIM = REAL(NLI - NLI0)/REAL(NNID)
      RCFN = REAL(NCFN - NCFN0)/REAL(NSTD)
      RCFL = REAL(NCFL - NCFL0)/REAL(NNID)
      LAVD = AVDIM .GT. (MAXL - 0.05E0)
      LCFN = RCFN .GT. 0.9E0
      LCFL = RCFL .GT. 0.9E0
      LWARN = LAVD .OR. LCFN .OR. LCFL
      IF (.NOT.LWARN) GO TO 255
      NWARN = NWARN + 1
      IF (NWARN .GT. 10) GO TO 255
      IF (LAVD) THEN
      MSG='SLSODKR- Warning. Poor iterative algorithm performance seen '
      CALL XERRWV (MSG, 60, 111, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (LAVD) THEN
      MSG='      at T = R1 by average no. of linear iterations = R2    '
      CALL XERRWV (MSG, 60, 111, 0, 0, 0, 0, 2, TN, AVDIM)
      ENDIF
      IF (LCFN) THEN
      MSG='SLSODKR- Warning. Poor iterative algorithm performance seen '
      CALL XERRWV (MSG, 60, 112, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (LCFN) THEN
      MSG='      at T = R1 by nonlinear convergence failure rate = R2  '
      CALL XERRWV (MSG, 60, 112, 0, 0, 0, 0, 2, TN, RCFN)
      ENDIF
      IF (LCFL) THEN
      MSG='SLSODKR- Warning. Poor iterative algorithm performance seen '
      CALL XERRWV (MSG, 60, 113, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (LCFL) THEN
      MSG='      at T = R1 by linear convergence failure rate = R2     '
      CALL XERRWV (MSG, 60, 113, 0, 0, 0, 0, 2, TN, RCFL)
      ENDIF
 255  CONTINUE
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODKR-  Warning.. Internal T(=R1) and H(=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODKR-  Above warning has been issued I1 times. '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL SSTOKA(NEQ,Y,YH,NYH,YH,EWT,SAVF,SAVX,ACOR,WM,IWM,F,JAC,PSOL)
C-----------------------------------------------------------------------
      CALL SSTOKA (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LSAVX), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), F, JAC, PSOL)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 550), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).
C Call SRCHEK to check for a root within the last step.
C Then, if no root was found, check for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
C
      IF (NGC .EQ. 0) GO TO 315
      CALL SRCHEK (3, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT)
      IF (IRT .NE. 1) GO TO 315
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 315  CONTINUE
C
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODKR.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
 425  CONTINUE
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNI
      IWORK(20) = NLI
      IWORK(21) = NPS
      IWORK(22) = NCFN
      IWORK(23) = NCFL
      IWORK(24) = NSFI
      IWORK(25) = NJEV
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODKR-  At current T (=R1), MXSTEP (=I1) steps  '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODKR-  At T(=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODKR-  At T (=R1), too much accuracy requested '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODKR- At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODKR-  At T (=R1) and step size H (=R2), the   '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 580
C KFLAG = -3.  Unrecoverable error from PSOL. --------------------------
 550  MSG = 'SLSODKR-  At T (=R1) an unrecoverable error return'
      CALL XERRWV (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      was made from Subroutine PSOL     '
      CALL XERRWV (MSG, 40, 206, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0E0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNI
      IWORK(20) = NLI
      IWORK(21) = NPS
      IWORK(22) = NCFN
      IWORK(23) = NCFL
      IWORK(24) = NSFI
      IWORK(25) = NJEV
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODKR-  ISTATE(=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODKR-  ITASK (=I1) illegal.'
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODKR- ISTATE.gt.1 but SLSODKR not initialized. '
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODKR-  NEQ (=I1) .lt. 1    '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODKR-  ISTATE = 3 and NEQ increased (I1 to I2).'
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODKR-  ITOL (=I1) illegal. '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODKR-  IOPT (=I1) illegal. '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODKR-  MF (=I1) illegal.   '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODKR-  MAXORD (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODKR-  MXSTEP (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODKR-  MXHNIL (=I1) .lt. 0 '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODKR-  TOUT (=R1) behind T (=R2)     '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODKR-  HMAX (=R1) .lt. 0.0 '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODKR-  HMIN (=R1) .lt. 0.0 '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG='SLSODKR-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) '
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG='SLSODKR-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) '
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODKR-  RTOL(I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODKR-  ATOL(I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODKR-  EWT(I1) is R1 .le. 0.0        '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODKR- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODKR-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2) '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODKR-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODKR-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODKR-  At start of problem, too much accuracy  '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODKR-  Trouble in SINTDY. ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
      GO TO 700
 630  MSG = 'SLSODKR-  NG (=I1) .lt. 0     '
      CALL XERRWV (MSG, 30, 30, 0, 1, NG, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 631  MSG = 'SLSODKR-  NG changed (from I1 to I2) illegally,   '
      CALL XERRWV (MSG, 50, 31, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      i.e. not immediately after a root was found.'
      CALL XERRWV (MSG, 50, 31, 0, 2, NGC, NG, 0, 0.0E0, 0.0E0)
      GO TO 700
 632  MSG = 'SLSODKR-  One or more components of g has a root  '
      CALL XERRWV (MSG, 50, 32, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      too near to the initial point.    '
      CALL XERRWV (MSG, 40, 32, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODKR-  Run aborted.. apparent infinite loop.   '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODKR ---------------------
      END
*DECK SLSODI
      SUBROUTINE SLSODI (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 18 November 2003 version of
C SLSODI: Livermore Solver for Ordinary Differential Equations
C         (Implicit form).
C
C This version is in single precision.
C
C SLSODI solves the initial value problem for linearly implicit
C systems of first order ODEs,
C     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
C or, in component form,
C     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
C        i,1      1                     i,NEQ      NEQ
C
C      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
C           i      1   2       NEQ
C
C If A is singular, this is a differential-algebraic system.
C
C SLSODI is a variant version of the SLSODE package.
C-----------------------------------------------------------------------
C Reference:
C     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C-----------------------------------------------------------------------
C Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODI package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First, provide a subroutine of the form:
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DIMENSION Y(*), S(*), R(*)
C which computes the residual function
C     r = g(t,y)  -  A(t,y) * s ,
C as a function of t and the vectors y and s.  (s is an internally
C generated approximation to dy/dt.)  The arrays Y and S are inputs
C to the RES routine and should not be altered.  The residual
C vector is to be stored in the array R.  The argument IRES should be
C ignored for casual use of SLSODI.  (For uses of IRES, see the
C paragraph on RES in the full description below.)
C
C B. Next, decide whether full or banded form is more economical
C for the storage of matrices.  SLSODI must deal internally with the
C matrices A and dr/dy, where r is the residual function defined above.
C SLSODI generates a linear combination of these two matrices, and
C this is treated in either full or banded form.
C     The matrix structure is communicated by a method flag MF,
C which is 21 or 22 for the full case, and 24 or 25 in the band case.
C     In the banded case, SLSODI requires two half-bandwidth
C parameters ML and MU.  These are, respectively, the widths of the
C lower and upper parts of the band, excluding the main diagonal.
C Thus the band consists of the locations (i,j) with
C i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
C Note that the band must accommodate the nonzero elements of
C A(t,y), dg/dy, and d(A*s)/dy (s fixed).  Alternatively, one
C can define a band that encloses only the elements that are relatively
C large in magnitude, and gain some economy in storage and possibly
C also efficiency, although the appropriate threshhold for
C retaining matrix elements is highly problem-dependent.
C
C C. You must also provide a subroutine of the form:
C               SUBROUTINE ADDA (NEQ, T, Y, ML, MU, P, NROWP)
C               DIMENSION Y(*), P(NROWP,*)
C which adds the matrix A = A(t,y) to the contents of the array P.
C T and the Y array are input and should not be altered.
C     In the full matrix case, this routine should add elements of
C to P in the usual order.  I.e., add A(i,j) to P(i,j).  (Ignore the
C ML and MU arguments in this case.)
C     In the band matrix case, this routine should add element A(i,j)
C to P(i-j+MU+1,j).  I.e., add the diagonal lines of A to the rows of
C P from the top down (the top line of A added to the first row of P).
C
C D. For the sake of efficiency, you are encouraged to supply the
C Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
C (s = a fixed vector) as above.  If dr/dy is being supplied,
C use MF = 21 or 24, and provide a subroutine of the form:
C               SUBROUTINE JAC (NEQ, T, Y, S, ML, MU, P, NROWP)
C               DIMENSION Y(*), S(*), P(NROWP,*)
C which computes dr/dy as a function of t, y, and s.  Here T, Y, and
C S are inputs, and the routine is to load dr/dy into P as follows:
C     In the full matrix case (MF = 21), load P(i,j) with dr(i)/dy(j),
C the partial derivative of r(i) with respect to y(j).  (Ignore the
C ML and MU arguments in this case.)
C     In the band matrix case (MF = 24), load P(i-j+mu+1,j) with
C dr(i)/dy(j), i.e. load the diagonal lines of dr/dy into the rows of
C P from the top down.
C     In either case, only nonzero elements need be loaded, and the
C indexing of P is the same as in the ADDA routine.
C     Note that if A is independent of y (or this dependence
C is weak enough to be ignored) then JAC is to compute dg/dy.
C     If it is not feasible to provide a JAC routine, use
C MF = 22 or 25, and SLSODI will compute an approximate Jacobian
C internally by difference quotients.
C
C E. Next decide whether or not to provide the initial value of the
C derivative vector dy/dt.  If the initial value of A(t,y) is
C nonsingular (and not too ill-conditioned), you may let SLSODI compute
C this vector (ISTATE = 0).  (SLSODI will solve the system A*s = g for
C s, with initial values of A and g.)  If A(t,y) is initially
C singular, then the system is a differential-algebraic system, and
C you must make use of the particular form of the system to compute the
C initial values of y and dy/dt.  In that case, use ISTATE = 1 and
C load the initial value of dy/dt into the array YDOTI.
C The input array YDOTI and the initial Y array must be consistent with
C the equations A*dy/dt = g.  This implies that the initial residual
C r = g(t,y) - A(t,y)*YDOTI  must be approximately zero.
C
C F. Write a main program which calls Subroutine SLSODI once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by SLSODI.  On the first call to SLSODI, supply arguments as follows:
C RES    = name of user subroutine for residual function r.
C ADDA   = name of user subroutine for computing and adding A(t,y).
C JAC    = name of user subroutine for Jacobian matrix dr/dy
C          (MF = 21 or 24).  If not used, pass a dummy name.
C Note: the names for the RES and ADDA routines and (if used) the
C        JAC routine must be declared External in the calling program.
C NEQ    = number of scalar equations in the system.
C Y      = array of initial values, of length NEQ.
C YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
C          initial dy/dt is supplied, and 0 otherwise.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
C             22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least 20 + NEQ.
C          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
C          and upper half-bandwidths ML,MU.
C LIW    = declared length of IWORK (in user's dimension).
C MF     = method flag.  Standard values are:
C          21 for a user-supplied full Jacobian.
C          22 for an internally generated full Jacobian.
C          24 for a user-supplied banded Jacobian.
C          25 for an internally generated banded Jacobian.
C          for other choices of MF, see the paragraph on MF in
C          the full description below.
C Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
C and possibly ATOL.
C
C G. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if SLSODI was successful, negative otherwise.
C          -1 means excess work done on this call (check all inputs).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 cannot occur in casual use.
C          -8 means SLSODI was unable to compute the initial dy/dt.
C             In casual use, this means A(t,y) is initially singular.
C             Supply YDOTI and use ISTATE = 1 on the first call.
C
C  If SLSODI returns ISTATE = -1, -4, or -5, then the output of
C  SLSODI also includes YDOTI = array containing residual vector
C  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
C
C H. To continue the integration after a successful return, simply
C reset TOUT and call SLSODI again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is a simple example problem, with the coding
C needed for its solution by SLSODI.  The problem is from chemical
C kinetics, and consists of the following three equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C       0.   = y1 + y2 + y3 - 1.
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.
C
C The following coding solves this problem with SLSODI, using MF = 21
C and printing results at t = .4, 4., ..., 4.e10.  It uses
C ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
C y2 has much smaller values.  dy/dt is supplied in YDOTI. We had
C obtained the initial value of dy3/dt by differentiating the
C third equation and evaluating the first two at t = 0.
C At the end of the run, statistical quantities of interest are
C printed (see optional outputs in the full description below).
C
C     EXTERNAL RESID, APLUSP, DGBYDY
C     REAL ATOL, RTOL, RWORK, T, TOUT, Y, YDOTI
C     DIMENSION Y(3), YDOTI(3), ATOL(3), RWORK(58), IWORK(23)
C     NEQ = 3
C     Y(1) = 1.
C     Y(2) = 0.
C     Y(3) = 0.
C     YDOTI(1) = -.04
C     YDOTI(2) =  .04
C     YDOTI(3) =  0.
C     T = 0.
C     TOUT = .4
C     ITOL = 2
C     RTOL = 1.E-4
C     ATOL(1) = 1.E-6
C     ATOL(2) = 1.E-10
C     ATOL(3) = 1.E-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 58
C     LIW = 23
C     MF = 21
C     DO 40  IOUT = 1,12
C       CALL SLSODI(RESID, APLUSP, DGBYDY, NEQ, Y, YDOTI, T, TOUT, ITOL,
C    1     RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF)
C       WRITE (6,20)  T, Y(1), Y(2), Y(3)
C  20   FORMAT(' At t =',E12.4,'   Y =',3E14.6)
C       IF (ISTATE .LT. 0 )  GO TO 80
C  40   TOUT = TOUT*10.
C     WRITE (6,60)  IWORK(11), IWORK(12), IWORK(13)
C  60 FORMAT(/' No. steps =',I4,'  No. r-s =',I4,'  No. J-s =',I4)
C     STOP
C  80 WRITE (6,90)  ISTATE
C  90 FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE RESID(NEQ, T, Y, S, R, IRES)
C     REAL T, Y, S, R
C     DIMENSION Y(3), S(3), R(3)
C     R(1) = -.04*Y(1) + 1.E4*Y(2)*Y(3) - S(1)
C     R(2) = .04*Y(1) - 1.E4*Y(2)*Y(3) - 3.E7*Y(2)*Y(2) - S(2)
C     R(3) = Y(1) + Y(2) + Y(3) - 1.
C     RETURN
C     END
C
C     SUBROUTINE APLUSP(NEQ, T, Y, ML, MU, P, NROWP)
C     REAL T, Y, P
C     DIMENSION Y(3), P(NROWP,3)
C     P(1,1) = P(1,1) + 1.
C     P(2,2) = P(2,2) + 1.
C     RETURN
C     END
C
C     SUBROUTINE DGBYDY(NEQ, T, Y, S, ML, MU, P, NROWP)
C     REAL T, Y, S, P
C     DIMENSION Y(3), S(3), P(NROWP,3)
C     P(1,1) = -.04
C     P(1,2) = 1.E4*Y(3)
C     P(1,3) = 1.E4*Y(2)
C     P(2,1) = .04
C     P(2,2) = -1.E4*Y(3) - 6.E7*Y(2)
C     P(2,3) = -1.E4*Y(2)
C     P(3,1) = 1.
C     P(3,2) = 1.
C     P(3,3) = 1.
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C   At t =  4.0000e-01   Y =  9.851726e-01  3.386406e-05  1.479357e-02
C   At t =  4.0000e+00   Y =  9.055142e-01  2.240418e-05  9.446344e-02
C   At t =  4.0000e+01   Y =  7.158050e-01  9.184616e-06  2.841858e-01
C   At t =  4.0000e+02   Y =  4.504846e-01  3.222434e-06  5.495122e-01
C   At t =  4.0000e+03   Y =  1.831701e-01  8.940379e-07  8.168290e-01
C   At t =  4.0000e+04   Y =  3.897016e-02  1.621193e-07  9.610297e-01
C   At t =  4.0000e+05   Y =  4.935213e-03  1.983756e-08  9.950648e-01
C   At t =  4.0000e+06   Y =  5.159269e-04  2.064759e-09  9.994841e-01
C   At t =  4.0000e+07   Y =  5.306413e-05  2.122677e-10  9.999469e-01
C   At t =  4.0000e+08   Y =  5.494532e-06  2.197826e-11  9.999945e-01
C   At t =  4.0000e+09   Y =  5.129457e-07  2.051784e-12  9.999995e-01
C   At t =  4.0000e+10   Y = -7.170472e-08 -2.868188e-13  1.000000e+00
C
C   No. steps = 330  No. r-s = 404  No. J-s =  69
C
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSODI.
C
C The user interface to SLSODI consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODI, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODI package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the SLSODI package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
C     IOPT, LRW, LIW, MF,
C and those used for both input and output are
C     Y, T, ISTATE, YDOTI.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODI to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C RES    = the name of the user-supplied subroutine which supplies
C          the residual vector for the ODE system, defined by
C            r = g(t,y) - A(t,y) * s
C          as a function of the scalar t and the vectors
C          s and y (s approximates dy/dt).  This subroutine
C          is to have the form
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DIMENSION Y(*), S(*), R(*)
C          where NEQ, T, Y, S, and IRES are input, and R and
C          IRES are output.  Y, S, and R are arrays of length NEQ.
C             On input, IRES indicates how SLSODI will use the
C          returned array R, as follows:
C             IRES = 1  means that SLSODI needs the full residual,
C                       r = g - A*s, exactly.
C             IRES = -1 means that SLSODI is using R only to compute
C                       the Jacobian dr/dy by difference quotients.
C          The RES routine can ignore IRES, or it can omit some terms
C          if IRES = -1.  If A does not depend on y, then RES can
C          just return R = g when IRES = -1.  If g - A*s contains other
C          additive terms that are independent of y, these can also be
C          dropped, if done consistently, when IRES = -1.
C             The subroutine should set the flag IRES if it
C          encounters a halt condition or illegal input.
C          Otherwise, it should not reset IRES.  On output,
C             IRES = 1 or -1 represents a normal return, and
C          SLSODI continues integrating the ODE.  Leave IRES
C          unchanged from its input value.
C             IRES = 2 tells SLSODI to immediately return control
C          to the calling program, with ISTATE = 3.  This lets
C          the calling program change parameters of the problem,
C          if necessary.
C             IRES = 3 represents an error condition (for example, an
C          illegal value of y).  SLSODI tries to integrate the system
C          without getting IRES = 3 from RES.  If it cannot, SLSODI
C          returns with ISTATE = -7 or -1.
C             On an SLSODI return with ISTATE = 3, -1, or -7, the values
C          of T and Y returned correspond to the last point reached
C          successfully without getting the flag IRES = 2 or 3.
C             The flag values IRES = 2 and 3 should not be used to
C          handle switches or root-stop conditions.  This is better
C          done by calling SLSODI in a one-step mode and checking the
C          stopping function for a sign change at each step.
C             If quantities computed in the RES routine are needed
C          externally to SLSODI, an extra call to RES should be made
C          for this purpose, for consistent and accurate results.
C          To get the current dy/dt for the S argument, use SINTDY.
C             RES must be declared External in the calling
C          program.  See note below for more about RES.
C
C ADDA   = the name of the user-supplied subroutine which adds the
C          matrix A = A(t,y) to another matrix stored in the same form
C          as A.  The storage form is determined by MITER (see MF).
C          This subroutine is to have the form
C               SUBROUTINE ADDA (NEQ, T, Y, ML, MU, P, NROWP)
C               DIMENSION Y(*), P(NROWP,*)
C          where NEQ, T, Y, ML, MU, and NROWP are input and P is
C          output.  Y is an array of length NEQ, and the matrix P is
C          stored in an NROWP by NEQ array.
C             In the full matrix case ( MITER = 1 or 2) ADDA should
C          add  A    to P(i,j).  ML and MU are ignored.
C                i,j
C             In the band matrix case ( MITER = 4 or 5) ADDA should
C          add  A    to  P(i-j+MU+1,j).
C                i,j
C          See JAC for details on this band storage form.
C             ADDA must be declared External in the calling program.
C          See note below for more information about ADDA.
C
C JAC    = the name of the user-supplied subroutine which supplies the
C          Jacobian matrix, dr/dy, where r = g - A*s.  The form of the
C          Jacobian matrix is determined by MITER.  JAC is required
C          if MITER = 1 or 4 -- otherwise a dummy name can be
C          passed.  This subroutine is to have the form
C               SUBROUTINE JAC ( NEQ, T, Y, S, ML, MU, P, NROWP )
C               DIMENSION Y(*), S(*), P(NROWP,*)
C          where NEQ, T, Y, S, ML, MU, and NROWP are input and P
C          is output.  Y and S are arrays of length NEQ, and the
C          matrix P is stored in an NROWP by NEQ array.
C          P is to be loaded with partial derivatives (elements
C          of the Jacobian matrix) on output.
C             In the full matrix case (MITER = 1), ML and MU
C          are ignored and the Jacobian is to be loaded into P
C          by columns-- i.e., dr(i)/dy(j) is loaded into P(i,j).
C             In the band matrix case (MITER = 4), the elements
C          within the band are to be loaded into P by columns,
C          with diagonal lines of dr/dy loaded into the
C          rows of P.  Thus dr(i)/dy(j) is to be loaded
C          into P(i-j+MU+1,j).  The locations in P in the two
C          triangular areas which correspond to nonexistent matrix
C          elements can be ignored or loaded arbitrarily, as they
C          they are overwritten by SLSODI.  ML and MU are the
C          half-bandwidth parameters (see IWORK).
C               In either case, P is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to RES with the same
C          arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by RES and not recomputed by JAC
C          if desired.  Also, JAC may alter the Y array, if desired.
C               JAC need not provide dr/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               JAC must be declared External in the calling program.
C               See note below for more about JAC.
C
C    Note on RES, ADDA, and JAC:
C          These subroutines may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in the subroutines) and/or Y has length
C          exceeding NEQ(1).  However, these routines should not alter
C          NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
C          See the descriptions of NEQ and Y below.
C
C NEQ    = the size of the system (number of first order ordinary
C          differential equations or scalar algebraic equations).
C          Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in RES, ADDA, or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODI package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to RES, ADDA, and JAC.  Hence, if it is an array,
C          locations NEQ(2),... may be used to store other integer data
C          and pass it to RES, ADDA, or JAC.  Each such subroutine
C          must include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 0 or 1), and only for output on other
C          calls.  On the first call, Y must contain the vector of
C          initial values.  On output, Y contains the computed solution
C          vector, evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to RES,
C          ADDA, and JAC.  Hence its length may exceed NEQ,
C          and locations Y(NEQ+1),... may be used to store other real
C          data and pass it to RES, ADDA, or JAC.  (The SLSODI
C          package accesses only Y(1),...,Y(NEQ). )
C
C YDOTI  = a real array for the initial value of the vector
C          dy/dt and for work space, of dimension at least NEQ.
C
C          On input:
C            If ISTATE = 0, then SLSODI will compute the initial value
C          of dy/dt, if A is nonsingular.  Thus YDOTI will
C          serve only as work space and may have any value.
C            If ISTATE = 1, then YDOTI must contain the initial value
C          of dy/dt.
C            If ISTATE = 2 or 3 (continuation calls), then YDOTI
C          may have any value.
C            Note: If the initial value of A is singular, then
C          SLSODI cannot compute the initial value of dy/dt, so
C          it must be provided in YDOTI, with ISTATE = 1.
C
C          On output, when SLSODI terminates abnormally with ISTATE =
C          -1, -4, or -5, YDOTI will contain the residual
C          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
C          its initial value, and YDOTI is supplied with ISTATE = 1,
C          then there may have been an incorrect input value of
C          YDOTI = dy/dt, or the problem (as given to SLSODI)
C          may not have a solution.
C
C          If desired, the YDOTI array may be used for other
C          purposes between calls to the solver.
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          on an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 0 or 1), TOUT may be
C          equal to T for one call, then should .ne. T for the next
C          call.  For the initial T, an input value of TOUT .ne. T is
C          used in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          0  means this is the first call for the problem, and
C             SLSODI is to compute the initial value of dy/dt
C             (while doing other initializations).  See note below.
C          1  means this is the first call for the problem, and
C             the initial value of dy/dt has been supplied in
C             YDOTI (SLSODI will do other initializations).  See note
C             below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
C             and any of the optional inputs except H0.
C             (See IWORK description for ML and MU.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 0 or 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           0 or 1  means nothing was done; TOUT = t and
C              ISTATE = 0 or 1 on input.
C           2  means that the integration was performed successfully.
C           3  means that the user-supplied Subroutine RES signalled
C              SLSODI to halt the integration and return (IRES = 2).
C              Integration as far as T was achieved with no occurrence
C              of IRES = 2, but this flag was set on attempting the
C              next step.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix.
C          -6  means EWT(i) became zero for some i during the
C              integration.  pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              the integration was successful as far as T.
C          -7  means that the user-supplied Subroutine RES set
C              its error flag (IRES = 3) despite repeated tries by
C              SLSODI to avoid that condition.
C          -8  means that ISTATE was 0 on input but SLSODI was unable
C              to compute the initial value of dy/dt.  See the
C              printed message for details.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Similarly, ISTATE (= 3) need not be reset if RES told
C          SLSODI to return because the calling program must change
C          the parameters of the problem.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real working array (single precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENWM   = NEQ**2 + 2    if MITER is 1 or 2, and
C          LENWM   = (2*ML+MU+1)*NEQ + 2 if MITER is 4 or 5.
C          (See MF description for the definition of METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is
C             22 + 16*NEQ + NEQ**2         for MF = 11 or 12,
C             22 + 17*NEQ + (2*ML+MU)*NEQ  for MF = 14 or 15,
C             22 +  9*NEQ + NEQ**2         for MF = 21 or 22,
C             22 + 10*NEQ + (2*ML+MU)*NEQ  for MF = 24 or 25.
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C          20 + NEQ .  The first few words of IWORK are used for
C          conditional and optional inputs and optional outputs.
C
C          The following 2 words in IWORK are conditional inputs:
C            IWORK(1) = ML     These are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       matrices in the problem-- the Jacobian dr/dy
C                       and the left-hand side matrix A. These
C                       half-bandwidths exclude the main diagonal,
C                       so the total bandwidth is ML + MU + 1 .
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if MITER is 4 or 5, and
C                       ignored otherwise.
C                       ML and MU may in fact be the band parameters
C                       for matrices to which dr/dy and A are only
C                       approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to SLSODI
C for the same problem, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODI between calls, if
C desired (but not for use by RES, ADDA, or JAC).
C
C MF     = the method flag.  Used only for input.  The legal values of
C          MF are 11, 12, 14, 15, 21, 22, 24, and 25.
C          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
C            METH indicates the basic linear multistep method:
C              METH = 1 means the implicit Adams method.
C              METH = 2 means the method based on Backward
C                       Differentiation Formulas (BDFs).
C                The BDF method is strongly preferred for stiff
C              problems, while the Adams method is preferred when
C              the problem is not stiff.  If the matrix A(t,y) is
C              nonsingular, stiffness here can be taken to mean that of
C              the explicit ODE system dy/dt = A-inverse * g.  If A is
C              singular, the concept of stiffness is not well defined.
C                If you do not know whether the problem is stiff, we
C              recommend using METH = 2.  If it is stiff, the advantage
C              of METH = 2 over METH = 1 will be great, while if it is
C              not stiff, the advantage of METH = 1 will be slight.
C              If maximum efficiency is important, some experimentation
C              with METH may be necessary.
C            MITER indicates the corrector iteration method:
C              MITER = 1 means chord iteration with a user-supplied
C                        full (NEQ by NEQ) Jacobian.
C              MITER = 2 means chord iteration with an internally
C                        generated (difference quotient) full Jacobian.
C                        This uses NEQ+1 extra calls to RES per dr/dy
C                        evaluation.
C              MITER = 4 means chord iteration with a user-supplied
C                        banded Jacobian.
C              MITER = 5 means chord iteration with an internally
C                        generated banded Jacobian (using ML+MU+2
C                        extra calls to RES per dr/dy evaluation).
C              If MITER = 1 or 4, the user must supply a Subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For other values of MITER, a dummy argument can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C the use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODI, the variables listed
C below are quantities related to the performance of SLSODI
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODI, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
C input) or -8, they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NRE     IWORK(12) the number of residual evaluations (RES calls)
C                   for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (each involving
C                   an evaluation of A and dr/dy) for the problem so
C                   far.  This equals the number of calls to ADDA and
C                   (if MITER = 1 or 4) JAC, and the number of matrix
C                   LU decompositions.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output to
C                        represent the estimated local error in y on the
C                        last step. This is the vector E in the descrip-
C                        tion of the error control.  It is defined only
C                        on a return from SLSODI with ISTATE = 2.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODI.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSODI, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSODI.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCOM(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODI (see Part 3 below).
C                             RSAV must be a real array of length 218
C                             or more, and ISAV must be an integer
C                             array of length 37 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCOM is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODI.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODI.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   CALL SINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODI).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by SLSODI directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODI is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODI, and
C   (2) the internal Common block
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C
C If SLSODI is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common block in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODI is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODI call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODI call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCOM (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the SLSODI package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODI call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the SVNORM
C routine (see below), and also used by SLSODI in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) SVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = SVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine SEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of SVNORM suitable for use in the error control in SLSODI.
C None of the arguments should be altered by SVNORM.
C For example, a user-supplied SVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19800424  DATE WRITTEN
C 19800519  Corrected access of YH on forced order reduction;
C           numerous corrections to prologues and other comments.
C 19800617  In main driver, added loading of SQRT(UROUND) in RWORK;
C           minor corrections to main prologue.
C 19800903  Corrected ISTATE logic; minor changes in prologue.
C 19800923  Added zero initialization of HU and NQU.
C 19801028  Reorganized RES calls in AINVG, STODI, and PREPJI;
C           in LSODI, corrected NRE increment and reset LDY0 at 580;
C           numerous corrections to main prologue.
C 19801218  Revised XERRWV routine; minor corrections to main prologue.
C 19810330  Added Common block /LSI001/; use LSODE's INTDY and SOLSY;
C           minor corrections to XERRWV and error message at 604;
C           minor corrections to declarations; corrections to prologues.
C 19810818  Numerous revisions: replaced EWT by 1/EWT; used flags
C           JCUR, ICF, IERPJ, IERSL between STODI and subordinates;
C           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF;
C           reorganized returns from STODI; reorganized type decls.;
C           fixed message length in XERRWV; changed default LUNIT to 6;
C           changed Common lengths; changed comments throughout.
C 19820906  Corrected use of ABS(H) in STODI; minor comment fixes.
C 19830510  Numerous revisions: revised diff. quotient increment;
C           eliminated block /LSI001/, using IERPJ flag;
C           revised STODI logic after PJAC return;
C           revised tuning of H change and step attempts in STODI;
C           corrections to main prologue and internal comments.
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODI;
C           in STODI, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           converted arithmetic IF statements to logical IF statements;
C           numerous corrections to prologues and internal comments.
C 20020502  Minor revisions to main prologue.
C 20031105  Restored 'own' variables to Common block, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C 20031117  Changed internal names NRE, LSAVR to NFE, LSAVF resp.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODI package.
C
C In addition to Subroutine SLSODI, the SLSODI package includes the
C following subroutines and function routines:
C  SAINVG   computes the initial value of the vector
C             dy/dt = A-inverse * g
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODI   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPREPJI  computes and preprocesses the Jacobian matrix
C           and the Newton iteration matrix P.
C  SSOLSY   manages solution of linear system in chord iteration.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted RMS-norm of a vector.
C  SSRCOM   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SGEFA and SGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  SGBFA and SGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SVNORM, RUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL SPREPJI, SSOLSY
      REAL RUMACH, SVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IER, IFLAG, IMXER, IRES, KGO,
     1   LENIW, LENRW, LENWM, LP, LYD0, ML, MORD, MU, MXHNL0, MXSTP0
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following internal Common block contains
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODI, SINTDY, SSTODI,
C SPREPJI, and SSOLSY.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 0 or 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 0 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .LE. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 0 or 1)
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .LE. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LE. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .EQ. 3)  GO TO 608
      IF (MITER .LT. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .LE. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .GT. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .LE. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .LE. 2) LENWM = N*N + 2
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 70     CONTINUE
      IF (ISTATE .LE. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to SSTODI. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
      DO 80 I = 1,N
 80     YDOTI(I) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
 90   RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 0 or 1).
C It contains all remaining initializations, the call to SAINVG
C (if ISTATE = 1), and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 105
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 105  JSTART = 0
      RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NFE = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Compute initial dy/dt, if necessary, and load it and initial Y into YH
      LYD0 = LYH + NYH
      LP = LWM + 1
      IF (ISTATE .EQ. 1) GO TO 120
C SLSODI must compute initial dy/dt (LYD0 points to YH(*,2)). ----------
         CALL SAINVG( RES, ADDA, NEQ, T, Y, RWORK(LYD0), MITER,
     1                ML, MU, RWORK(LP), IWORK(21), IER )
         NFE = NFE + 1
         IF (IER .LT. 0) GO TO 560
         IF (IER .GT. 0) GO TO 565
         DO 115 I = 1,N
 115        RWORK(I+LYH-1) = Y(I)
         GO TO 130
C Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
 120     DO 125 I = 1,N
            RWORK(I+LYH-1) = Y(I)
 125        RWORK(I+LYD0-1) = YDOTI(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
 130  CONTINUE
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 135 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 135    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
C                                       1
C where   w0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(i) = i-th component of initial value of dy/dt,
C         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 145
      DO 140 I = 1,N
 140    TOL = MAX(TOL,RTOL(I))
 145  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SVNORM (N, RWORK(LYD0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LYD0-1) = H0*RWORK(I+LYD0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODI.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODI-  Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODI-  Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL SSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
C                 ADDA,JAC,SPREPJI,SSOLSY)
C Note: SAVF in SSTODI occupies the same space as YDOTI in SLSODI.
C-----------------------------------------------------------------------
      CALL SSTODI (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   YDOTI, RWORK(LSAVF), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), RES, ADDA, JAC, SPREPJI, SSOLSY )
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 400, 550), KGO
C
C KGO = 1:success; 2:error test failure; 3:convergence failure;
C       4:RES ordered return. 5:RES returned error.
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  see if TOUT or TCRIT was reached.  adjust h if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODI.
C if ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      IF (KFLAG .EQ. -3) ISTATE = 3
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODI-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODI-  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 590
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODI-  At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 590
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODI-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 570
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODI-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 570
C IRES = 3 returned by RES, despite retries by SSTODI. -----------------
 550  MSG = 'SLSODI-  At T (=R1) residual routine returned     '
      CALL XERRWV (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      error IRES = 3 repeatedly.        '
      CALL XERRWV (MSG, 40, 206, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 590
C SAINVG failed because matrix A was singular. -------------------------
 560  IER = -IER
      MSG='SLSODI- Attempt to initialize dy/dt failed:  Matrix A is    '
      CALL XERRWV (MSG, 60, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      singular.  SGEFA or SGBFA returned INFO = I1'
      CALL XERRWV (MSG, 50, 207, 0, 1, IER, 0, 0, 0.0E0, 0.0E0)
      ISTATE = -8
      RETURN
C SAINVG failed because RES set IRES to 2 or 3. ------------------------
 565  MSG = 'SLSODI-  Attempt to initialize dy/dt failed       '
      CALL XERRWV (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      because residual routine set its error flag '
      CALL XERRWV (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      to IRES = (I1)'
      CALL XERRWV (MSG, 20, 208, 0, 1, IER, 0, 0, 0.0E0, 0.0E0)
      ISTATE = -8
      RETURN
C Compute IMXER if relevant. -------------------------------------------
 570  BIG = 0.0E0
      IMXER = 1
      DO 575 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 575
        BIG = SIZE
        IMXER = I
 575    CONTINUE
      IWORK(16) = IMXER
C Compute residual if relevant. ----------------------------------------
 580  LYD0 = LYH + NYH
      DO 585  I = 1,N
         RWORK(I+LSAVF-1) = RWORK(I+LYD0-1)/H
 585     Y(I) = RWORK(I+LYH-1)
      IRES = 1
      CALL RES (NEQ, TN, Y, RWORK(LSAVF), YDOTI, IRES)
      NFE = NFE + 1
      IF (IRES .LE. 1) GO TO 595
      MSG = 'SLSODI-  Residual routine set its flag IRES       '
      CALL XERRWV (MSG, 50, 210, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      to (I1) when called for final output.       '
      CALL XERRWV (MSG, 50, 210, 0, 1, IRES, 0, 0, 0.0E0, 0.0E0)
      GO TO 595
C Set Y vector, T, and optional outputs. -------------------------------
 590  DO 592 I = 1,N
 592    Y(I) = RWORK(I+LYH-1)
 595  T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODI-  ISTATE (=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODI-  ITASK (=I1) illegal. '
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODI-  ISTATE .gt. 1 but SLSODI not initialized.'
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODI-  NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODI-  ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODI-  ITOL (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODI-  IOPT (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODI-  MF (=I1) illegal.    '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 609  MSG = 'SLSODI-  ML(=I1) illegal: .lt. 0 or .ge. NEQ(=I2) '
      CALL XERRWV (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 610  MSG = 'SLSODI-  MU(=I1) illegal: .lt. 0 or .ge. NEQ(=I2) '
      CALL XERRWV (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODI-  MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODI-  MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODI-  MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODI-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODI-  HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODI-  HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG='SLSODI-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG='SLSODI-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODI-  RTOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODI-  ATOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODI-  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODI-  TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODI-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODI-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODI-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODI-  At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODI-  Trouble in SINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODI-  Run aborted.. apparent infinite loop.    '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODI ----------------------
      END
*DECK SLSOIBT
      SUBROUTINE SLSOIBT (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 18 November 2003 version of
C SLSOIBT: Livermore Solver for Ordinary differential equations given
C          in Implicit form, with Block-Tridiagonal Jacobian treatment.
C
C This version is in single precision.
C
C SLSOIBT solves the initial value problem for linearly implicit
C systems of first order ODEs,
C     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
C or, in component form,
C     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
C        i,1      1                     i,NEQ      NEQ
C
C      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
C           i      1   2       NEQ
C
C If A is singular, this is a differential-algebraic system.
C
C SLSOIBT is a variant version of the SLSODI package, for the case where
C the matrices A, dg/dy, and d(A*s)/dy are all block-tridiagonal.
C-----------------------------------------------------------------------
C Reference:
C     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing,  R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C-----------------------------------------------------------------------
C Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Charles S. Kenney
C formerly at:   Naval Weapons Center
C                China Lake, CA 93555
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSOIBT package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First, provide a subroutine of the form:
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DIMENSION Y(*), S(*), R(*)
C which computes the residual function
C     r = g(t,y)  -  A(t,y) * s ,
C as a function of t and the vectors y and s.  (s is an internally
C generated approximation to dy/dt.)  The arrays Y and S are inputs
C to the RES routine and should not be altered.  The residual
C vector is to be stored in the array R.  The argument IRES should be
C ignored for casual use of SLSOIBT.  (For uses of IRES, see the
C paragraph on RES in the full description below.)
C
C B. Next, identify the block structure of the matrices A = A(t,y) and
C dr/dy.  SLSOIBT must deal internally with a linear combination, P, of
C these two matrices.  The matrix P (hence both A and dr/dy) must have
C a block-tridiagonal form with fixed structure parameters
C     MB = block size, MB .ge. 1, and
C     NB = number of blocks in each direction, NB .ge. 4,
C with MB*NB = NEQ.  In each of the NB block-rows of the matrix P
C (each consisting of MB consecutive rows), the nonzero elements are
C to lie in three consecutive MB by MB blocks.  In block-rows
C 2 through NB - 1, these are centered about the main diagonal.
C in block-rows 1 and NB, they are the diagonal blocks and the two
C blocks adjacent to the diagonal block.  (Thus block positions (1,3)
C and (NB,NB-2) can be nonzero.)
C Alternatively, P (hence A and dr/dy) may be only approximately
C equal to matrices with this form, and SLSOIBT should still succeed.
C The block-tridiagonal matrix P is described by three arrays,
C each of size MB by MB by NB:
C     PA = array of diagonal blocks,
C     PB = array of superdiagonal (and one subdiagonal) blocks, and
C     PC = array of subdiagonal (and one superdiagonal) blocks.
C Specifically, the three MB by MB blocks in the k-th block-row of P
C are stored in (reading across):
C     PC(*,*,k) = block to the left of the diagonal block,
C     PA(*,*,k) = diagonal block, and
C     PB(*,*,k) = block to the right of the diagonal block,
C except for k = 1, where the three blocks (reading across) are
C     PA(*,*,1) (= diagonal block), PB(*,*,1), and PC(*,*,1),
C and k = NB, where they are
C     PB(*,*,NB), PC(*,*,NB), and PA(*,*,NB) (= diagonal block).
C (Each asterisk * stands for an index that ranges from 1 to MB.)
C
C C. You must also provide a subroutine of the form:
C     SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
C     DIMENSION Y(*), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C which adds the nonzero blocks of the matrix A = A(t,y) to the
C contents of the arrays PA, PB, and PC, following the structure
C description in Paragraph B above.
C T and the Y array are input and should not be altered.
C Thus the affect of ADDA should be the following:
C     DO 30 K = 1,NB
C       DO 20 J = 1,MB
C         DO 10 I = 1,MB
C           PA(I,J,K) = PA(I,J,K) +
C             ( (I,J) element of K-th diagonal block of A)
C           PB(I,J,K) = PB(I,J,K) +
C             ( (I,J) element of block in block position (K,K+1) of A,
C             or in block position (NB,NB-2) if K = NB)
C           PC(I,J,K) = PC(I,J,K) +
C             ( (I,J) element of block in block position (K,K-1) of A,
C             or in block position (1,3) if K = 1)
C 10        CONTINUE
C 20      CONTINUE
C 30    CONTINUE
C
C D. For the sake of efficiency, you are encouraged to supply the
C Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
C (s = a fixed vector) as above.  If dr/dy is being supplied,
C use MF = 21, and provide a subroutine of the form:
C     SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
C     DIMENSION Y(*), S(*), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C which computes dr/dy as a function of t, y, and s.  Here T, Y, and
C S are inputs, and the routine is to load dr/dy into PA, PB, PC,
C according to the structure description in Paragraph B above.
C That is, load the diagonal blocks into PA, the superdiagonal blocks
C (and block (NB,NB-2) ) into PB, and the subdiagonal blocks (and
C block (1,3) ) into PC.  The blocks in block-row k of dr/dy are to
C be loaded into PA(*,*,k), PB(*,*,k), and PC(*,*,k).
C     Only nonzero elements need be loaded, and the indexing
C of PA, PB, and PC is the same as in the ADDA routine.
C     Note that if A is independent of Y (or this dependence
C is weak enough to be ignored) then JAC is to compute dg/dy.
C     If it is not feasible to provide a JAC routine, use
C MF = 22, and SLSOIBT will compute an approximate Jacobian
C internally by difference quotients.
C
C E. Next decide whether or not to provide the initial value of the
C derivative vector dy/dt.  If the initial value of A(t,y) is
C nonsingular (and not too ill-conditioned), you may let SLSOIBT compute
C this vector (ISTATE = 0).  (SLSOIBT will solve the system A*s = g for
C s, with initial values of A and g.)  If A(t,y) is initially
C singular, then the system is a differential-algebraic system, and
C you must make use of the particular form of the system to compute the
C initial values of y and dy/dt.  In that case, use ISTATE = 1 and
C load the initial value of dy/dt into the array YDOTI.
C The input array YDOTI and the initial Y array must be consistent with
C the equations A*dy/dt = g.  This implies that the initial residual
C r = g(t,y) - A(t,y)*YDOTI  must be approximately zero.
C
C F. Write a main program which calls Subroutine SLSOIBT once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages by
C SLSOIBT.  on the first call to SLSOIBT, supply arguments as follows:
C RES    = name of user subroutine for residual function r.
C ADDA   = name of user subroutine for computing and adding A(t,y).
C JAC    = name of user subroutine for Jacobian matrix dr/dy
C          (MF = 21).  If not used, pass a dummy name.
C Note: the names for the RES and ADDA routines and (if used) the
C        JAC routine must be declared External in the calling program.
C NEQ    = number of scalar equations in the system.
C Y      = array of initial values, of length NEQ.
C YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
C          initial dy/dt is supplied, and 0 otherwise.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 + 9*NEQ + 3*MB*MB*NB        for MF = 21 or 22.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least 20 + NEQ.
C          Input in IWORK(1) the block size MB and in IWORK(2) the
C          number NB of blocks in each direction along the matrix A.
C          These must satisfy  MB .ge. 1, NB .ge. 4, and MB*NB = NEQ.
C LIW    = declared length of IWORK (in user's dimension).
C MF     = method flag.  Standard values are:
C          21 for a user-supplied Jacobian.
C          22 for an internally generated Jacobian.
C          For other choices of MF, see the paragraph on MF in
C          the full description below.
C Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
C and possibly ATOL.
C
C G. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if SLSOIBT was successful, negative otherwise.
C          -1 means excess work done on this call (check all inputs).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 cannot occur in casual use.
C          -8 means SLSOIBT was unable to compute the initial dy/dt.
C             In casual use, this means A(t,y) is initially singular.
C             Supply YDOTI and use ISTATE = 1 on the first call.
C
C  If SLSOIBT returns ISTATE = -1, -4, or -5, then the output of
C  SLSOIBT also includes YDOTI = array containing residual vector
C  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
C
C H. To continue the integration after a successful return, simply
C reset TOUT and call SLSOIBT again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is an example problem, with the coding needed
C for its solution by SLSOIBT.  The problem comes from the partial
C differential equation (the Burgers equation)
C   du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
C on -1 .le. x .le. 1.  The boundary conditions are
C   du/dx = 0  at x = -1 and at x = 1.
C The initial profile is a square wave,
C   u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.
C The PDE is discretized in x by a simplified Galerkin method,
C using piecewise linear basis functions, on a grid of 40 intervals.
C The equations at x = -1 and 1 use a 3-point difference approximation
C for the right-hand side.  The result is a system A * dy/dt = g(y),
C of size NEQ = 41, where y(i) is the approximation to u at x = x(i),
C with x(i) = -1 + (i-1)*delx, delx = 2/(NEQ-1) = .05.  The individual
C equations in the system are
C   dy(1)/dt = ( y(3) - 2*y(2) + y(1) ) * eta / delx**2,
C   dy(NEQ)/dt = ( y(NEQ-2) - 2*y(NEQ-1) + y(NEQ) ) * eta / delx**2,
C and for i = 2, 3, ..., NEQ-1,
C   (1/6) dy(i-1)/dt + (4/6) dy(i)/dt + (1/6) dy(i+1)/dt
C       = ( y(i-1)**2 - y(i+1)**2 ) / (4*delx)
C         + ( y(i+1) - 2*y(i) + y(i-1) ) * eta / delx**2.
C The following coding solves the problem with MF = 21, with output
C of solution statistics at t = .1, .2, .3, and .4, and of the
C solution vector at t = .4.  Here the block size is just MB = 1.
C
C     EXTERNAL RESID, ADDABT, JACBT
C     REAL ATOL, RTOL, RWORK, T, TOUT, Y, YDOTI
C     DIMENSION Y(41), YDOTI(41), RWORK(514), IWORK(61)
C     NEQ = 41
C     DO 10 I = 1,NEQ
C  10   Y(I) = 0.0
C     Y(11) = 0.5
C     DO 20 I = 12,30
C  20   Y(I) = 1.0
C     Y(31) = 0.5
C     T = 0.0
C     TOUT = 0.1
C     ITOL = 1
C     RTOL = 1.0E-4
C     ATOL = 1.0E-5
C     ITASK = 1
C     ISTATE = 0
C     IOPT = 0
C     LRW = 514
C     LIW = 61
C     IWORK(1) = 1
C     IWORK(2) = NEQ
C     MF = 21
C     DO 40 IO = 1,4
C       CALL SLSOIBT (RESID, ADDABT, JACBT, NEQ, Y, YDOTI, T, TOUT,
C    1     ITOL,RTOL,ATOL, ITASK, ISTATE, IOPT, RWORK,LRW,IWORK,LIW, MF)
C       WRITE (6,30) T, IWORK(11), IWORK(12), IWORK(13)
C  30   FORMAT(' At t =',F5.2,'   No. steps =',I4,'  No. r-s =',I4,
C    1         '  No. J-s =',I3)
C       IF (ISTATE .NE. 2)  GO TO 90
C       TOUT = TOUT + 0.1
C  40   CONTINUE
C     WRITE(6,50) (Y(I),I=1,NEQ)
C  50 FORMAT(/' Final solution values..'/9(5E12.4/))
C     STOP
C  90 WRITE(6,95) ISTATE
C  95 FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE RESID (N, T, Y, S, R, IRES)
C     REAL T, Y, S, R, ETA, DELX, EODSQ
C     DIMENSION Y(N), S(N), R(N)
C     DATA ETA/0.05/, DELX/0.05/
C     EODSQ = ETA/DELX**2
C     R(1) = EODSQ*(Y(3) - 2.0*Y(2) + Y(1)) - S(1)
C     NM1 = N - 1
C     DO 10 I = 2,NM1
C       R(I) = (Y(I-1)**2 - Y(I+1)**2)/(4.0*DELX)
C    1        + EODSQ*(Y(I+1) - 2.0*Y(I) + Y(I-1))
C    2        - (S(I-1) + 4.0*S(I) + S(I+1))/6.0
C  10   CONTINUE
C     R(N) = EODSQ*(Y(N-2) - 2.0*Y(NM1) + Y(N)) - S(N)
C     RETURN
C     END
C
C     SUBROUTINE ADDABT (N, T, Y, MB, NB, PA, PB, PC)
C     REAL T, Y, PA, PB, PC
C     DIMENSION Y(N), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C     PA(1,1,1) = PA(1,1,1) + 1.0
C     NM1 = N - 1
C     DO 10 K = 2,NM1
C       PA(1,1,K) = PA(1,1,K) + (4.0/6.0)
C       PB(1,1,K) = PB(1,1,K) + (1.0/6.0)
C       PC(1,1,K) = PC(1,1,K) + (1.0/6.0)
C  10   CONTINUE
C     PA(1,1,N) = PA(1,1,N) + 1.0
C     RETURN
C     END
C
C     SUBROUTINE JACBT (N, T, Y, S, MB, NB, PA, PB, PC)
C     REAL T, Y, S, PA, PB, PC, ETA, DELX, EODSQ
C     DIMENSION Y(N), S(N), PA(MB,MB,NB),PB(MB,MB,NB),PC(MB,MB,NB)
C     DATA ETA/0.05/, DELX/0.05/
C     EODSQ = ETA/DELX**2
C     PA(1,1,1) = EODSQ
C     PB(1,1,1) = -2.0*EODSQ
C     PC(1,1,1) = EODSQ
C     DO 10 K = 2,N
C       PA(1,1,K) = -2.0*EODSQ
C       PB(1,1,K) = -Y(K+1)*(0.5/DELX) + EODSQ
C       PC(1,1,K) = Y(K-1)*(0.5/DELX) + EODSQ
C  10   CONTINUE
C     PB(1,1,N) = EODSQ
C     PC(1,1,N) = -2.0*EODSQ
C     PA(1,1,N) = EODSQ
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C At t = 0.10   No. steps =  35  No. r-s =  45  No. J-s =  9
C At t = 0.20   No. steps =  43  No. r-s =  54  No. J-s = 10
C At t = 0.30   No. steps =  48  No. r-s =  60  No. J-s = 11
C At t = 0.40   No. steps =  51  No. r-s =  64  No. J-s = 12
C
C Final solution values..
C  1.2747e-02  1.1997e-02  1.5560e-02  2.3767e-02  3.7224e-02
C  5.6646e-02  8.2645e-02  1.1557e-01  1.5541e-01  2.0177e-01
C  2.5397e-01  3.1104e-01  3.7189e-01  4.3530e-01  5.0000e-01
C  5.6472e-01  6.2816e-01  6.8903e-01  7.4612e-01  7.9829e-01
C  8.4460e-01  8.8438e-01  9.1727e-01  9.4330e-01  9.6281e-01
C  9.7632e-01  9.8426e-01  9.8648e-01  9.8162e-01  9.6617e-01
C  9.3374e-01  8.7535e-01  7.8236e-01  6.5321e-01  5.0003e-01
C  3.4709e-01  2.1876e-01  1.2771e-01  7.3671e-02  5.0642e-02
C  5.4496e-02
C
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSOIBT.
C
C The user interface to SLSOIBT consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSOIBT, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSOIBT package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the SLSOIBT package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
C     IOPT, LRW, LIW, MF,
C and those used for both input and output are
C     Y, T, ISTATE, YDOTI.
C The work arrays RWORK and IWORK are also used for additional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSOIBT to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C RES    = the name of the user-supplied subroutine which supplies
C          the residual vector for the ODE system, defined by
C            r = g(t,y) - A(t,y) * s
C          as a function of the scalar t and the vectors
C          s and y (s approximates dy/dt).  This subroutine
C          is to have the form
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DIMENSION Y(*), S(*), R(*)
C          where NEQ, T, Y, S, and IRES are input, and R and
C          IRES are output. Y, S, and R are arrays of length NEQ.
C             On input, IRES indicates how SLSOIBT will use the
C          returned array R, as follows:
C             IRES = 1  means that SLSOIBT needs the full residual,
C                       r = g - A*s, exactly.
C             IRES = -1 means that SLSOIBT is using R only to compute
C                       the Jacobian dr/dy by difference quotients.
C          The RES routine can ignore IRES, or it can omit some terms
C          if IRES = -1.  If A does not depend on y, then RES can
C          just return R = g when IRES = -1.  If g - A*s contains other
C          additive terms that are independent of y, these can also be
C          dropped, if done consistently, when IRES = -1.
C             The subroutine should set the flag IRES if it
C          encounters a halt condition or illegal input.
C          Otherwise, it should not reset IRES.  On output,
C             IRES = 1 or -1 represents a normal return, and
C          SLSOIBT continues integrating the ODE.  Leave IRES
C          unchanged from its input value.
C             IRES = 2 tells SLSOIBT to immediately return control
C          to the calling program, with ISTATE = 3.  This lets
C          the calling program change parameters of the problem
C          if necessary.
C             IRES = 3 represents an error condition (for example, an
C          illegal value of y).  SLSOIBT tries to integrate the system
C          without getting IRES = 3 from RES.  If it cannot, SLSOIBT
C          returns with ISTATE = -7 or -1.
C             On an SLSOIBT return with ISTATE = 3, -1, or -7, the
C          values of T and Y returned correspond to the last point
C          reached successfully without getting the flag IRES = 2 or 3.
C             The flag values IRES = 2 and 3 should not be used to
C          handle switches or root-stop conditions.  This is better
C          done by calling SLSOIBT in a one-step mode and checking the
C          stopping function for a sign change at each step.
C             If quantities computed in the RES routine are needed
C          externally to SLSOIBT, an extra call to RES should be made
C          for this purpose, for consistent and accurate results.
C          To get the current dy/dt for the S argument, use SINTDY.
C             RES must be declared External in the calling
C          program. See note below for more about RES.
C
C ADDA   = the name of the user-supplied subroutine which adds the
C          matrix A = A(t,y) to another matrix, P, stored in
C          block-tridiagonal form.  This routine is to have the form
C               SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
C               DIMENSION Y(*), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C          where NEQ, T, Y, MB, NB, and the arrays PA, PB, and PC
C          are input, and the arrays PA, PB, and PC are output.
C          Y is an array of length NEQ, and the arrays PA, PB, PC
C          are all MB by MB by NB.
C             Here a block-tridiagonal structure is assumed for A(t,y),
C          and also for the matrix P to which A is added here,
C          as described in Paragraph B of the Summary of Usage above.
C          Thus the affect of ADDA should be the following:
C               DO 30 K = 1,NB
C                 DO 20 J = 1,MB
C                   DO 10 I = 1,MB
C                     PA(I,J,K) = PA(I,J,K) +
C                       ( (I,J) element of K-th diagonal block of A)
C                     PB(I,J,K) = PB(I,J,K) +
C                       ( (I,J) element of block (K,K+1) of A,
C                       or block (NB,NB-2) if K = NB)
C                     PC(I,J,K) = PC(I,J,K) +
C                       ( (I,J) element of block (K,K-1) of A,
C                       or block (1,3) if K = 1)
C           10        CONTINUE
C           20      CONTINUE
C           30    CONTINUE
C             ADDA must be declared External in the calling program.
C          See note below for more information about ADDA.
C
C JAC    = the name of the user-supplied subroutine which supplies
C          the Jacobian matrix, dr/dy, where r = g - A*s.  JAC is
C          required if MITER = 1.  Otherwise a dummy name can be
C          passed.  This subroutine is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
C               DIMENSION Y(*), S(*),
C              1    PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C          where NEQ, T, Y, S, MB, NB, and the arrays PA, PB, and PC
C          are input, and the arrays PA, PB, and PC are output.
C          Y and S are arrays of length NEQ, and the arrays PA, PB, PC
C          are all MB by MB by NB.
C          PA, PB, and PC are to be loaded with partial derivatives
C          (elements of the Jacobian matrix) on output, in terms of the
C          block-tridiagonal structure assumed, as described
C          in Paragraph B of the Summary of Usage above.
C          That is, load the diagonal blocks into PA, the
C          superdiagonal blocks (and block (NB,NB-2) ) into PB, and
C          the subdiagonal blocks (and block (1,3) ) into PC.
C          The blocks in block-row k of dr/dy are to be loaded into
C          PA(*,*,k), PB(*,*,k), and PC(*,*,k).
C          Thus the affect of JAC should be the following:
C               DO 30 K = 1,NB
C                 DO 20 J = 1,MB
C                   DO 10 I = 1,MB
C                     PA(I,J,K) = ( (I,J) element of
C                       K-th diagonal block of dr/dy)
C                     PB(I,J,K) = ( (I,J) element of block (K,K+1)
C                       of dr/dy, or block (NB,NB-2) if K = NB)
C                     PC(I,J,K) = ( (I,J) element of block (K,K-1)
C                       of dr/dy, or block (1,3) if K = 1)
C           10        CONTINUE
C           20      CONTINUE
C           30    CONTINUE
C               PA, PB, and PC are preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to RES with the same
C          arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by RES and not recomputed by JAC
C          if desired.  Also, JAC may alter the Y array, if desired.
C               JAC need not provide dr/dy exactly.  A crude
C          approximation will do, so that SLSOIBT may be used when
C          A and dr/dy are not really block-tridiagonal, but are close
C          to matrices that are.
C               JAC must be declared External in the calling program.
C               See note below for more about JAC.
C
C    Note on RES, ADDA, and JAC:
C          These subroutines may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in the subroutines) and/or Y has length
C          exceeding NEQ(1).  However, these routines should not alter
C          NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
C          See the descriptions of NEQ and Y below.
C
C NEQ    = the size of the system (number of first order ordinary
C          differential equations or scalar algebraic equations).
C          Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in RES, ADDA, or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSOIBT package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to RES, ADDA, and JAC.  Hence, if it is an array,
C          locations NEQ(2),... may be used to store other integer data
C          and pass it to RES, ADDA, or JAC.  Each such subroutine
C          must include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 0 or 1), and only for output on other
C          calls.  On the first call, Y must contain the vector of
C          initial values.  On output, Y contains the computed solution
C          vector, evaluated at t.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to RES,
C          ADDA, and JAC.  Hence its length may exceed NEQ,
C          and locations Y(NEQ+1),... may be used to store other real
C          data and pass it to RES, ADDA, or JAC.  (The SLSOIBT
C          package accesses only Y(1),...,Y(NEQ). )
C
C YDOTI  = a real array for the initial value of the vector
C          dy/dt and for work space, of dimension at least NEQ.
C
C          On input:
C            If ISTATE = 0 then SLSOIBT will compute the initial value
C          of dy/dt, if A is nonsingular.  Thus YDOTI will
C          serve only as work space and may have any value.
C            If ISTATE = 1 then YDOTI must contain the initial value
C          of dy/dt.
C            If ISTATE = 2 or 3 (continuation calls) then YDOTI
C          may have any value.
C            Note: If the initial value of A is singular, then
C          SLSOIBT cannot compute the initial value of dy/dt, so
C          it must be provided in YDOTI, with ISTATE = 1.
C
C          On output, when SLSOIBT terminates abnormally with ISTATE =
C          -1, -4, or -5, YDOTI will contain the residual
C          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
C          its initial value, and YDOTI is supplied with ISTATE = 1,
C          there may have been an incorrect input value of
C          YDOTI = dy/dt, or the problem (as given to SLSOIBT)
C          may not have a solution.
C
C          If desired, the YDOTI array may be used for other
C          purposes between calls to the solver.
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 0 or 1), TOUT may be
C          equal to T for one call, then should .ne. T for the next
C          call.  For the initial T, an input value of TOUT .ne. T is
C          used in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          0  means this is the first call for the problem, and
C             SLSOIBT is to compute the initial value of dy/dt
C             (while doing other initializations).  See note below.
C          1  means this is the first call for the problem, and
C             the initial value of dy/dt has been supplied in
C             YDOTI (SLSOIBT will do other initializations).
C             See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, MB, NB,
C             and any of the optional inputs except H0.
C             (See IWORK description for MB and NB.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 0 or 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           0 or 1  means nothing was done; TOUT = t and
C              ISTATE = 0 or 1 on input.
C           2  means that the integration was performed successfully.
C           3  means that the user-supplied Subroutine RES signalled
C              SLSOIBT to halt the integration and return (IRES = 2).
C              Integration as far as T was achieved with no occurrence
C              of IRES = 2, but this flag was set on attempting the
C              next step.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i) = 0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means that the user-supplied Subroutine RES set
C              its error flag (IRES = 3) despite repeated tries by
C              SLSOIBT to avoid that condition.
C          -8  means that ISTATE was 0 on input but SLSOIBT was unable
C              to compute the initial value of dy/dt.  See the
C              printed message for details.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Similarly, ISTATE (= 3) need not be reset if RES told
C          SLSOIBT to return because the calling program must change
C          the parameters of the problem.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real working array (single precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENWM  = 3*MB*MB*NB + 2.
C          (See MF description for the definition of METH.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is
C             22 + 16*NEQ + 3*MB*MB*NB     for MF = 11 or 12,
C             22 + 9*NEQ + 3*MB*MB*NB      for MF = 21 or 22.
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C          20 + NEQ .  The first few words of IWORK are used for
C          additional and optional inputs and optional outputs.
C
C          The following 2 words in IWORK are additional required
C          inputs to SLSOIBT:
C            IWORK(1) = MB = block size
C            IWORK(2) = NB = number of blocks in the main diagonal
C          These must satisfy  MB .ge. 1, NB .ge. 4, and MB*NB = NEQ.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to SLSOIBT
C for the same problem, except possibly for the additional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSOIBT between calls, if
C desired (but not for use by RES, ADDA, or JAC).
C
C MF     = the method flag.  used only for input.  The legal values of
C          MF are 11, 12, 21, and 22.
C          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
C            METH indicates the basic linear multistep method:
C              METH = 1 means the implicit Adams method.
C              METH = 2 means the method based on Backward
C                       Differentiation Formulas (BDFS).
C                The BDF method is strongly preferred for stiff
C              problems, while the Adams method is preferred when the
C              problem is not stiff.  If the matrix A(t,y) is
C              nonsingular, stiffness here can be taken to mean that of
C              the explicit ODE system dy/dt = A-inverse * g.  If A is
C              singular, the concept of stiffness is not well defined.
C                If you do not know whether the problem is stiff, we
C              recommend using METH = 2.  If it is stiff, the advantage
C              of METH = 2 over METH = 1 will be great, while if it is
C              not stiff, the advantage of METH = 1 will be slight.
C              If maximum efficiency is important, some experimentation
C              with METH may be necessary.
C            MITER indicates the corrector iteration method:
C              MITER = 1 means chord iteration with a user-supplied
C                        block-tridiagonal Jacobian.
C              MITER = 2 means chord iteration with an internally
C                        generated (difference quotient) block-
C                        tridiagonal Jacobian approximation, using
C                        3*MB+1 extra calls to RES per dr/dy evaluation.
C              If MITER = 1, the user must supply a Subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For MITER = 2, a dummy argument can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSOIBT, the variables listed
C below are quantities related to the performance of SLSOIBT
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSOIBT, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
C input) or -8, they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NRE     IWORK(12) the number of residual evaluations (RES calls)
C                   for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (each involving
C                   an evaluation of a and dr/dy) for the problem so
C                   far.  This equals the number of calls to ADDA and
C                   (if MITER = 1) to JAC, and the number of matrix
C                   LU decompositions.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output to
C                        represent the estimated local error in y on
C                        the last step.  This is the vector E in the
C                        description of the error control.  It is
C                        defined only on a return from SLSOIBT with
C                        ISTATE = 2.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSOIBT.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSOIBT, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSOIBT.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCOM(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSOIBT (see Part 3 below).
C                             RSAV must be a real array of length 218
C                             or more, and ISAV must be an integer
C                             array of length 37 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCOM is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSOIBT.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSOIBT.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   CALL SINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the t last returned by SLSOIBT).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(t), is already provided
C             by SLSOIBT directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSOIBT is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSOIBT, and
C   (2) the internal Common block
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C
C If SLSOIBT is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common block in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSOIBT is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSOIBT call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSOIBT call for that problem.  To save and restore the Common
C blocks, use Subroutine SSRCOM (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the SLSOIBT package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSOIBT call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the SVNORM
C routine (see below), and also used by SLSOIBT in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) SVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = SVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine SEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of SVNORM suitable for use in the error control in SLSOIBT.
C None of the arguments should be altered by SVNORM.
C For example, a user-supplied SVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19840625  DATE WRITTEN
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODI;
C           in STODI, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           converted arithmetic IF statements to logical IF statements;
C           numerous corrections to prologues and internal comments.
C 20020502  Minor revisions to main prologue.
C 20031105  Restored 'own' variables to Common block, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C 20031117  Changed internal names NRE, LSAVR to NFE, LSAVF resp.
C
C-----------------------------------------------------------------------
C Other routines in the SLSOIBT package.
C
C In addition to Subroutine SLSOIBT, the SLSOIBT package includes the
C following subroutines and function routines:
C  SAIGBT   computes the initial value of the vector
C             dy/dt = A-inverse * g
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODI   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted RMS-norm of a vector.
C  SSRCOM   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SPJIBT   computes and preprocesses the Jacobian matrix
C           and the Newton iteration matrix P.
C  SSLSBT   manages solution of linear system in chord iteration.
C  SDECBT and SSOLBT   are routines for solving block-tridiagonal
C           systems of linear algebraic equations.
C  SGEFA and SGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  SDOT     is one of the basic linear algebra modules (BLAS).
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SVNORM, SDOT, RUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL SPJIBT, SSLSBT
      REAL RUMACH, SVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IER, IFLAG, IMXER, IRES, KGO,
     1   LENIW, LENRW, LENWM, LP, LYD0, MB, MORD, MXHNL0, MXSTP0, NB
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following internal Common block contains
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSOIBT, SINTDY, SSTODI,
C SPJIBT, and SSLSBT.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 0 or 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 0 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .LE. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 0 or 1)
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, MB, and NB.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .LE. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 1 .OR. MITER .GT. 2) GO TO 608
      MB = IWORK(1)
      NB = IWORK(2)
      IF (MB .LT. 1 .OR. MB .GT. N) GO TO 609
      IF (NB .LT. 4) GO TO 610
      IF (MB*NB .NE. N) GO TO 609
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .LE. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .GT. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .LE. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      LENWM = 3*MB*MB*NB + 2
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 70     CONTINUE
      IF (ISTATE .LE. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to SSTODI. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
      DO 80 I = 1,N
 80     YDOTI(I) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
 90   RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 0 or 1).
C It contains all remaining initializations, the call to SAIGBT
C (if ISTATE = 1), and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = RUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 105
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
 105  JSTART = 0
      RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NFE = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Compute initial dy/dt, if necessary, and load it and initial Y into YH
      LYD0 = LYH + NYH
      LP = LWM + 1
      IF ( ISTATE .EQ. 1 )  GO TO 120
C SLSOIBT must compute initial dy/dt (LYD0 points to YH(*,2)). ---------
         CALL SAIGBT( RES, ADDA, NEQ, T, Y, RWORK(LYD0),
     1               MB, NB, RWORK(LP), IWORK(21), IER )
         NFE = NFE + 1
         IF (IER .LT. 0) GO TO 560
         IF (IER .GT. 0) GO TO 565
         DO 115  I = 1,N
  115       RWORK(I+LYH-1) = Y(I)
         GO TO 130
C Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
  120    DO 125  I = 1,N
            RWORK(I+LYH-1) = Y(I)
  125       RWORK(I+LYD0-1) = YDOTI(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
  130 CONTINUE
      NQ = 1
      H = 1.0E0
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 135 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 135    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
C                                       1
C where   w0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(i) = i-th component of initial value of dy/dt,
C         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 145
      DO 140 I = 1,N
 140    TOL = MAX(TOL,RTOL(I))
 145  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SVNORM (N, RWORK(LYD0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LYD0-1) = H0*RWORK(I+LYD0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODI.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSOIBT- Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSOIBT- Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL SSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
C                 ADDA,JAC,SPJIBT,SSLSBT)
C Note: SAVF in SSTODI occupies the same space as YDOTI in SLSOIBT.
C-----------------------------------------------------------------------
      CALL SSTODI (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   YDOTI, RWORK(LSAVF), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), RES, ADDA, JAC, SPJIBT, SSLSBT )
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 400, 550), KGO
C
C KGO = 1:success; 2:error test failure; 3:convergence failure;
C       4:RES ordered return; 5:RES returned error.
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  see if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSOIBT.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
  420 ISTATE = 2
      IF ( KFLAG .EQ. -3 )  ISTATE = 3
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSOIBT- At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSOIBT- At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 590
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSOIBT- At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 590
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSOIBT- At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = 'error test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 570
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSOIBT- At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 570
C IRES = 3 returned by RES, despite retries by SSTODI.------------------
 550  MSG = 'SLSOIBT- At T (=R1) residual routine returned     '
      CALL XERRWV (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      error IRES = 3 repeatedly.        '
      CALL XERRWV (MSG, 40, 206, 0, 0, 0, 0, 1, TN, 0.0E0)
      ISTATE = -7
      GO TO 590
C SAIGBT failed because a diagonal block of A matrix was singular. -----
 560  IER = -IER
      MSG='SLSOIBT- Attempt to initialize dy/dt failed:  Matrix A has a'
      CALL XERRWV (MSG, 60, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      singular diagonal block, block no. = (I1)   '
      CALL XERRWV (MSG, 50, 207, 0, 1, IER, 0, 0, 0.0E0, 0.0E0)
      ISTATE = -8
      RETURN
C SAIGBT failed because RES set IRES to 2 or 3. ------------------------
 565  MSG = 'SLSOIBT- Attempt to initialize dy/dt failed       '
      CALL XERRWV (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      because residual routine set its error flag '
      CALL XERRWV (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      to IRES = (I1)'
      CALL XERRWV (MSG, 20, 208, 0, 1, IER, 0, 0, 0.0E0, 0.0E0)
      ISTATE = -8
      RETURN
C Compute IMXER if relevant. -------------------------------------------
 570  BIG = 0.0E0
      IMXER = 1
      DO 575 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 575
        BIG = SIZE
        IMXER = I
 575    CONTINUE
      IWORK(16) = IMXER
C Compute residual if relevant. ----------------------------------------
 580  LYD0 = LYH + NYH
      DO 585 I = 1,N
         RWORK(I+LSAVF-1) = RWORK(I+LYD0-1)/H
 585     Y(I) = RWORK(I+LYH-1)
      IRES = 1
      CALL RES (NEQ, TN, Y, RWORK(LSAVF), YDOTI, IRES)
      NFE = NFE + 1
      IF (IRES .LE. 1)  GO TO 595
      MSG = 'SLSOIBT- Residual routine set its flag IRES       '
      CALL XERRWV (MSG, 50, 210, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      to (I1) when called for final output.       '
      CALL XERRWV (MSG, 50, 210, 0, 1, IRES, 0, 0, 0.0E0, 0.0E0)
      GO TO 595
C Set Y vector, T, and optional outputs. -------------------------------
 590  DO 592 I = 1,N
 592    Y(I) = RWORK(I+LYH-1)
 595  T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSOIBT- ISTATE (=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSOIBT- ITASK (=I1) illegal. '
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSOIBT- ISTATE.gt.1 but SLSOIBT not initialized. '
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSOIBT- NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSOIBT- ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSOIBT- ITOL (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSOIBT- IOPT (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSOIBT- MF (=I1) illegal.    '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 609  MSG = 'SLSOIBT- MB (=I1) or NB (=I2) illegal.  '
      CALL XERRWV (MSG, 40, 9, 0, 2, MB, NB, 0, 0.0E0, 0.0E0)
      GO TO 700
 610  MSG = 'SLSOIBT- NB (=I1) .lt. 4 illegal.       '
      CALL XERRWV (MSG, 40, 10, 0, 1, NB, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSOIBT- MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSOIBT- MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSOIBT- MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSOIBT- TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSOIBT- HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSOIBT- HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG='SLSOIBT- RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG='SLSOIBT- IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSOIBT- RTOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSOIBT- ATOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSOIBT- EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSOIBT- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSOIBT- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSOIBT- At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSOIBT- Trouble in SINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSOIBT- Run aborted.. apparent infinite loop.    '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSOIBT ---------------------
      END
*DECK SLSODIS
      SUBROUTINE SLSODIS (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      REAL Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 18 November 2003 version of
C SLSODIS: Livermore Solver for Ordinary Differential equations
C          (Implicit form) with general Sparse Jacobian matrices.
C
C This version is in single precision.
C
C SLSODIS solves the initial value problem for linearly implicit
C systems of first order ODEs,
C     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
C or, in component form,
C     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
C        i,1      1                     i,NEQ      NEQ
C
C      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
C           i      1   2       NEQ
C
C If A is singular, this is a differential-algebraic system.
C
C SLSODIS is a variant version of the SLSODI package, and is intended
C for stiff problems in which the matrix A and the Jacobian matrix
C d(g - A*s)/dy have arbitrary sparse structures.
C
C Authors:       Alan C. Hindmarsh
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Sheila Balsdon
C                Zycor, Inc.
C                Austin, TX 78741
C-----------------------------------------------------------------------
C References:
C 1.  M. K. Seager and S. Balsdon,  LSODIS, A Sparse Implicit
C     ODE Solver, in Proceedings of the IMACS 10th World Congress,
C     Montreal, August 8-13, 1982.
C
C 2.  Alan C. Hindmarsh,  LSODE and LSODI, Two New Initial Value
C     Ordinary Differential Equation Solvers,
C     ACM-SIGNUM Newsletter, vol. 15, no. 4 (1980), pp. 10-11.
C
C 3.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
C     Yale Sparse Matrix Package: I. The Symmetric Codes,
C     Int. J. Num. Meth. Eng., vol. 18 (1982), pp. 1145-1151.
C
C 4.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
C     Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
C     Research Report No. 114, Dept. of Computer Sciences, Yale
C     University, 1977.
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the SLSODIS package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First, provide a subroutine of the form:
C                SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C                DIMENSION Y(*), S(*), R(*)
C which computes the residual function
C      r = g(t,y)  -  A(t,y) * s ,
C as a function of t and the vectors y and s.  (s is an internally
C generated approximation to dy/dt.)  The arrays Y and S are inputs
C to the RES routine and should not be altered.  The residual
C vector is to be stored in the array R.  The argument IRES should be
C ignored for casual use of SLSODIS.  (For uses of IRES, see the
C paragraph on RES in the full description below.)
C
C B. SLSODIS must deal internally with the matrices A and dr/dy, where
C r is the residual function defined above.  SLSODIS generates a linear
C combination of these two matrices in sparse form.
C      The matrix structure is communicated by a method flag, MF:
C         MF =  21 or  22     when the user provides the structures of
C                             matrix A and dr/dy,
C         MF = 121 or 222     when the user does not provide structure
C                             information, and
C         MF = 321 or 422     when the user provides the structure
C                             of matrix A.
C
C C. You must also provide a subroutine of the form:
C                SUBROUTINE ADDA (NEQ, T, Y, J, IAN, JAN, P)
C                DIMENSION Y(*), IAN(*), JAN(*), P(*)
C which adds the matrix A = A(t,y) to the contents of the array P.
C NEQ, T, Y, and J are input arguments and should not be altered.
C This routine should add the J-th column of matrix A to the array
C P (of length NEQ).  I.e. add A(i,J) to P(i) for all relevant
C values of i.  The arguments IAN and JAN should be ignored for normal
C situations.  SLSODIS will call the ADDA routine with J = 1,2,...,NEQ.
C
C D. For the sake of efficiency, you are encouraged to supply the
C Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
C (s = a fixed vector) as above.  If dr/dy is being supplied,
C use MF = 21, 121, or 321, and provide a subroutine of the form:
C               SUBROUTINE JAC (NEQ, T, Y, S, J, IAN, JAN, PDJ)
C               DIMENSION Y(*), S(*), IAN(*), JAN(*), PDJ(*)
C which computes dr/dy as a function of t, y, and s.  Here NEQ, T, Y, S,
C and J are input arguments, and the JAC routine is to load the array
C PDJ (of length NEQ) with the J-th column of dr/dy.  I.e. load PDJ(i)
C with dr(i)/dy(J) for all relevant values of i.  The arguments IAN and
C JAN should be ignored for normal situations.  SLSODIS will call the
C JAC routine with J = 1,2,...,NEQ.
C      Only nonzero elements need be loaded.  A crude approximation
C to dr/dy, possibly with fewer nonzero elememts, will suffice.
C Note that if A is independent of y (or this dependence
C is weak enough to be ignored) then JAC is to compute dg/dy.
C      If it is not feasible to provide a JAC routine, use
C MF = 22, 222, or 422 and SLSODIS will compute an approximate
C Jacobian internally by difference quotients.
C
C E. Next decide whether or not to provide the initial value of the
C derivative vector dy/dt.  If the initial value of A(t,y) is
C nonsingular (and not too ill-conditioned), you may let SLSODIS compute
C this vector (ISTATE = 0).  (SLSODIS will solve the system A*s = g for
C s, with initial values of A and g.)  If A(t,y) is initially
C singular, then the system is a differential-algebraic system, and
C you must make use of the particular form of the system to compute the
C initial values of y and dy/dt.  In that case, use ISTATE = 1 and
C load the initial value of dy/dt into the array YDOTI.
C The input array YDOTI and the initial Y array must be consistent with
C the equations A*dy/dt = g.  This implies that the initial residual
C r = g(t,y) - A(t,y)*YDOTI   must be approximately zero.
C
C F. Write a main program which calls Subroutine SLSODIS once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages by
C SLSODIS.  On the first call to SLSODIS, supply arguments as follows:
C RES    = name of user subroutine for residual function r.
C ADDA   = name of user subroutine for computing and adding A(t,y).
C JAC    = name of user subroutine for Jacobian matrix dr/dy
C          (MF = 121).  If not used, pass a dummy name.
C Note: The names for the RES and ADDA routines and (if used) the
C        JAC routine must be declared External in the calling program.
C NEQ    = number of scalar equations in the system.
C Y      = array of initial values, of length NEQ.
C YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          The estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
C          initial dy/dt is supplied, and 0 otherwise.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
C          where:
C          NNZ    = the number of nonzero elements in the sparse
C                   iteration matrix  P = A - con*dr/dy (con = scalar)
C                   (If NNZ is unknown, use an estimate of it.)
C          LENRAT = the real to integer wordlength ratio (usually 1 in
C                   single precision and 2 in double precision).
C          In any case, the required size of RWORK cannot generally
C          be predicted in advance for any value of MF, and the
C          value above is a rough estimate of a crude lower bound.
C          Some experimentation with this size may be necessary.
C          (When known, the correct required length is an optional
C          output, available in IWORK(17).)
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least 30.
C LIW    = declared length of IWORK (in user's dimension).
C MF     = method flag.  Standard values are:
C          121 for a user-supplied sparse Jacobian.
C          222 for an internally generated sparse Jacobian.
C          For other choices of MF, see the paragraph on MF in
C          the full description below.
C Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
C and possibly ATOL.
C
C G. The output from the first call, or any call, is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE =  2  if SLSODIS was successful, negative otherwise.
C          -1 means excess work done on this call (check all inputs).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 cannot occur in casual use.
C          -8 means SLSODIS was unable to compute the initial dy/dt.
C             in casual use, this means A(t,y) is initially singular.
C             Supply YDOTI and use ISTATE = 1 on the first call.
C          -9 means a fatal error return flag came from sparse solver
C             SCDRV by way of SPRJIS or SSOLSS.  Should never happen.
C
C          A return with ISTATE = -1, -4, or -5, may result from using
C          an inappropriate sparsity structure, one that is quite
C          different from the initial structure.  Consider calling
C          SLSODIS again with ISTATE = 3 to force the structure to be
C          reevaluated.  See the full description of ISTATE below.
C
C  If SLSODIS returns ISTATE = -1, -4  or -5, then the output of
C  SLSODIS also includes YDOTI = array containing residual vector
C  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
C
C H. To continue the integration after a successful return, simply
C reset TOUT and call SLSODIS again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is an example problem, with the coding needed
C for its solution by SLSODIS.  The problem comes from the partial
C differential equation (the Burgers equation)
C   du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
C on -1 .le. x .le. 1.  The boundary conditions are periodic:
C   u(-1,t) = u(1,t)  and  du/dx(-1,t) = du/dx(1,t)
C The initial profile is a square wave,
C   u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.
C The PDE is discretized in x by a simplified Galerkin method,
C using piecewise linear basis functions, on a grid of 40 intervals.
C The result is a system A * dy/dt = g(y), of size NEQ = 40,
C where y(i) is the approximation to u at x = x(i), with
C x(i) = -1 + (i-1)*delx, delx = 2/NEQ = .05.
C The individual equations in the system are (in order):
C  (1/6)dy(NEQ)/dt+(4/6)dy(1)/dt+(1/6)dy(2)/dt
C       = r4d*(y(NEQ)**2-y(2)**2)+eodsq*(y(2)-2*y(1)+y(NEQ))
C for i = 2,3,...,nm1,
C  (1/6)dy(i-1)/dt+(4/6)dy(i)/dt+(1/6)dy(i+1)/dt
C       = r4d*(y(i-1)**2-y(i+1)**2)+eodsq*(y(i+1)-2*y(i)+y(i-1))
C and finally
C  (1/6)dy(nm1)/dt+(4/6)dy(NEQ)/dt+(1/6)dy(1)/dt
C       = r4d*(y(nm1)**2-y(1)**2)+eodsq*(y(1)-2*y(NEQ)+y(nm1))
C where r4d = 1/(4*delx), eodsq = eta/delx**2 and nm1 = NEQ-1.
C The following coding solves the problem with MF = 121, with output
C of solution statistics at t = .1, .2, .3, and .4, and of the
C solution vector at t = .4.  Optional outputs (run statistics) are
C also printed.
C
C     EXTERNAL RESID, ADDASP, JACSP
C     REAL ATOL, RTOL, RW, T, TOUT, Y, YDOTI, R4D, EODSQ, DELX
C     DIMENSION Y(40), YDOTI(40), RW(1409), IW(30)
C     COMMON /TEST1/ R4D, EODSQ, NM1
C     DATA ITOL/1/, RTOL/1.0E-3/, ATOL/1.0E-3/, ITASK/1/, IOPT/0/
C     DATA NEQ/40/, LRW/1409/, LIW/30/, MF/121/
C
C     DELX = 2.0/NEQ
C     R4D = 0.25/DELX
C     EODSQ = 0.05/DELX**2
C     NM1 = NEQ - 1
C     DO 10 I = 1,NEQ
C 10    Y(I) = 0.0
C     Y(11) = 0.5
C     DO 15 I = 12,30
C 15    Y(I) = 1.0
C     Y(31) = 0.5
C     T = 0.0
C     TOUT = 0.1
C     ISTATE = 0
C     DO 30 IO = 1,4
C       CALL SLSODIS (RESID, ADDASP, JACSP, NEQ, Y, YDOTI, T, TOUT,
C    1    ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RW, LRW, IW, LIW, MF)
C       WRITE(6,20) T,IW(11),RW(11)
C 20    FORMAT(' At t =',F5.2,'   No. steps =',I4,
C    1    '    Last step =',E12.4)
C       IF (ISTATE .NE. 2) GO TO 90
C       TOUT = TOUT + 0.1
C 30  CONTINUE
C     WRITE (6,40) (Y(I),I=1,NEQ)
C 40  FORMAT(/' Final solution values..'/8(5E12.4/))
C     WRITE(6,50) IW(17),IW(18),IW(11),IW(12),IW(13)
C     NNZLU = IW(25) + IW(26) + NEQ
C     WRITE(6,60) IW(19),NNZLU
C 50  FORMAT(/' Required RW size =',I5,'   IW size =',I4/
C    1  ' No. steps =',I4,'   No. r-s =',I4,'   No. J-s =',i4)
C 60  FORMAT(' No. of nonzeros in P matrix =',I4,
C    1  '   No. of nonzeros in LU =',I4)
C     STOP
C 90  WRITE (6,95) ISTATE
C 95  FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE GFUN (N, T, Y, G)
C     REAL T, Y, G, R4D, EODSQ
C     DIMENSION G(N), Y(N)
C     COMMON /TEST1/ R4D, EODSQ, NM1
C     G(1) = R4D*(Y(N)**2-Y(2)**2) + EODSQ*(Y(2)-2.0*Y(1)+Y(N))
C     DO 10 I = 2,NM1
C       G(I) = R4D*(Y(I-1)**2 - Y(I+1)**2)
C    1        + EODSQ*(Y(I+1) - 2.0*Y(I) + Y(I-1))
C 10    CONTINUE
C     G(N) = R4D*(Y(NM1)**2-Y(1)**2) + EODSQ*(Y(1)-2.0*Y(N)+Y(NM1))
C     RETURN
C     END
C
C     SUBROUTINE RESID (N, T, Y, S, R, IRES)
C     REAL T, Y, S, R, R4D, EODSQ
C     DIMENSION Y(N), S(N), R(N)
C     COMMON /TEST1/ R4D, EODSQ, NM1
C     CALL GFUN (N, T, Y, R)
C     R(1) = R(1) - (S(N) + 4.0*S(1) + S(2))/6.0
C     DO 10 I = 2,NM1
C 10    R(I) = R(I) - (S(I-1) + 4.0*S(I) + S(I+1))/6.0
C     R(N) = R(N) - (S(NM1) + 4.0*S(N) + S(1))/6.0
C     RETURN
C     END
C
C     SUBROUTINE ADDASP (N, T, Y, J, IP, JP, P)
C     REAL T, Y, P
C     DIMENSION Y(N), IP(*), JP(*), P(N)
C     JM1 = J - 1
C     JP1 = J + 1
C     IF (J .EQ. N) JP1 = 1
C     IF (J .EQ. 1) JM1 = N
C     P(J) = P(J) + (2.0/3.0)
C     P(JP1) = P(JP1) + (1.0/6.0)
C     P(JM1) = P(JM1) + (1.0/6.0)
C     RETURN
C     END
C
C     SUBROUTINE JACSP (N, T, Y, S, J, IP, JP, PDJ)
C     REAL T, Y, S, PDJ, R4D, EODSQ
C     DIMENSION Y(N), S(N), IP(*), JP(*), PDJ(N)
C     COMMON /TEST1/ R4D, EODSQ, NM1
C     JM1 = J - 1
C     JP1 = J + 1
C     IF (J .EQ. 1) JM1 = N
C     IF (J .EQ. N) JP1 = 1
C     PDJ(JM1) = -2.0*R4D*Y(J) + EODSQ
C     PDJ(J) = -2.0*EODSQ
C     PDJ(JP1) = 2.0*R4D*Y(J) + EODSQ
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C At t = 0.10   No. steps =  15    Last step =  1.6863e-02
C At t = 0.20   No. steps =  19    Last step =  2.4101e-02
C At t = 0.30   No. steps =  22    Last step =  4.3143e-02
C At t = 0.40   No. steps =  24    Last step =  5.7819e-02
C
C Final solution values..
C  1.8371e-02  1.3578e-02  1.5864e-02  2.3805e-02  3.7245e-02
C  5.6630e-02  8.2538e-02  1.1538e-01  1.5522e-01  2.0172e-01
C  2.5414e-01  3.1150e-01  3.7259e-01  4.3608e-01  5.0060e-01
C  5.6482e-01  6.2751e-01  6.8758e-01  7.4415e-01  7.9646e-01
C  8.4363e-01  8.8462e-01  9.1853e-01  9.4500e-01  9.6433e-01
C  9.7730e-01  9.8464e-01  9.8645e-01  9.8138e-01  9.6584e-01
C  9.3336e-01  8.7497e-01  7.8213e-01  6.5315e-01  4.9997e-01
C  3.4672e-01  2.1758e-01  1.2461e-01  6.6208e-02  3.3784e-02
C
C Required RW size = 1409   IW size =  30
C No. steps =  24   No. r-s =  33   No. J-s =   8
C No. of nonzeros in P matrix = 120   No. of nonzeros in LU = 194
C
C-----------------------------------------------------------------------
C Full Description of User Interface to SLSODIS.
C
C The user interface to SLSODIS consists of the following parts.
C
C 1.   The call sequence to Subroutine SLSODIS, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the SLSODIS package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the SLSODIS package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
C     IOPT, LRW, LIW, MF,
C and those used for both input and output are
C     Y, T, ISTATE, YDOTI.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine SLSODIS to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C RES    = the name of the user-supplied subroutine which supplies
C          the residual vector for the ODE system, defined by
C            r = g(t,y) - A(t,y) * s
C          as a function of the scalar t and the vectors
C          s and y (s approximates dy/dt).  This subroutine
C          is to have the form
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DIMENSION Y(*), S(*), R(*)
C          where NEQ, T, Y, S, and IRES are input, and R and
C          IRES are output.  Y, S, and R are arrays of length NEQ.
C             On input, IRES indicates how SLSODIS will use the
C          returned array R, as follows:
C             IRES = 1  means that SLSODIS needs the full residual,
C                       r = g - A*s, exactly.
C             IRES = -1 means that SLSODIS is using R only to compute
C                       the Jacobian dr/dy by difference quotients.
C          The RES routine can ignore IRES, or it can omit some terms
C          if IRES = -1.  If A does not depend on y, then RES can
C          just return R = g when IRES = -1.  If g - A*s contains other
C          additive terms that are independent of y, these can also be
C          dropped, if done consistently, when IRES = -1.
C             The subroutine should set the flag IRES if it
C          encounters a halt condition or illegal input.
C          Otherwise, it should not reset IRES.  On output,
C             IRES = 1 or -1 represents a normal return, and
C          SLSODIS continues integrating the ODE.  Leave IRES
C          unchanged from its input value.
C             IRES = 2 tells SLSODIS to immediately return control
C          to the calling program, with ISTATE = 3.  This lets
C          the calling program change parameters of the problem
C          if necessary.
C             IRES = 3 represents an error condition (for example, an
C          illegal value of y).  SLSODIS tries to integrate the system
C          without getting IRES = 3 from RES.  If it cannot, SLSODIS
C          returns with ISTATE = -7 or -1.
C             On a return with ISTATE = 3, -1, or -7, the values
C          of T and Y returned correspond to the last point reached
C          successfully without getting the flag IRES = 2 or 3.
C             The flag values IRES = 2 and 3 should not be used to
C          handle switches or root-stop conditions.  This is better
C          done by calling SLSODIS in a one-step mode and checking the
C          stopping function for a sign change at each step.
C             If quantities computed in the RES routine are needed
C          externally to SLSODIS, an extra call to RES should be made
C          for this purpose, for consistent and accurate results.
C          To get the current dy/dt for the S argument, use SINTDY.
C             RES must be declared External in the calling
C          program.  See note below for more about RES.
C
C ADDA   = the name of the user-supplied subroutine which adds the
C          matrix A = A(t,y) to another matrix stored in sparse form.
C          This subroutine is to have the form
C               SUBROUTINE ADDA (NEQ, T, Y, J, IAN, JAN, P)
C               DIMENSION Y(*), IAN(*), JAN(*), P(*)
C          where NEQ, T, Y, J, IAN, JAN, and P  are input.  This routine
C          should add the J-th column of matrix A to the array P, of
C          length NEQ.  Thus a(i,J) is to be added to P(i) for all
C          relevant values of i.  Here T and Y have the same meaning as
C          in Subroutine RES, and J is a column index (1 to NEQ).
C          IAN and JAN are undefined in calls to ADDA for structure
C          determination (MOSS .ne. 0).  Otherwise, IAN and JAN are
C          structure descriptors, as defined under optional outputs
C          below, and so can be used to determine the relevant row
C          indices i, if desired.
C               Calls to ADDA are made with J = 1,...,NEQ, in that
C          order.  ADDA must not alter its input arguments.
C               ADDA must be declared External in the calling program.
C          See note below for more information about ADDA.
C
C JAC    = the name of the user-supplied subroutine which supplies
C          the Jacobian matrix, dr/dy, where r = g - A*s.  JAC is
C          required if MITER = 1, or MOSS = 1 or 3.  Otherwise a dummy
C          name can be passed.  This subroutine is to have the form
C               SUBROUTINE JAC ( NEQ, T, Y, S, J, IAN, JAN, PDJ)
C               DIMENSION Y(*), S(*), IAN(*), JAN(*), PDJ(*)
C         where NEQ, T, Y, S, J, IAN, and JAN are input.  The
C         array PDJ, of length NEQ, is to be loaded with column J
C         of the Jacobian on output.  Thus dr(i)/dy(J) is to be
C         loaded into PDJ(i) for all relevant values of i.
C         Here T, Y, and S have the same meaning as in Subroutine RES,
C         and J is a column index (1 to NEQ).  IAN and JAN
C         are undefined in calls to JAC for structure determination
C         (MOSS .ne. 0).  Otherwise, IAN and JAN are structure
C         descriptors, as defined under optional outputs below, and
C         so can be used to determine the relevant row indices i, if
C         desired.
C              JAC need not provide dr/dy exactly.  A crude
C         approximation (possibly with greater sparsity) will do.
C              In any case, PDJ is preset to zero by the solver,
C         so that only the nonzero elements need be loaded by JAC.
C         Calls to JAC are made with J = 1,...,NEQ, in that order, and
C         each such set of calls is preceded by a call to RES with the
C         same arguments NEQ, T, Y, S, and IRES.  Thus to gain some
C         efficiency intermediate quantities shared by both calculations
C         may be saved in a user Common block by RES and not recomputed
C         by JAC, if desired.  JAC must not alter its input arguments.
C              JAC must be declared External in the calling program.
C              See note below for more about JAC.
C
C    Note on RES, ADDA, and JAC:
C          These subroutines may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in the subroutines) and/or Y has length
C          exceeding NEQ(1).  However, these subroutines should not
C          alter NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
C          See the descriptions of NEQ and Y below.
C
C NEQ    = the size of the system (number of first order ordinary
C          differential equations or scalar algebraic equations).
C          Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in RES, ADDA, or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The SLSODIS package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to RES, ADDA, and JAC.  Hence, if it is an array,
C          locations NEQ(2),... may be used to store other integer data
C          and pass it to RES, ADDA, or JAC.  Each such subroutine
C          must include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 0 or 1), and only for output on other
C          calls.  On the first call, Y must contain the vector of
C          initial values.  On output, Y contains the computed solution
C          vector, evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to RES,
C          ADDA, and JAC.  Hence its length may exceed NEQ,
C          and locations Y(NEQ+1),... may be used to store other real
C          data and pass it to RES, ADDA, or JAC.  (The SLSODIS
C          package accesses only Y(1),...,Y(NEQ). )
C
C YDOTI  = a real array for the initial value of the vector
C          dy/dt and for work space, of dimension at least NEQ.
C
C          On input:
C            If ISTATE = 0 then SLSODIS will compute the initial value
C          of dy/dt, if A is nonsingular.  Thus YDOTI will
C          serve only as work space and may have any value.
C            If ISTATE = 1 then YDOTI must contain the initial value
C          of dy/dt.
C            If ISTATE = 2 or 3 (continuation calls) then YDOTI
C          may have any value.
C            Note: If the initial value of A is singular, then
C          SLSODIS cannot compute the initial value of dy/dt, so
C          it must be provided in YDOTI, with ISTATE = 1.
C
C          On output, when SLSODIS terminates abnormally with ISTATE =
C          -1, -4, or -5, YDOTI will contain the residual
C          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
C          its initial value, and YDOTI is supplied with ISTATE = 1,
C          there may have been an incorrect input value of
C          YDOTI = dy/dt, or the problem (as given to SLSODIS)
C          may not have a solution.
C
C          If desired, the YDOTI array may be used for other
C          purposes between calls to the solver.
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 0 or 1), TOUT may be
C          equal to T for one call, then should .ne. T for the next
C          call.  For the initial T, an input value of TOUT .ne. T is
C          used in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          0  means this is the first call for the problem, and
C             SLSODIS is to compute the initial value of dy/dt
C             (while doing other initializations).  See note below.
C          1  means this is the first call for the problem, and
C             the initial value of dy/dt has been supplied in
C             YDOTI (SLSODIS will do other initializations).
C             See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C             the conditional inputs IA, JA, IC, and JC,
C             and any of the optional inputs except H0.
C             A call with ISTATE = 3 will cause the sparsity
C             structure of the problem to be recomputed.
C             (Structure information is reread from IA and JA if
C             MOSS = 0, 3, or 4 and from IC and JC if MOSS = 0).
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 0 or 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           0 or 1  means nothing was done; TOUT = T and
C              ISTATE = 0 or 1 on input.
C           2  means that the integration was performed successfully.
C           3  means that the user-supplied Subroutine RES signalled
C              SLSODIS to halt the integration and return (IRES = 2).
C              Integration as far as T was achieved with no occurrence
C              of IRES = 2, but this flag was set on attempting the
C              next step.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i) = 0.0)
C              was requested on a variable which has now vanished.
C              the integration was successful as far as T.
C          -7  means that the user-supplied Subroutine RES set
C              its error flag (IRES = 3) despite repeated tries by
C              SLSODIS to avoid that condition.
C          -8  means that ISTATE was 0 on input but SLSODIS was unable
C              to compute the initial value of dy/dt.  See the
C              printed message for details.
C          -9  means a fatal error return flag came from the sparse
C              solver SCDRV by way of SPRJIS or SSOLSS (numerical
C              factorization or backsolve).  This should never happen.
C              The integration was successful as far as T.
C
C          Note: An error return with ISTATE = -1, -4, or -5
C          may mean that the sparsity structure of the
C          problem has changed significantly since it was last
C          determined (or input).  In that case, one can attempt to
C          complete the integration by setting ISTATE = 3 on the next
C          call, so that a new structure determination is done.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          similarly, ISTATE (= 3) need not be reset if RES told
C          SLSODIS to return because the calling program must change
C          the parameters of the problem.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a work array used for a mixture of real (single precision)
C          and integer work space.
C          The length of RWORK (in real words) must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LWM = 2*NNZ + 2*NEQ + (NNZ+9*NEQ)/LENRAT   if MITER = 1,
C          LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT  if MITER = 2.
C          in the above formulas,
C          NNZ    = number of nonzero elements in the iteration matrix
C                   P = A - con*J  (con is a constant and J is the
C                   Jacobian matrix dr/dy).
C          LENRAT = the real to integer wordlength ratio (usually 1 in
C                   single precision and 2 in double precision).
C          (See the MF description for METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          the minimum length of RWORK is:
C             20 + 16*NEQ + LWM  for MF = 11, 111, 311, 12, 212, 412,
C             20 +  9*NEQ + LWM  for MF = 21, 121, 321, 22, 222, 422.
C          The above formula for LWM is only a crude lower bound.
C          The required length of RWORK cannot be readily predicted
C          in general, as it depends on the sparsity structure
C          of the problem.  Some experimentation may be necessary.
C
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C            32 + 2*NEQ + NZA + NZC   for MOSS = 0,
C            30                       for MOSS = 1 or 2,
C            31 + NEQ + NZA           for MOSS = 3 or 4.
C          (NZA is the number of nonzero elements in matrix A, and
C          NZC is the number of nonzero elements in dr/dy.)
C
C          In SLSODIS, IWORK is used for conditional and
C          optional inputs and optional outputs.
C
C          The following two blocks of words in IWORK are conditional
C          inputs, required if MOSS = 0, 3, or 4, but not otherwise
C          (see the description of MF for MOSS).
C            IWORK(30+j) = IA(j)     (j=1,...,NEQ+1)
C            IWORK(31+NEQ+k) = JA(k) (k=1,...,NZA)
C          The two arrays IA and JA describe the sparsity structure
C          to be assumed for the matrix A.  JA contains the row
C          indices where nonzero elements occur, reading in columnwise
C          order, and IA contains the starting locations in JA of the
C          descriptions of columns 1,...,NEQ, in that order, with
C          IA(1) = 1.  Thus, for each column index j = 1,...,NEQ, the
C          values of the row index i in column j where a nonzero
C          element may occur are given by
C            i = JA(k),  where   IA(j) .le. k .lt. IA(j+1).
C          If NZA is the total number of nonzero locations assumed,
C          then the length of the JA array is NZA, and IA(NEQ+1) must
C          be NZA + 1.  Duplicate entries are not allowed.
C          The following additional blocks of words are required
C          if MOSS = 0, but not otherwise.  If LC = 31 + NEQ + NZA, then
C            IWORK(LC+j) = IC(j)       (j=1,...,NEQ+1), and
C            IWORK(LC+NEQ+1+k) = JC(k) (k=1,...,NZC)
C          The two arrays IC and JC describe the sparsity
C          structure to be assumed for the Jacobian matrix dr/dy.
C          They are used in the same manner as the above IA and JA
C          arrays.  If NZC is the number of nonzero locations
C          assumed, then the length of the JC array is NZC, and
C          IC(NEQ+1) must be NZC + 1.  Duplicate entries are not
C          allowed.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to SLSODIS
C for the same problem, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside SLSODIS between calls, if
C desired (but not for use by RES, ADDA, or JAC).
C
C MF     = the method flag.  Used only for input.
C          MF has three decimal digits-- MOSS, METH, and MITER.
C          For standard options:
C             MF = 100*MOSS + 10*METH + MITER.
C          MOSS indicates the method to be used to obtain the sparsity
C          structure of the Jacobian matrix:
C            MOSS = 0 means the user has supplied IA, JA, IC, and JC
C                     (see descriptions under IWORK above).
C            MOSS = 1 means the user has supplied JAC (see below) and
C                     the structure will be obtained from NEQ initial
C                     calls to JAC and NEQ initial calls to ADDA.
C            MOSS = 2 means the structure will be obtained from NEQ+1
C                     initial calls to RES and NEQ initial calls to ADDA
C            MOSS = 3 like MOSS = 1, except user has supplied IA and JA.
C            MOSS = 4 like MOSS = 2, except user has supplied IA and JA.
C          METH indicates the basic linear multistep method:
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on Backward
C                     Differentiation Formulas (BDFs).
C              The BDF method is strongly preferred for stiff problems,
C            while the Adams method is preferred when the problem is
C            not stiff.  If the matrix A(t,y) is nonsingular,
C            stiffness here can be taken to mean that of the explicit
C            ODE system dy/dt = A-inverse * g.  If A is singular,
C            the concept of stiffness is not well defined.
C              If you do not know whether the problem is stiff, we
C            recommend using METH = 2.  If it is stiff, the advantage
C            of METH = 2 over METH = 1 will be great, while if it is
C            not stiff, the advantage of METH = 1 will be slight.
C            If maximum efficiency is important, some experimentation
C            with METH may be necessary.
C          MITER indicates the corrector iteration method:
C            MITER = 1 means chord iteration with a user-supplied
C                      sparse Jacobian, given by Subroutine JAC.
C            MITER = 2 means chord iteration with an internally
C                      generated (difference quotient) sparse
C                      Jacobian (using NGP extra calls to RES per
C                      dr/dy value, where NGP is an optional
C                      output described below.)
C            If MITER = 1 or MOSS = 1 or 3 the user must supply a
C            Subroutine JAC (the name is arbitrary) as described above
C            under JAC.  Otherwise, a dummy argument can be used.
C
C          The standard choices for MF are:
C            MF = 21 or 22 for a stiff problem with IA/JA and IC/JC
C                 supplied,
C            MF = 121 for a stiff problem with JAC supplied, but not
C                 IA/JA or IC/JC,
C            MF = 222 for a stiff problem with neither IA/JA, IC/JC/,
C                 nor JAC supplied,
C            MF = 321 for a stiff problem with IA/JA and JAC supplied,
C                 but not IC/JC,
C            MF = 422 for a stiff problem with IA/JA supplied, but not
C                 IC/JC or JAC.
C
C          The sparseness structure can be changed during the problem
C          by making a call to SLSODIS with ISTATE = 3.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from SLSODIS, the variables listed
C below are quantities related to the performance of SLSODIS
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from SLSODIS, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
C input) or -8, they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NRE     IWORK(12) the number of residual evaluations (RES calls)
C                   for the problem so far, excluding those for
C                   structure determination (MOSS = 2 or 4).
C
C NJE     IWORK(13) the number of Jacobian evaluations (each involving
C                   an evaluation of A and dr/dy) for the problem so
C                   far, excluding those for structure determination
C                   (MOSS = 1 or 3).  This equals the number of calls
C                   to ADDA and (if MITER = 1) JAC.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NNZ     IWORK(19) the number of nonzero elements in the iteration
C                   matrix  P = A - con*J  (con is a constant and
C                   J is the Jacobian matrix dr/dy).
C
C NGP     IWORK(20) the number of groups of column indices, used in
C                   difference quotient Jacobian aproximations if
C                   MITER = 2.  This is also the number of extra RES
C                   evaluations needed for each Jacobian evaluation.
C
C NLU     IWORK(21) the number of sparse LU decompositions for the
C                   problem so far. (Excludes the LU decomposition
C                   necessary when ISTATE = 0.)
C
C LYH     IWORK(22) the base address in RWORK of the history array YH,
C                   described below in this list.
C
C IPIAN   IWORK(23) the base address of the structure descriptor array
C                   IAN, described below in this list.
C
C IPJAN   IWORK(24) the base address of the structure descriptor array
C                   JAN, described below in this list.
C
C NZL     IWORK(25) the number of nonzero elements in the strict lower
C                   triangle of the LU factorization used in the chord
C                   iteration.
C
C NZU     IWORK(26) the number of nonzero elements in the strict upper
C                   triangle of the LU factorization used in the chord
C                   iteration.  The total number of nonzeros in the
C                   factorization is therefore NZL + NZU + NEQ.
C
C The following four arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address, and its description.
C For YH and ACOR, the base addresses are in RWORK (a real array).
C The integer arrays IAN and JAN are to be obtained by declaring an
C integer array IWK and identifying IWK(1) with RWORK(21), using either
C an equivalence statement or a subroutine call.  Then the base
C addresses IPIAN (of IAN) and IPJAN (of JAN) in IWK are to be obtained
C as optional outputs IWORK(23) and IWORK(24), respectively.
C Thus IAN(1) is IWK(ipian), etc.
C
C Name    Base Address      Description
C
C IAN    IPIAN (in IWK)  structure descriptor array of size NEQ + 1.
C JAN    IPJAN (in IWK)  structure descriptor array of size NNZ.
C         (see above)    IAN and JAN together describe the sparsity
C                        structure of the iteration matrix
C                          P = A - con*J,  as used by SLSODIS.
C                        JAN contains the row indices of the nonzero
C                        locations, reading in columnwise order, and
C                        IAN contains the starting locations in JAN of
C                        the descriptions of columns 1,...,NEQ, in
C                        that order, with IAN(1) = 1.  Thus for each
C                        j = 1,...,NEQ, the row indices i of the
C                        nonzero locations in column j are
C                        i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).
C                        Note that IAN(NEQ+1) = NNZ + 1.
C YH      LYH            the Nordsieck history array, of size NYH by
C          (optional     (NQCUR + 1), where NYH is the initial value
C           output)      of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.  The base address LYH
C                        is another optional output, listed above.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output to
C                        represent the estimated local error in y on the
C                        last step.  This is the vector E in the
C                        description of the error control. It is defined
C                        only on a return from SLSODIS with ISTATE = 2.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with SLSODIS.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from SLSODIS, if
C                             The default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by SLSODIS.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL SSRCMS(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             SLSODIS (see Part 3 below).
C                             RSAV must be a real array of length 224
C                             or more, and ISAV must be an integer
C                             array of length 71 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                SSRCMS is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with SLSODIS.
C
C   CALL SINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from SLSODIS.
C
C The detailed instructions for using SINTDY are as follows.
C The form of the call is:
C
C   LYH = IWORK(22)
C   CALL SINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by SLSODIS).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(t), is already provided
C             by SLSODIS directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with SINTDY.
C LYH       = the base address of the history array YH, obtained
C             as an optional output as shown above.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If SLSODIS is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to SLSODIS, and
C   (2) the two internal Common blocks
C         /SLS001/  of length  255  (218 single precision words
C                      followed by 37 integer words),
C         /SLSS01/  of length  40  (6 single precision words
C                      followed by 34 integer words).
C
C If SLSODIS is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by SLSODIS is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last SLSODIS call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next SLSODIS call for that problem.  To save and restore the Common
C blocks, use Subroutines SSRCMS (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the SLSODIS package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) SEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the SLSODIS call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by SEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by SEWSET is passed to the SVNORM
C routine (see below), and also used by SLSODIS in the computation
C of the optional output IMXER, and the increments for difference
C quotient Jacobians.
C
C In the user-supplied version of SEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In SEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in SEWSET the statements:
C     REAL RLS
C     COMMON /SLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) SVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = SVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine SEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of SVNORM suitable for use in the error control in SLSODIS.
C None of the arguments should be altered by SVNORM.
C For example, a user-supplied SVNORM routine might:
C   -substitute a max-norm of (V(i)*w(I)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19820714  DATE WRITTEN
C 19830812  Major update, based on recent LSODI and LSODES revisions:
C           Upgraded MDI in ODRV package: operates on M + M-transpose.
C           Numerous revisions in use of work arrays;
C           use wordlength ratio LENRAT; added IPISP & LRAT to Common;
C           added optional outputs IPIAN/IPJAN;
C           Added routine CNTNZU; added NZL and NZU to /LSS001/;
C           changed ADJLR call logic; added optional outputs NZL & NZU;
C           revised counter initializations; revised PREPI stmt. nos.;
C           revised difference quotient increment;
C           eliminated block /LSI001/, using IERPJ flag;
C           revised STODI logic after PJAC return;
C           revised tuning of H change and step attempts in STODI;
C           corrections to main prologue and comments throughout.
C 19870320  Corrected jump on test of umax in CDRV routine.
C 20010125  Numerous revisions: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODI;
C           in STODI, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           converted arithmetic IF statements to logical IF statements;
C           numerous corrections to prologues and internal comments.
C 20020502  Minor revisions to main prologue.
C 20031021  Fixed address offset bugs in Subroutine SPREPI.
C 20031027  Changed 0. to 0.0E0 in Subroutine SPREPI.
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C 20031117  Changed internal names NRE, LSAVR to NFE, LSAVF resp.
C
C-----------------------------------------------------------------------
C Other routines in the SLSODIS package.
C
C In addition to Subroutine SLSODIS, the SLSODIS package includes the
C following subroutines and function routines:
C  SIPREPI  acts as an interface between SLSODIS and SPREPI, and also
C           does adjusting of work space pointers and work arrays.
C  SPREPI   is called by SIPREPI to compute sparsity and do sparse
C           matrix preprocessing.
C  SAINVGS  computes the initial value of the vector
C             dy/dt = A-inverse * g
C  ADJLR    adjusts the length of required sparse matrix work space.
C           It is called by SPREPI.
C  CNTNZU   is called by SPREPI and counts the nonzero elements in the
C           strict upper triangle of P + P-transpose.
C  JGROUP   is called by SPREPI to compute groups of Jacobian column
C           indices for use when MITER = 2.
C  SINTDY   computes an interpolated value of the y vector at t = TOUT.
C  SSTODI   is the core integrator, which does one step of the
C           integration and the associated error control.
C  SCFODE   sets all method coefficients and test constants.
C  SPRJIS   computes and preprocesses the Jacobian matrix J = dr/dy
C           and the Newton iteration matrix P = A - h*l0*J.
C  SSOLSS   manages solution of linear system in chord iteration.
C  SEWSET   sets the error weight vector EWT before each step.
C  SVNORM   computes the weighted RMS-norm of a vector.
C  SSRCMS   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  SODRV    constructs a reordering of the rows and columns of
C           a matrix by the minimum degree algorithm.  SODRV is a
C           driver routine which calls Subroutines MD, MDI, MDM,
C           MDP, MDU, and SSRO.  See Ref. 3 for details.  (The SODRV
C           module has been modified since Ref. 3, however.)
C  SCDRV    performs reordering, symbolic factorization, numerical
C           factorization, or linear system solution operations,
C           depending on a path argument IPATH.  SCDRV is a
C           driver routine which calls Subroutines SNROC, SNSFC,
C           SNNFC, SNNSC, and SNNTC.  See Ref. 4 for details.
C           SLSODIS uses SCDRV to solve linear systems in which the
C           coefficient matrix is  P = A - con*J, where A is the
C           matrix for the linear system A(t,y)*dy/dt = g(t,y),
C           con is a scalar, and J is an approximation to
C           the Jacobian dr/dy.  Because SCDRV deals with rowwise
C           sparsity descriptions, SCDRV works with P-transpose, not P.
C           SLSODIS also uses SCDRV to solve the linear system
C             A(t,y)*dy/dt = g(t,y)  for dy/dt when ISTATE = 0.
C           (For this, SCDRV works with A-transpose, not A.)
C  RUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWV, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWV is machine-dependent.
C Note:  SVNORM, RUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL SPRJIS, SSOLSS
      REAL RUMACH, SVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, I1, I2, IER, IGO, IFLAG, IMAX, IMUL, IMXER, IPFLAG,
     1   IPGO, IREM, IRES, J, KGO, LENRAT, LENYHT, LENIW, LENRW,
     2   LIA, LIC, LJA, LJC, LRTEM, LWTEM, LYD0, LYHD, LYHN, MF1,
     3   MORD, MXHNL0, MXSTP0, NCOLM
      REAL ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      REAL CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
      REAL ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE LENRAT, MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following two internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block SLS001 is declared in subroutines SLSODIS, SIPREPI,SPREPI,
C SINTDY, SSTODI, SPRJIS, and SSOLSS.  
C The block SLSS01 is declared in subroutines SLSODIS, SAINVGS,
C SIPREPI, SPREPI, SPRJIS, and SSOLSS.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /SLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /SLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH,
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C In the Data statement below, set LENRAT equal to the ratio of
C the wordlength for a real number to that for an integer.  Usually,
C LENRAT = 1 for single precision and 2 for double precision.  If the
C true ratio is not an integer, use the next smaller integer (.ge. 1),
C-----------------------------------------------------------------------
      DATA LENRAT/1/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropirately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 0 or 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 0 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .LE. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 0 or 1)
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C If ISTATE = 0 or 1, the final setting of work space pointers, the
C matrix preprocessing, and other initializations are done in Block C.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT, and
C MF.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .LE. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      MOSS = MF/100
      MF1 = MF - 100*MOSS
      METH = MF1/10
      MITER = MF1 - 10*METH
      IF (MOSS .LT. 0 .OR. MOSS .GT. 4) GO TO 608
      IF (MITER .EQ. 2 .AND. MOSS .EQ. 1) MOSS = MOSS + 1
      IF (MITER .EQ. 2 .AND. MOSS .EQ. 3) MOSS = MOSS + 1
      IF (MITER .EQ. 1 .AND. MOSS .EQ. 2) MOSS = MOSS - 1
      IF (MITER .EQ. 1 .AND. MOSS .EQ. 4) MOSS = MOSS - 1
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 1 .OR. MITER .GT. 2) GO TO 608
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .LE. 1) H0 = 0.0E0
      HMXI = 0.0E0
      HMIN = 0.0E0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .GT. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0E0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0E0) GO TO 615
      HMXI = 0.0E0
      IF (HMAX .GT. 0.0E0) HMXI = 1.0E0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0E0) GO TO 616
C Check RTOL and ATOL for legality. ------------------------------------
 60   RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 65 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0E0) GO TO 619
        IF (ATOLI .LT. 0.0E0) GO TO 620
 65     CONTINUE
C-----------------------------------------------------------------------
C Compute required work array lengths, as far as possible, and test
C these against LRW and LIW.  Then set tentative pointers for work
C arrays.  Pointers to RWORK/IWORK segments are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  WM, YH, SAVR, EWT, ACOR.
C The required length of the matrix work space WM is not yet known,
C and so a crude minimum value is used for the initial tests of LRW
C and LIW, and YH is temporarily stored as far to the right in RWORK
C as possible, to leave the maximum amount of space for WM for matrix
C preprocessing.  Thus if MOSS .ne. 2 or 4, some of the segments of
C RWORK are temporarily omitted, as they are not needed in the
C preprocessing.  These omitted segments are: ACOR if ISTATE = 1,
C EWT and ACOR if ISTATE = 3 and MOSS = 1, and SAVR, EWT, and ACOR if
C ISTATE = 3 and MOSS = 0.
C-----------------------------------------------------------------------
      LRAT = LENRAT
      IF (ISTATE .LE. 1) NYH = N
      IF (MITER .EQ. 1) LWMIN = 4*N + 10*N/LRAT
      IF (MITER .EQ. 2) LWMIN = 4*N + 11*N/LRAT
      LENYH = (MAXORD+1)*NYH
      LREST = LENYH + 3*N
      LENRW = 20 + LWMIN + LREST
      IWORK(17) = LENRW
      LENIW = 30
      IF (MOSS .NE. 1 .AND. MOSS .NE. 2) LENIW = LENIW + N + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
      LIA = 31
      IF (MOSS .NE. 1 .AND. MOSS .NE. 2)
     1   LENIW = LENIW + IWORK(LIA+N) - 1
      IWORK(18) = LENIW
      IF (LENIW .GT. LIW) GO TO 618
      LJA = LIA + N + 1
      LIA = MIN(LIA,LIW)
      LJA = MIN(LJA,LIW)
      LIC = LENIW + 1
      IF (MOSS .EQ. 0) LENIW = LENIW + N + 1
      IWORK(18) = LENIW
      IF (LENIW .GT. LIW) GO TO 618
      IF (MOSS .EQ. 0) LENIW =  LENIW + IWORK(LIC+N) - 1
      IWORK(18) = LENIW
      IF (LENIW .GT. LIW) GO TO 618
      LJC = LIC + N + 1
      LIC = MIN(LIC,LIW)
      LJC = MIN(LJC,LIW)
      LWM = 21
      IF (ISTATE .LE. 1) NQ = ISTATE
      NCOLM = MIN(NQ+1,MAXORD+2)
      LENYHM = NCOLM*NYH
      LENYHT = LENYHM
      IMUL = 2
      IF (ISTATE .EQ. 3) IMUL = MOSS
      IF (ISTATE .EQ. 3 .AND. MOSS .EQ. 3) IMUL = 1
      IF (MOSS .EQ. 2 .OR. MOSS .EQ. 4) IMUL = 3
      LRTEM = LENYHT + IMUL*N
      LWTEM = LRW - 20 - LRTEM
      LENWK = LWTEM
      LYHN = LWM + LWTEM
      LSAVF = LYHN + LENYHT
      LEWT = LSAVF + N
      LACOR = LEWT + N
      ISTATC = ISTATE
      IF (ISTATE .LE. 1) GO TO 100
C-----------------------------------------------------------------------
C ISTATE = 3.  Move YH to its new location.
C Note that only the part of YH needed for the next step, namely
C MIN(NQ+1,MAXORD+2) columns, is actually moved.
C A temporary error weight array EWT is loaded if MOSS = 2 or 4.
C Sparse matrix processing is done in SIPREPI/SPREPI.
C If MAXORD was reduced below NQ, then the pointers are finally set
C so that SAVR is identical to (YH*,MAXORD+2)
C-----------------------------------------------------------------------
      LYHD = LYH - LYHN
      IMAX = LYHN - 1 + LENYHM
C Move YH.  Move right if LYHD < 0; move left if LYHD > 0. -------------
      IF (LYHD .LT. 0) THEN
        DO 72 I = LYHN,IMAX
          J = IMAX + LYHN - I
 72       RWORK(J) = RWORK(J+LYHD)
      ENDIF
      IF (LYHD .GT. 0) THEN
        DO 76 I = LYHN,IMAX
 76       RWORK(I) = RWORK(I+LYHD)
      ENDIF
 80   LYH = LYHN
      IWORK(22) = LYH
      IF (MOSS .NE. 2 .AND. MOSS .NE. 4) GO TO 85
C Temporarily load EWT if MOSS = 2 or 4.
      CALL SEWSET (N,ITOL,RTOL,ATOL,RWORK(LYH),RWORK(LEWT))
      DO 82 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 82     RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 85     CONTINUE
C SIPREPI and SPREPI do sparse matrix preprocessing. -------------------
      LSAVF = MIN(LSAVF,LRW)
      LEWT = MIN(LEWT,LRW)
      LACOR = MIN(LACOR,LRW)
      CALL SIPREPI (NEQ, Y, YDOTI, RWORK, IWORK(LIA), IWORK(LJA),
     1   IWORK(LIC), IWORK(LJC), IPFLAG, RES, JAC, ADDA)
      LENRW = LWM - 1 + LENWK + LREST
      IWORK(17) = LENRW
      IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
      IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
      IPGO = -IPFLAG + 1
      GO TO (90, 628, 629, 630, 631, 632, 633, 634, 634), IPGO
 90   IWORK(22) = LYH
      LYD0 = LYH + N
      IF (LENRW .GT. LRW) GO TO 617
C Set flag to signal changes to SSTODI.---------------------------------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 94
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI. --------
      DO 92 I = 1,N
 92     YDOTI(I) = RWORK(I+LSAVF-1)
 94   IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0E0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 0 or 1).
C It contains all remaining initializations, the call to SAINVGS
C (if ISTATE = 0), the sparse matrix preprocessing, and the
C calculation if the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  CONTINUE
      LYH = LYHN
      IWORK(22) = LYH
      TN = T
      NST = 0
      NFE = 0
      H = 1.0E0
      NNZ = 0
      NGP = 0
      NZL = 0
      NZU = 0
C Load the initial value vector in YH.----------------------------------
      DO 105 I = 1,N
 105    RWORK(I+LYH-1) = Y(I)
      IF (ISTATE .NE. 1) GO TO 108
C Initial dy/dt was supplied.  Load it into YH (LYD0 points to YH(*,2).)
      LYD0 = LYH + NYH
      DO 106 I = 1,N
 106    RWORK(I+LYD0-1) = YDOTI(I)
 108  CONTINUE
C Load and invert the EWT array.  (H is temporarily set to 1.0.)--------
      CALL SEWSET (N,ITOL,RTOL,ATOL,RWORK(LYH),RWORK(LEWT))
      DO 110 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 621
 110    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
C Call SIPREPI and SPREPI to do sparse matrix preprocessing.------------
      LACOR = MIN(LACOR,LRW)
      CALL SIPREPI (NEQ, Y, YDOTI, RWORK, IWORK(LIA), IWORK(LJA),
     1   IWORK(LIC), IWORK(LJC), IPFLAG, RES, JAC, ADDA)
      LENRW = LWM - 1 + LENWK + LREST
      IWORK(17) = LENRW
      IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
      IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
      IPGO = -IPFLAG + 1
      GO TO (115, 628, 629, 630, 631, 632, 633, 634, 634), IPGO
 115  IWORK(22) = LYH
      IF (LENRW .GT. LRW) GO TO 617
C Compute initial dy/dt, if necessary, and load it into YH.-------------
      LYD0 = LYH + N
      IF (ISTATE .NE. 0) GO TO 120
      CALL SAINVGS (NEQ, T, Y, RWORK(LWM), RWORK(LWM), RWORK(LACOR),
     1              RWORK(LYD0), IER, RES, ADDA)
      NFE = NFE + 1
      IGO = IER + 1
      GO TO (120, 565, 560, 560), IGO
C Check TCRIT for legality (ITASK = 4 or 5). ---------------------------
 120  CONTINUE
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 125
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0E0) GO TO 625
      IF (H0 .NE. 0.0E0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0E0)
     1   H0 = TCRIT - T
C Initialize all remaining parameters. ---------------------------------
 125  UROUND = RUMACH()
      JSTART = 0
      RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NJE = 0
      NLU = 0
      NSLAST = 0
      HU = 0.0E0
      NQU = 0
      CCMAX = 0.3E0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
C                                       1
C where   w0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(i) = i-th component of initial value of dy/dt,
C         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0E0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0E0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 145
      DO 140 I = 1,N
 140    TOL = MAX(TOL,RTOL(I))
 145  IF (TOL .GT. 0.0E0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0E0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0E0*UROUND)
      TOL = MIN(TOL,0.001E0)
      SUM = SVNORM (N, RWORK(LYD0), RWORK(LEWT))
      SUM = 1.0E0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0E0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0E0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LYD0-1) = H0*RWORK(I+LYD0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0E0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0E0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 245
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0E0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator SSTODI.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0E0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0E0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0E0) GO TO 280
      TOLSF = TOLSF*2.0E0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SLSODIS- Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWV (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SLSODIS- Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWV (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL SSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,WM,RES,
C                 ADDA,JAC,SPRJIS,SSOLSS)
C Note: SAVF in SSTODI occupies the same space as YDOTI in SLSODIS.
C-----------------------------------------------------------------------
      CALL SSTODI (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   YDOTI, RWORK(LSAVF), RWORK(LACOR), RWORK(LWM),
     2   RWORK(LWM), RES, ADDA, JAC, SPRJIS, SSOLSS )
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 400, 550, 555), KGO
C
C KGO = 1:success; 2:error test failure; 3:convergence failure;
C       4:RES ordered return; 5:RES returned error;
C       6:fatal error from SCDRV via SPRJIS or SSOLSS.
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  iF ((TN - TOUT)*H .LT. 0.0E0) GO TO 250
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0E0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0E0) GO TO 345
      CALL SINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0E0 + 4.0E0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0E0) GO TO 250
      H = (TCRIT - TN)*(1.0E0 - 4.0E0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0E0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from SLSODIS.
C if ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      IF ( KFLAG .EQ. -3 )  ISTATE = 3
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNZ
      IWORK(20) = NGP
      IWORK(21) = NLU
      IWORK(25) = NZL
      IWORK(26) = NZU
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SLSODIS- At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0E0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODIS- At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 590
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SLSODIS- At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 590
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SLSODIS- At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='     error test failed repeatedly or with ABS(H) = HMIN     '
      CALL XERRWV (MSG, 60, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 570
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'SLSODIS- At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 570
C IRES = 3 returned by RES, despite retries by SSTODI. -----------------
 550  MSG = 'SLSODIS- At T (=R1) residual routine returned     '
      CALL XERRWV (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '    error IRES = 3 repeatedly.'
      CALL XERRWV (MSG, 30, 206, 1, 0, 0, 0, 0, TN, 0.0E0)
      ISTATE = -7
      GO TO 590
C KFLAG = -5.  Fatal error flag returned by SPRJIS or SSOLSS (SCDRV). --
 555  MSG = 'SLSODIS- At T (=R1) and step size H (=R2), a fatal'
      CALL XERRWV (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      error flag was returned by SCDRV (by way of '
      CALL XERRWV (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      Subroutine SPRJIS or SSOLSS)      '
      CALL XERRWV (MSG, 40, 207, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -9
      GO TO 580
C SAINVGS failed because matrix A was singular. ------------------------
 560  MSG='SLSODIS- Attempt to initialize dy/dt failed because matrix A'
      CALL XERRWV (MSG, 60, 208, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='     was singular.  SCDRV returned zero pivot error flag.   '
      CALL XERRWV (MSG, 60, 208, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = 'SAINVGS set its error flag to IER = (I1)'
      CALL XERRWV (MSG, 40, 208, 0, 1, IER, 0, 0, 0.0E0, 0.0E0)
      ISTATE = -8
      RETURN
C SAINVGS failed because RES set IRES to 2 or 3. -----------------------
 565  MSG = 'SLSODIS- Attempt to initialize dy/dt failed       '
      CALL XERRWV (MSG, 50, 209, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      because residual routine set its error flag '
      CALL XERRWV (MSG, 50, 209, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      to IRES = (I1)'
      CALL XERRWV (MSG, 20, 209, 0, 1, IER, 0, 0, 0.0E0, 0.0E0)
      ISTATE = -8
      RETURN
C Compute IMXER if relevant. -------------------------------------------
 570  BIG = 0.0E0
      IMXER = 1
      DO 575 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 575
        BIG = SIZE
        IMXER = I
 575    CONTINUE
      IWORK(16) = IMXER
C Compute residual if relevant. ----------------------------------------
 580  LYD0 = LYH + NYH
      DO 585  I = 1, N
         RWORK(I+LSAVF-1) = RWORK(I+LYD0-1) / H
 585     Y(I) = RWORK(I+LYH-1)
      IRES = 1
      CALL RES (NEQ, TN, Y, RWORK(LSAVF), YDOTI, IRES)
      NFE = NFE + 1
      IF ( IRES .LE. 1 )  GO TO 595
      MSG = 'SLSODIS- Residual routine set its flag IRES       '
      CALL XERRWV (MSG, 50, 210, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG = '      to (I1) when called for final output.       '
      CALL XERRWV (MSG, 50, 210, 0, 1, IRES, 0, 0, 0.0E0, 0.0E0)
      GO TO 595
C set y vector, t, and optional outputs. -------------------------------
 590  DO 592 I = 1,N
 592    Y(I) = RWORK(I+LYH-1)
 595  T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNZ
      IWORK(20) = NGP
      IWORK(21) = NLU
      IWORK(25) = NZL
      IWORK(26) = NZU
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'SLSODIS- ISTATE (=I1) illegal.'
      CALL XERRWV (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0E0, 0.0E0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SLSODIS- ITASK (=I1) illegal. '
      CALL XERRWV (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 603  MSG = 'SLSODIS-ISTATE .gt. 1 but SLSODIS not initialized.'
      CALL XERRWV (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 604  MSG = 'SLSODIS- NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 605  MSG = 'SLSODIS- ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWV (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0E0, 0.0E0)
      GO TO 700
 606  MSG = 'SLSODIS- ITOL (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 607  MSG = 'SLSODIS- IOPT (=I1) illegal.  '
      CALL XERRWV (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 608  MSG = 'SLSODIS- MF (=I1) illegal.    '
      CALL XERRWV (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 611  MSG = 'SLSODIS- MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 612  MSG = 'SLSODIS- MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 613  MSG = 'SLSODIS- MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0E0, 0.0E0)
      GO TO 700
 614  MSG = 'SLSODIS- TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0E0)
      GO TO 700
 615  MSG = 'SLSODIS- HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0E0)
      GO TO 700
 616  MSG = 'SLSODIS- HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0E0)
      GO TO 700
 617  MSG = 'SLSODIS- RWORK length is insufficient to proceed. '
      CALL XERRWV (MSG, 50, 17, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 618  MSG = 'SLSODIS- IWORK length is insufficient to proceed. '
      CALL XERRWV (MSG, 50, 18, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0E0, 0.0E0)
      GO TO 700
 619  MSG = 'SLSODIS- RTOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0E0)
      GO TO 700
 620  MSG = 'SLSODIS- ATOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWV (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0E0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SLSODIS- EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0E0)
      GO TO 700
 622  MSG='SLSODIS- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWV (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='SLSODIS- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='SLSODIS- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='SLSODIS- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SLSODIS- At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0E0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'SLSODIS- Trouble in SINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0E0)
      GO TO 700
 628  MSG='SLSODIS- RWORK length insufficient (for Subroutine SPREPI). '
      CALL XERRWV (MSG, 60, 28, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 28, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 629  MSG='SLSODIS- RWORK length insufficient (for Subroutine JGROUP). '
      CALL XERRWV (MSG, 60, 29, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 29, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 630  MSG='SLSODIS- RWORK length insufficient (for Subroutine SODRV).  '
      CALL XERRWV (MSG, 60, 30, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 30, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 631  MSG='SLSODIS- Error from SODRV in Yale Sparse Matrix Package.    '
      CALL XERRWV (MSG, 60, 31, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      IMUL = (IYS - 1)/N
      IREM = IYS - IMUL*N
      MSG='      At T (=R1), SODRV returned error flag = I1*NEQ + I2.  '
      CALL XERRWV (MSG, 60, 31, 0, 2, IMUL, IREM, 1, TN, 0.0E0)
      GO TO 700
 632  MSG='SLSODIS- RWORK length insufficient (for Subroutine SCDRV).  '
      CALL XERRWV (MSG, 60, 32, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 32, 0, 2, LENRW, LRW, 0, 0.0E0, 0.0E0)
      GO TO 700
 633  MSG='SLSODIS- Error from SCDRV in Yale Sparse Matrix Package.    '
      CALL XERRWV (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      IMUL = (IYS - 1)/N
      IREM = IYS - IMUL*N
      MSG='      At T (=R1), SCDRV returned error flag = I1*NEQ + I2.  '
      CALL XERRWV (MSG, 60, 33, 0, 2, IMUL, IREM, 1, TN, 0.0E0)
      IF (IMUL .EQ. 2) THEN
      MSG='        Duplicate entry in sparsity structure descriptors.  '
      CALL XERRWV (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      IF (IMUL .EQ. 3 .OR. IMUL .EQ. 6) THEN
      MSG='        Insufficient storage for SNSFC (called by SCDRV).   '
      CALL XERRWV (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      ENDIF
      GO TO 700
 634  MSG='SLSODIS- At T (=R1) residual routine (called by SPREPI)     '
      CALL XERRWV (MSG, 60, 34, 0, 0, 0, 0, 0, 0.0E0, 0.0E0)
      IER = -IPFLAG - 5
      MSG = '     returned error IRES (=I1)'
      CALL XERRWV (MSG, 30, 34, 0, 1, IER, 0, 1, TN, 0.0E0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'SLSODIS- Run aborted.. apparent infinite loop.    '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0E0, 0.0E0)
      RETURN
C----------------------- End of Subroutine SLSODIS ---------------------
      END
