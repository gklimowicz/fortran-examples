!     To integrate a function with logarithmic singularity using Gaussian formulas
!
!     RINT : (output) Calculated value of the integral
!     A : (input) The upper limit
!     AEPS : (input) The required absolute accuracy
!     REPS : (input) The required relative accuracy
!		The estimated error should be less than MAX(AEPS,REPS*ABS(RINT))
!     DIF : (output) estimated (absolute) error achieved by the subroutine
!     F : (input) Name of the function routine to calculate the
!		integrand (divided by LOG(A/X))
!     NPT : (output) Number of function evaluations used by the subroutine
!     IER : (output) Error parameter, IER=0 implies successful execution
!     	IER=30 implies specified accuracy was not achieved
!     		DIF will contain the estimated accuracy
!
!     Function F(X) must be supplied by the user
!     Note that subroutine calculates integral of F(X)*LOG(A/X)
!
!	Required routines : F

      SUBROUTINE GAULOG(RINT,A,AEPS,REPS,DIF,F,NPT,IER)
!      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(30),X(30)
 
!     Weights and abscissas for Gaussian formula with logarithmic singularity
!     W(N-1),...,W(2N-2), are the weights for N-point rule and
!     X(N-1),...,X(2N-2), the corresponding abscissas
!     Weights and abscissas are available for N=2,4,8,16
 
      DATA (X(I),I=1,30)/1.120088061669761830D-1,6.022769081187381028D-1
     *,       4.144848019938322080D-2, 2.452749143206022519D-1,
     *        5.561654535602758372D-1, 8.489823945329851746D-1,
     *        1.332024416089246501D-2, 7.975042901389493841D-2,
     *        1.978710293261880538D-1, 3.541539943519094197D-1,
     *        5.294585752349172777D-1, 7.018145299390999638D-1,
     *        8.493793204411066760D-1, 9.533264500563597888D-1,
     *        3.897834487115909095D-3, 2.302894561687320045D-2,
     *        5.828039830624031972D-2, 1.086783650910538817D-1,
     *        1.726094549098437244D-1, 2.479370544705782363D-1,
     *        3.320945491299168705D-1, 4.221839105819483085D-1,
     *        5.150824733814623250D-1, 6.075561204477284747D-1,
     *        6.963756532282138523D-1, 7.784325658732652431D-1,
     *        8.508502697153909688D-1, 9.110868572222718348D-1,
     *        9.570255717035421226D-1, 9.870478002479844660D-1/

      DATA (W(I),I=1,30)/7.185393190303844407D-1,2.814606809696155593D-1
     *,       3.834640681451351249D-1, 3.868753177747626273D-1,
     *        1.904351269501424154D-1, 3.922548712995983245D-2,
     *        1.644166047280028868D-1, 2.375256100233060205D-1,
     *        2.268419844319191264D-1, 1.757540790060702450D-1,
     *        1.129240302467590519D-1, 5.787221071778207240D-2,
     *        2.097907374213297804D-2, 3.686407104027619013D-3,
     *        6.079171004359114509D-2, 1.029156775175820228D-1,
     *        1.223556620460090919D-1, 1.275692469370159323D-1,
     *        1.230135746000709083D-1, 1.118472448554857552D-1,
     *        9.659638515212439849D-2, 7.935666435147320573D-2,
     *        6.185049458196527197D-2, 4.543524650772672381D-2,
     *        3.109897475158184829D-2, 1.945976592736087029D-2,
     *        1.077625496320554213D-2, 4.972542890087649610D-3,
     *        1.678201110051197249D-3, 2.823537646684367889D-4/
 
      IER=0
!     The 2-point formula
      R1=(F(A*X(1))*W(1)+F(A*X(2))*W(2))*A
      NPT=2
      N=2
 
!     Use higher order formula until convergence
      DO 2000 J=2,4
        N=N*2
        R2=0.0
        DO 1000 I=N-1,2*N-2
1000    R2=R2+F(X(I)*A)*W(I)
        R2=R2*A
 
        NPT=NPT+N
        DIF=R2-R1
        RINT=R2
        IF(ABS(DIF).LT.MAX(AEPS,REPS*ABS(RINT))) RETURN
        R1=R2
2000  CONTINUE
 
!     Integral fails to converge
      IER=30
      RETURN
      END
