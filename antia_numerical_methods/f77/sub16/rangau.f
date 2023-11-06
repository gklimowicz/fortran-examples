!	To generate random numbers with Gaussian probability distribution
!	It generates random numbers with zero mean and variance of 1.
!	
!	SEED : (input/output) real seed, it should be positive and
!		less than AM. It is updated by the routine and should
!		not be modified between two calls, unless a fresh
!		sequence is required
!
!       THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
!       VERSION AS THE SEED IS NOW REAL*16 INSTEAD OF INTEGER.
!
!	Required routines : None

      FUNCTION RANGAU(SEED)
      IMPLICIT REAL*16(A-H,O-Z)
      PARAMETER(AM=2147483648Q0,A=45875Q0,AC=453816693Q0,
     1                  AN=2147483647Q0)
      PARAMETER(PI=3.14159265358979323846264338327950288Q0)

      R1=1+MOD(A*SEED+AC,AM)
      IF(SEED.EQ.0.0) SEED=0.1Q0
      RANGAU=SQRT(2.Q0*LOG(AN/SEED))*COS(2.0*PI*R1/AN)
      SEED=R1
      END
