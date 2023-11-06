!	To generate uniformly distributed random numbers in interval (0,1)
!
!	SEED : (input/output) is a real value used as the seed
!		It should be positive during initial call and
!		should not be modified between different calls.
!
!	Required routines : None

      FUNCTION RAN1(SEED)
      IMPLICIT REAL*16(A-H,O-Z)
!	Retain the following declaration even for REAL*4 version
!	otherwise AM and AC will be rounded
      PARAMETER(AM=2147483648Q0,A=45875,AC=453816693Q0,AN=2147483647Q0)

      SEED=MOD(SEED*A+AC,AM)
      RAN1=SEED/AN
      END
