!	To generate uniformly distributed random numbers in interval (0,1)
!
!	ISEED : (input/output) is an integer value used as the seed
!		It should be initialised to negative value before first call
!		and should not be modified between successive calls.
!
!	Required routines : None

      FUNCTION RANF(ISEED)
!      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(M1=714025,IA1=1366,IC1=150889)
      PARAMETER(M2=214326,IA2=3613,IC2=45289)
      PARAMETER(M3=139968,IA3=3877,IC3=29573,ISH=43)
      DIMENSION RAN(ISH)
      SAVE
      DATA IFLG/0/

!	Initialise on first call or when ISEED<0
      IF(ISEED.LT.0.OR.IFLG.EQ.0) THEN
        IFLG=1
        RM1=1.D0/M1
        RM2=1.D0/M2

!	Seeds for the three random number generators
        IS1=MOD(-ISEED,M1)
        IS2=MOD(IA1*IS1+IC1,M1)
        IS3=MOD(IA2*IS2+IC2,M2)
        ISEED=1

!	Store ISH random numbers in the array RAN
        DO 1000 J=1,ISH
          IS1=MOD(IA1*IS1+IC1,M1)
          IS2=MOD(IA2*IS2+IC2,M2)
          RAN(J)=(FLOAT(IS1)+FLOAT(IS2)*RM2)*RM1
1000    CONTINUE
      ENDIF

      IS1=MOD(IA1*IS1+IC1,M1)
      IS2=MOD(IA2*IS2+IC2,M2)
      IS3=MOD(IA3*IS3+IC3,M3)
!	Select a random entry from RAN and store a new number in its place
      I=1+(ISH*IS3)/M3
      RANF=RAN(I)
      RAN(I)=(FLOAT(IS1)+FLOAT(IS2)*RM2)*RM1
      END
