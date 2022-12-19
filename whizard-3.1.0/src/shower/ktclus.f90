! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module ktclus

  use kinds, only: default

  public :: ktclur
  public :: ktreco

contains

!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
!C     KTCLUS: written by Mike Seymour, July 1992.
!C     Last modified November 2000.
!C     Please send comments or suggestions to Mike.Seymour@rl.ac.uk
!C
!C     This is a general-purpose kt clustering package.
!C     It can handle ee, ep and pp collisions.
!C     It is loosely based on the program of Siggi Bethke.
!C
!C     The time taken (on a 10MIP machine) is (0.2microsec)*N**3
!C     where N is the number of particles.
!C     Over 90 percent of this time is used in subroutine KTPMIN, which
!C     simply finds the minimum member of a one-dimensional array.
!C     It is well worth thinking about optimization: on the SPARCstation
!C     a factor of two increase was obtained simply by increasing the
!C     optimization level from its default value.
!C
!C     The approach is to separate the different stages of analysis.
!C     KTCLUS does all the clustering and records a merging history.
!C     It returns a simple list of the y values at which each merging
!C     occured. Then the following routines can be called to give extra
!C     information on the most recently analysed event.
!C     KTCLUR is identical but includes an R parameter, see below.
!C     KTYCUT gives the number of jets at each given YCUT value.
!C     KTYSUB gives the number of sub-jets at each given YCUT value.
!C     KTBEAM gives same info as KTCLUS but only for merges with the beam
!C     KTJOIN gives same info as KTCLUS but for merges of sub-jets.
!C     KTRECO reconstructs the jet momenta at a given value of YCUT.
!C     It also gives information on which jets at scale YCUT belong to
!C     which macro-jets at scale YMAC, for studying sub-jet properties.
!C     KTINCL reconstructs the jet momenta according to the inclusive jet
!C     definition of Ellis and Soper.
!C     KTISUB, KTIJOI and KTIREC are like KTYSUB, KTJOIN and KTRECO,
!C     except that they only apply to one inclusive jet at a time,
!C     with the pt of that jet automatically used for ECUT.
!C     KTWICH gives a list of which particles ended up in which jets.
!C     KTWCHS gives the same thing, but only for subjets.
!C     Note that the numbering of jets used by these two routines is
!C     guaranteed to be the same as that used by KTRECO.
!C
!C     The collision type and analysis type are indicated by the first
!C     argument of KTCLUS. IMODE=<TYPE><ANGLE><MONO><RECOM> where
!C     TYPE:  1=>ee, 2=>ep with p in -z direction, 3=>pe, 4=>pp
!C     ANGLE: 1=>angular kt def., 2=>DeltaR, 3=>f(DeltaEta,DeltaPhi)
!C            where f()=2(cosh(eta)-cos(phi)) is the QCD emission metric
!C     MONO:  1=>derive relative pseudoparticle angles from jets
!C            2=>monotonic definitions of relative angles
!C     RECOM: 1=>E recombination scheme, 2=>pt scheme, 3=>pt**2 scheme
!C
!C     There are also abbreviated forms for the most common combinations:
!C     IMODE=1 => E scheme in e+e-                              (=1111)
!C           2 => E scheme in ep                                (=2111)
!C           3 => E scheme in pe                                (=3111)
!C           4 => E scheme in pp                                (=4111)
!C           5 => covariant E scheme in pp                      (=4211)
!C           6 => covariant pt-scheme in pp                     (=4212)
!C           7 => covariant monotonic pt**2-scheme in pp        (=4223)
!C
!C     KTRECO no longer needs to reconstruct the momenta according to the
!C     same recombination scheme in which they were clustered. Its first
!C     argument gives the scheme, taking the same values as RECOM above.
!C
!C     Note that unlike previous versions, all variables which hold y
!C     values have been named in a consistent way:
!C     Y()  is the output scale at which jets were merged,
!C     YCUT is the input scale at which jets should be counted, and
!C          jet-momenta reconstructed etc,
!C     YMAC is the input macro-jet scale, used in determining whether
!C          or not each jet is a sub-jet.
!C     The original scheme defined in our papers is equivalent to always
!C     setting YMAC=1.
!C     Whenever a YCUT or YMAC variable is used, it is rounded down
!C     infinitesimally, so that for example, setting YCUT=Y(2) refers
!C     to the scale where the event is 2-jet, even if rounding errors
!C     have shifted its value slightly.
!C
!C     An R parameter can be used in hadron-hadron collisions by
!C     calling KTCLUR instead of KTCLUS.  This is as suggested by
!C     Ellis and Soper, but implemented slightly differently,
!C     as in M.H. Seymour, LU TP 94/2 (submitted to Nucl. Phys. B.).
!C     R**2 multiplies the single Kt everywhere it is used.
!C     Calling KTCLUR with R=1 is identical to calling KTCLUS.
!C     R plays a similar role to the jet radius in a cone-type algorithm,
!C     but is scaled up by about 40% (ie R=0.7 in a cone algorithm is
!C     similar to this algorithm with R=1).
!C     Note that R.EQ.1 must be used for the e+e- and ep versions,
!C     and is strongly recommended for the hadron-hadron version.
!C     However, R values smaller than 1 have been found to be useful for
!C     certain applications, particularly the mass reconstruction of
!C     highly-boosted colour-singlets such as high-pt hadronic Ws,
!C     as in M.H. Seymour, LU TP 93/8 (to appear in Z. Phys. C.).
!C     Situations in which R<1 is useful are likely to also be those in
!C     which the inclusive reconstruction method is more useful.
!C
!C     Also included is a set of routines for doing Lorentz boosts:
!C     KTLBST finds the boost matrix to/from the cm frame of a 4-vector
!C     KTRROT finds the rotation matrix from one vector to another
!C     KTMMUL multiplies together two matrices
!C     KTVMUL multiplies a vector by a matrix
!C     KTINVT inverts a transformation matrix (nb NOT a general 4 by 4)
!C     KTFRAM boosts a list of vectors between two arbitrary frames
!C     KTBREI boosts a list of vectors between the lab and Breit frames
!C     KTHADR boosts a list of vectors between the lab and hadronic cmf
!C       The last two need the momenta in the +z direction of the lepton
!C       and hadron beams, and the 4-momentum of the outgoing lepton.
!C
!C     The main reference is:
!C       S. Catani, Yu.L. Dokshitzer, M.H. Seymour and B.R. Webber,
!C         Nucl.Phys.B406(1993)187.
!C     The ep version was proposed in:
!C       S. Catani, Yu.L. Dokshitzer and B.R. Webber,
!C         Phys.Lett.285B(1992)291.
!C     The inclusive reconstruction method was proposed in:
!C       S.D. Ellis and D.E. Soper,
!C         Phys.Rev.D48(1993)3160.
!C
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
  SUBROUTINE KTCLUR(IMODE,PP,NN,R,ECUT,Y,*)
    use io_units
    IMPLICIT NONE
!C---DO CLUSTER ANALYSIS OF PARTICLES IN PP
!C
!C   IMODE   = INPUT  : DESCRIBED ABOVE
!C   PP(I,J) = INPUT  : 4-MOMENTUM OF Jth PARTICLE: I=1,4 => PX,PY,PZ,E
!C   NN      = INPUT  : NUMBER OF PARTICLES
!C   R       = INPUT  : ELLIS AND SOPER'S R PARAMETER, SEE ABOVE.
!C   ECUT    = INPUT  : DENOMINATOR OF KT MEASURE. IF ZERO, ETOT IS USED
!C   Y(J)    = OUTPUT : VALUE OF Y FOR WHICH EVENT CHANGES FROM BEING
!C                        J JET TO J-1 JET
!C   LAST ARGUMENT IS LABEL TO JUMP TO IF FOR ANY REASON THE EVENT
!C   COULD NOT BE PROCESSED (MOST LIKELY DUE TO TOO MANY PARTICLES)
!C
!C   NOTE THAT THE MOMENTA ARE DECLARED DOUBLE PRECISION,
!C   AND ALL OTHER FLOATING POINT VARIABLES ARE DECLARED DOUBLE PRECISION
!C
    INTEGER NMAX,IM,IMODE,TYPE,ANGL,MONO,RECO,N,I,J,NN, &
         IMIN,JMIN,KMIN,NUM,HIST,INJET,IABBR,NABBR
    PARAMETER (NMAX=512,NABBR=7)
    DOUBLE PRECISION PP(4,*)
    integer :: u
!CHANGE    DOUBLE PRECISION R,ECUT,Y(*),P,KT,ETOT,RSQ,KTP,KTS,KTPAIR,KTSING, &
    DOUBLE PRECISION R,ECUT,Y(*),P,KT,ETOT,RSQ,KTP,KTS, &
         KTMIN,ETSQ,KTLAST,KTMAX,KTTMP
    LOGICAL FIRST
    CHARACTER TITLE(4,4)*10
!C---KT RECORDS THE KT**2 OF EACH MERGING.
!C---KTLAST RECORDS FOR EACH MERGING, THE HIGHEST ECUT**2 FOR WHICH THE
!C   RESULT IS NOT MERGED WITH THE BEAM (COULD BE LARGER THAN THE
!C   KT**2 AT WHICH IT WAS MERGED IF THE KT VALUES ARE NOT MONOTONIC).
!C   THIS MAY SOUND POINTLESS, BUT ITS USEFUL FOR DETERMINING WHETHER
!C   SUB-JETS SURVIVED TO SCALE Y=YMAC OR NOT.
!C---HIST RECORDS MERGING HISTORY:
!C   N=>DELETED TRACK N, M*NMAX+N=>MERGED TRACKS M AND N (M<N).
    COMMON /KTCOMM/ETOT,RSQ,P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX), &
         KT(NMAX),KTLAST(NMAX),HIST(NMAX),NUM
    DIMENSION INJET(NMAX),IABBR(NABBR)
    DATA FIRST,TITLE,IABBR/.TRUE., &
         'e+e-      ','ep        ','pe        ','pp        ', &
         'angle     ','DeltaR    ','f(DeltaR) ','**********', &
         'no        ','yes       ','**********','**********', &
         'E         ','Pt        ','Pt**2     ','**********', &
         1111,2111,3111,4111,4211,4212,4223/
!C---CHECK INPUT
    IM=IMODE
    IF (IM.GE.1.AND.IM.LE.NABBR) IM=IABBR(IM)
    TYPE=MOD(IM/1000,10)
    ANGL=MOD(IM/100 ,10)
    MONO=MOD(IM/10  ,10)
    RECO=MOD(IM     ,10)
    IF (NN.GT.NMAX) CALL KTWARN('KT-MAX',100,*999)
    IF (NN.LT.1) CALL KTWARN('KT-LT1',100,*999)
    IF (NN.LT.2.AND.TYPE.EQ.1) CALL KTWARN('KT-LT2',100,*999)
    IF (TYPE.LT.1.OR.TYPE.GT.4.OR.ANGL.LT.1.OR.ANGL.GT.4.OR. &
         MONO.LT.1.OR.MONO.GT.2.OR.RECO.LT.1.OR.RECO.GT.3) CALL KTWARN('KTCLUS',101,*999)
    u = given_output_unit ()
    IF (FIRST) THEN
       WRITE (u,'(/,1X,54(''*'')/A)') &
            ' KTCLUS: written by Mike Seymour, July 1992.'
       WRITE (u,'(A)') &
            ' Last modified November 2000.'
       WRITE (u,'(A)') &
            ' Please send comments or suggestions to Mike.Seymour@rl.ac.uk'
       WRITE (u,'(/A,I2,2A)') &
            '       Collision type =',TYPE,' = ',TITLE(TYPE,1)
       WRITE (u,'(A,I2,2A)') &
            '     Angular variable =',ANGL,' = ',TITLE(ANGL,2)
       WRITE (u,'(A,I2,2A)') &
            ' Monotonic definition =',MONO,' = ',TITLE(MONO,3)
       WRITE (u,'(A,I2,2A)') &
            ' Recombination scheme =',RECO,' = ',TITLE(RECO,4)
       IF (R.NE.1) THEN
          WRITE (u,'(A,F5.2)') &
               '     Radius parameter =',R
          IF (TYPE.NE.4) WRITE (u,'(A)') &
               ' R.NE.1 is strongly discouraged for this collision type!'
       ENDIF
       WRITE (u,'(1X,54(''*'')/)')
       FIRST=.FALSE.
    ENDIF
!C---COPY PP TO P
    N=NN
    NUM=NN
    CALL KTCOPY(PP,N,P,(RECO.NE.1))
    ETOT=0
    DO I=1,N
       ETOT=ETOT+P(4,I)
    END DO
    IF (ETOT.EQ.0) CALL KTWARN('KTCLUS',102,*999)
    IF (ECUT.EQ.0) THEN
       ETSQ=1/ETOT**2
    ELSE
       ETSQ=1/ECUT**2
    ENDIF
    RSQ=R**2
!C---CALCULATE ALL PAIR KT's
    DO I=1,N-1
       DO J=I+1,N
          KTP(J,I)=-1
          KTP(I,J)=KTPAIR(ANGL,P(1,I),P(1,J),KTP(J,I))
       END DO
    END DO
!C---CALCULATE ALL SINGLE KT's
    DO I=1,N
       KTS(I)=KTSING(ANGL,TYPE,P(1,I))
    END DO
    KTMAX=0
!C---MAIN LOOP
300 CONTINUE
!C---FIND MINIMUM MEMBER OF KTP
    CALL KTPMIN(KTP,NMAX,N,IMIN,JMIN)
!C---FIND MINIMUM MEMBER OF KTS
    CALL KTSMIN(KTS,NMAX,N,KMIN)
!C---STORE Y VALUE OF TRANSITION FROM N TO N-1 JETS
    KTMIN=KTP(IMIN,JMIN)
    KTTMP=RSQ*KTS(KMIN)
    IF ((TYPE.GE.2.AND.TYPE.LE.4).AND. &
         (KTTMP.LE.KTMIN.OR.N.EQ.1)) &
         KTMIN=KTTMP
    KT(N)=KTMIN
    Y(N)=KT(N)*ETSQ
!C---IF MONO.GT.1, SEQUENCE IS SUPPOSED TO BE MONOTONIC, IF NOT, WARN
    IF (KTMIN.LT.KTMAX.AND.MONO.GT.1) CALL KTWARN('KTCLUS',1,*999)
    IF (KTMIN.GE.KTMAX) KTMAX=KTMIN
!C---IF LOWEST KT IS TO A BEAM, THROW IT AWAY AND MOVE LAST ENTRY UP
    IF (KTMIN.EQ.KTTMP) THEN
       CALL KTMOVE(P,KTP,KTS,NMAX,N,KMIN,1)
!C---UPDATE HISTORY AND CROSS-REFERENCES
       HIST(N)=KMIN
       INJET(N)=KMIN
       DO I=N,NN
          IF (INJET(I).EQ.KMIN) THEN
             KTLAST(I)=KTMAX
             INJET(I)=0
          ELSEIF (INJET(I).EQ.N) THEN
             INJET(I)=KMIN
          ENDIF
       END DO
!C---OTHERWISE MERGE JETS IMIN AND JMIN AND MOVE LAST ENTRY UP
    ELSE
       CALL KTMERG(P,KTP,KTS,NMAX,IMIN,JMIN,N,TYPE,ANGL,MONO,RECO)
       CALL KTMOVE(P,KTP,KTS,NMAX,N,JMIN,1)
!C---UPDATE HISTORY AND CROSS-REFERENCES
       HIST(N)=IMIN*NMAX+JMIN
       INJET(N)=IMIN
       DO I=N,NN
          IF (INJET(I).EQ.JMIN) THEN
             INJET(I)=IMIN
          ELSEIF (INJET(I).EQ.N) THEN
             INJET(I)=JMIN
          ENDIF
       END DO
    ENDIF
!C---THATS ALL THERE IS TO IT
    N=N-1
    IF (N.GT.1 .OR. N.GT.0.AND.(TYPE.GE.2.AND.TYPE.LE.4)) GOTO 300
    IF (N.EQ.1) THEN
       KT(N)=1D20
       Y(N)=KT(N)*ETSQ
    ENDIF
    RETURN
999 RETURN 1
  END SUBROUTINE KTCLUR
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
    SUBROUTINE KTRECO(RECO,PP,NN,ECUT,YCUT,YMAC,PJET,JET,NJET,NSUB,*)
      IMPLICIT NONE
!C---RECONSTRUCT KINEMATICS OF JET SYSTEM, WHICH HAS ALREADY BEEN
!C   ANALYSED BY KTCLUS. NOTE THAT NO CONSISTENCY CHECK IS MADE: USER
!C   IS TRUSTED TO USE THE SAME PP VALUES AS FOR KTCLUS
!C
!C   RECO     = INPUT : RECOMBINATION SCHEME (NEED NOT BE SAME AS KTCLUS)
!C   PP(I,J)  = INPUT : 4-MOMENTUM OF Jth PARTICLE: I=1,4 => PX,PY,PZ,E
!C   NN       = INPUT : NUMBER OF PARTICLES
!C   ECUT     = INPUT : DENOMINATOR OF KT MEASURE. IF ZERO, ETOT IS USED
!C   YCUT     = INPUT : Y VALUE AT WHICH TO RECONSTRUCT JET MOMENTA
!C   YMAC     = INPUT : Y VALUE USED TO DEFINE MACRO-JETS, TO DETERMINE
!C                        WHICH JETS ARE SUB-JETS
!C   PJET(I,J)=OUTPUT : 4-MOMENTUM OF Jth JET AT SCALE YCUT
!C   JET(J)   =OUTPUT : THE MACRO-JET WHICH CONTAINS THE Jth JET,
!C                        SET TO ZERO IF JET IS NOT A SUB-JET
!C   NJET     =OUTPUT : THE NUMBER OF JETS
!C   NSUB     =OUTPUT : THE NUMBER OF SUB-JETS (EQUAL TO THE NUMBER OF
!C                        NON-ZERO ENTRIES IN JET())
!C   LAST ARGUMENT IS LABEL TO JUMP TO IF FOR ANY REASON THE EVENT
!C   COULD NOT BE PROCESSED
!C
!C   NOTE THAT THE MOMENTA ARE DECLARED DOUBLE PRECISION,
!C   AND ALL OTHER FLOATING POINT VARIABLES ARE DECLARED DOUBLE PRECISION
!C
      INTEGER NMAX,RECO,NUM,N,NN,NJET,NSUB,JET(*),HIST,IMIN,JMIN,I,J
      PARAMETER (NMAX=512)
      DOUBLE PRECISION PP(4,*),PJET(4,*)
      DOUBLE PRECISION ECUT,P,KT,KTP,KTS,ETOT,RSQ,ETSQ,YCUT,YMAC,KTLAST, &
           ROUND
      PARAMETER (ROUND=0.99999D0)
      COMMON /KTCOMM/ETOT,RSQ,P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX), &
           KT(NMAX),KTLAST(NMAX),HIST(NMAX),NUM
!C---CHECK INPUT
      IF (RECO.LT.1.OR.RECO.GT.3) THEN
         PRINT *,'RECO=',RECO
         CALL KTWARN('KTRECO',100,*999)
      ENDIF
!C---COPY PP TO P
      N=NN
      IF (NUM.NE.NN) CALL KTWARN('KTRECO',101,*999)
      CALL KTCOPY(PP,N,P,(RECO.NE.1))
      IF (ECUT.EQ.0) THEN
         ETSQ=1/ETOT**2
      ELSE
         ETSQ=1/ECUT**2
      ENDIF
!C---KEEP MERGING UNTIL YCUT
 100  IF (ETSQ*KT(N).LT.ROUND*YCUT) THEN
         IF (HIST(N).LE.NMAX) THEN
            CALL KTMOVE(P,KTP,KTS,NMAX,N,HIST(N),0)
         ELSE
            IMIN=HIST(N)/NMAX
            JMIN=HIST(N)-IMIN*NMAX
            CALL KTMERG(P,KTP,KTS,NMAX,IMIN,JMIN,N,0,0,0,RECO)
            CALL KTMOVE(P,KTP,KTS,NMAX,N,JMIN,0)
         ENDIF
         N=N-1
         IF (N.GT.0) GOTO 100
      ENDIF
!C---IF YCUT IS TOO LARGE THERE ARE NO JETS
      NJET=N
      NSUB=N
      IF (N.EQ.0) RETURN
!C---SET UP OUTPUT MOMENTA
      DO I=1,NJET
         IF (RECO.EQ.1) THEN
            DO J=1,4
               PJET(J,I)=P(J,I)
            END DO
         ELSE
            PJET(1,I)=P(6,I)*COS(P(8,I))
            PJET(2,I)=P(6,I)*SIN(P(8,I))
            PJET(3,I)=P(6,I)*SINH(P(7,I))
            PJET(4,I)=P(6,I)*COSH(P(7,I))
         ENDIF
         JET(I)=I
      END DO
!C---KEEP MERGING UNTIL YMAC TO FIND THE FATE OF EACH JET
300   IF (ETSQ*KT(N).LT.ROUND*YMAC) THEN
         IF (HIST(N).LE.NMAX) THEN
            IMIN=0
            JMIN=HIST(N)
            NSUB=NSUB-1
         ELSE
            IMIN=HIST(N)/NMAX
            JMIN=HIST(N)-IMIN*NMAX
            IF (ETSQ*KTLAST(N).LT.ROUND*YMAC) NSUB=NSUB-1
         ENDIF
         DO I=1,NJET
            IF (JET(I).EQ.JMIN) JET(I)=IMIN
            IF (JET(I).EQ.N) JET(I)=JMIN
         END DO
         N=N-1
         IF (N.GT.0) GOTO 300
      ENDIF
      RETURN
999   RETURN 1
    END SUBROUTINE KTRECO
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
    FUNCTION KTPAIR(ANGL,P,Q,ANGLE)
      IMPLICIT NONE
!C---CALCULATE LOCAL KT OF PAIR, USING ANGULAR SCHEME:
!C   1=>ANGULAR, 2=>DeltaR, 3=>f(DeltaEta,DeltaPhi)
!C   WHERE f(eta,phi)=2(COSH(eta)-COS(phi)) IS THE QCD EMISSION METRIC
!C---IF ANGLE<0, IT IS SET TO THE ANGULAR PART OF THE LOCAL KT ON RETURN
!C   IF ANGLE>0, IT IS USED INSTEAD OF THE ANGULAR PART OF THE LOCAL KT
      INTEGER ANGL
! CHANGE      DOUBLE PRECISION P(9),Q(9),KTPAIR,R,KTMDPI,ANGLE,ETA,PHI,ESQ
      DOUBLE PRECISION P(9),Q(9),KTPAIR,R,ANGLE,ETA,PHI,ESQ
!C---COMPONENTS OF MOMENTA ARE PX,PY,PZ,E,1/P,PT,ETA,PHI,PT**2
      R=ANGLE
      IF (ANGL.EQ.1) THEN
         IF (R.LE.0) R=2*(1-(P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3))*(P(5)*Q(5)))
         ESQ=MIN(P(4),Q(4))**2
      ELSEIF (ANGL.EQ.2.OR.ANGL.EQ.3) THEN
         IF (R.LE.0) THEN
            ETA=P(7)-Q(7)
            PHI=KTMDPI(P(8)-Q(8))
            IF (ANGL.EQ.2) THEN
               R=ETA**2+PHI**2
            ELSE
               R=2*(COSH(ETA)-COS(PHI))
            ENDIF
         ENDIF
         ESQ=MIN(P(9),Q(9))
      ELSEIF (ANGL.EQ.4) THEN
         ESQ=(1d0/(P(5)*Q(5))-P(1)*Q(1)-P(2)*Q(2)- &
              P(3)*Q(3))*2D0/(P(5)*Q(5))/(0.0001D0+1d0/P(5)+1d0/Q(5))**2
         R=1d0
      ELSE
         CALL KTWARN('KTPAIR',200,*999)
         STOP
      ENDIF
      KTPAIR=ESQ*R
      IF (ANGLE.LT.0) ANGLE=R
999 END FUNCTION KTPAIR
!C-----------------------------------------------------------------------
    FUNCTION KTSING(ANGL,TYPE,P)
      IMPLICIT NONE
!C---CALCULATE KT OF PARTICLE, USING ANGULAR SCHEME:
!C   1=>ANGULAR, 2=>DeltaR, 3=>f(DeltaEta,DeltaPhi)
!C---TYPE=1 FOR E+E-, 2 FOR EP, 3 FOR PE, 4 FOR PP
!C   FOR EP, PROTON DIRECTION IS DEFINED AS -Z
!C   FOR PE, PROTON DIRECTION IS DEFINED AS +Z
      INTEGER ANGL,TYPE
      DOUBLE PRECISION P(9),KTSING,COSTH,R,SMALL
      DATA SMALL/1D-4/
      IF (ANGL.EQ.1.OR.ANGL.EQ.4) THEN
         COSTH=P(3)*P(5)
         IF (TYPE.EQ.2) THEN
            COSTH=-COSTH
         ELSEIF (TYPE.EQ.4) THEN
            COSTH=ABS(COSTH)
         ELSEIF (TYPE.NE.1.AND.TYPE.NE.3) THEN
            CALL KTWARN('KTSING',200,*999)
            STOP
         ENDIF
         R=2*(1-COSTH)
!C---IF CLOSE TO BEAM, USE APPROX 2*(1-COS(THETA))=SIN**2(THETA)
         IF (R.LT.SMALL) R=(P(1)**2+P(2)**2)*P(5)**2
         KTSING=P(4)**2*R
      ELSEIF (ANGL.EQ.2.OR.ANGL.EQ.3) THEN
         KTSING=P(9)
      ELSE
         CALL KTWARN('KTSING',201,*999)
         STOP
      ENDIF
999 END FUNCTION KTSING
!C-----------------------------------------------------------------------
    SUBROUTINE KTPMIN(A,NMAX,N,IMIN,JMIN)
      IMPLICIT NONE
!C---FIND THE MINIMUM MEMBER OF A(NMAX,NMAX) WITH IMIN < JMIN <= N
      INTEGER NMAX,N,IMIN,JMIN,KMIN,I,J,K
!C---REMEMBER THAT A(X+(Y-1)*NMAX)=A(X,Y)
!C   THESE LOOPING VARIABLES ARE J=Y-2, I=X+(Y-1)*NMAX
      DOUBLE PRECISION A(*),AMIN
      K=1+NMAX
      KMIN=K
      AMIN=A(KMIN)
      DO J=0,N-2
         DO I=K,K+J
            IF (A(I).LT.AMIN) THEN
               KMIN=I
               AMIN=A(KMIN)
            ENDIF
         END DO
         K=K+NMAX
      END DO
      JMIN=KMIN/NMAX+1
      IMIN=KMIN-(JMIN-1)*NMAX
    END SUBROUTINE KTPMIN
!C-----------------------------------------------------------------------
    SUBROUTINE KTSMIN(A,NMAX,N,IMIN)
      IMPLICIT NONE
!C---FIND THE MINIMUM MEMBER OF A
      INTEGER N,NMAX,IMIN,I
      DOUBLE PRECISION A(NMAX)
      IMIN=1
      DO I=1,N
         IF (A(I).LT.A(IMIN)) IMIN=I
      END DO
    END SUBROUTINE KTSMIN
!C-----------------------------------------------------------------------
    SUBROUTINE KTCOPY(A,N,B,ONSHLL)
      IMPLICIT NONE
!C---COPY FROM A TO B. 5TH=1/(3-MTM), 6TH=PT, 7TH=ETA, 8TH=PHI, 9TH=PT**2
!C   IF ONSHLL IS .TRUE. PARTICLE ENTRIES ARE PUT ON-SHELL BY SETTING E=P
      INTEGER I,N
      DOUBLE PRECISION A(4,N)
      LOGICAL ONSHLL
      DOUBLE PRECISION B(9,N),ETAMAX,SINMIN,EPS
      DATA ETAMAX,SINMIN,EPS/10,0,1D-6/
!C---SINMIN GETS CALCULATED ON FIRST CALL
      IF (SINMIN.EQ.0) SINMIN=1/COSH(ETAMAX)
      DO I=1,N
         B(1,I)=A(1,I)
         B(2,I)=A(2,I)
         B(3,I)=A(3,I)
         B(4,I)=A(4,I)
         B(5,I)=SQRT(A(1,I)**2+A(2,I)**2+A(3,I)**2)
         IF (ONSHLL) B(4,I)=B(5,I)
         IF (B(5,I).EQ.0) B(5,I)=1D-10
         B(5,I)=1/B(5,I)
         B(9,I)=A(1,I)**2+A(2,I)**2
         B(6,I)=SQRT(B(9,I))
         B(7,I)=B(6,I)*B(5,I)
         IF (B(7,I).GT.SINMIN) THEN
            B(7,I)=A(4,I)**2-A(3,I)**2
            IF (B(7,I).LE.EPS*B(4,I)**2.OR.ONSHLL) B(7,I)=B(9,I)
            B(7,I)=LOG((B(4,I)+ABS(B(3,I)))**2/B(7,I))/2
         ELSE
            B(7,I)=ETAMAX+2
         ENDIF
         B(7,I)=SIGN(B(7,I),B(3,I))
         IF (A(1,I).EQ.0 .AND. A(2,I).EQ.0) THEN
            B(8,I)=0
         ELSE
            B(8,I)=ATAN2(A(2,I),A(1,I))
         ENDIF
      END DO
    END SUBROUTINE KTCOPY
!C-----------------------------------------------------------------------
    SUBROUTINE KTMERG(P,KTP,KTS,NMAX,I,J,N,TYPE,ANGL,MONO,RECO)
      IMPLICIT NONE
!C---MERGE THE Jth PARTICLE IN P INTO THE Ith PARTICLE
!C   J IS ASSUMED GREATER THAN I. P CONTAINS N PARTICLES BEFORE MERGING.
!C---ALSO RECALCULATING THE CORRESPONDING KTP AND KTS VALUES IF MONO.GT.0
!C   FROM THE RECOMBINED ANGULAR MEASURES IF MONO.GT.1
!C---NOTE THAT IF MONO.LE.0, TYPE AND ANGL ARE NOT USED
      INTEGER ANGL,RECO,TYPE,I,J,K,N,NMAX,MONO
      DOUBLE PRECISION P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX),PT,PTT, &
! CHANGE           KTMDPI,KTUP,PI,PJ,ANG,KTPAIR,KTSING,ETAMAX,EPS
           KTUP,PI,PJ,ANG,ETAMAX,EPS
      KTUP(I,J)=KTP(MAX(I,J),MIN(I,J))
      DATA ETAMAX,EPS/10,1D-6/
      IF (J.LE.I) CALL KTWARN('KTMERG',200,*999)
!C---COMBINE ANGULAR MEASURES IF NECESSARY
      IF (MONO.GT.1) THEN
         DO K=1,N
            IF (K.NE.I.AND.K.NE.J) THEN
               IF (RECO.EQ.1) THEN
                  PI=P(4,I)
                  PJ=P(4,J)
               ELSEIF (RECO.EQ.2) THEN
                  PI=P(6,I)
                  PJ=P(6,J)
               ELSEIF (RECO.EQ.3) THEN
                  PI=P(9,I)
                  PJ=P(9,J)
               ELSE
                  CALL KTWARN('KTMERG',201,*999)
                  STOP
               ENDIF
               IF (PI.EQ.0.AND.PJ.EQ.0) THEN
                  PI=1
                  PJ=1
               ENDIF
               KTP(MAX(I,K),MIN(I,K))= &
                    (PI*KTUP(I,K)+PJ*KTUP(J,K))/(PI+PJ)
            ENDIF
         END DO
      ENDIF
      IF (RECO.EQ.1) THEN
!C---VECTOR ADDITION
         P(1,I)=P(1,I)+P(1,J)
         P(2,I)=P(2,I)+P(2,J)
         P(3,I)=P(3,I)+P(3,J)
!c         P(4,I)=P(4,I)+P(4,J) ! JA
         P(5,I)=SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)
         P(4,I)=P(5,I) ! JA (Massless scheme)
         IF (P(5,I).EQ.0) THEN
            P(5,I)=1
         ELSE
            P(5,I)=1/P(5,I)
         ENDIF
      ELSEIF (RECO.EQ.2) THEN
!C---PT WEIGHTED ETA-PHI ADDITION
         PT=P(6,I)+P(6,J)
         IF (PT.EQ.0) THEN
            PTT=1
         ELSE
            PTT=1/PT
         ENDIF
         P(7,I)=(P(6,I)*P(7,I)+P(6,J)*P(7,J))*PTT
         P(8,I)=KTMDPI(P(8,I)+P(6,J)*PTT*KTMDPI(P(8,J)-P(8,I)))
         P(6,I)=PT
         P(9,I)=PT**2
      ELSEIF (RECO.EQ.3) THEN
!C---PT**2 WEIGHTED ETA-PHI ADDITION
         PT=P(9,I)+P(9,J)
         IF (PT.EQ.0) THEN
            PTT=1
         ELSE
            PTT=1/PT
         ENDIF
         P(7,I)=(P(9,I)*P(7,I)+P(9,J)*P(7,J))*PTT
         P(8,I)=KTMDPI(P(8,I)+P(9,J)*PTT*KTMDPI(P(8,J)-P(8,I)))
         P(6,I)=P(6,I)+P(6,J)
         P(9,I)=P(6,I)**2
      ELSE
         CALL KTWARN('KTMERG',202,*999)
         STOP
      ENDIF
!C---IF MONO.GT.0 CALCULATE NEW KT MEASURES. IF MONO.GT.1 USE ANGULAR ONES.
      IF (MONO.LE.0) RETURN
!C---CONVERTING BETWEEN 4-MTM AND PT,ETA,PHI IF NECESSARY
      IF (ANGL.NE.1.AND.RECO.EQ.1) THEN
         P(9,I)=P(1,I)**2+P(2,I)**2
         P(7,I)=P(4,I)**2-P(3,I)**2
         IF (P(7,I).LE.EPS*P(4,I)**2) P(7,I)=P(9,I)
         IF (P(7,I).GT.0) THEN
            P(7,I)=LOG((P(4,I)+ABS(P(3,I)))**2/P(7,I))/2
            IF (P(7,I).GT.ETAMAX) P(7,I)=ETAMAX+2
         ELSE
            P(7,I)=ETAMAX+2
         ENDIF
         P(7,I)=SIGN(P(7,I),P(3,I))
         IF (P(1,I).NE.0.AND.P(2,I).NE.0) THEN
            P(8,I)=ATAN2(P(2,I),P(1,I))
         ELSE
            P(8,I)=0
         ENDIF
      ELSEIF (ANGL.EQ.1.AND.RECO.NE.1) THEN
         P(1,I)=P(6,I)*COS(P(8,I))
         P(2,I)=P(6,I)*SIN(P(8,I))
         P(3,I)=P(6,I)*SINH(P(7,I))
         P(4,I)=P(6,I)*COSH(P(7,I))
         IF (P(4,I).NE.0) THEN
            P(5,I)=1/P(4,I)
         ELSE
            P(5,I)=1
         ENDIF
      ENDIF
      ANG=0
      DO K=1,N
         IF (K.NE.I.AND.K.NE.J) THEN
            IF (MONO.GT.1) ANG=KTUP(I,K)
            KTP(MIN(I,K),MAX(I,K))= &
                 KTPAIR(ANGL,P(1,I),P(1,K),ANG)
         ENDIF
      END DO
      KTS(I)=KTSING(ANGL,TYPE,P(1,I))
999 END SUBROUTINE KTMERG
!C-----------------------------------------------------------------------
    SUBROUTINE KTMOVE(P,KTP,KTS,NMAX,N,J,IOPT)
      IMPLICIT NONE
!C---MOVE THE Nth PARTICLE IN P TO THE Jth POSITION
!C---ALSO MOVING KTP AND KTS IF IOPT.GT.0
      INTEGER I,J,N,NMAX,IOPT
      DOUBLE PRECISION P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX)
      DO I=1,9
         P(I,J)=P(I,N)
      END DO
      IF (IOPT.LE.0) RETURN
      DO I=1,J-1
         KTP(I,J)=KTP(I,N)
         KTP(J,I)=KTP(N,I)
      END DO
      DO I=J+1,N-1
         KTP(J,I)=KTP(I,N)
         KTP(I,J)=KTP(N,I)
      END DO
      KTS(J)=KTS(N)
    END SUBROUTINE KTMOVE
!C-----------------------------------------------------------------------
    FUNCTION KTMDPI(PHI)
      IMPLICIT NONE
!C---RETURNS PHI, MOVED ONTO THE RANGE [-PI,PI)
      DOUBLE PRECISION KTMDPI,PHI,PI,TWOPI,THRPI,EPS
      PARAMETER (PI=3.14159265358979324D0,TWOPI=6.28318530717958648D0, &
           THRPI=9.42477796076937972D0)
      PARAMETER (EPS=1D-15)
      KTMDPI=PHI
      IF (KTMDPI.LE.PI) THEN
        IF (KTMDPI.GT.-PI) THEN
          GOTO 100
        ELSEIF (KTMDPI.GT.-THRPI) THEN
          KTMDPI=KTMDPI+TWOPI
        ELSE
          KTMDPI=-MOD(PI-KTMDPI,TWOPI)+PI
        ENDIF
      ELSEIF (KTMDPI.LE.THRPI) THEN
        KTMDPI=KTMDPI-TWOPI
      ELSE
        KTMDPI=MOD(PI+KTMDPI,TWOPI)-PI
      ENDIF
 100  IF (ABS(KTMDPI).LT.EPS) KTMDPI=0
    END FUNCTION KTMDPI
!C-----------------------------------------------------------------------
    SUBROUTINE KTWARN(SUBRTN,ICODE,*)
!C     DEALS WITH ERRORS DURING EXECUTION
!C     SUBRTN = NAME OF CALLING SUBROUTINE
!C     ICODE  = ERROR CODE:    - 99 PRINT WARNING & CONTINUE
!C                          100-199 PRINT WARNING & JUMP
!C                          200-    PRINT WARNING & STOP DEAD
!C-----------------------------------------------------------------------
      INTEGER ICODE
      CHARACTER(len=6) SUBRTN
      WRITE (6,10) SUBRTN,ICODE
10    FORMAT(/' KTWARN CALLED FROM SUBPROGRAM ',A6,': CODE =',I4/)
      IF (ICODE.LT.100) RETURN
      IF (ICODE.LT.200) RETURN 1
      STOP
    END SUBROUTINE KTWARN
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------


end module ktclus
