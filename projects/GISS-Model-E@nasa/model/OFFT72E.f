      Module OFFT72
C****
!@sum   OFFT72 calculates the Fast Fourier Transform for KM = 72
!@auth  Gary L. Russell
!@ver   2008/12/29, KM=72
!@param KM length of input array; to change KM => rewrite module !!!
C****
      Implicit  None
      Integer*4,Parameter ::  KM = 72
      Real*8,   Parameter :: zKM = 1d0/KM, zKMH = 2d0/KM
      Real*8 C,S, C2,S2, C24Q0,C24Q1,S24Q1, C8QN,S8QN, C4QN,S4QN,
     *       C20N,C21N,S20N,S21N
      Common /OFFTCOM/
     *   C(KM), S(KM),  !   C(K) = cos(2*Pi*K/KM)
     *  C2(KM),S2(KM),  !  C2(K) = cos(2*Pi*K/KM) * 2
     *  C24Q0(0:23),C24Q1(0:23),S24Q1(0:23),
     *  C8QN(0:7,0:4),S8QN(0:7,4),
     *  C4QN(0:3,0:9),S4QN(0:3,9),
     *  C20N(0:18),C21N(0:17),S20N(17),S21N(18)
      EndModule OFFT72

      Subroutine OFFT0 (IM)
C****
!@sum   OFFT0 initializes sines and cosines used by OFFT routines
!@auth  Gary L. Russell
!@ver   2008/12/24, KM=72
C****
      Use OFFT72
      Implicit  None
      Real*8,   Parameter  :: TWOPI = 6.283185307179586477d0
      Integer*4,Intent(In) :: IM  !  @var IM size of arrays (must=KM)
      Integer*4 N
C****
      If (IM /= KM)  GoTo 800
      Do N=1,KM/4-1
         C(N) = Cos (TWOPI*N/KM)
         S(N) = Sin (TWOPI*N/KM)
         C(KM/2-N) = -C(N)
         C(KM/2+N) = -C(N)
         C(KM  -N) =  C(N)
         S(KM/2-N) =  S(N)
         S(KM/2+N) = -S(N)
         S(KM  -N) = -S(N)  ;  EndDo
         C(KM  /4) =  0
         C(KM*2/4) = -1
         C(KM*3/4) =  0
         C(KM)     =  1
         S(KM  /4) =  1
         S(KM*2/4) =  0
         S(KM*3/4) = -1
         S(KM)     =  0
C****
      Do N=1,KM
         C2(N) = C(N) * 2
         S2(N) = S(N) * 2  ;  EndDo
      Return
C****
  800 Write (6,*) 'This version of OFFT is for ',KM,'. IM =',IM
      Call STOP_MODEL ('stopped in OFFT72.f',255)
      EndSubroutine OFFT0

      Subroutine OFFT (F,A,B)
C****
!@sum   OFFT calculates fast fourier transform of an input array F
!@auth  Gary L. Russell
!@ver   2008/12/29, KM=72
C**** Input:  F   = grid point values
C**** Output: A,B = spectral coefficients
C****
C**** OFFT calculates a fast fourier transform of the input array F,
C**** producing the cosine and sine coefficients in the output
C**** arrays A and B.  F is dimensioned by 72 = KM, and A and B are
C**** dimensioned by 0:36.  Upon entering OFFT, the total energy is:
C****   .5*sum[F(K)^2]
C**** with the sum being taken over all K from 1 to KM.
C**** The Fourier coefficients are defined by:
C****   A(N) = sum [F(K)*cos(-2*Pii*N*K/KM)] / KMH
C****   B(N) = sum [F(K)*sin(-2*Pii*N*K/KM)] / KMH
C**** with the sum being taken over all K from 1 to KM.
C**** KMH = KM/2 when 0 < N < KM/2, and KMH = KM when N = 0 or KM/2.
C**** In the program's notation, CPQN and SPQN mean:
C****   CPQN = sum [F(K)*cos(2*Pi*N*K/KM)]
C****   SPQN = sum [F(K)*sin(2*Pi*N*K/KM)]
C**** with the sum being taken over all K = Q mod(P).
C**** The same total energy can be calculated from the spectral
C**** coefficients as:  .5 * sum [A(N)^2 + B(N)^2]*KMH
C**** with the sum being taken over all wave numbers N from 0 to KM/2.
C****
      Use OFFT72
      Implicit None
      Real*8,Intent(In)  :: F(KM)      !  @var F input gridpoint array
      Real*8,Intent(Out) :: A(0:KM/2)  !  @var A fourier coeffs. (cos)
      Real*8,Intent(Out) :: B(0:KM/2)  !  @var B fourier coeffs. (sin)
      Integer*4 N,Q
C****
C**** Factor of 3: from (P=72, Q=0:71, N=0) to (P=24, Q=0:23, N=0:1)
C**** In each section, the number of variables defined should be KM.
C****
      C24Q0(0) = (F(24)+F(48))       + F(72)
      C24Q1(0) = (F(24)+F(48))*C(24) + F(72)
      S24Q1(0) = (F(24)-F(48))*S(24)
      Do 10 Q=1,23
      C24Q0(Q) = F(Q)      + F(Q+24)         + F(Q+48)
      C24Q1(Q) = F(Q)*C(Q) + F(Q+24)*C(Q+24) + F(Q+48)*C(Q+48)
   10 S24Q1(Q) = F(Q)*S(Q) + F(Q+24)*S(Q+24) + F(Q+48)*S(Q+48)
C****
C**** Factor of 3: from (P=24, Q=0:23, N=0:1) to (P=8, Q=0:7, N=0:4)
C****
      C8QN(0,0) = C24Q0(0) +  C24Q0(8) + C24Q0(16)
      C8QN(0,1) = C24Q1(0) +  C24Q1(8) + C24Q1(16)
      S8QN(0,1) = S24Q1(0) +  S24Q1(8) + S24Q1(16)
      C8QN(0,3) = C24Q0(0) + (C24Q0(8) + C24Q0(16))*C(24)
      S8QN(0,3) =            (C24Q0(8) - C24Q0(16))*S(24)
      C8QN(0,2) = C24Q1(0) + (C24Q1(8) + C24Q1(16))*C(24)
     +                     + (S24Q1(8) - S24Q1(16))*S(24)
      C8QN(0,4) = C24Q1(0) + (C24Q1(8) + C24Q1(16))*C(24)
     -                     - (S24Q1(8) - S24Q1(16))*S(24)
      S8QN(0,2) =            (C24Q1(8) - C24Q1(16))*S(24)
     -          - S24Q1(0) - (S24Q1(8) + S24Q1(16))*C(24)
      S8QN(0,4) =            (C24Q1(8) - C24Q1(16))*S(24)
     +          + S24Q1(0) + (S24Q1(8) + S24Q1(16))*C(24)
      Do 20 Q=1,7
      C8QN(Q,0) = C24Q0(Q) + C24Q0(Q+8) + C24Q0(Q+16)
      C8QN(Q,1) = C24Q1(Q) + C24Q1(Q+8) + C24Q1(Q+16)
      S8QN(Q,1) = S24Q1(Q) + S24Q1(Q+8) + S24Q1(Q+16)
      C8QN(Q,3) = C24Q0(Q)   *C(Q*3)    + C24Q0(Q+8) *C(Q*3+24)
     +          + C24Q0(Q+16)*C(Q*3+48)
      S8QN(Q,3) = C24Q0(Q)   *S(Q*3)    + C24Q0(Q+8) *S(Q*3+24)
     +          + C24Q0(Q+16)*S(Q*3+48)
      C8QN(Q,2) = C24Q1(Q)   *C(Q*3)    + S24Q1(Q)   *S(Q*3)
     +          + C24Q1(Q+8) *C(Q*3+24) + S24Q1(Q+8) *S(Q*3+24)
     +          + C24Q1(Q+16)*C(Q*3+48) + S24Q1(Q+16)*S(Q*3+48)
      C8QN(Q,4) = C24Q1(Q)   *C(Q*3)    - S24Q1(Q)   *S(Q*3)
     +          + C24Q1(Q+8) *C(Q*3+24) - S24Q1(Q+8) *S(Q*3+24)
     +          + C24Q1(Q+16)*C(Q*3+48) - S24Q1(Q+16)*S(Q*3+48)
      S8QN(Q,2) = C24Q1(Q)   *S(Q*3)    - S24Q1(Q)   *C(Q*3)
     +          + C24Q1(Q+8) *S(Q*3+24) - S24Q1(Q+8) *C(Q*3+24)
     +          + C24Q1(Q+16)*S(Q*3+48) - S24Q1(Q+16)*C(Q*3+48)
   20 S8QN(Q,4) = C24Q1(Q)   *S(Q*3)    + S24Q1(Q)   *C(Q*3)
     +          + C24Q1(Q+8) *S(Q*3+24) + S24Q1(Q+8) *C(Q*3+24)
     +          + C24Q1(Q+16)*S(Q*3+48) + S24Q1(Q+16)*C(Q*3+48)
C****
C**** Factor of 2: from (P=8, Q=0:7, N=0:4) to (P=4, Q=0:3, N=0:9)
C****
      C4QN(0,0) = C8QN(0,0) + C8QN(4,0)
      C4QN(0,9) = C8QN(0,0) - C8QN(4,0)
C     S4QN(0,9) = 0
      Do 30 N=1,4
      C4QN(0,N)   = C8QN(0,N) + C8QN(4,N)
      C4QN(0,9-N) = C8QN(0,N) - C8QN(4,N)
      S4QN(0,N)   = S8QN(0,N) + S8QN(4,N)
   30 S4QN(0,9-N) = S8QN(4,N) - S8QN(0,N)
      Do 31 Q=1,3
      C4QN(Q,0) =  C8QN(Q,0)+C8QN(Q+4,0)
      C4QN(Q,9) = (C8QN(Q,0)-C8QN(Q+4,0))*C(Q*9)
      S4QN(Q,9) = (C8QN(Q,0)-C8QN(Q+4,0))*S(Q*9)  !  ?????
      Do 31 N=1,4
      C4QN(Q,N)   =  C8QN(Q,N)+C8QN(Q+4,N)
      C4QN(Q,9-N) = (C8QN(Q,N)-C8QN(Q+4,N))*C(Q*9)
     +            + (S8QN(Q,N)-S8QN(Q+4,N))*S(Q*9)
      S4QN(Q,N)   =  S8QN(Q,N)+S8QN(Q+4,N)
   31 S4QN(Q,9-N) = (C8QN(Q,N)-C8QN(Q+4,N))*S(Q*9)
     -            - (S8QN(Q,N)-S8QN(Q+4,N))*C(Q*9)
C****
C**** Factor of 2: from (P=4, Q=0:3, N=0:9) to (P=2, Q=0:1, N=0:18)
C****
      C20N(0)  = C4QN(0,0) + C4QN(2,0)
      C20N(18) = C4QN(0,0) - C4QN(2,0)
      C21N(0)  = C4QN(1,0) + C4QN(3,0)
C     C21N(18) = 0
C     S20N(18) = 0
      S21N(18) = C4QN(1,0) - C4QN(3,0)
      C20N(9)  = C4QN(0,9)              !  C4QN(2,9) = 0
      S20N(9)  = S4QN(2,9)              !  S4QN(0,9) = 0
      C21N(9)  = C4QN(1,9) + C4QN(3,9)
      S21N(9)  = C4QN(1,9) - C4QN(3,9)  !  = S4QN(1,9) + S4QN(3,9)
      Do 40 N=1,8
      C20N(N)    = C4QN(0,N) + C4QN(2,N)
      C20N(18-N) = C4QN(0,N) - C4QN(2,N)
      C21N(N)    = C4QN(1,N) + C4QN(3,N)
      S21N(18-N) = C4QN(1,N) - C4QN(3,N)
      S20N(N)    = S4QN(2,N) + S4QN(0,N)
      S20N(18-N) = S4QN(2,N) - S4QN(0,N)
      S21N(N)    = S4QN(1,N) + S4QN(3,N)
   40 C21N(18-N) = S4QN(1,N) - S4QN(3,N)
C****
C**** Factor of 2: from (P=2, Q=0:1, N=0:18) to (P=1, Q=0, N=0:36)
C****
      A(0)  = (C20N(0)+C21N(0))*zKM
      A(36) = (C20N(0)-C21N(0))*zKM
      B(0)  = 0
      B(36) = 0
      A(18) = C20N(18)*zKMH  !  C21N(18) = 0
      B(18) = S21N(18)*zKMH  !  S20N(18) = 0
      Do 50 N=1,17
      A(N)    = (C20N(N)+C21N(N))*zKMH
      A(36-N) = (C20N(N)-C21N(N))*zKMH
      B(N)    = (S20N(N)+S21N(N))*zKMH
   50 B(36-N) = (S21N(N)-S20N(N))*zKMH
      Return
      EndSubroutine OFFT

      Subroutine OFFTI (A,B,F)
C****
!@sum   OFFTI performs an inverse fast fourier transform
!@auth  Gary L. Russell
!@ver   2008/12/24, KM=72
C**** Input: A,B = spectral coefficients
C**** Output:  F = grid point values
C****
      Use OFFT72
      Implicit None
      Real*8,Intent(Out) :: F(KM)      !  @var F output gridpoint array
      Real*8,Intent(In)  :: A(0:KM/2)  !  @var A fourier coeffs. (cos)
      Real*8,Intent(In)  :: B(0:KM/2)  !  @var B fourier coeffs. (sin)
      Integer*4 N,Q
C****
C**** Factor of 2: from (P=1, Q=0, N=0:36) to (P=2, Q=0:1, N=0:18)
C**** New variables should be multiplied by: 1/2*zKMH = KM/4
C****
      C20N(0) = (A(0) + A(36))*2
      C21N(0) = (A(0) - A(36))*2
      C20N(18) = A(18)*2
      S21N(18) = B(18)*2
C     C21N(18) = 0
C     S20N(18) = 0
      Do 10 N=1,17
      S20N(N) = B(N) - B(36-N)
      S21N(N) = B(N) + B(36-N)
      C20N(N) = A(N) + A(36-N)
   10 C21N(N) = A(N) - A(36-N)
C****
C**** Factor of 2: from (P=2, Q=0:1, N=0:18) to (P=4, Q=0:3, N=0:9)
C**** New variables should be multiplied by: 1/4*zKMH = KM/8
C****
      C4QN(0,0)  = C20N(0) + C20N(18)
      C4QN(2,0)  = C20N(0) - C20N(18)
      C4QN(1,0)  = C21N(0) + S21N(18)
      C4QN(3,0)  = C21N(0) - S21N(18)
      C4QN(0,9) = C20N(9)*2
C     C4QN(2,9) = 0
      C4QN(1,9) = C21N(9) + S21N(9)
      C4QN(3,9) = C21N(9) - S21N(9)
C     S4QN(0,9) = 0
      S4QN(2,9) = S20N(9)*2
C     S4QN(1,9) = S21N(9) + C21N(9)  !  = C4QN(1,9)
C     S4QN(3,9) = S21N(9) - C21N(9)  !  = - C4QN(3,9)
      Do 20 N=1,8
      C4QN(0,N) = C20N(N) + C20N(18-N)
      C4QN(2,N) = C20N(N) - C20N(18-N)
      C4QN(1,N) = C21N(N) + S21N(18-N)
      C4QN(3,N) = C21N(N) - S21N(18-N)
      S4QN(0,N) = S20N(N) - S20N(18-N)
      S4QN(2,N) = S20N(N) + S20N(18-N)
      S4QN(1,N) = S21N(N) + C21N(18-N)
   20 S4QN(3,N) = S21N(N) - C21N(18-N)
C****
C**** Factor of 2: from (P=4, Q=0:3, N=0:9) to (P=8, Q=0:7, N=0:4)
C**** New variables should be multiplied by: 1/8*zKMH = KM/16
C****
      C8QN(0,0) = C4QN(0,0) + C4QN(0,9)
      C8QN(4,0) = C4QN(0,0) - C4QN(0,9)
      C8QN(2,0) = C4QN(2,0) + S4QN(2,9)
      C8QN(6,0) = C4QN(2,0) - S4QN(2,9)
      C8QN(1,0) = C4QN(1,0) + C4QN(1,9)*S2(9)
      C8QN(5,0) = C4QN(1,0) - C4QN(1,9)*S2(9)
      C8QN(3,0) = C4QN(3,0) - C4QN(3,9)*S2(9)
      C8QN(7,0) = C4QN(3,0) + C4QN(3,9)*S2(9)
      Do 30 N=1,4
      C8QN(0,N) = C4QN(0,N) +  C4QN(0,9-N)
      C8QN(4,N) = C4QN(0,N) -  C4QN(0,9-N)
      C8QN(2,N) = C4QN(2,N) +  S4QN(2,9-N)
      C8QN(6,N) = C4QN(2,N) -  S4QN(2,9-N)
      C8QN(1,N) = C4QN(1,N) + (C4QN(1,9-N)+S4QN(1,9-N))*S(9)
      C8QN(5,N) = C4QN(1,N) - (C4QN(1,9-N)+S4QN(1,9-N))*S(9)
      C8QN(3,N) = C4QN(3,N) - (C4QN(3,9-N)-S4QN(3,9-N))*S(9)
      C8QN(7,N) = C4QN(3,N) + (C4QN(3,9-N)-S4QN(3,9-N))*S(9)
      S8QN(0,N) = S4QN(0,N) -  S4QN(0,9-N)
      S8QN(4,N) = S4QN(0,N) +  S4QN(0,9-N)
      S8QN(2,N) = S4QN(2,N) +  C4QN(2,9-N)
      S8QN(6,N) = S4QN(2,N) -  C4QN(2,9-N)
      S8QN(1,N) = S4QN(1,N) + (C4QN(1,9-N)-S4QN(1,9-N))*S(9)
      S8QN(5,N) = S4QN(1,N) - (C4QN(1,9-N)-S4QN(1,9-N))*S(9)
      S8QN(3,N) = S4QN(3,N) + (C4QN(3,9-N)+S4QN(3,9-N))*S(9)
   30 S8QN(7,N) = S4QN(3,N) - (C4QN(3,9-N)+S4QN(3,9-N))*S(9)
C****
C**** Factor of 3: from (P=8, Q=0:7, N=0:4) to (P=24, Q=0:23, N=0:1)
C**** New variables should be multiplied by: 1/24*zKMH = KM/48
C****
      C24Q0(0)  = C8QN(0,0) +  C8QN(0,3)*2
      C24Q0(8)  = C8QN(0,0) +  C8QN(0,3)*C2(24) + S8QN(0,3)*S2(24)
      C24Q0(16) = C8QN(0,0) +  C8QN(0,3)*C2(24) - S8QN(0,3)*S2(24)
      C24Q1(0)  = C8QN(0,1) + (C8QN(0,2)+C8QN(0,4))
      S24Q1(0)  = S8QN(0,1) - (S8QN(0,2)-S8QN(0,4))
      C24Q1(8)  = C8QN(0,1) + (C8QN(0,2)+C8QN(0,4))*C(24)
     +                      + (S8QN(0,2)+S8QN(0,4))*S(24)
      C24Q1(16) = C8QN(0,1) + (C8QN(0,2)+C8QN(0,4))*C(24)
     +                      - (S8QN(0,2)+S8QN(0,4))*S(24)
      S24Q1(8)  = S8QN(0,1) + (C8QN(0,2)-C8QN(0,4))*S(24)
     +                      - (S8QN(0,2)-S8QN(0,4))*C(24)
      S24Q1(16) = S8QN(0,1) - (C8QN(0,2)-C8QN(0,4))*S(24)
     +                      - (S8QN(0,2)-S8QN(0,4))*C(24)
      Do 40 Q=1,7
      C24Q0(Q)    = C8QN(Q,0) +  C8QN(Q,3)*C2(Q*3)
     +                        +  S8QN(Q,3)*S2(Q*3)
      C24Q0(Q+8)  = C8QN(Q,0) +  C8QN(Q,3)*C2(Q*3+24)
     +                        +  S8QN(Q,3)*S2(Q*3+24)
      C24Q0(Q+16) = C8QN(Q,0) +  C8QN(Q,3)*C2(Q*3+48)
     +                        +  S8QN(Q,3)*S2(Q*3+48)
      C24Q1(Q)    = C8QN(Q,1) + (C8QN(Q,2)+C8QN(Q,4))*C(Q*3)
     +                        + (S8QN(Q,2)+S8QN(Q,4))*S(Q*3)
      C24Q1(Q+8)  = C8QN(Q,1) + (C8QN(Q,2)+C8QN(Q,4))*C(Q*3+24)
     +                        + (S8QN(Q,2)+S8QN(Q,4))*S(Q*3+24)
      C24Q1(Q+16) = C8QN(Q,1) + (C8QN(Q,2)+C8QN(Q,4))*C(Q*3+48)
     +                        + (S8QN(Q,2)+S8QN(Q,4))*S(Q*3+48)
      S24Q1(Q)    = S8QN(Q,1) + (C8QN(Q,2)-C8QN(Q,4))*S(Q*3)
     +                        - (S8QN(Q,2)-S8QN(Q,4))*C(Q*3)
      S24Q1(Q+8)  = S8QN(Q,1) + (C8QN(Q,2)-C8QN(Q,4))*S(Q*3+24)
     +                        - (S8QN(Q,2)-S8QN(Q,4))*C(Q*3+24)
   40 S24Q1(Q+16) = S8QN(Q,1) + (C8QN(Q,2)-C8QN(Q,4))*S(Q*3+48)
     +                        - (S8QN(Q,2)-S8QN(Q,4))*C(Q*3+48)
C****
C**** Factor of 3: from (P=24, Q=0:23, N=0:1) to (P=72, Q=0:71, N=0)
C**** New variables should be multiplied by: 2/72*zKMH = KM/72 = 1
C****
      F(72) = .5*C24Q0(0) + C24Q1(0)
      F(24) = .5*C24Q0(0) + C24Q1(0)*C(24) + S24Q1(0)*S(24)
      F(48) = .5*C24Q0(0) + C24Q1(0)*C(48) + S24Q1(0)*S(48)
      Do 50 Q=1,23
      F(Q)    = .5*C24Q0(Q) + C24Q1(Q)*C(Q)    + S24Q1(Q)*S(Q)
      F(Q+24) = .5*C24Q0(Q) + C24Q1(Q)*C(Q+24) + S24Q1(Q)*S(Q+24)
   50 F(Q+48) = .5*C24Q0(Q) + C24Q1(Q)*C(Q+48) + S24Q1(Q)*S(Q+48)
      Return
      EndSubroutine OFFTI

      Subroutine OFFTE (F,E)
C****
!@sum   OFFTE calcs. spectral energy E from input gridpoint values F
!@auth  Gary L. Russell
!@ver   2008/12/24, KM=72
C****
      Use OFFT72
      Implicit None
      Real*8,Intent(In)  :: F(KM)      !  @var F input gridpoint array
      Real*8,Intent(Out) :: E(0:KM/2)  !  @var E spectral energy
      Integer*4 N
      Real*8    A(0:KM/2),B(0:KM/2)
C****
      Call OFFT (F,A,B)
      E(0) = .5*KM * A(0)**2
      Do 10 N=1,KM/2-1
   10 E(N) = .25*KM * (A(N)**2 + B(N)**2)
      E(KM/2) = .5*KM * A(KM/2)**2
      Return
      EndSubroutine OFFTE
