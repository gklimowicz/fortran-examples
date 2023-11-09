      Module OFFT144
C****                                                  2008/12/11
!@sum   OFFT144 calculates the Fast Fourier Transform for KM = 144
!@auth  Gary L. Russell
!@ver   2.0
!@param KM length of input array; to change KM => rewrite module !!!
C****
      Implicit  None
      Integer*4,Parameter ::  KM = 144
      Real*8,   Parameter :: zKM = 1d0/KM, zKMH = 2d0/KM
      Real*8 C,S, C2,S2, C48Q0,C48Q1,S48Q1, C16QN,S16QN, C8QN,S8QN,
     *  C4QN,S4QN, C20N,C21N,S20N,S21N
      Common /OFFTCOM/
     *   C(KM), S(KM),  !   C(K) = cos(2*Pi*K/KM)
     *  C2(KM),S2(KM),  !  C2(K) = cos(2*Pi*K/KM) * 2
     *  C48Q0(0:47),C48Q1(0:47),S48Q1(0:47),
     *  C16QN(0:15,0:4),S16QN(0:15,4),
     *  C8QN(0:7,0:9),S8QN(0:7,9),
     *  C4QN(0:3,0:18),S4QN(0:3,18),
     *  C20N(0:36),C21N(0:35),S20N(35),S21N(36)
      EndModule OFFT144

      Subroutine OFFT0 (IM)
C****
!@sum   OFFT0 initializes sines and cosines used by OFFT routines
!@auth  Gary L. Russell
!@ver   2.0
C****
      Use OFFT144
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
      Call STOP_MODEL ('stopped in OFFT144.f',255)
      EndSubroutine OFFT0

      Subroutine OFFT (F,A,B)
C****
!@sum   OFFT calculates fast fourier transform of an input array F
!@auth  Gary L. Russell
!@ver   2.0
C**** Input:  F   = grid point values
C**** Output: A,B = spectral coefficients
C****
C**** OFFT calculates a fast fourier transform of the input array F,
C**** producing the cosine and sine coefficients in the output
C**** arrays A and B.  F is dimensioned by 144 = KM, and A and B are
C**** dimensioned by 0:72.  Upon entering OFFT, the total energy is:
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
      Use OFFT144
      Implicit None
      Real*8,Intent(In)  :: F(KM)      !  @var F input gridpoint array
      Real*8,Intent(Out) :: A(0:KM/2)  !  @var A fourier coeffs. (cos)
      Real*8,Intent(Out) :: B(0:KM/2)  !  @var B fourier coeffs. (sin)
      Integer*4 N,Q
C****
C**** Factor of 3: from (P=144, Q=0:143, N=0) to (P=48, Q=0:47, N=0:1)
C**** In each section, the number of variables defined should be KM.
C****
      C48Q0(0) = (F(48)+F(96))       + F(144)
      C48Q1(0) = (F(48)+F(96))*C(48) + F(144)
      S48Q1(0) = (F(48)-F(96))*S(48)
      Do 10 Q=1,47
      C48Q0(Q) = F(Q)      + F(Q+48)         + F(Q+96)
      C48Q1(Q) = F(Q)*C(Q) + F(Q+48)*C(Q+48) + F(Q+96)*C(Q+96)
   10 S48Q1(Q) = F(Q)*S(Q) + F(Q+48)*S(Q+48) + F(Q+96)*S(Q+96)
C****
C**** Factor of 3: from (P=48, Q=0:47, N=0:1) to (P=16, Q=0:15, N=0:4)
C****
      C16QN(0,0) = C48Q0(0) + C48Q0(16) + C48Q0(32)
      C16QN(0,1) = C48Q1(0) + C48Q1(16) + C48Q1(32)
      S16QN(0,1) = S48Q1(0) + S48Q1(16) + S48Q1(32)
      C16QN(0,3) = C48Q0(0) + (C48Q0(16)+C48Q0(32))*C(48)
      S16QN(0,3) =            (C48Q0(16)-C48Q0(32))*S(48)
      C16QN(0,2) = C48Q1(0) + (C48Q1(16)+C48Q1(32))*C(48)
     +                      + (S48Q1(16)-S48Q1(32))*S(48)
      C16QN(0,4) = C48Q1(0) + (C48Q1(16)+C48Q1(32))*C(48)
     -                      - (S48Q1(16)-S48Q1(32))*S(48)
      S16QN(0,2) =            (C48Q1(16)-C48Q1(32))*S(48)
     -           - S48Q1(0) - (S48Q1(16)+S48Q1(32))*C(48)
      S16QN(0,4) =            (C48Q1(16)-C48Q1(32))*S(48)
     +           + S48Q1(0) + (S48Q1(16)+S48Q1(32))*C(48)
      Do 20 Q=1,15
      C16QN(Q,0) = C48Q0(Q) + C48Q0(Q+16) + C48Q0(Q+32)
      C16QN(Q,1) = C48Q1(Q) + C48Q1(Q+16) + C48Q1(Q+32)
      S16QN(Q,1) = S48Q1(Q) + S48Q1(Q+16) + S48Q1(Q+32)
      C16QN(Q,3) = C48Q0(Q)   *C(Q*3)    + C48Q0(Q+16)*C(Q*3+48)
     +           + C48Q0(Q+32)*C(Q*3+96)
      S16QN(Q,3) = C48Q0(Q)   *S(Q*3)    + C48Q0(Q+16)*S(Q*3+48)
     +           + C48Q0(Q+32)*S(Q*3+96)
      C16QN(Q,2) = C48Q1(Q)   *C(Q*3)    + S48Q1(Q)   *S(Q*3)
     +           + C48Q1(Q+16)*C(Q*3+48) + S48Q1(Q+16)*S(Q*3+48)
     +           + C48Q1(Q+32)*C(Q*3+96) + S48Q1(Q+32)*S(Q*3+96)
      C16QN(Q,4) = C48Q1(Q)   *C(Q*3)    - S48Q1(Q)   *S(Q*3)
     +           + C48Q1(Q+16)*C(Q*3+48) - S48Q1(Q+16)*S(Q*3+48)
     +           + C48Q1(Q+32)*C(Q*3+96) - S48Q1(Q+32)*S(Q*3+96)
      S16QN(Q,2) = C48Q1(Q)   *S(Q*3)    - S48Q1(Q)   *C(Q*3)
     +           + C48Q1(Q+16)*S(Q*3+48) - S48Q1(Q+16)*C(Q*3+48)
     +           + C48Q1(Q+32)*S(Q*3+96) - S48Q1(Q+32)*C(Q*3+96)
   20 S16QN(Q,4) = C48Q1(Q)   *S(Q*3)    + S48Q1(Q)   *C(Q*3)
     +           + C48Q1(Q+16)*S(Q*3+48) + S48Q1(Q+16)*C(Q*3+48)
     +           + C48Q1(Q+32)*S(Q*3+96) + S48Q1(Q+32)*C(Q*3+96)
C****
C**** Factor of 2: from (P=16, Q=0:15, N=0:4) to (P=8, Q=0:7, N=0:9)
C****
      C8QN(0,0) = C16QN(0,0) + C16QN(8,0)
      C8QN(0,9) = C16QN(0,0) - C16QN(8,0)
      Do 30 N=1,4
      C8QN(0,N)   = C16QN(0,N) + C16QN(8,N)
      C8QN(0,9-N) = C16QN(0,N) - C16QN(8,N)
      S8QN(0,N)   = S16QN(0,N) + S16QN(8,N)
   30 S8QN(0,9-N) = S16QN(8,N) - S16QN(0,N)
      Do 31 Q=1,7
      C8QN(Q,0) =  C16QN(Q,0)+C16QN(Q+8,0)
      C8QN(Q,9) = (C16QN(Q,0)-C16QN(Q+8,0))*C(Q*9)
      S8QN(Q,9) = (C16QN(Q,0)-C16QN(Q+8,0))*S(Q*9)  !  ?????
      Do 31 N=1,4
      C8QN(Q,N)   =  C16QN(Q,N)+C16QN(Q+8,N)
      C8QN(Q,9-N) = (C16QN(Q,N)-C16QN(Q+8,N))*C(Q*9)
     +            + (S16QN(Q,N)-S16QN(Q+8,N))*S(Q*9)
      S8QN(Q,N)   =  S16QN(Q,N)+S16QN(Q+8,N)
   31 S8QN(Q,9-N) = (C16QN(Q,N)-C16QN(Q+8,N))*S(Q*9)
     -            - (S16QN(Q,N)-S16QN(Q+8,N))*C(Q*9)
C****
C**** Factor of 2: from (P=8, Q=0:7, N=0:9) to (P=4, Q=0:3, N=0:18)
C****
      C4QN(0,0)  =  C8QN(0,0) + C8QN(4,0)
      C4QN(0,18) =  C8QN(0,0) - C8QN(4,0)
      C4QN(1,0)  =  C8QN(1,0) + C8QN(5,0)
      C4QN(1,18) = (C8QN(1,0) - C8QN(5,0))*C(18)
      C4QN(2,0)  =  C8QN(2,0) + C8QN(6,0)
C     C4QN(2,18) = 0
C     S4QN(0,18) = 0
C     S4QN(1,18) = (C8QN(1,0) - C8QN(5,0))*C(18)  !  = C4QN(1,18)
      S4QN(2,18) =  C8QN(2,0) - C8QN(6,0)
C     S4QN(3,18) = (C8QN(3,0) - C8QN(7,0))*C(18)  !  = - C4QN(3,18)
      C4QN(3,0)  =  C8QN(7,0) + C8QN(3,0)
      C4QN(3,18) = (C8QN(7,0) - C8QN(3,0))*C(18)
      Do 40 Q=0,3
      C4QN(Q,9) = C8QN(Q,9) + C8QN(Q+4,9)
   40 S4QN(Q,9) = S8QN(Q,9) + S8QN(Q+4,9)  !  ?????
      Do 41 N=1,8
      C4QN(0,N)    =  C8QN(0,N)+C8QN(4,N)
      C4QN(0,18-N) =  C8QN(0,N)-C8QN(4,N)
      C4QN(2,N)    =  C8QN(2,N)+C8QN(6,N)
      S4QN(2,18-N) =  C8QN(2,N)-C8QN(6,N)
      S4QN(0,N)    =  S8QN(4,N)+S8QN(0,N)
      S4QN(0,18-N) =  S8QN(4,N)-S8QN(0,N)
      S4QN(2,N)    =  S8QN(2,N)+S8QN(6,N)
      C4QN(2,18-N) =  S8QN(2,N)-S8QN(6,N)
      C4QN(1,N)    =  C8QN(1,N)+C8QN(5,N)
      S4QN(1,N)    =                       S8QN(1,N)+S8QN(5,N)
      C4QN(1,18-N) = (C8QN(1,N)-C8QN(5,N)+(S8QN(1,N)-S8QN(5,N)))*C(18)
      S4QN(1,18-N) = (C8QN(1,N)-C8QN(5,N)-(S8QN(1,N)-S8QN(5,N)))*C(18)
      S4QN(3,N)    =  S8QN(3,N)+S8QN(7,N)
      C4QN(3,N)    =                       C8QN(3,N)+C8QN(7,N)
      C4QN(3,18-N) = (S8QN(3,N)-S8QN(7,N)-(C8QN(3,N)-C8QN(7,N)))*C(18)
   41 S4QN(3,18-N) = (S8QN(3,N)-S8QN(7,N)+(C8QN(3,N)-C8QN(7,N)))*C(18)
C****
C**** Factor of 2: from (P=4, Q=0:3, N=0:18) to (P=2, Q=0:1, N=0:36)
C****
      C20N(0)  = C4QN(0,0) + C4QN(2,0)
      C20N(36) = C4QN(0,0) - C4QN(2,0)
      C21N(0)  = C4QN(1,0) + C4QN(3,0)
C     C21N(36) = 0
C     S20N(36) = 0
      S21N(36) = C4QN(1,0) - C4QN(3,0)
      C20N(18) = C4QN(0,18)               !  C4QN(2,18) = 0
      S20N(18) = S4QN(2,18)               !  S4QN(0,18) = 0
      C21N(18) = C4QN(1,18) + C4QN(3,18)
      S21N(18) = C4QN(1,18) - C4QN(3,18)  !  = S4QN(1,18) + S4QN(3,18)
      Do 50 N=1,17
      C20N(N)    = C4QN(0,N) + C4QN(2,N)
      C20N(36-N) = C4QN(0,N) - C4QN(2,N)
      C21N(N)    = C4QN(1,N) + C4QN(3,N)
      S21N(36-N) = C4QN(1,N) - C4QN(3,N)
      S20N(N)    = S4QN(2,N) + S4QN(0,N)
      S20N(36-N) = S4QN(2,N) - S4QN(0,N)
      S21N(N)    = S4QN(1,N) + S4QN(3,N)
   50 C21N(36-N) = S4QN(1,N) - S4QN(3,N)
C****
C**** Factor of 2: from (P=2, Q=0:1, N=0:36) to (P=1, Q=0, N=0:72)
C****
      A(0)  = (C20N(0)+C21N(0))*zKM
      A(72) = (C20N(0)-C21N(0))*zKM
      B(0)  = 0
      B(72) = 0
      A(36) = C20N(36)*zKMH  !  C21N(36) = 0
      B(36) = S21N(36)*zKMH  !  S20N(36) = 0
      Do 60 N=1,35
      A(N)    = (C20N(N)+C21N(N))*zKMH
      A(72-N) = (C20N(N)-C21N(N))*zKMH
      B(N)    = (S20N(N)+S21N(N))*zKMH
   60 B(72-N) = (S21N(N)-S20N(N))*zKMH
      Return
      EndSubroutine OFFT

      Subroutine OFFTI (A,B,F)
C****
!@sum   OFFTI performs an inverse fast fourier transform
!@auth  Gary L. Russell
!@ver   2.0 (KM=144)
C**** Input: A,B = spectral coefficients
C**** Output:  F = grid point values
C****
      Use OFFT144
      Implicit None
      Real*8,Intent(Out) :: F(KM)      !  @var F output gridpoint array
      Real*8,Intent(In)  :: A(0:KM/2)  !  @var A fourier coeffs. (cos)
      Real*8,Intent(In)  :: B(0:KM/2)  !  @var B fourier coeffs. (sin)
      Integer*4 N,Q
C****
C**** Factor of 2: from (P=1, Q=0, N=0:72) to (P=2, Q=0:1, N=0:36)
C**** New variables should be multiplied by: 1/2*zKMH = KM/4
C****
      C20N(0) = (A(0) + A(72))*2
      C21N(0) = (A(0) - A(72))*2
      C20N(36) = A(36)*2
      S21N(36) = B(36)*2
C     C21N(36) = 0
C     S20N(36) = 0
      Do 10 N=1,35
      S20N(N) = B(N) - B(72-N)
      S21N(N) = B(N) + B(72-N)
      C20N(N) = A(N) + A(72-N)
   10 C21N(N) = A(N) - A(72-N)
C****
C**** Factor of 2: from (P=2, Q=0:1, N=0:36) to (P=4, Q=0:3, N=0:18)
C**** New variables should be multiplied by: 1/4*zKMH = KM/8
C****
      C4QN(0,0)  = C20N(0) + C20N(36)
      C4QN(2,0)  = C20N(0) - C20N(36)
      C4QN(1,0)  = C21N(0) + S21N(36)
      C4QN(3,0)  = C21N(0) - S21N(36)
      C4QN(0,18) = C20N(18)*2
C     C4QN(2,18) = 0
      C4QN(1,18) = C21N(18) + S21N(18)
      C4QN(3,18) = C21N(18) - S21N(18)
C     S4QN(0,18) = 0
      S4QN(2,18) = S20N(18)*2
C     S4QN(1,18) = S21N(18) + C21N(18)  !  = C4QN(1,18)
C     S4QN(3,18) = S21N(18) - C21N(18)  !  = - C4QN(3,18)
      Do 20 N=1,17
      C4QN(0,N) = C20N(N) + C20N(36-N)
      C4QN(2,N) = C20N(N) - C20N(36-N)
      C4QN(1,N) = C21N(N) + S21N(36-N)
      C4QN(3,N) = C21N(N) - S21N(36-N)
      S4QN(0,N) = S20N(N) - S20N(36-N)
      S4QN(2,N) = S20N(N) + S20N(36-N)
      S4QN(1,N) = S21N(N) + C21N(36-N)
   20 S4QN(3,N) = S21N(N) - C21N(36-N)
C****
C**** Factor of 2: from (P=4, Q=0:3, N=0:18) to (P=8, Q=0:7, N=0:9)
C**** New variables should be multiplied by: 1/8*zKMH = KM/16
C****
      C8QN(0,0) = C4QN(0,0) + C4QN(0,18)
      C8QN(4,0) = C4QN(0,0) - C4QN(0,18)
      C8QN(2,0) = C4QN(2,0) + S4QN(2,18)
      C8QN(6,0) = C4QN(2,0) - S4QN(2,18)
      C8QN(1,0) = C4QN(1,0) + C4QN(1,18)*S2(18)
      C8QN(5,0) = C4QN(1,0) - C4QN(1,18)*S2(18)
      C8QN(3,0) = C4QN(3,0) - C4QN(3,18)*S2(18)
      C8QN(7,0) = C4QN(3,0) + C4QN(3,18)*S2(18)
      C8QN(0,9) = C4QN(0,9)*2
C     C8QN(4,9) = 0
      C8QN(2,9) = C4QN(2,9) +  S4QN(2,9)
      C8QN(6,9) = C4QN(2,9) -  S4QN(2,9)
      C8QN(1,9) = C4QN(1,9) + (C4QN(1,9)+S4QN(1,9))*S(18)
      C8QN(5,9) = C4QN(1,9) - (C4QN(1,9)+S4QN(1,9))*S(18)
      C8QN(3,9) = C4QN(3,9) - (C4QN(3,9)-S4QN(3,9))*S(18)
      C8QN(7,9) = C4QN(3,9) + (C4QN(3,9)-S4QN(3,9))*S(18)
C     S8QN(0,9) = 0
      S8QN(4,9) = S4QN(0,9)*2
      S8QN(2,9) = S4QN(2,9) +  C4QN(2,9)
      S8QN(6,9) = S4QN(2,9) -  C4QN(2,9)
      S8QN(1,9) = S4QN(1,9) + (C4QN(1,9)-S4QN(1,9))*S(18)
      S8QN(5,9) = S4QN(1,9) - (C4QN(1,9)-S4QN(1,9))*S(18)
      S8QN(3,9) = S4QN(3,9) + (C4QN(3,9)+S4QN(3,9))*S(18)
      S8QN(7,9) = S4QN(3,9) - (C4QN(3,9)+S4QN(3,9))*S(18)
      Do 30 N=1,8
      C8QN(0,N) = C4QN(0,N) +  C4QN(0,18-N)
      C8QN(4,N) = C4QN(0,N) -  C4QN(0,18-N)
      C8QN(2,N) = C4QN(2,N) +  S4QN(2,18-N)
      C8QN(6,N) = C4QN(2,N) -  S4QN(2,18-N)
      C8QN(1,N) = C4QN(1,N) + (C4QN(1,18-N)+S4QN(1,18-N))*S(18)
      C8QN(5,N) = C4QN(1,N) - (C4QN(1,18-N)+S4QN(1,18-N))*S(18)
      C8QN(3,N) = C4QN(3,N) - (C4QN(3,18-N)-S4QN(3,18-N))*S(18)
      C8QN(7,N) = C4QN(3,N) + (C4QN(3,18-N)-S4QN(3,18-N))*S(18)
      S8QN(0,N) = S4QN(0,N) -  S4QN(0,18-N)
      S8QN(4,N) = S4QN(0,N) +  S4QN(0,18-N)
      S8QN(2,N) = S4QN(2,N) +  C4QN(2,18-N)
      S8QN(6,N) = S4QN(2,N) -  C4QN(2,18-N)
      S8QN(1,N) = S4QN(1,N) + (C4QN(1,18-N)-S4QN(1,18-N))*S(18)
      S8QN(5,N) = S4QN(1,N) - (C4QN(1,18-N)-S4QN(1,18-N))*S(18)
      S8QN(3,N) = S4QN(3,N) + (C4QN(3,18-N)+S4QN(3,18-N))*S(18)
   30 S8QN(7,N) = S4QN(3,N) - (C4QN(3,18-N)+S4QN(3,18-N))*S(18)
C****
C**** Factor of 2: from (P=8, Q=0:7, N=0:9) to (P=16, Q=0:15, N=0:4)
C**** New variables should be multiplied by: 1/16*zKMH = KM/32
C****
      C16QN( 0,0) = C8QN(0,0) +  C8QN(0,9)
      C16QN( 8,0) = C8QN(0,0) -  C8QN(0,9)
      C16QN( 4,0) = C8QN(4,0) +  S8QN(4,9)
      C16QN(12,0) = C8QN(4,0) -  S8QN(4,9)
      C16QN( 2,0) = C8QN(2,0) + (S8QN(2,9)+C8QN(2,9))*S(18)
      C16QN(10,0) = C8QN(2,0) - (S8QN(2,9)+C8QN(2,9))*S(18)
      C16QN( 6,0) = C8QN(6,0) + (S8QN(6,9)-C8QN(6,9))*S(18)
      C16QN(14,0) = C8QN(6,0) - (S8QN(6,9)-C8QN(6,9))*S(18)
      Do 40 N=1,4
      C16QN( 0,N) = C8QN(0,N) +  C8QN(0,9-N)
      C16QN( 8,N) = C8QN(0,N) -  C8QN(0,9-N)
      S16QN( 0,N) = S8QN(0,N) -  S8QN(0,9-N)
      S16QN( 8,N) = S8QN(0,N) +  S8QN(0,9-N)
      C16QN( 4,N) = C8QN(4,N) +  S8QN(4,9-N)
      C16QN(12,N) = C8QN(4,N) -  S8QN(4,9-N)
      S16QN( 4,N) = S8QN(4,N) +  C8QN(4,9-N)
      S16QN(12,N) = S8QN(4,N) -  C8QN(4,9-N)
      C16QN( 2,N) = C8QN(2,N) + (C8QN(2,9-N)+S8QN(2,9-N))*S(18)
      C16QN(10,N) = C8QN(2,N) - (C8QN(2,9-N)+S8QN(2,9-N))*S(18)
      S16QN( 2,N) = S8QN(2,N) + (C8QN(2,9-N)-S8QN(2,9-N))*S(18)
      S16QN(10,N) = S8QN(2,N) - (C8QN(2,9-N)-S8QN(2,9-N))*S(18)
      C16QN( 6,N) = C8QN(6,N) - (C8QN(6,9-N)-S8QN(6,9-N))*S(18)
      C16QN(14,N) = C8QN(6,N) + (C8QN(6,9-N)-S8QN(6,9-N))*S(18)
      S16QN( 6,N) = S8QN(6,N) + (C8QN(6,9-N)+S8QN(6,9-N))*S(18)
   40 S16QN(14,N) = S8QN(6,N) - (C8QN(6,9-N)+S8QN(6,9-N))*S(18)
      Do 41 Q=1,7,2
      C16QN(Q  ,0) = C8QN(Q,0) + C8QN(Q,9)*C(Q*9) + S8QN(Q,9)*S(Q*9)
      C16QN(Q+8,0) = C8QN(Q,0) - C8QN(Q,9)*C(Q*9) - S8QN(Q,9)*S(Q*9)
      Do 41 N=1,4
      C16QN(Q  ,N) = C8QN(Q,N) + C8QN(Q,9-N)*C(Q*9) + S8QN(Q,9-N)*S(Q*9)
      C16QN(Q+8,N) = C8QN(Q,N) - C8QN(Q,9-N)*C(Q*9) - S8QN(Q,9-N)*S(Q*9)
      S16QN(Q  ,N) = S8QN(Q,N) + C8QN(Q,9-N)*S(Q*9) - S8QN(Q,9-N)*C(Q*9)
   41 S16QN(Q+8,N) = S8QN(Q,N) - C8QN(Q,9-N)*S(Q*9) + S8QN(Q,9-N)*C(Q*9)
C****
C**** Factor of 3: from (P=16, Q=0:15, N=0:4) to (P=48, Q=0:47, N=0:1)
C**** New variables should be multiplied by: 1/48*zKMH = KM/96
C****
      C48Q0(0)  = C16QN(0,0) +  C16QN(0,3)*2
      C48Q0(16) = C16QN(0,0) +  C16QN(0,3)*C2(48) + S16QN(0,3)*S2(48)
      C48Q0(32) = C16QN(0,0) +  C16QN(0,3)*C2(48) - S16QN(0,3)*S2(48)
      C48Q1(0)  = C16QN(0,1) + (C16QN(0,2)+C16QN(0,4))
      S48Q1(0)  = S16QN(0,1) - (S16QN(0,2)-S16QN(0,4))
      C48Q1(16) = C16QN(0,1) + (C16QN(0,2)+C16QN(0,4))*C(48)
     +                       + (S16QN(0,2)+S16QN(0,4))*S(48)
      C48Q1(32) = C16QN(0,1) + (C16QN(0,2)+C16QN(0,4))*C(48)
     +                       - (S16QN(0,2)+S16QN(0,4))*S(48)
      S48Q1(16) = S16QN(0,1) + (C16QN(0,2)-C16QN(0,4))*S(48)
     +                       - (S16QN(0,2)-S16QN(0,4))*C(48)
      S48Q1(32) = S16QN(0,1) - (C16QN(0,2)-C16QN(0,4))*S(48)
     +                       - (S16QN(0,2)-S16QN(0,4))*C(48)
      Do 50 Q=1,15
      C48Q0(Q)    = C16QN(Q,0) +  C16QN(Q,3)*C2(Q*3)
     +                         +  S16QN(Q,3)*S2(Q*3)
      C48Q0(Q+16) = C16QN(Q,0) +  C16QN(Q,3)*C2(Q*3+48)
     +                         +  S16QN(Q,3)*S2(Q*3+48)
      C48Q0(Q+32) = C16QN(Q,0) +  C16QN(Q,3)*C2(Q*3+96)
     +                         +  S16QN(Q,3)*S2(Q*3+96)
      C48Q1(Q)    = C16QN(Q,1) + (C16QN(Q,2)+C16QN(Q,4))*C(Q*3)
     +                         + (S16QN(Q,2)+S16QN(Q,4))*S(Q*3)
      C48Q1(Q+16) = C16QN(Q,1) + (C16QN(Q,2)+C16QN(Q,4))*C(Q*3+48)
     +                         + (S16QN(Q,2)+S16QN(Q,4))*S(Q*3+48)
      C48Q1(Q+32) = C16QN(Q,1) + (C16QN(Q,2)+C16QN(Q,4))*C(Q*3+96)
     +                         + (S16QN(Q,2)+S16QN(Q,4))*S(Q*3+96)
      S48Q1(Q)    = S16QN(Q,1) + (C16QN(Q,2)-C16QN(Q,4))*S(Q*3)
     +                         - (S16QN(Q,2)-S16QN(Q,4))*C(Q*3)
      S48Q1(Q+16) = S16QN(Q,1) + (C16QN(Q,2)-C16QN(Q,4))*S(Q*3+48)
     +                         - (S16QN(Q,2)-S16QN(Q,4))*C(Q*3+48)
   50 S48Q1(Q+32) = S16QN(Q,1) + (C16QN(Q,2)-C16QN(Q,4))*S(Q*3+96)
     +                         - (S16QN(Q,2)-S16QN(Q,4))*C(Q*3+96)
C****
C**** Factor of 3: from (P=48, Q=0:47, N=0:1) to (P=144, Q=0:143, N=0)
C**** New variables should be multiplied by: 2/144*zKMH = KM/144 = 1
C****
      F(144) = .5*C48Q0(0) + C48Q1(0)
      F(48)  = .5*C48Q0(0) + C48Q1(0)*C(48) + S48Q1(0)*S(48)
      F(96)  = .5*C48Q0(0) + C48Q1(0)*C(96) + S48Q1(0)*S(96)
      Do 60 Q=1,47
      F(Q)    = .5*C48Q0(Q) + C48Q1(Q)*C(Q)    + S48Q1(Q)*S(Q)
      F(Q+48) = .5*C48Q0(Q) + C48Q1(Q)*C(Q+48) + S48Q1(Q)*S(Q+48)
   60 F(Q+96) = .5*C48Q0(Q) + C48Q1(Q)*C(Q+96) + S48Q1(Q)*S(Q+96)
      Return
      EndSubroutine OFFTI

      Subroutine OFFTE (F,E)
C****
!@sum   OFFTE calcs. spectral energy E from input gridpoint values F
!@auth  Gary L. Russell
!@ver   2.0
C****
      Use OFFT144
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
