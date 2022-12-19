! $Id: ct10pdf.f90 7177 2015-08-26 16:02:41Z jr_reuter $

!
! ---- FORTRAN 90 ----
! 
! The contents of this file are identical to CT10Pdf.f with a
! small hack to facilitate custom PDF table search paths. In addition,
! the code has been wrapped in a module to avoid symbol collisions (in
! particular with LHApdf). For the same reason, the common blocks have
! been postfixed with _ct10_WO_. ct10pdf has been renamed to getct10pdf.
! A trim has been added when complaining about missing table files.
!
!============================================================================
!                CTEQ-TEA Parton Distribution Functions: version 2010
!                             April 23, 2010, v1.00
!                             July 13, 2010, v1.01
!                             July 26, 2010, v1.02
!                             August 3, 2010, v1.03
!                             August 9, 2010, v1.04
!
!   Ref[1]: "New parton distributions for collider physics"
!      By : H.L. Lai, M. Guzzi, J. Huston, Z. Li, P. Nadolsky, J. Pumplin, C.-P. Yuan
!           arXiv:1007.2241 (hep-ph)
!
!   This package contains
!   (1) 1+52 sets of CT10 PDF's (4/23);
!   (2) 1+52 sets of CT10W PDF's (4/23);
!   (3) 10 sets of CT10 alternative alpha_s PDF's (7/13, 7/26);
!   (4) 10 sets of CT10W alternative alpha_s PDF's (7/13, 7/26);
!   (5) 4 sets of CT10 & CT10W in Fixed Flavor Scheme PDF's (8/3).

!  Details about the calling convention are:
! --------------------------------------------------------------------------------
!  Iset   PDF-set     Description           Alpha_s(Mz) Table_File     Ref
! ================================================================================
!  100    CT10.00     Central CT10           0.118      ct10.00.pds    [1]
!  1xx    CT10.xx     +/- sets               0.118      ct10.xx.pds    [1]
!        where xx = 01-52: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
! --------------------------
!  200    CT10W.00    Central CT10W          0.118      ct10w00.pds    [1]
!  2xx    CT10W.xx    +/- sets               0.118      ct10wxx.pds    [1]
!        where xx = 01-52: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
! ------------------------------------------------------------------------
!  Alpha_s range recommended to estimate the uncertainty of PDF+alpha_s
!   10    CT10.AS0    CT10 alpha_s series    0.116      ct10.as0.pds   [1]
!   11    CT10.AS1    CT10 alpha_s series    0.117      ct10.as1.pds   [1]
!   12    CT10.AS2    CT10 alpha_s series    0.119      ct10.as2.pds   [1]
!   13    CT10.AS3    CT10 alpha_s series    0.120      ct10.as3.pds   [1]
! --------------------------
!   20    CT10W.AS0   CT10W alpha_s series   0.116      ct10was0.pds   [1]
!   21    CT10W.AS1   CT10W alpha_s series   0.117      ct10was1.pds   [1]
!   22    CT10W.AS2   CT10W alpha_s series   0.119      ct10was2.pds   [1]
!   23    CT10W.AS3   CT10W alpha_s series   0.120      ct10was3.pds   [1]
! ------------------------------------------------------------------------
!  Extended alpha_s range for further exploration
!   14    CT10.AS4    CT10 alpha_s series    0.113      ct10.as4.pds   [1]
!   15    CT10.AS5    CT10 alpha_s series    0.114      ct10.as5.pds   [1]
!   16    CT10.AS6    CT10 alpha_s series    0.115      ct10.as6.pds   [1]
!   17    CT10.AS7    CT10 alpha_s series    0.121      ct10.as7.pds   [1]
!   18    CT10.AS8    CT10 alpha_s series    0.122      ct10.as8.pds   [1]
!   19    CT10.AS9    CT10 alpha_s series    0.123      ct10.as9.pds   [1]
! --------------------------
!   24    CT10W.AS4   CT10W alpha_s series   0.113      ct10was4.pds   [1]
!   25    CT10W.AS5   CT10W alpha_s series   0.114      ct10was5.pds   [1]
!   26    CT10W.AS6   CT10W alpha_s series   0.115      ct10was6.pds   [1]
!   27    CT10W.AS7   CT10W alpha_s series   0.121      ct10was7.pds   [1]
!   28    CT10W.AS8   CT10W alpha_s series   0.122      ct10was8.pds   [1]
!   29    CT10W.AS9   CT10W alpha_s series   0.123      ct10was9.pds   [1]
! ------------------------------------------------------------------------
!  Fixed Flavor Scheme
!   30    CT10.3F     CT10 3-flavor          0.118      ct10.3f.pds    [1]
!   31    CT10.4F     CT10 4-flavor          0.118      ct10.4f.pds    [1]
!   32    CT10W.3F    CT10W 3-flavor         0.118      ct10w3f.pds    [1]
!   33    CT10W.4F    CT10W 4-flavor         0.118      ct10w4f.pds    [1]
! ===========================================================================
!   ** ALL fits are obtained by using the same coupling strength
!   \alpha_s(Mz)=0.118 and the NLO (HOPPET) running \alpha_s formula.
!
!   The table grids are generated for 
!    *  10^-8 < x < 1 and 1.3 < Q < 10^5 (GeV) for CT10 & CT10W series;
!
!   PDF values outside of the above range are returned using extrapolation.
!
!   The Table_Files are assumed to be in the working directory.
!
!   Before using the PDF, it is necessary to do the initialization by
!       Call SetCT10(Iset)
!   where Iset is the desired PDF specified in the above table.
!
!   The function CT10Pdf (Iparton, X, Q)
!   returns the parton distribution inside the proton for parton [Iparton]
!   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
!   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
!                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
!
!   For detailed information on the parameters used, e.q. quark masses,
!   AlphaS, ... etc.,  see info lines at the beginning of the
!   Table_Files.
!
!   These programs, as provided, are in double precision.  By removing the
!   "Implicit Double Precision" lines, they can also be run in single
!   precision.
!
!   If you have detailed questions concerning these CT10 distributions,
!   or if you find problems/bugs using this package, direct inquires to
!   hllai@tmue.edu.tw, pumplin@pa.msu.edu or nadolsky@physics.smu.edu.
!
!===========================================================================

module ct10pdf

contains

  Function getCT10Pdf (Iparton, X, Q)
    Implicit Double Precision (A-H,O-Z)
    Logical Warn
    Common / CtqPar2_ct10_WO_ / Nx, Nt, NfMx, MxVal
    Common / QCDtbl_ct10_WO_ /  AlfaQ, Qalfa, Ipk, Iorder, Nfl !for external use

    Data Warn /.true./
    Data Qsml /.3d0/
    save Warn

    If (X .lt. 0d0 .or. X .gt. 1D0) Then
       Print *, 'X out of range in getCT10Pdf: ', X
       getCt10Pdf = 0D0
       Return
    End if

    If (Q .lt. Qsml) Then
       Print *, 'Q out of range in getCT10Pdf: ', Q
       Stop
    End if

    If (abs(Iparton).gt. NfMx) Then
       If (Warn) Then
          !        print a warning for calling extra flavor
          Warn = .false.
          Print *, 'Warning: Iparton out of range in getCT10Pdf! '
          Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
       End if
       getCT10Pdf = 0D0
    else

       getCT10Pdf = PartonX10 (Iparton, X, Q)
       if (getCT10Pdf.lt.0D0) getCT10Pdf = 0D0
    end if                     !if (abs(Iparton...

    Return

    !                             ********************
  end function getCT10Pdf

  Subroutine SetCT10 (prefix, Iset)
    Implicit Double Precision (A-H,O-Z)
    Parameter (Isetmax0=2)
    Character Flnm(Isetmax0)*5, nn*3, Tablefile*40
    character(*) prefix
    Data (Flnm(I), I=1,Isetmax0) / 'ct10.', 'ct10w' /
    Data Isetold, Isetmin1, Isetmax1 /-987,100,152/
    Data Isetmin2,Isetmax2 /200,252/
    Data IsetASmn1,IsetASmx1 /10,19/
    Data IsetASmn2,IsetASmx2 /20,29/
    Data IsetFFSmn,IsetFFSmx /30,33/
    Common /Setchange_ct10_WO_/ Isetch
    Common /CT10Jset_ct10_WO_/ Jset
    save

    Jset=Iset
    ! If data file not initialized, do so.
    If(Iset.ne.Isetold) then

       If (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
          !                                                     101 - 152
          write(nn,'(I3)') Iset
          Tablefile=Flnm(1)//nn(2:3)//'.pds'
       Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
          !                                                     200 - 252
          write(nn,'(I3)') Iset
          Tablefile=Flnm(2)//nn(2:3)//'.pds'
       Elseif (Iset.ge.IsetASmn1 .and. Iset.le.IsetASmx1) Then
          !                                                       10 - 19
          write(nn,'(I2)') Iset
          Tablefile=Flnm(1)//'as'//nn(2:2)//'.pds'
       Elseif (Iset.ge.IsetASmn2 .and. Iset.le.IsetASmx2) Then
          !                                                       20 - 29
          write(nn,'(I2)') Iset
          Tablefile=Flnm(2)//'as'//nn(2:2)//'.pds'
       Elseif (Iset.ge.IsetFFSmn .and. Iset.le.IsetFFSmx) Then
          !                                                       30 - 33
          Js=(Iset-28)/2
          write(nn,'(I1)') Iset-25-Js*2
          Tablefile=Flnm(Js)//nn(1:1)//'f.pds'
       Else
          Print *, 'Invalid Iset number in SetCT10 :', Iset
          Stop
       End if
       IU= NextUn()
       Open(IU, File=prefix // "/" // Tablefile, Status='OLD', Err=100)
21     Call Readpds0 (IU)
       Close (IU)
       Isetold=Iset
       Isetch=1
    End if
    Return

100 Print *, ' Data file ', prefix // "/" // trim (Tablefile), &
         ' cannot be opened in SetCT10!!'
    Stop
    ! ********************
  End subroutine SetCT10

  Subroutine Readpds0 (Nu)
    Implicit Double Precision (A-H,O-Z)
    Character Line*80
    PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
    PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
    Common / CtqPar1_ct10_WO_ / qBase, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
    Common / CtqPar2_ct10_WO_ / Nx, Nt, NfMx, MxVal
    Common / XQrange_ct10_WO_ / Qini, Qmax, Xmin
    Common / Masstbl_ct10_WO_ / Amass(6)
    Common / QCDtbl_ct10_WO_ /  AlfaQ, Qalfa, Ipk, Iorder, Nfl !for external use

    Read  (Nu, '(A)') Line
    Read  (Nu, '(A)') Line
    Read  (Nu, *) Dr, Fl, qBase, (Amass(I),I=1,6)
    Iorder = Nint(Dr)
    Nfl = Nint(Fl)

    Read  (Nu, '(A)') Line
    Read  (Nu, *) Ipk,AlfaQ,Qalfa, NfMx, MxVal, N0

    Read  (Nu, '(A)') Line
    Read  (Nu, *) NX,  NT, N0, NG, N0

    Read  (Nu, '(A)') (Line,I=1,NG+2)
    Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

    Read  (Nu, '(A)') Line
    Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
    XV(0)=0D0

    Nblk = (NX+1) * (NT+1)
    Npts =  Nblk  * (NfMx+1+MxVal)
    Read  (Nu, '(A)') Line
    Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

    Return
    !  ****************************
  end subroutine Readpds0

  Function PartonX10 (IPRTN, XX, QQ)
!  Given the parton distribution function in the array U in
!  COMMON / PEVLDT / , this routine interpolates to find
!  the parton distribution at an arbitray point in x and q.
    Implicit Double Precision (A-H,O-Z)
    
    PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
    PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

    Common / CtqPar1_ct10_WO_ / qBase, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
    Common / CtqPar2_ct10_WO_ / Nx, Nt, NfMx, MxVal
    Common / XQrange_ct10_WO_ / Qini, Qmax, Xmin
    Common /Setchange_ct10_WO_/ Isetch

    Dimension fvec(4), fij(4)
    Dimension xvpow(0:mxx)
    Data OneP / 1.00001 /
    Data xpow / 0.3d0 /       !**** choice of interpolation variable
    Data nqvec / 4 /
    Data ientry / 0 /
    Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
    Save xvpow
    Save X, Q, JX, JQ, JLX, JLQ
    Save ss, const1, const2, const3, const4, const5, const6
    Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
    Save tmp1, tmp2, tdet

! store the powers used for interpolation on first call...
    if(Isetch .eq. 1) then
       Isetch = 0

       xvpow(0) = 0D0
       do i = 1, nx
          xvpow(i) = xv(i)**xpow
       end do
    else if((XX.eq.X).and.(QQ.eq.Q)) then
       goto 99
    endif

    X = XX
    Q = QQ
    tt = log(log(Q/qBase))

!      -------------    find lower end of interval containing x, i.e.,
!                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
    JLx = -1
    JU = Nx+1
11  If (JU-JLx .GT. 1) Then
       JM = (JU+JLx) / 2
       If (X .Ge. XV(JM)) Then
          JLx = JM
       Else
          JU = JM
       Endif
       Goto 11
    End if
!                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
!                           |---|---|---|...|---|-x-|---|...|---|---|
!                     x     0  Xmin               x                 1
!
    If     (JLx .LE. -1) Then
       Print '(A,1pE12.4)','Severe error: x <= 0 in PartonX10! x = ',x
       Stop
    ElseIf (JLx .Eq. 0) Then
       Jx = 0
    Elseif (JLx .LE. Nx-2) Then

! For interrior points, keep x in the middle, as shown above
       Jx = JLx - 1
    Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

!                  We tolerate a slight over-shoot of one (OneP=1.00001),
!              perhaps due to roundoff or whatever, but not more than that.
!                                      Keep at least 4 points >= Jx
       Jx = JLx - 2
    Else
       Print '(A,1pE12.4)','Severe error: x > 1 in PartonX10! x = ',x
       Stop
    Endif
!          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.
!                       This is the variable to be interpolated in
    ss = x**xpow
    
    If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

!     initiation work for "interior bins": store the lattice points in s...
       svec1 = xvpow(jx)
       svec2 = xvpow(jx+1)
       svec3 = xvpow(jx+2)
       svec4 = xvpow(jx+3)

       s12 = svec1 - svec2
       s13 = svec1 - svec3
       s23 = svec2 - svec3
       s24 = svec2 - svec4
       s34 = svec3 - svec4

       sy2 = ss - svec2
       sy3 = ss - svec3

! constants needed for interpolating in s at fixed t lattice points...
       const1 = s13/s23
       const2 = s12/s23
       const3 = s34/s23
       const4 = s24/s23
       s1213 = s12 + s13
       s2434 = s24 + s34
       sdet = s12*s34 - s1213*s2434
       tmp = sy2*sy3/sdet
       const5 = (s34*sy2-s2434*sy3)*tmp/s12
       const6 = (s1213*sy2-s12*sy3)*tmp/s34

    End If

!         --------------Now find lower end of interval containing Q, i.e.,
!                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
    JLq = -1
    JU = NT+1
12  If (JU-JLq .GT. 1) Then
       JM = (JU+JLq) / 2
       If (tt .GE. TV(JM)) Then
          JLq = JM
       Else
          JU = JM
       End if
       Goto 12
    End if

    If (JLq .LE. 0) Then
       Jq = 0
    Else if (JLq .LE. Nt-2) Then
! keep q in the middle, as shown above
       Jq = JLq - 1
    Else
! JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
       Jq = Nt - 3
       
    End if
! This is the interpolation variable in Q
    
    If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
!                                        store the lattice points in t...
       tvec1 = Tv(jq)
       tvec2 = Tv(jq+1)
       tvec3 = Tv(jq+2)
       tvec4 = Tv(jq+3)

       t12 = tvec1 - tvec2
       t13 = tvec1 - tvec3
       t23 = tvec2 - tvec3
       t24 = tvec2 - tvec4
       t34 = tvec3 - tvec4

       ty2 = tt - tvec2
       ty3 = tt - tvec3

       tmp1 = t12 + t13
       tmp2 = t24 + t34

       tdet = t12*t34 - tmp1*tmp2

    End If


! get the pdf function values at the lattice points...

99  If (Iprtn .Gt. MxVal) Then
       Ip = - Iprtn
    Else
       Ip = Iprtn
    End If
    jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

    Do it = 1, nqvec

       J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
!                          For the first 4 x points, interpolate x^2*f(x,Q)
!                           This applies to the two lowest bins JLx = 0, 1
!            We can not put the JLx.eq.1 bin into the "interrior" section
!                           (as we do for q), since Upd(J1) is undefined.
          fij(1) = 0
          fij(2) = Upd(J1+1) * XV(1)**2
          fij(3) = Upd(J1+2) * XV(2)**2
          fij(4) = Upd(J1+3) * XV(3)**2
!
!                 Use Polint which allows x to be anywhere w.r.t. the grid

          Call Polint4F (XVpow(0), Fij(1), ss, Fx)
          
          If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
!                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
!                                                This is the highest x bin:

          Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

          Fvec(it) = Fx
          
       Else
!                       for all interior points, use Jon's in-line function
!                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
          sf2 = Upd(J1+1)
          sf3 = Upd(J1+2)

          g1 =  sf2*const1 - sf3*const2
          g4 = -sf2*const3 + sf3*const4

          Fvec(it) = (const5*(Upd(J1)-g1) &
               + const6*(Upd(J1+3)-g4) &
               + sf2*sy3 - sf3*sy2) / s23

       End if

    end do
!                                   We now have the four values Fvec(1:4)
!     interpolate in t...

    If (JLq .LE. 0) Then
!                         1st Q-bin, as well as extrapolation to lower Q
       Call Polint4F (TV(0), Fvec(1), tt, ff)

    Else If (JLq .GE. Nt-1) Then
!                         Last Q-bin, as well as extrapolation to higher Q
       Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
    Else
!                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
!       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
!                         the full range QV(0:Nt)  (in contrast to XV)
       tf2 = fvec(2)
       tf3 = fvec(3)

       g1 = ( tf2*t13 - tf3*t12) / t23
       g4 = (-tf2*t34 + tf3*t24) / t23

       h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 &
            +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

       ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
    End If

    PartonX10 = ff

    Return
    ! ********************
  End function PartonX10

  SUBROUTINE POLINT4F (XA,YA,X,Y)

    IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
!  but assuming N=4, and ignoring the error estimation.
!  suggested by Z. Sullivan.
    DIMENSION XA(*),YA(*)

    H1=XA(1)-X
    H2=XA(2)-X
    H3=XA(3)-X
    H4=XA(4)-X

    W=YA(2)-YA(1)
    DEN=W/(H1-H2)
    D1=H2*DEN
    C1=H1*DEN
    
    W=YA(3)-YA(2)
    DEN=W/(H2-H3)
    D2=H3*DEN
    C2=H2*DEN

    W=YA(4)-YA(3)
    DEN=W/(H3-H4)
    D3=H4*DEN
    C3=H3*DEN

    W=C2-D1
    DEN=W/(H1-H3)
    CD1=H3*DEN
    CC1=H1*DEN

    W=C3-D2
    DEN=W/(H2-H4)
    CD2=H4*DEN
    CC2=H2*DEN

    W=CC2-CD1
    DEN=W/(H1-H4)
    DD1=H4*DEN
    DC1=H1*DEN

    If((H3+H4).lt.0D0) Then
       Y=YA(4)+D3+CD2+DD1
    Else if((H2+H3).lt.0D0) Then
       Y=YA(3)+D2+CD1+DC1
    Else if((H1+H2).lt.0D0) Then
       Y=YA(2)+C2+CD1+DC1
    ELSE
       Y=YA(1)+C1+CC1+DC1
    END IF

    RETURN
    ! *************************
  END subroutine POLINT4F

  Function NextUn()
    ! Returns an unallocated FORTRAN i/o unit.
    Logical ::  EX

    Do N = 10, 300
       INQUIRE (UNIT=N, OPENED=EX)
       If (.NOT. EX) then
          NextUn = N
          Return
       End if
    end Do
    Stop ' There is no available I/O unit. '
    ! *************************
  End function NextUn

end module ct10pdf

