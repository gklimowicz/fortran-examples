!  This is the set of user routines which can be used to read the CJ PDF
!  tables and extract values of the CJ PDFs at specified values of
!  momentum fraction x and factorization scale Q.
!
!
!  Patterned after the CTEQ6 PDF routines
!
!  J.F. Owens   June 2012 - December 2012
!  A. Accardi   July 2015
!
!  CJ12:
!  Three sets of PDFs are available: CJ12_min, CJ12_mid, and CJ12_max
!  corresponding to the minimum, midddle, and maximum nuclear corrections.
!
!  CJ15:
!  Two sets available: CJ15_LO and CJ15_NLO
!
!  The following numbering scheme is used for the variable ISET
!
!  ISET             PDF
!
!  100      | CJ12_min central PDF
!  101-138  | CJ12_min error PDF sets
!           |
!  200      | CJ12_mid central PDF
!  201-238  | CJ12_mid error PDF sets
!           |
!  300      | CJ12_max central PDF
!  301-338  | CJ12_max error PDF sets
!
!  400      | CJ15lo central PDF
!  401-448  | CJ15lo error PDF sets
!
!  500      | CJ15nlo central PDF
!  501-548  | CJ15nlo error PDF sets
!
!  The tables cover the x range 10^-6 < x < 1. and Q range 1.3 < Q < 10^5 GeV.
!  Values outside these ranges must be considered as extrapolations.
!
!  The user should initalize the PDF set by first calling SetCJ(ISET)
!
!  A function call to CJpdf(Iptn,x,Q) returns the number density PDF of
!  flavor Iptn at momentum fraction x and factorization scale Q
!
!  A subroutine call to CJpdf_all(x,Q,bb,cb,sb,db,ub,g,u,d,s,c,b) returns all
!  the PDFs with one call
!
!  Flavor index: Iptn = -5:5 for bb, cb, sb, db, ub, g, u, d, s, c, b,
!
!  For the CJ12 PDFs b=bb, c=cb, and s=sb
!
!  An example of how to use the CJ PDFs is included in the program tst_CJpdf.f
!  This will produce a table of PDFs for a specified value of ISET
!

module CJ_pdf
  use io_units

contains

  Function CJpdf(Iptn, x, Q)
    Implicit double precision (A-H,O-Z)
    Common /CJPAR2/ Nx, Nq, Nfmx
    Common /QCDtable/ Alam, Nfl, Iord

    If(x.lt.0.d0.or.x.gt.1.d0)then
       Print*,'x out of range in CJPDF:',x
       Stop
    Endif

    If(Q.lt.Alam)then
       Print*,'Q out of range in CJpdf:',Q
       Stop
    Endif

    If(Iptn.lt.-Nfmx.or.Iptn.gt.Nfmx)then
       Print*,'Iptn out of range in CJpdf:',Iptn
       Print*,'Maximum number of flavors is:',Nfmx
       Stop
    Endif
    CJpdf=Getpdf(Iptn,x,Q)
    ! if(CJpdf.lt.0.d0)CJPDF=0.d0

    Return
  End function CJpdf

  subroutine CJpdf_all(x,Q,bb,cb,sb,db,ub,g,u,d,s,c,b)
    Implicit Double Precision (a-h,o-z)
    bb=CJpdf(-5,x,Q)
    cb=CJpdf(-4,x,Q)
    sb=CJpdf(-3,x,Q)
    db=CJpdf(-2,x,Q)
    ub=CJpdf(-1,x,Q)
    g=CJpdf(0,x,Q)
    u=CJpdf(1,x,Q)
    d=CJpdf(2,x,Q)
    s=sb
    c=cb
    b=bb
    return
  end subroutine CJpdf_all

  Subroutine SetCJ (prefix, Iset)
    Implicit Double Precision (A-H,O-Z)
    Character Flnm(5)*9, nn*3, Tablefile*100
    character (*) prefix
    Data (Flnm(I), I=1,5) &
         / 'CJ12_min_', 'CJ12_mid_', &
           'CJ12_max_', 'CJ15LO_', 'CJ15NLO_' /
    Data Isetold/-1/
    save

    ! If data file not initialized, do so.
    If(Iset.ne.Isetold) then

       IU= NextUn()
       If (Iset.ge.100 .and. Iset.le.140) Then
          write(nn,'(I3)') Iset
          Tablefile=trim(Flnm(1))//nn(2:3)//'.tbl'
       Elseif (Iset.ge.200 .and. Iset.le.240) Then
          write(nn,'(I3)') Iset
          Tablefile=trim(Flnm(2))//nn(2:3)//'.tbl'
       Elseif (Iset.ge.300 .and. Iset.le.340) Then
          write(nn,'(I3)') Iset
          Tablefile=trim(Flnm(3))//nn(2:3)//'.tbl'
       Elseif (Iset.ge.400 .and. Iset.le.442) Then
          write(nn,'(I3)') Iset
          Tablefile=trim(Flnm(4))//nn(2:3)//'.tbl'
       Elseif (Iset.ge.500 .and. Iset.le.542) Then
          write(nn,'(I3)') Iset
          Tablefile=trim(Flnm(5))//nn(2:3)//'.tbl'
       Else
          Print *, 'Invalid Iset number in SetCJ :', Iset
          Stop
       Endif
       ! print*,'Opening ',Tablefile
       Open(IU, File=prefix // "/" // Tablefile, Status='OLD', Err=100)
21     Call ReadTbl (IU,Iset)
       Close (IU)
       Isetold=Iset
    Endif
    Return

100 Print *, ' Data file ', prefix // "/" // Tablefile, &
    ' cannot be opened in SetCJ!!'
    Stop
    !                             ********************
  End subroutine SetCJ

  Subroutine ReadTbl (Nu,iset)
    Implicit Double Precision (A-H,O-Z)
    Character Line*80
    PARAMETER (MXX = 105, MXQ = 26, MXF = 6)
    PARAMETER (MXPQX=(MXF+3)*MXQ*MXX)
    Common / CJPAR1 / Al, XV(MXX), QV(MXQ), SV(MXQ), UPD(MXPQX)
    Common / CJPAR2 / Nx, NQ, NfMx
    Common / XQrange / Qini, Qmax, Xmin
    Common / QCDtable /  Alam, Nfl, Iorder
    Common / Masstbl / Amass(6)

    Read  (Nu, '(A)') Line
    Read  (Nu, '(A)') Line
    Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
    Iorder = Nint(Dr)
    Nfl = Nint(Fl)
    Alam = Al

    Read  (Nu, '(A)') Line
    Read  (Nu, *) NX,  NQ, NfMx
    Read  (Nu, '(A)') Line
    Read  (Nu, *) QINI, QMAX, (QV(I), I =1, NQ)
    Read  (Nu, '(A)') Line
    Read  (Nu, *) XMIN, (XV(I), I =1, NX)
    do Iq = 1, NQ
       SV(Iq) = Log(Log (QV(Iq) /Alam)/Log(Qini/Alam))
    end do

    Nblk=Nx*NQ
    Npts=Nblk*(NfMx+3)
    Read(Nu,'(A)') Line
    if (iset.le.399) then
       Read(Nu,25,err=999) (UPD(I), I=1,Npts) ! CJ12 format
    else
       Read(Nu,26,err=999) (UPD(I), I=1,Npts) ! CJ15 and higher format
    end if
25  format(5E14.6)
26  format(5E12.5)
    Return
999 print*,i,UPD(i-1),UPD(i)
    stop
  End subroutine ReadTbl

  Function NextUn()
    ! Returns an unallocated FORTRAN i/o unit.
    Logical :: EX

    do N = 10, 300
       INQUIRE (UNIT=N, OPENED=EX)
       If (.NOT. EX) then
          NextUn = N
          Return
       Endif
    end do
    Stop ' There is no available I/O unit. '
    ! *************************
  End function NextUn

  Function GetPDF(Iptn,x,Q)
    IMPLICIT REAL*8 (A-H,O-Z)
    PARAMETER (MXX = 105, MXQ = 26, MXF = 6)
    PARAMETER (MXPQX=(MXF+3)*MXQ*MXX)
    Common / CJPAR1 / Al, XV(MXX), QV(MXQ), SV(MXQ), UPD(MXPQX)
    Common / CJPAR2 / Nx, NQ, NfMx
    Common / XQrange / Qini, Qmax, Xmin
    Common / QCDtable /  Alam, Nfl, Iorder
    Common / Masstbl / Amass(6)

!
!  Linear interpolation in s
!

    s=log(log(Q/Alam)/log(Qini/Alam))

    do j=1,NQ
       if(s.lt.SV(j))then
          J2=J
          if(J2.eq.1)J2=2
          J1=J2-1
          S2=SV(j2)
          S1=SV(j2-1)
          goto 1
       endif
    end do
1   continue
!
!  In CJ12 s=sbar, c=cbar, and b=bbar
!
    if(Iptn.gt.2)then
       Iptn_temp=-Iptn
    else
       Iptn_temp=Iptn
    endif
    A1=GetFV(Iptn_temp,x,J1)
    A2=GetFV(Iptn_temp,x,J2)
    ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)
    if(ans.lt.0.d0)ans=0.d0
    GetPDF=ANS
    RETURN
  END function GetPDF

  FUNCTION GetFV(Iptn,x,J)
    IMPLICIT REAL*8 (A-H,O-Z)
    PARAMETER (MXX = 105, MXQ = 26, MXF = 6)
    PARAMETER (MXPQX=(MXF+3)*MXQ*MXX)
    Dimension xx(4),fx(4)
    Common / CJPAR1 / Al, XV(MXX), QV(MXQ), SV(MXQ), UPD(MXPQX)
    Common / CJPAR2 / Nx, NQ, NfMx
    Common / XQrange / Qini, Qmax, Xmin
    Common / QCDtable /  Alam, Nfl, Iorder
    Common / Masstbl / Amass(6)
!
!  Interpolate in x using the form A*x**alpha*(1-x)**beta
!  which is locally valid for each PDF
!
    DO 1 I=1,Nx
       IF(X.LT.XV(I)) GO TO 2
1      CONTINUE
       !    2 I=I-1
2      I=I-2
    If(I.le.0.d0) I=2
    if(i.gt.(nx-3))i=nx-3
    In1=I+Nx*(J-1)+Nx*NQ*(Iptn+5)
    xx(1)=xv(I)
    xx(2)=xv(i+1)
    xx(3)=xv(I+2)
    xx(4)=xv(I+3)
    fx(1)=UPD(In1)
    fx(2)=UPD(In1+1)
    fx(3)=UPD(In1+2)
    fx(4)=UPD(In1+3)
    !      call polint(xx,fx,4,x,ans,dy)
    call polint4(xx,fx,x,ans)
    getfv=ans
  END function GetFV

  SUBROUTINE POLINT4 (XA,YA,X,Y)

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
    Elseif((H2+H3).lt.0D0) Then
       Y=YA(3)+D2+CD1+DC1
    Elseif((H1+H2).lt.0D0) Then
       Y=YA(2)+C2+CD1+DC1
    ELSE
       Y=YA(1)+C1+CC1+DC1
    ENDIF

    RETURN
  END subroutine polint4

  SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
    IMPLICIT REAL*8 (A-H,O-Z)
    PARAMETER (NMAX=25)
    DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
    NS=1
    DIF=ABS(X-XA(1))
    DO I=1,N
       DIFT=ABS(X-XA(I))
       IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
       ENDIF
       C(I)=YA(I)
       D(I)=YA(I)
    end DO
    Y=YA(NS)
    NS=NS-1
    DO M=1,N-1
       DO I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)then
             print*,xa(i),xa(i+m),x
             print*, 'Enter <CR> to continue'
             read(*,*)
          endif
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
       end DO
       IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
       ELSE
          DY=D(NS)
          NS=NS-1
       ENDIF
       Y=Y+DY
    end DO
    RETURN
  END subroutine polint

end module CJ_pdf
