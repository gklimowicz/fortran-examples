!------------------------------------------------------------------------------
module GaussianQuadrature_mod
!------------------------------------------------------------------------------
!@sum Module that contains methods to perform integration using Gaussian
!@+ quadrature
!@auth SSSO ASTG
  implicit none
  private
  public OrthoPolyCoeffs
  public GaussianCoeffs

  integer, parameter :: DP = selected_real_kind(14)
  integer, parameter :: NMAX=30

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  SUBROUTINE OrthoPolyCoeffs(N,ANU,A,B)
!@sum Computes coefficients A, B given modified moments ANU
    INTEGER, intent(in) :: N
    REAL(kind=DP), intent(in) :: ANU(2*N)
    REAL(kind=DP), intent(out) :: A(N), B(N)
!
    INTEGER :: k
    REAL(kind=DP) :: sig(2*NMAX+1,2*NMAX+1)

    sig(1,3:2*n) = 0.0
    sig(2,2:2*n+1) = anu(1:2*n)
    a(1) = anu(2)/anu(1)
    b(1) = 0.0d0
    DO k = 3,n+1
      sig(k,k:2*n-k+3) = sig(k-1,k+1:2*n-k+4) - a(k-2)*sig(k-1,k:2*n-k+3) &
        - b(k-2)*sig(k-2,k:2*n-k+3)
      a(k-1) = sig(k,k+1)/sig(k,k) - sig(k-1,k)/sig(k-1,k-1)
      b(k-1) = sig(k,k)/sig(k-1,k-1)
    end do

  END SUBROUTINE OrthoPolyCoeffs


  SUBROUTINE GaussianCoeffs(N,A,B,AMU0,X,W,QLerror)
!@sum Computes abscissas and weights for Gaussian quadrature formula
    INTEGER, intent(in ) :: N
    INTEGER, intent(out) :: QLerror
    REAL(kind=DP), intent(inout) :: A(N),B(N)
    REAL(kind=DP), intent(in ) :: AMU0
    REAL(kind=DP), intent(out) :: W(N),X(N)
!
    REAL(kind=DP) :: Z(NMAX,NMAX)
    INTEGER :: I,J

    QLerror = 0
    b(2:n) = sqrt(b(2:n))
    z = 0.0d0
    do i = 1,n
      z(i,i) = 1.D0
    end do
    CALL TriDiagQL(A,B,Z,QLerror)
    IF(QLerror > 0) RETURN

    ! Ordering of the abscissas is usually not needed.
    !CALL sortEigen(A,Z,N,NMAX)
    x = a
    w = amu0*z(1,1:n)**2

  END SUBROUTINE GaussianCoeffs


  SUBROUTINE TriDiagQL(triMat, subDiag, eigen, ierr)
!@sum QL algorithm to determine the eigenvalues and eigenvectors of a real,
!@+ symmetric, tridiagonal matrix (triMat)
    REAL(kind=DP), DIMENSION(:), INTENT(INOUT) :: triMat,subDiag 
    REAL(kind=DP), DIMENSION(:,:), optional, INTENT(INOUT) :: eigen
    INTEGER, intent(out) :: ierr
    INTEGER :: i, k, l, m, matSize, iter
    REAL(kind=DP) :: s, c, p, b, f
    REAL(kind=DP) :: absM, shift, norm

    matSize = size(triMat)
    subDiag(:) = eoshift(subDiag(:),1)
    subDiag(matSize) = 0.D0

    do L = 1, matSize
      ITER = 0
      iterate: do
        do M = L, matSize-1
          absM = ABS(triMat(M)) + ABS(triMat(M+1))
          IF (ABS(subDiag(M))+absM == absM) exit
        end do
        if (m == l) exit iterate
        if (iter == 300) then
          ierr = 1
          RETURN
        end if
        ITER=ITER+1
        shift = (triMat(L+1)-triMat(L))/(2.D+00*subDiag(L))
        norm = EuclideanNorm(shift, 1.0D0)
        shift = triMat(M)-triMat(L)+subDiag(L)/(shift+SIGN(norm,shift))
        S = 1.D0
        C = 1.D0
        P = 0.D0
        do I = M-1,L,-1
          F = S*subDiag(I)
          B = C*subDiag(I)
          norm = EuclideanNorm(F,shift)
          subDiag(I+1) = norm
          if (norm == 0.D0) then
            triMat(I+1) = triMat(I+1)-P
            subDiag(M) = 0.D0
            cycle iterate
          end if
          S = F/norm
          C = shift/norm
          shift = triMat(I+1)-P
          norm = (triMat(I)-shift)*S+2.D0*C*B
          P = S*norm
          triMat(I+1) = shift+P
          shift = C*norm-B
          if (present(eigen)) then !  form eigenvectors
            do K = 1,matSize
              F = eigen(K,I+1)
              eigen(K,I+1) = S*eigen(K,I)+C*F
              eigen(K,I) = C*eigen(K,I)-S*F
            end do
          end if
        end do ! I
        triMat(L) = triMat(L)-P
        subDiag(L) = shift
        subDiag(M) = 0.D0
      end do iterate
    end do ! L

  END SUBROUTINE TriDiagQL


  REAL(kind=DP) FUNCTION EuclideanNorm(A,B)
!@sum Computes (a**2+b**2)**(1/2)
    REAL(kind=DP) :: A,B
!
    REAL(kind=DP) :: ABSA,ABSB

    ABSA = ABS(A)
    ABSB = ABS(B)
    IF (ABSA > ABSB) THEN
      EuclideanNorm = ABSA*SQRT(1.D+00+(ABSB/ABSA)**2)
    ELSE
      IF (ABSB == 0.D+00) THEN
        EuclideanNorm = 0.D+00
      ELSE
        EuclideanNorm = ABSB*SQRT(1.D+00+(ABSA/ABSB)**2)
      END IF
    ENDIF

  END FUNCTION EuclideanNorm


  SUBROUTINE sortEigen(eVal,eVec,N,NP)
!@sum Given the eigenvalues and eigenvectors, sort them in descending
!@+ order using insertion method
    INTEGER , intent(in) :: N,NP
    REAL(kind=DP), intent(out) :: eVal(NP),eVec(NP,NP)
!
    INTEGER :: K,J,I
    REAL(kind=DP) :: P

    DO I = 1,N-1
      K = I
      P = eVal(I)
      DO J = I+1,N
        IF(eVal(J) >= P)THEN
          K = J
          P = eVal(J)
        ENDIF
      end do ! J
      IF(K /= I)THEN
        eVal(K) = eVal(I)
        eVal(I) = P
        DO J = 1,N
          P = eVec(J,I)
          eVec(J,I) = eVec(J,K)
          eVec(J,K) = P
        end do ! J
      ENDIF
    end do !I

  END SUBROUTINE sortEigen

end module GaussianQuadrature_mod
