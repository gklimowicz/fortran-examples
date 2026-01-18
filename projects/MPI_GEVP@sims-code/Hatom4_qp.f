      program Hatom4_qp
      implicit real(16)(a-h,o-z), integer(4) (i-n)
!
! This quadruple precision program generates H and S matrices for use 
! as a test of the eigenvalue routine which solves the generalized 
! symmetric indefinite eigenvalue problem H x = L S x
!
! The GEVP (generalized eigenvalue program) computes the S-state 
! energies for the four non-interacting H-atom system using the 
! non-orthogonal basis
!         r**(i-1) exp(-alpha*r), i=1,2,3,...n.
! This program calculates the S- and H-matrices for a single H-atom and
! then uses these matrices to build up H and S for the 4 atom case.
! H and S are stored in canonical upper triangle by column format
! and written to unit 7 (H) and unit 8 (S) in that format.

      INTERFACE
         SUBROUTINE Hatom4_MEBuild_qp(m, alpha, List, HH, SS)
         implicit real(16)(a-h,o-z), integer(4) (i-n)
         integer(4), intent (IN) :: m                !maximum s orbital q.n.
         real(16), intent (IN)   :: alpha            !STO exponent
         logical,  intent (IN)   :: List             !T = print one electron integrals 
         real(16), intent (OUT)  :: HH(*)            !hamiltonian matrix
         real(16), intent (OUT)  :: SS(*)            !overlap matrix
         END
      END INTERFACE

      parameter (MXN = 20)
      real(16), allocatable :: H(:),S(:)
      character(60) :: Printstr
      real(16) :: alpha 
      logical :: List, Already_open
      real(4) :: tin, tout

      print*,' Hatom4_qp computes H and S matrices to be used'
      print*,' in calculating S-state energies for the four'
      print*,' non-interacting H-atom system using the'
      print*,' non-orthogonal basis'
      print*,'       r**(i-1) exp(-alpha*r), i=1,2,3,...n.'

   10 print*,' Enter n, alpha:'
      read*,m, alpha
      if(m.le.0) go to 10
      if(m.gt.MXN) then
         print*,'maximum value of n (=20) exceeded...try again'
         go to 10
      endif
      Nterms = m**4

      print*,' Nterms = ',Nterms
      print*,' Nterms = n**4 term H4 expansion. '
      print*
      print*,' The matrices will be written to unit 7 (H)'
      print*,' and unit 8 (S) in canonical upper triangle'
      print*,' by column format.'

      List = .false.
      kH = 7
      kS = 8
      open(unit=kH,form='unformatted',status='REPLACE')
      open(unit=kS,form='unformatted',status='REPLACE')
      rewind kH; rewind kS
      MESize = Nterms*(Nterms+1)/2
      allocate(H(MESize),S(MESize))
      print*,' Calling Hatom4_MEBuild_qp ...'
      call CPU_time(tin)
      call Hatom4_MEBuild_qp(m, alpha, List, H, S)
      call CPU_time(tout)
      print '(f8.3,a)', tout-tin,
     &      ' seconds required to build H and S partitions.'

! Now write out MEs in canonical upper triangle format
!
      kkk = 0
      do k = 1, Nterms
         write(kH) H(kkk+1:kkk+k)
         write(kS) S(kkk+1:kkk+k)        
         kkk = kkk + k
      enddo   !on k


      inquire(kH, OPENED = Already_open)
      if(Already_open) close(unit=kH)
      inquire(kS, OPENED = Already_open)
      if(Already_open) close(unit=kS)
      if(allocated(H)) deallocate(H)
      if(allocated(S)) deallocate(S)

      END PROGRAM Hatom4_qp

      subroutine Hatom4_MEBuild_qp(m, alpha, List, HH, SS)
      implicit real(16)(a-h,o-z), integer(4) (i-n)
      integer(4), intent (IN) :: m                !maximum s orbital q.n.
      real(16), intent (IN)   :: alpha            !STO exponent
      logical, intent (IN)    :: List             !T = print one electron integrals 
      real(16), intent (OUT)  :: HH(*)            !Natom hamiltonian matrix
      real(16), intent (OUT)  :: SS(*)            !Natom overlap matrix
*
* This subroutine builds the S and H matrices for the the 4 non-interacting 
* H-atom system (singlet S states * only). The one-electron basis used is the 
* normalized non-orthogonal basis
*
*       N(i,alpha) * r**(i-1) exp(-alpha*r), i=1,2,3,...n.
*
* where N(i,alpha) is the normalization factor. The program first calculates 
* the S- and H-matrices for a single H-atom and then uses these matrices to 
* build up the S and H matrices for the 4 atom case.
*
      dimension H(m,m), S(m,m)
      dimension fac(0:60),t(0:60)
      real(16) :: Xnormfac(20)
      dimension ijprod(m**4,4)
*
      one = 1
      MM = m**4
      fac(0)=1
      t(0)=1
      do k=1,2*m+4
         fac(k)=fac(k-1)*k
         t(k)=t(k-1)*2*alpha
      enddo
*
      do j=1,m
         do i=1,m
            S(i,j)=fac(i+j)/t(i+j+1)
            A0=fac(i+j-2)/t(i+j-1)
            A1=fac(i+j-1)/t(i+j)
            A2=fac(i+j)/t(i+j+1)
            H(i,j)=(-one/2*( j*(j-1)*A0 - 2*(j*alpha-1)*A1
     &              + alpha**2*A2 ) )
         enddo   !on i
      enddo   !on j
* At this point we have H and S for one H-atom...using these matrix 
* elements we can build H and S for a system of 4 non-interacting 
* H-atoms. Before doing this, we normalize the STOs (scaling H and S)
* to get a better behaved 4-atom system in so far as the Inverse 
*  Iteration step is concerned.
      do i = 1,m
         Xnormfac(i) = sqrt(1/S(i,i))
      enddo
      do j=1,m
         do i=1,m
            S(i,j) = S(i,j)*Xnormfac(i)*Xnormfac(j)
            H(i,j) = H(i,j)*Xnormfac(i)*Xnormfac(j)
         enddo
      enddo
      if(List) then
         print*,'One-electron integrals over normalized STOs'
         print*,'          i         S(i,i)                H(i,i)'
         do i=1,m
            print*,i,S(i,i),H(i,i)
         enddo
      endif
* Normalization of the STOs is done...

! Now build the full sized H and S
! for a 4 H-atom system...
      kk=0
      do i=1,m
         do j=1,m
            do k=1,m
               do l=1,m
                  kk=kk+1
                  if(kk.gt.MM)then
                     print'(a,i6)',' Too many terms...limit is', MM
                  endif
                  ijprod(kk,1) = i
                  ijprod(kk,2) = j
                  ijprod(kk,3) = k
                  ijprod(kk,4) = l
               enddo   !on l
            enddo   !on k
         enddo   !on j
      enddo   !on i
!     allocate(HH(kk*(kk+1)/2))    !,HHH(kk,kk))
!     allocate(SS(kk*(kk+1)/2))    !,SSS(kk,kk))
      kkk = 0
      do k2 = 1, kk
         j1 = ijprod(k2,1)
         j2 = ijprod(k2,2)
         j3 = ijprod(k2,3)
         j4 = ijprod(k2,4)
         do k1=1,k2
            i1 = ijprod(k1,1)
            i2 = ijprod(k1,2)
            i3 = ijprod(k1,3)
            i4 = ijprod(k1,4)
            kkk = kkk+1
            HH(kkk) = H(i1,j1)*S(i2,j2)*S(i3,j3)*S(i4,j4)
     *              + S(i1,j1)*H(i2,j2)*S(i3,j3)*S(i4,j4)
     *              + S(i1,j1)*s(i2,j2)*H(i3,j3)*S(i4,j4)
     *              + S(i1,j1)*S(i2,j2)*S(i3,j3)*H(i4,j4)
            SS(kkk) = S(i1,j1)*S(i2,j2)*S(i3,j3)*S(i4,j4)
         enddo   !on k1
      enddo   !on k2

      RETURN
      END
