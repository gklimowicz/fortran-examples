
! Copyright (C) 2002-2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: nonlinopt
! !INTERFACE:
subroutine nonlinopt
! !USES:
use modmain
use modmpi
use modomp
use modtest
! !DESCRIPTION:
!   Calculates the second-order response tensor
!   $\chi^{abc}(-2\omega;\omega,\omega)$, where $a$, $b$ and $c$ label Cartesian
!   directions. This tensor is used for determining the optical second-harmonic
!   generation of materials. We follow the convention of Sipe and Ghahramani in
!   {\it Phys. Rev. B} {\bf 48}, 11705 (1993); and Hughes and Sipe in
!   {\it Phys. Rev. B} {\bf 53}, 10751 (1996). The individual contributions
!   $\chi_{II}^{abc}(-2\omega;\omega,\omega)$,
!   $\eta_{II}^{abc}(-2\omega;\omega,\omega)$ and
!   $\frac{i}{2\omega}\sigma_{II}^{abc}(-2\omega;\omega,\omega)$ are also
!   written separately to file.
!
! !REVISION HISTORY:
!   Rewrote earlier version, June 2010 (Sharma)
!   Improved parallelism, January 2020 (R. Cohen)
!   Rewrote, thanks to corrections from X. Gonze, March 2022 (JKD)
!EOP
!BOC
implicit none
! local variables
logical tssr
integer ik,jk,l,m,n,i
integer iw,ioc,a,b,c
integer nthd
real(8) t0,t1
complex(8) eta,z1,z2
character(64) fname
! allocatable arrays
real(8), allocatable :: w(:),e(:,:),d(:,:,:),f(:,:)
complex(8), allocatable :: r(:,:,:),zv(:)
complex(8), allocatable :: cc1(:,:),cc2(:,:)
complex(8), allocatable :: ce1(:,:),ce2(:,:),cs1(:,:)
complex(8), allocatable :: chi2w(:),eta2w(:),sigma2w(:)
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
! read the eigenvalues and occupancies from file
call readevalsv
call readoccsv
! set flag for scissor operator
if (abs(scissor).gt.1.d-8) then
  tssr=.true.
else
  tssr=.false.
end if
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
! generate energy grid (starting from zero)
allocate(w(nwplot))
t1=wplot(2)/dble(nwplot)
do iw=1,nwplot
  w(iw)=t1*dble(iw-1)
end do
allocate(chi2w(nwplot),eta2w(nwplot),sigma2w(nwplot))
t0=wkptnr/omega
! begin loop over components
do ioc=1,noptcomp
  a=optcomp(1,ioc)
  b=optcomp(2,ioc)
  c=optcomp(3,ioc)
  chi2w(:)=0.d0
  eta2w(:)=0.d0
  sigma2w(:)=0.d0
! parallel loop over non-reduced k-points
  call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(e,f,d,r,zv) &
!$OMP PRIVATE(cc1,cc2,ce1,ce2,cs1) &
!$OMP PRIVATE(jk,n,m,i,t1,z1,z2,l) &
!$OMP REDUCTION(+:chi2w,eta2w,sigma2w) &
!$OMP NUM_THREADS(nthd)
  allocate(e(nstsv,nstsv),f(nstsv,nstsv),d(nstsv,nstsv,3))
  allocate(r(nstsv,nstsv,3),zv(nwplot))
  allocate(cc1(nstsv,nstsv),cc2(nstsv,nstsv))
  allocate(ce1(nstsv,nstsv),ce2(nstsv,nstsv))
  allocate(cs1(nstsv,nstsv))
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(nonlinopt_)
    write(*,'("Info(nonlinopt): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(nonlinopt_)
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! calculate differences in eigenvalues and occupation numbers
    do n=1,nstsv
      do m=1,nstsv
        e(m,n)=evalsv(m,jk)-evalsv(n,jk)
        f(m,n)=occsv(m,jk)-occsv(n,jk)
      end do
    end do
! read momentum matrix elements from file
    call getpmat(vkl(:,ik),r)
! compute the Delta matrix elements
    do i=1,3
      do n=1,nstsv
        do m=1,nstsv
          d(m,n,i)=dble(r(m,m,i))-dble(r(n,n,i))
        end do
      end do
    end do
! compute the matrix elements of the position operator
    do i=1,3
      do n=1,nstsv
        do m=1,nstsv
          t1=e(m,n)
          if (abs(t1).gt.swidth) then
            z1=r(m,n,i)/t1
            r(m,n,i)=cmplx(aimag(z1),-dble(z1),8)
          else
            r(m,n,i)=0.d0
          end if
        end do
      end do
    end do
! apply the scissor operator to the eigenvalue differences
    if (tssr) then
      do m=1,nstsv
        if (evalsv(m,jk).gt.efermi) then
          e(m,:)=e(m,:)+scissor
          e(:,m)=e(:,m)-scissor
        end if
      end do
    end if
! zero the coefficients for χ_II, η_II and i/2ω σ_II
    cc1(:,:)=0.d0; cc2(:,:)=0.d0
    ce1(:,:)=0.d0; ce2(:,:)=0.d0
    cs1(:,:)=0.d0
! sum over states
    do n=1,nstsv
      do m=1,nstsv
        do l=1,nstsv
! terms involving a triple summation
          z1=0.5d0*r(n,m,a)*(r(m,l,b)*r(l,n,c)+r(m,l,c)*r(l,n,b))
! χ_II(-2ω;ω,ω) terms
          t1=e(l,n)-e(m,l)
          if (abs(t1).gt.swidth) then
! Eq. (B4)
            z2=z1/t1
            if (abs(f(n,m)).gt.epsocc) then
              cc2(m,n)=cc2(m,n)+2.d0*f(n,m)*z2
            end if
            if (abs(f(m,l)).gt.epsocc) then
              cc1(m,l)=cc1(m,l)+f(m,l)*z2
            end if
            if (abs(f(l,n)).gt.epsocc) then
              cc1(l,n)=cc1(l,n)+f(l,n)*z2
            end if
          end if
! η_II(-2ω;ω,ω) terms
          z2=z1*e(m,n)
! Eq. (B13b)
          if (abs(f(n,l)).gt.epsocc) then
            t1=e(l,n)
            if (abs(t1).gt.swidth) then
              ce1(l,n)=ce1(l,n)+f(n,l)*z2/t1**2
            end if
          end if
          if (abs(f(l,m)).gt.epsocc) then
            t1=e(m,l)
            if (abs(t1).gt.swidth) then
              ce1(m,l)=ce1(m,l)-f(l,m)*z2/t1**2
            end if
          end if
          if (abs(f(n,m)).gt.epsocc) then
! Eq. (B13a)
            t1=e(m,n)
            if (abs(t1).gt.swidth) then
              t1=1.d0/t1**2
              z1=2.d0*f(n,m)*(e(m,l)-e(l,n))*t1*z1
              ce2(m,n)=ce2(m,n)+z1
! i/2ω σ_II(-2ω;ω,ω) term
! Eq. (B17)
              z1=e(n,l)*r(l,m,a)*(r(m,n,b)*r(n,l,c)+r(m,n,c)*r(n,l,b)) &
                -e(l,m)*r(n,l,a)*(r(l,m,b)*r(m,n,c)+r(l,m,c)*r(m,n,b))
              z1=0.25d0*f(n,m)*t1*z1
              cs1(m,n)=cs1(m,n)+z1
            end if
          end if
        end do
! terms involving a double summation
        if (abs(f(n,m)).gt.epsocc) then
! Eq. (B12a)
          t1=e(m,n)
          if (abs(t1).gt.swidth) then
            t1=1.d0/t1**2
            z1=r(n,m,a)*(d(m,n,b)*r(m,n,c)+d(m,n,c)*r(m,n,b))
            z1=cmplx(aimag(z1),-dble(z1),8)
            z1=4.d0*f(n,m)*t1*z1
            ce2(m,n)=ce2(m,n)+z1
! Eq. (B16b)
            z1=r(n,m,a)*(r(m,n,b)*d(m,n,c)+r(m,n,c)*d(m,n,b))
            z1=cmplx(-aimag(z1),dble(z1),8)
            z1=0.25d0*f(n,m)*t1*z1
            cs1(m,n)=cs1(m,n)+z1
          end if
        end if
      end do
    end do
    do n=1,nstsv
      do m=1,nstsv
        zv(:)=1.d0/(e(m,n)-w(:)+eta)
        chi2w(:)=chi2w(:)+cc1(m,n)*zv(:)
        eta2w(:)=eta2w(:)+ce1(m,n)*zv(:)
        sigma2w(:)=sigma2w(:)+cs1(m,n)*zv(:)
        zv(:)=1.d0/(e(m,n)-2.d0*(w(:)-eta))
        chi2w(:)=chi2w(:)+cc2(m,n)*zv(:)
        eta2w(:)=eta2w(:)+ce2(m,n)*zv(:)
      end do
    end do
  end do
!$OMP END DO
  deallocate(e,f,d,r,zv)
  deallocate(cc1,cc2,ce1,ce2,cs1)
!$OMP END PARALLEL
  call freethd(nthd)
! multiply response functions by prefactor
  chi2w(:)=t0*chi2w(:)
  eta2w(:)=t0*eta2w(:)
  sigma2w(:)=t0*sigma2w(:)
! add response functions from each process and redistribute
  if (np_mpi.gt.1) then
    call mpi_allreduce(mpi_in_place,chi2w,nwplot,mpi_double_complex,mpi_sum, &
     mpicom,ierror)
    call mpi_allreduce(mpi_in_place,eta2w,nwplot,mpi_double_complex,mpi_sum, &
     mpicom,ierror)
    call mpi_allreduce(mpi_in_place,sigma2w,nwplot,mpi_double_complex,mpi_sum, &
     mpicom,ierror)
  end if
! write χ_II(-2ω;ω,ω), η_II(-2ω;ω,ω) and i/2ω σ_II(-2ω;ω,ω) to file
  if (mp_mpi) then
    write(fname,'("CHI_II_2WWW_",3I1,".OUT")') a,b,c
    open(50,file=trim(fname),form='FORMATTED')
    write(fname,'("ETA_II_2WWW_",3I1,".OUT")') a,b,c
    open(51,file=trim(fname),form='FORMATTED')
    write(fname,'("SIGMA_II_2WWW_",3I1,".OUT")') a,b,c
    open(52,file=trim(fname),form='FORMATTED')
    do iw=1,nwplot
      t1=dble(w(iw))
      write(50,'(2G18.10)') t1,dble(chi2w(iw))
      write(51,'(2G18.10)') t1,dble(eta2w(iw))
      write(52,'(2G18.10)') t1,dble(sigma2w(iw))
    end do
    write(50,*)
    write(51,*)
    write(52,*)
    do iw=1,nwplot
      t1=dble(w(iw))
      write(50,'(2G18.10)') t1,aimag(chi2w(iw))
      write(51,'(2G18.10)') t1,aimag(eta2w(iw))
      write(52,'(2G18.10)') t1,aimag(sigma2w(iw))
    end do
    close(50)
    close(51)
    close(52)
! write χ(-2ω;ω,ω) to file
    chi2w(:)=chi2w(:)+eta2w(:)+sigma2w(:)
    write(fname,'("CHI_2WWW_",3I1,".OUT")') a,b,c
    open(50,file=trim(fname),form='FORMATTED')
    do iw=1,nwplot
      t1=dble(w(iw))
      write(50,'(2G18.10)') t1,dble(chi2w(iw))
    end do
    write(50,*)
    do iw=1,nwplot
      t1=dble(w(iw))
      write(50,'(2G18.10)') t1,aimag(chi2w(iw))
    end do
    close(50)
  end if
! end loop over components
end do
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(nonlinopt):")')
  write(*,'(" Following the convention in Phys. Rev. B 48, 11705 (1993) and")')
  write(*,'(" Phys. Rev. B 53, 10751 (1996), the second-order response")')
  write(*,'(" functions χ_II(-2ω;ω,ω), η_II(-2ω;ω,ω) and i/2ω σ_II(-2ω;ω,ω)")')
  write(*,'(" were written to the files CHI_II_2WWW_abc.OUT,")')
  write(*,'(" ETA_II_2WWW_abc.OUT and SIGMA_II_2WWW_abc.OUT, respectively")')
  write(*,*)
  write(*,'(" The total second-order response function χ(-2ω;ω,ω) was")')
  write(*,'(" written to the file CHI_2WWW_abc.OUT")')
  write(*,*)
  write(*,'(" This was done for Cartesian components :")')
  do ioc=1,noptcomp
    write(*,'("  a = ",I1,", b = ",I1,", c = ",I1)') optcomp(1:3,ioc)
  end do
end if
! write chi2w to test file if required
call writetest(125,'non-linear susceptibility',nv=nwplot,tol=1.d-2,zva=chi2w)
deallocate(w,chi2w,eta2w,sigma2w)
end subroutine
!EOC

