
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dielectric_tdrt
use modmain
use modtddft
implicit none
! local variables
integer its,iw,i,j
real(8) w1,w2,a0,t0,t1,t2
complex(8) eta,z1,z2
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:),wt(:),jt(:,:)
real(8), allocatable :: f1(:),f2(:),g(:)
complex(8), allocatable :: ew(:,:),jw(:,:),eps(:)
! initialise global variables
call init0
call init1
! generate energy grid (always non-negative)
allocate(w(nwplot))
w1=max(wplot(1),0.d0)
w2=max(wplot(2),w1)
t1=(w2-w1)/dble(nwplot)
do iw=1,nwplot
  w(iw)=w1+t1*dble(iw-1)
end do
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
! determine the weights for the spline integration
allocate(wt(ntimes))
call wsplint(ntimes,times,wt)
! compute the electric field from E = -1/c dA/dt and Fourier transform
allocate(f1(ntimes),f2(ntimes),g(ntimes))
allocate(ew(nwplot,3))
t0=-1.d0/solsc
do i=1,3
  f1(:)=afieldt(i,:)
  call fderiv(1,ntimes,times,f1,g)
! constant term corresponding to instantaneous A-field at t=0
  a0=f1(1)
  if (task.eq.480) then
! Fourier transform E(t) numerically
    do iw=1,nwplot
      do its=1,ntimes
        t1=g(its)
        t2=w(iw)*times(its)
        f1(its)=t1*cos(t2)
        f2(its)=t1*sin(t2)
      end do
      t1=dot_product(wt(:),f1(:))
      t2=dot_product(wt(:),f2(:))
      ew(iw,i)=t0*cmplx(t1+a0,t2,8)
    end do
  else
! analytic Fourier transform of E(t) assumed to be a delta function at t=0
    t1=t0*a0
    ew(1:nwplot,i)=t1
  end if
end do
! read in the total current from file
allocate(jt(3,ntimes))
call readjtot(jt)
! divide by the unit cell volume
jt(:,:)=jt(:,:)/omega
! filter the high-frequency components from the current
do its=1,ntimes
  t1=exp(-swidth*times(its))
  jt(:,its)=t1*jt(:,its)
end do
! Fourier transform the current
allocate(jw(nwplot,3))
do i=1,3
  do iw=1,nwplot
    do its=1,ntimes
      t1=jt(i,its)
      t2=w(iw)*times(its)
      f1(its)=t1*cos(t2)
      f2(its)=t1*sin(t2)
    end do
    t1=dot_product(wt(:),f1(:))
    t2=dot_product(wt(:),f2(:))
    jw(iw,i)=cmplx(t1,t2,8)
  end do
end do
deallocate(wt,f1,f2,g,jt)
! compute the dielectric function and write to file
allocate(eps(nwplot))
do i=1,3
  do j=1,3
    do iw=1,nwplot
      z1=jw(iw,i)
      z2=ew(iw,j)
      t1=abs(dble(z2))+abs(aimag(z2))
      if (t1.gt.1.d-8) then
        z1=z1/z2
      else
        z1=0.d0
      end if
      z1=fourpi*cmplx(-aimag(z1),dble(z1),8)
      z1=z1/(w(iw)+eta)
      if (i.eq.j) z1=z1+1.d0
      eps(iw)=z1
    end do
    write(fname,'("EPSILON_TDRT_",2I1,".OUT")') i,j
    open(50,file=trim(fname),form='FORMATTED')
    do iw=1,nwplot
      write(50,'(2G18.10)') w(iw),dble(eps(iw))
    end do
    write(50,*)
    do iw=1,nwplot
      write(50,'(2G18.10)') w(iw),aimag(eps(iw))
    end do
    close(50)
  end do
end do
write(*,*)
write(*,'("Info(dielectric_tdrt):")')
write(*,'(" dielectric tensor determined from real-time evolution")')
write(*,'(" written to EPSILON_TDRT_ij.OUT for components i,j = 1,2,3")')
write(*,*)
write(*,'("(Note that only those components which are not orthogonal to the")')
write(*,'(" applied A-field will be calculated correctly)")')
deallocate(w,ew,jw,eps)
end subroutine

