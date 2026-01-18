
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfplot(np,vpl,zfmt,zfir,fp)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: np
real(8), intent(in) :: vpl(3,np)
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
complex(8), intent(out) :: fp(np)
! local variables
integer ias,is,ip,nthd
! allocatable arrays
complex(8), allocatable :: zfmt1(:,:,:),zfft(:)
! unpack the muffin-tin function
allocate(zfmt1(lmmaxo,nrcmtmax,natmtot))
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call zfmtpack(.false.,nrcmt(is),nrcmti(is),zfmt(:,ias),zfmt1(:,:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! Fourier transform rfir to G-space
allocate(zfft(ngtot))
zfft(:)=zfir(:)
call zfftifc(3,ngridg,-1,zfft)
! begin loop over all points
call holdthd(np,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ip=1,np
  call zfip(ip)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(zfmt1,zfft)
return

contains

subroutine zfip(ip)
implicit none
! arguments
integer, intent(in) :: ip
! local variables
integer is,ia,ias,nrc,nrci
integer irc0,irc,lmax,l,m,lm
integer ig,ifg,i1,i2,i3,i,j
real(8) rmt2,r,ya1(4),ya2(4),t1,t2
real(8) v1(3),v2(3),v3(3),v4(3),v5(3)
complex(8) z1
! automatic arrays
complex(8) ylm(lmmaxo)
v2(:)=vpl(:,ip)
call r3frac(epslat,v2)
! convert point to Cartesian coordinates
call r3mv(avec,v2,v1)
! check if point is in a muffin-tin
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    v2(:)=v1(:)-atposc(:,ia,is)
    do i1=-1,1
      v3(:)=v2(:)+dble(i1)*avec(:,1)
      do i2=-1,1
        v4(:)=v3(:)+dble(i2)*avec(:,2)
        do i3=-1,1
          v5(:)=v4(:)+dble(i3)*avec(:,3)
          t1=v5(1)**2+v5(2)**2+v5(3)**2
          if (t1.lt.rmt2) then
            r=sqrt(t1)
            call genylmv(lmaxo,v5,ylm)
            do irc=1,nrc
              if (rcmt(irc,is).ge.r) then
                if (irc.le.3) then
                  irc0=1
                else if (irc.gt.nrc-2) then
                  irc0=nrc-3
                else
                  irc0=irc-2
                end if
                r=max(r,rcmt(1,is))
                if (irc0.le.nrci) then
                  lmax=lmaxi
                else
                  lmax=lmaxo
                end if
                z1=0.d0
                lm=0
                do l=0,lmax
                  do m=-l,l
                    lm=lm+1
                    do j=1,4
                      i=irc0+j-1
                      ya1(j)=dble(zfmt1(lm,i,ias))
                      ya2(j)=aimag(zfmt1(lm,i,ias))
                    end do
                    t1=poly4(rcmt(irc0,is),ya1,r)
                    t2=poly4(rcmt(irc0,is),ya2,r)
                    z1=z1+cmplx(t1,t2,8)*ylm(lm)
                  end do
                end do
                goto 10
              end if
            end do
          end if
        end do
      end do
    end do
  end do
end do
! otherwise use direct Fourier transform of interstitial function
z1=0.d0
do ig=1,ngvec
  ifg=igfft(ig)
  t1=vgc(1,ig)*v1(1)+vgc(2,ig)*v1(2)+vgc(3,ig)*v1(3)
  z1=z1+zfft(ifg)*cmplx(cos(t1),sin(t1),8)
end do
10 continue
fp(ip)=z1
end subroutine

pure real(8) function poly4(xa,ya,x)
implicit none
! arguments
real(8), intent(in) :: xa(4),ya(4),x
! local variables
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
! evaluate the polynomial coefficients
x0=xa(1)
x1=xa(2)-x0; x2=xa(3)-x0; x3=xa(4)-x0
t4=x1-x2; t5=x1-x3; t6=x2-x3
y0=ya(1)
y1=ya(2)-y0; y2=ya(3)-y0; y3=ya(4)-y0
t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
t0=1.d0/(x2*t3*t4*t5*t6)
t3=t3*y2
c3=t1*t4+t2*t6-t3*t5
t4=x1**2; t5=x2**2; t6=x3**2
c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
t1=x-x0
! evaluate the polynomial
poly4=y0+t0*t1*(c1+t1*(c2+c3*t1))
end function

end subroutine

