
! Copyright (C) 2002-2011 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wavefcr(tsh,lrstp,is,ia,ist,m,ld,wfcr)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: lrstp,is,ia,ist
! pass in m-1/2
integer, intent(in) :: m
integer, intent(in) :: ld
complex(4), intent(out) :: wfcr(ld,2)
! local variables
integer ias,nr,nri,ir
integer k,l,lm,lm1,lm2
integer i,i1,i2
real(8) c1,c2,t0,t1,t2
l=lsp(ist,is)
k=ksp(ist,is)
if (((k.ne.l+1).and.(k.ne.l)).or.(m.lt.-k).or.(m.gt.k-1)) then
  write(*,*)
  write(*,'("Error(wavefcr): mismatched l, k or m : ",3I4)') l,k,m
  write(*,'(" for species ",I4)') is
  write(*,'(" atom ",I4)') ia
  write(*,'(" and state ",I6)') ist
  write(*,*)
  stop
end if
if (l.gt.lmaxo) then
  wfcr(:,:)=0.e0
  return
end if
ias=idxas(ia,is)
! calculate the Clebsch-Gordon coefficients
t1=sqrt(dble(l+m+1)/dble(2*l+1))
t2=sqrt(dble(l-m)/dble(2*l+1))
if (k.eq.l+1) then
  c1=t1
  c2=t2
else
  c1=t2
  c2=-t1
end if
if (abs(m).le.l) then
  lm1=l*(l+1)+m+1
else
  lm1=0
end if
if (abs(m+1).le.l) then
  lm2=l*(l+1)+m+2
else
  lm2=0
end if
nr=nrmt(is)
nri=nrmti(is)
! zero the wavefunction
if (lrstp.eq.1) then
  wfcr(1:npmt(is),:)=0.e0
else if (lrstp.eq.lradstp) then
  wfcr(1:npcmt(is),:)=0.e0
else
  write(*,*)
  write(*,'("Error(wavefcr): invalid lrstp : ",I8)') lrstp
  write(*,*)
  stop
end if
!----------------------------------!
!     inner part of muffin-tin     !
!----------------------------------!
if (l.gt.lmaxi) goto 10
if (tsh) then
  i1=lm1
  i2=lm2
else
  i=0
end if
do ir=1,nri,lrstp
! major component of radial wavefunction
  t0=rwfcr(ir,1,ist,ias)*rlmt(ir,-1,is)
  if (tsh) then
    if (lm1.gt.0) wfcr(i1,1)=t0*c1
    if (lm2.gt.0) wfcr(i2,2)=t0*c2
    i1=i1+lmmaxi
    i2=i2+lmmaxi
  else
    t1=t0*c1
    t2=t0*c2
    if (lm1.gt.0) then
      do lm=1,lmmaxi
        wfcr(i+lm,1)=t1*zbshti(lm,lm1)
      end do
    end if
    if (lm2.gt.0) then
      do lm=1,lmmaxi
        wfcr(i+lm,2)=t2*zbshti(lm,lm2)
      end do
    end if
    i=i+lmmaxi
  end if
end do
!----------------------------------!
!     outer part of muffin-tin     !
!----------------------------------!
10 continue
if (lrstp.eq.1) then
  i=lmmaxi*nrmti(is)
else
  i=lmmaxi*nrcmti(is)
end if
if (tsh) then
  i1=i+lm1
  i2=i+lm2
end if
do ir=nri+lrstp,nr,lrstp
  t0=rwfcr(ir,1,ist,ias)*rlmt(ir,-1,is)
  if (tsh) then
    if (lm1.gt.0) wfcr(i1,1)=t0*c1
    if (lm2.gt.0) wfcr(i2,2)=t0*c2
    i1=i1+lmmaxo
    i2=i2+lmmaxo
  else
    t1=t0*c1
    t2=t0*c2
    if (lm1.gt.0) then
      do lm=1,lmmaxo
        wfcr(i+lm,1)=t1*zbshto(lm,lm1)
      end do
    end if
    if (lm2.gt.0) then
      do lm=1,lmmaxo
        wfcr(i+lm,2)=t2*zbshto(lm,lm2)
      end do
    end if
    i=i+lmmaxo
  end if
end do
end subroutine

