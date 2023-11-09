
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phlwidth
use modmain
use modphonon
implicit none
! local variables
integer i,j,iq,iv
real(8) gmin,gmax,t1
! allocatable arrays
real(8), allocatable :: wq(:),gq(:,:),gp(:,:)
complex(8), allocatable :: dynq(:,:,:),dynr(:,:,:),dynp(:,:)
complex(8), allocatable :: gmq(:,:,:),gmr(:,:,:),gmp(:,:)
complex(8), allocatable :: ev(:,:),b(:,:)
! initialise universal variables
call init0
call init2
allocate(wq(nbph),gq(nbph,nqpt),gp(nbph,npp1d))
allocate(dynq(nbph,nbph,nqpt))
allocate(dynr(nbph,nbph,nqptnr))
allocate(dynp(nbph,nbph))
allocate(gmq(nbph,nbph,nqpt))
allocate(gmr(nbph,nbph,nqptnr))
allocate(gmp(nbph,nbph))
allocate(ev(nbph,nbph),b(nbph,nbph))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! read in the phonon linewidths for each q-point
call readgamma(gq)
! loop over phonon q-points
do iq=1,nqpt
! find the eigenvalues and eigenvectors of the dynamical matrix
  call dynev(dynq(:,:,iq),wq,ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues squared are the phonon linewidths
  do i=1,nbph
    t1=sqrt(abs(gq(i,iq)))
    do j=1,nbph
      b(i,j)=t1*conjg(ev(j,i))
    end do
  end do
  call zgemm('N','N',nbph,nbph,nbph,zone,ev,nbph,b,nbph,zzero,gmq(:,:,iq),nbph)
end do
! Fourier transform the gamma matrices to real-space
call dynqtor(gmq,gmr)
! generate a set of q-point vectors along a path in the Brillouin zone
call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
gmin=1.d8
gmax=0.d0
! compute the linewidths along the path
do iq=1,npp1d
! compute the dynamical matrix at this particular q-point
  call dynrtoq(vplp1d(:,iq),dynr,dynp)
! find the phonon eigenvalues and eigenvectors
  call dynev(dynp,wq,ev)
! compute the gamma matrix at this particular q-point
  call dynrtoq(vplp1d(:,iq),gmr,gmp)
! diagonalise the gamma matrix simultaneously with the dynamical matrix
  call dynevs(ev,gmp,gp(:,iq))
! square the eigenvalues to recover the linewidths
  gp(:,iq)=gp(:,iq)**2
  do i=1,nbph
    t1=gp(i,iq)
    if (t1.lt.gmin) gmin=t1
    if (t1.gt.gmax) gmax=t1
  end do
end do
t1=(gmax-gmin)*0.5d0
gmin=gmin-t1
gmax=gmax+t1
! output the vertex location lines
open(50,file='PHLWLINES.OUT',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),gmin
  write(50,'(2G18.10)') dvp1d(iv),gmax
  write(50,*)
end do
close(50)
! output the phonon linewidth dispersion
open(50,file='PHLWIDTH.OUT',form='FORMATTED')
do i=1,nbph
  do iq=1,npp1d
    write(50,'(2G18.10)') dpp1d(iq),gp(i,iq)
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(phlwidth):")')
write(*,'(" phonon linewidth dispersion written to PHLWIDTH.OUT")')
write(*,'(" vertex location lines written to PHLWLINES.OUT")')
deallocate(wq,gq,gp)
deallocate(dynq,dynr,dynp)
deallocate(gmq,gmr,gmp)
deallocate(ev,b)
end subroutine

