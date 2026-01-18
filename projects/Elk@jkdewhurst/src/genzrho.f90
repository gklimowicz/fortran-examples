
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genzrho
! !INTERFACE:
subroutine genzrho(tsh,tspc,ngt,wfmt1,wfir1,wfmt2,wfir2,zrhomt,zrhoir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh    : .true. if the muffin-tin density is to be in spherical harmonics
!            (in,logical)
!   tspc   : .true. if the density should be contracted over spin (in,logical)
!   ngt    : total number of interstitial grid points (in,integer)
!   wfmt1  : muffin-tin part of wavefunction 1 in spherical coordinates
!            (in,complex(npcmtmax,natmtot,*))
!   wfir1  : interstitial wavefunction 1 (in,complex(ngt,*))
!   wfmt2  : muffin-tin part of wavefunction 2 in spherical coordinates
!            (in,complex(npcmtmax,natmtot,*))
!   wfir2  : interstitial wavefunction 2 (in,complex(ngt,*))
!   zrhomt : muffin-tin charge density in spherical harmonics/coordinates
!            (out,complex(npcmtmax,natmtot))
!   zrhoir : interstitial charge density (out,complex(ngt))
! !DESCRIPTION:
!   Calculates the complex overlap charge density from two input wavefunctions:
!   $$ \rho({\bf r})\equiv\Psi_1^{\dag}({\bf r})\cdot\Psi_2({\bf r}). $$
!   Note that the muffin-tin wavefunctions are provided in spherical coordinates
!   and the returned density is either in terms of spherical harmonic
!   coefficients or spherical coordinates when {\tt tsh} is {\tt .true.} or
!   {\tt .false.}, respectively.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!   Changed to single-precision wavefunctions, November 2022 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh,tspc
integer, intent(in) :: ngt
complex(4), intent(in) :: wfmt1(npcmtmax,natmtot,*),wfir1(ngt,*)
complex(4), intent(in) :: wfmt2(npcmtmax,natmtot,*),wfir2(ngt,*)
complex(8), intent(out) :: zrhomt(npcmtmax,natmtot),zrhoir(ngt)
! local variables
integer is,ias
! automatic arrays
complex(8) zfmt(npcmtmax)
! muffin-tin part
do ias=1,natmtot
  is=idxis(ias)
  if (tsh) then
    if (tspc.and.spinpol) then
! contract over spin
      call zrho2(npcmt(is),wfmt1(:,ias,1),wfmt1(:,ias,2),wfmt2(:,ias,1), &
       wfmt2(:,ias,2),zfmt)
    else
! no spin contraction
      call zrho1(npcmt(is),wfmt1(:,ias,1),wfmt2(:,ias,1),zfmt)
    end if
! convert to spherical harmonics
    call zfsht(nrcmt(is),nrcmti(is),zfmt,zrhomt(:,ias))
  else
    if (tspc.and.spinpol) then
      call zrho2(npcmt(is),wfmt1(:,ias,1),wfmt1(:,ias,2),wfmt2(:,ias,1), &
       wfmt2(:,ias,2),zrhomt(:,ias))
    else
      call zrho1(npcmt(is),wfmt1(:,ias,1),wfmt2(:,ias,1),zrhomt(:,ias))
    end if
  end if
end do
! interstitial part
if (tspc.and.spinpol) then
  call zrho2(ngt,wfir1,wfir1(:,2),wfir2,wfir2(:,2),zrhoir)
else
  call zrho1(ngt,wfir1,wfir2,zrhoir)
end if
return

contains

pure subroutine zrho1(n,wf1,wf2,zrho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf1(n),wf2(n)
complex(8), intent(out) :: zrho(n)
zrho(:)=conjg(wf1(:))*wf2(:)
end subroutine

pure subroutine zrho2(n,wf11,wf12,wf21,wf22,zrho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(8), intent(out) :: zrho(n)
zrho(:)=conjg(wf11(:))*wf21(:)+conjg(wf12(:))*wf22(:)
end subroutine

end subroutine
!EOC

