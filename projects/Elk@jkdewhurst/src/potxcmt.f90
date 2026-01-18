
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potxcmt(tsh,ias,xctype_,rhomt_,magmt_,taumt_,exmt_,ecmt_,vxcmt_, &
 bxcmt_,wxcmt_)
use modmain
use modxcifc
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: ias,xctype_(3)
real(8), intent(in) :: rhomt_(npmtmax,natmtot),magmt_(npmtmax,natmtot,ndmag)
real(8), intent(in) :: taumt_(npmtmax,natmtot,nspinor)
real(8), intent(out) :: exmt_(npmtmax,natmtot),ecmt_(npmtmax,natmtot)
real(8), intent(out) :: vxcmt_(npmtmax,natmtot),bxcmt_(npmtmax,natmtot,ndmag)
real(8), intent(out) :: wxcmt_(npmtmax,natmtot)
! local variables
integer ispn,idm,is
integer nr,nri,n,i
real(8) t0,t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: rho(:),rhoup(:),rhodn(:)
real(8), allocatable :: gvrho(:,:),gvup(:,:),gvdn(:,:)
real(8), allocatable :: grho(:),gup(:),gdn(:)
real(8), allocatable :: g2rho(:),g2up(:),g2dn(:)
real(8), allocatable :: g3rho(:),g3up(:),g3dn(:)
real(8), allocatable :: grho2(:),gup2(:),gdn2(:),gupdn(:)
real(8), allocatable :: ex(:),ec(:),vxc(:)
real(8), allocatable :: vx(:),vxup(:),vxdn(:)
real(8), allocatable :: vc(:),vcup(:),vcdn(:)
real(8), allocatable :: mag(:,:),bxc(:,:),tau(:,:)
real(8), allocatable :: dxdgr2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
real(8), allocatable :: dcdgr2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
real(8), allocatable :: dxdg2r(:),dxdg2u(:),dxdg2d(:)
real(8), allocatable :: dcdg2r(:),dcdg2u(:),dcdg2d(:)
real(8), allocatable :: dtdr(:),dtdru(:),dtdrd(:)
real(8), allocatable :: dtdgr2(:),dtdgu2(:),dtdgd2(:)
real(8), allocatable :: wx(:),wxup(:),wxdn(:)
real(8), allocatable :: wc(:),wcup(:),wcdn(:)
is=idxis(ias)
n=npmt(is)
! allocate local arrays
allocate(rho(n),ex(n),ec(n),vxc(n))
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) allocate(tau(n,nspinor))
if (spinpol) then
  allocate(mag(n,3),bxc(n,3))
end if
if (spinpol) then
  allocate(rhoup(n),rhodn(n))
  allocate(vxup(n),vxdn(n),vcup(n),vcdn(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),gup(n),gdn(n))
    allocate(g2up(n),g2dn(n))
    allocate(g3rho(n),g3up(n),g3dn(n))
  else if (xcgrad.eq.2) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(n,3),gvdn(n,3))
    allocate(gup2(n),gdn2(n),gupdn(n))
    allocate(dxdgu2(n),dxdgd2(n),dxdgud(n))
    allocate(dcdgu2(n),dcdgd2(n),dcdgud(n))
  else if (xcgrad.eq.3) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(n,3),gvdn(n,3))
    allocate(gup2(n),gdn2(n),gupdn(n))
  else if (xcgrad.eq.4) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(n,3),gvdn(n,3))
    allocate(gup2(n),gdn2(n),gupdn(n))
    allocate(dxdgu2(n),dxdgd2(n),dxdgud(n))
    allocate(dcdgu2(n),dcdgd2(n),dcdgud(n))
    allocate(dxdg2u(n),dxdg2d(n))
    allocate(dcdg2u(n),dcdg2d(n))
    allocate(dtdru(n),dtdrd(n))
    allocate(dtdgu2(n),dtdgd2(n))
    allocate(wxup(n),wxdn(n),wcup(n),wcdn(n))
  end if
else
  allocate(vx(n),vc(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),g2rho(n),g3rho(n))
  else if (xcgrad.eq.2) then
    allocate(g2rho(n),gvrho(n,3),grho2(n))
    allocate(dxdgr2(n),dcdgr2(n))
  else if (xcgrad.eq.3) then
    allocate(g2rho(n),gvrho(n,3),grho2(n))
  else if (xcgrad.eq.4) then
    allocate(g2rho(n),gvrho(n,3),grho2(n))
    allocate(dxdgr2(n),dcdgr2(n))
    allocate(dxdg2r(n),dcdg2r(n))
    allocate(dtdr(n),dtdgr2(n))
    allocate(wx(n),wc(n))
  end if
end if
nr=nrmt(is)
nri=nrmti(is)
if (tsh) then
! convert the density to spherical coordinates
  call rbsht(nr,nri,rhomt_(:,ias),rho)
else
  rho(1:n)=rhomt_(1:n,ias)
end if
! convert tau to spherical coordinates if required
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) then
  do ispn=1,nspinor
    if (tsh) then
      call rbsht(nr,nri,taumt_(:,ias,ispn),tau(:,ispn))
    else
      tau(1:n,ispn)=taumt_(1:n,ias,ispn)
    end if
  end do
end if
if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
! magnetisation in spherical coordinates
  do idm=1,ndmag
    if (tsh) then
      call rbsht(nr,nri,magmt_(:,ias,idm),mag(:,idm))
    else
      mag(1:n,idm)=magmt_(1:n,ias,idm)
    end if
  end do
! use scaled spin exchange-correlation if required
  if (tssxc) mag(:,1:ndmag)=mag(:,1:ndmag)*sxcscf
  if (ncmag) then
! non-collinear (use Kubler's trick)
    if (xcgrad.eq.0) then
! LSDA
      do i=1,n
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
        t0=rho(i)
        t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
        rhoup(i)=0.5d0*(t0+t1)
        rhodn(i)=0.5d0*(t0-t1)
      end do
    else
! functionals which require gradients
      do i=1,n
        t0=rho(i)
        t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2+dncgga)
        rhoup(i)=0.5d0*(t0+t1)
        rhodn(i)=0.5d0*(t0-t1)
      end do
    end if
  else
! collinear
    do i=1,n
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
      t0=rho(i)
      t1=mag(i,1)
      rhoup(i)=0.5d0*(t0+t1)
      rhodn(i)=0.5d0*(t0-t1)
    end do
  end if
! call the exchange-correlation interface routine
  if (xcgrad.le.0) then
    call xcifc(xctype_,n,tempa=swidth,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec, &
     vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.1) then
    call ggamt_sp_1(is,n,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    call xcifc(xctype_,n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup,gdn=gdn, &
     g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex,ec=ec,vxup=vxup,&
     vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.2) then
    call ggamt_sp_2a(is,n,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype_,n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
     gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
     dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
     dcdgud=dcdgud)
    call ggamt_sp_2b(is,n,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
     dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
  else if (xcgrad.eq.3) then
    call ggamt_sp_2a(is,n,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype_,n,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn,g2up=g2up, &
     g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tau(:,1),taudn=tau(:,2), &
     vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
     ex(:)=0.d0; ec(:)=0.d0
  else if (xcgrad.eq.4) then
    call ggamt_sp_2a(is,n,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
! enforce the von Weizsacker lower bound
    call k_vwlb(n,rhoup,gup2,tau(:,1))
    call k_vwlb(n,rhodn,gdn2,tau(:,2))
    call xcifc(xctype_,n,rhoup=rhoup,rhodn=rhodn,g2up=g2up,g2dn=g2dn,gup2=gup2,&
     gdn2=gdn2,gupdn=gupdn,tauup=tau(:,1),taudn=tau(:,2),ex=ex,ec=ec,vxup=vxup,&
     vxdn=vxdn,vcup=vcup,vcdn=vcdn,dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud, &
     dcdgu2=dcdgu2,dcdgd2=dcdgd2,dcdgud=dcdgud,dxdg2u=dxdg2u,dxdg2d=dxdg2d, &
     dcdg2u=dcdg2u,dcdg2d=dcdg2d,wxup=wxup,wxdn=wxdn,wcup=wcup,wcdn=wcdn)
    call ggamt_sp_2b(is,n,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
     dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
! determine δτ/δρ using an approximate kinetic energy functional
    call xcifc(ktype,n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2,dtdru=dtdru,&
     dtdrd=dtdrd,dtdgu2=dtdgu2,dtdgd2=dtdgd2)
    call ggamt_4(is,n,gvup,dtdru,dtdgu2,wxup,wcup,vxup,vcup)
    call ggamt_4(is,n,gvdn,dtdrd,dtdgd2,wxdn,wcdn,vxdn,vcdn)
    wxup(:)=0.5d0*(wxup(:)+wxdn(:)+wcup(:)+wcdn(:))
    if (tsh) then
      call rfsht(nr,nri,wxup,wxcmt_(:,ias))
    else
      wxcmt_(1:n,ias)=wxup(1:n)
    end if
  end if
! hybrid functionals
  if (hybrid) then
    t1=1.d0-hybridc
! scale exchange part of energy
    ex(:)=t1*ex(:)
! scale exchange part of potential
    vxup(:)=t1*vxup(:)
    vxdn(:)=t1*vxdn(:)
  end if
  if (ncmag) then
! non-collinear: locally spin rotate the exchange-correlation potential
    do i=1,n
      t1=vxup(i)+vcup(i)
      t2=vxdn(i)+vcdn(i)
      vxc(i)=0.5d0*(t1+t2)
! determine the exchange-correlation magnetic field
      t3=0.5d0*(t1-t2)
! |m| = rhoup - rhodn
      t4=rhoup(i)-rhodn(i)
      if (abs(t4).gt.1.d-8) t4=t3/t4
      bxc(i,1:3)=mag(i,1:3)*t4
    end do
  else
! collinear
    do i=1,n
      t1=vxup(i)+vcup(i)
      t2=vxdn(i)+vcdn(i)
      vxc(i)=0.5d0*(t1+t2)
      bxc(i,1)=0.5d0*(t1-t2)
    end do
  end if
! scale B_xc for if required
  if (tssxc) bxc(:,1:ndmag)=bxc(:,1:ndmag)*sxcscf
  do idm=1,ndmag
    if (tsh) then
! convert field to spherical harmonics
      call rfsht(nr,nri,bxc(:,idm),bxcmt_(:,ias,idm))
    else
      bxcmt_(1:n,ias,idm)=bxc(1:n,idm)
    end if
  end do
else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
  if (xcgrad.le.0) then
    call xcifc(xctype_,n,tempa=swidth,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
  else if (xcgrad.eq.1) then
    call ggamt_1(tsh,is,n,rhomt_(:,ias),grho,g2rho,g3rho)
    call xcifc(xctype_,n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex,ec=ec,&
     vx=vx,vc=vc)
  else if (xcgrad.eq.2) then
    call ggamt_2a(tsh,is,n,rhomt_(:,ias),g2rho,gvrho,grho2)
    call xcifc(xctype_,n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
     dxdgr2=dxdgr2,dcdgr2=dcdgr2)
    call ggamt_2b(is,n,g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
  else if (xcgrad.eq.3) then
    call ggamt_2a(tsh,is,n,rhomt_(:,ias),g2rho,gvrho,grho2)
    call xcifc(xctype_,n,c_tb09=c_tb09,rho=rho,g2rho=g2rho,grho2=grho2,tau=tau,&
     vx=vx,vc=vc)
    ex(:)=0.d0; ec(:)=0.d0
  else if (xcgrad.eq.4) then
    call ggamt_2a(tsh,is,n,rhomt_(:,ias),g2rho,gvrho,grho2)
! enforce the von Weizsacker lower bound
    call k_vwlb(n,rho,grho2,tau)
    call xcifc(xctype_,n,rho=rho,g2rho=g2rho,grho2=grho2,tau=tau,ex=ex,ec=ec, &
     vx=vx,vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2,dxdg2r=dxdg2r,dcdg2r=dcdg2r,wx=wx,&
     wc=wc)
    call ggamt_2b(is,n,g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
! determine δτ/δρ using an approximate kinetic energy functional
    call xcifc(ktype,n,rho=rho,grho2=grho2,dtdr=dtdr,dtdgr2=dtdgr2)
    call ggamt_4(is,n,gvrho,dtdr,dtdgr2,wx,wc,vx,vc)
    wx(:)=wx(:)+wc(:)
    if (tsh) then
      call rfsht(nr,nri,wx,wxcmt_(:,ias))
    else
      wxcmt_(1:n,ias)=wx(1:n)
    end if
  end if
! hybrid functionals
  if (hybrid) then
    t1=1.d0-hybridc
! scale exchange part of energy
    ex(:)=t1*ex(:)
! scale exchange part of potential
    vxc(:)=t1*vx(:)+vc(:)
  else
    vxc(:)=vx(:)+vc(:)
  end if
end if
if (tsh) then
! convert exchange and correlation energy densities to spherical harmonics
  call rfsht(nr,nri,ex,exmt_(:,ias))
  call rfsht(nr,nri,ec,ecmt_(:,ias))
! convert exchange-correlation potential to spherical harmonics
  call rfsht(nr,nri,vxc,vxcmt_(:,ias))
else
  exmt_(1:n,ias)=ex(1:n)
  ecmt_(1:n,ias)=ec(1:n)
  vxcmt_(1:n,ias)=vxc(1:n)
end if
deallocate(rho,ex,ec,vxc)
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) deallocate(tau)
if (spinpol) then
  deallocate(mag,bxc)
  deallocate(rhoup,rhodn,vxup,vxdn,vcup,vcdn)
  if (xcgrad.eq.1) then
    deallocate(grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
  else if (xcgrad.eq.2) then
    deallocate(g2up,g2dn)
    deallocate(gvup,gvdn)
    deallocate(gup2,gdn2,gupdn)
    deallocate(dxdgu2,dxdgd2,dxdgud)
    deallocate(dcdgu2,dcdgd2,dcdgud)
  else if (xcgrad.eq.3) then
    deallocate(g2up,g2dn)
    deallocate(gvup,gvdn)
    deallocate(gup2,gdn2,gupdn)
  else if (xcgrad.eq.4) then
    deallocate(g2up,g2dn)
    deallocate(gvup,gvdn)
    deallocate(gup2,gdn2,gupdn)
    deallocate(dxdgu2,dxdgd2,dxdgud)
    deallocate(dcdgu2,dcdgd2,dcdgud)
    deallocate(dxdg2u,dxdg2d)
    deallocate(dcdg2u,dcdg2d)
    deallocate(dtdru,dtdrd)
    deallocate(dtdgu2,dtdgd2)
    deallocate(wxup,wxdn,wcup,wcdn)
  end if
else
  deallocate(vx,vc)
  if (xcgrad.eq.1) then
    deallocate(grho,g2rho,g3rho)
  else if (xcgrad.eq.2) then
    deallocate(g2rho,gvrho,grho2)
    deallocate(dxdgr2,dcdgr2)
  else if (xcgrad.eq.3) then
    deallocate(g2rho,gvrho,grho2)
  else if (xcgrad.eq.4) then
    deallocate(g2rho,gvrho,grho2)
    deallocate(dxdgr2,dcdgr2,dxdg2r,dcdg2r)
    deallocate(dtdr,dtdgr2)
    deallocate(wx,wc)
  end if
end if
end subroutine

