
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potxcir(xctype_,rhoir_,magir_,tauir_,exir_,ecir_,vxcir_,bxcir_, &
 wxcir_)
use modmain
use modxcifc
implicit none
! arguments
integer, intent(in) :: xctype_(3)
real(8), intent(in) :: rhoir_(ngtot),magir_(ngtot,ndmag),tauir_(ngtot,nspinor)
real(8), intent(out) :: exir_(ngtot),ecir_(ngtot)
real(8), intent(out) :: vxcir_(ngtot),bxcir_(ngtot,ndmag),wxcir_(ngtot)
! local variables
integer n,i
real(8) t0,t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: rhoup(:),rhodn(:)
real(8), allocatable :: gvrho(:,:),gvup(:,:),gvdn(:,:)
real(8), allocatable :: grho(:),gup(:),gdn(:)
real(8), allocatable :: g2rho(:),g2up(:),g2dn(:)
real(8), allocatable :: g3rho(:),g3up(:),g3dn(:)
real(8), allocatable :: grho2(:),gup2(:),gdn2(:),gupdn(:)
real(8), allocatable :: vx(:),vxup(:),vxdn(:)
real(8), allocatable :: vc(:),vcup(:),vcdn(:)
real(8), allocatable :: dxdgr2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
real(8), allocatable :: dcdgr2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
real(8), allocatable :: dxdg2r(:),dxdg2u(:),dxdg2d(:)
real(8), allocatable :: dcdg2r(:),dcdg2u(:),dcdg2d(:)
real(8), allocatable :: dtdr(:),dtdru(:),dtdrd(:)
real(8), allocatable :: dtdgr2(:),dtdgu2(:),dtdgd2(:)
real(8), allocatable :: wx(:),wxup(:),wxdn(:)
real(8), allocatable :: wc(:),wcup(:),wcdn(:)
n=ngtot
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
if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
  if (ncmag) then
! non-collinear
    if (xcgrad.eq.0) then
! LSDA
      do i=1,n
        t0=rhoir_(i)
        t1=sqrt(magir_(i,1)**2+magir_(i,2)**2+magir_(i,3)**2)*sxcscf
        rhoup(i)=0.5d0*(t0+t1)
        rhodn(i)=0.5d0*(t0-t1)
      end do
    else
! functionals which require gradients
      do i=1,n
        t0=rhoir_(i)
        t1=sqrt(magir_(i,1)**2+magir_(i,2)**2+magir_(i,3)**2+dncgga)*sxcscf
        rhoup(i)=0.5d0*(t0+t1)
        rhodn(i)=0.5d0*(t0-t1)
      end do
    end if
  else
! collinear
    do i=1,n
      t0=rhoir_(i)
      t1=magir_(i,1)*sxcscf
      rhoup(i)=0.5d0*(t0+t1)
      rhodn(i)=0.5d0*(t0-t1)
    end do
  end if
  if (xcgrad.le.0) then
    call xcifc(xctype_,n,tempa=swidth,rhoup=rhoup,rhodn=rhodn,ex=exir_, &
     ec=ecir_,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.1) then
    call ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    call xcifc(xctype_,n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup,gdn=gdn, &
     g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir_,ec=ecir_, &
     vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.2) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype_,n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
     gupdn=gupdn,ex=exir_,ec=ecir_,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
     dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
     dcdgud=dcdgud)
    call ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
     dxdgud,dcdgu2,dcdgd2,dcdgud)
  else if (xcgrad.eq.3) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype_,n,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn,g2up=g2up, &
     g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tauir_(:,1), &
     taudn=tauir_(:,2),vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
    exir_(:)=0.d0; ecir_(:)=0.d0
  else if (xcgrad.eq.4) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
! enforce the von Weizsacker lower bound
    call k_vwlb(n,rhoup,gup2,tauir_(:,1))
    call k_vwlb(n,rhodn,gdn2,tauir_(:,2))
    call xcifc(xctype_,n,rhoup=rhoup,rhodn=rhodn,g2up=g2up,g2dn=g2dn, &
     gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tauir_(:,1),taudn=tauir_(:,2), &
     ex=exir_,ec=ecir_,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn,dxdgu2=dxdgu2, &
     dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2,dcdgud=dcdgud, &
     dxdg2u=dxdg2u,dxdg2d=dxdg2d,dcdg2u=dcdg2u,dcdg2d=dcdg2d,wxup=wxup, &
     wxdn=wxdn,wcup=wcup,wcdn=wcdn)
    call ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
     dxdgud,dcdgu2,dcdgd2,dcdgud)
! determine δτ/δρ using an approximate kinetic energy functional
    call xcifc(ktype,n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2,dtdru=dtdru,&
     dtdrd=dtdrd,dtdgu2=dtdgu2,dtdgd2=dtdgd2)
    call ggair_4(gvup,dtdru,dtdgu2,wxup,wcup,vxup,vcup)
    call ggair_4(gvdn,dtdrd,dtdgd2,wxdn,wcdn,vxdn,vcdn)
    wxcir_(:)=0.5d0*(wxup(:)+wxdn(:)+wcup(:)+wcdn(:))
  end if
! hybrid functionals
  if (hybrid) then
    t1=1.d0-hybridc
! scale exchange part of energy
    exir_(:)=t1*exir_(:)
! scale exchange part of potential
    vxup(:)=t1*vxup(:)
    vxdn(:)=t1*vxdn(:)
  end if
  if (ncmag) then
! non-collinear: spin rotate the local exchange potential
    do i=1,n
      t1=vxup(i)+vcup(i)
      t2=vxdn(i)+vcdn(i)
      vxcir_(i)=0.5d0*(t1+t2)
! determine the exchange-correlation magnetic field
      t3=0.5d0*(t1-t2)
      t4=rhoup(i)-rhodn(i)
      if (abs(t4).gt.1.d-8) t4=t3/t4
      bxcir_(i,:)=magir_(i,:)*t4
    end do
  else
! collinear
    do i=1,n
      t1=vxup(i)+vcup(i)
      t2=vxdn(i)+vcdn(i)
      vxcir_(i)=0.5d0*(t1+t2)
      bxcir_(i,1)=0.5d0*(t1-t2)
    end do
  end if
! scale field if required
  if (tssxc) bxcir_(:,1:ndmag)=bxcir_(:,1:ndmag)*sxcscf
else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
  if (xcgrad.le.0) then
    call xcifc(xctype_,n,tempa=swidth,rho=rhoir_,ex=exir_,ec=ecir_,vx=vx,vc=vc)
  else if (xcgrad.eq.1) then
    call ggair_1(rhoir_,grho,g2rho,g3rho)
    call xcifc(xctype_,n,rho=rhoir_,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=exir_,&
     ec=ecir_,vx=vx,vc=vc)
  else if (xcgrad.eq.2) then
    call ggair_2a(rhoir_,g2rho,gvrho,grho2)
    call xcifc(xctype_,n,rho=rhoir_,grho2=grho2,ex=exir_,ec=ecir_,vx=vx,vc=vc, &
     dxdgr2=dxdgr2,dcdgr2=dcdgr2)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
  else if (xcgrad.eq.3) then
    call ggair_2a(rhoir_,g2rho,gvrho,grho2)
    call xcifc(xctype_,n,c_tb09=c_tb09,rho=rhoir_,g2rho=g2rho,grho2=grho2, &
     tau=tauir_,vx=vx,vc=vc)
    exir_(:)=0.d0; ecir_(:)=0.d0
  else if (xcgrad.eq.4) then
    call ggair_2a(rhoir_,g2rho,gvrho,grho2)
! enforce the von Weizsacker lower bound
    call k_vwlb(n,rhoir_,grho2,tauir_)
    call xcifc(xctype_,n,rho=rhoir_,g2rho=g2rho,grho2=grho2,tau=tauir_, &
     ex=exir_,ec=ecir_,vx=vx,vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2,dxdg2r=dxdg2r, &
     dcdg2r=dcdg2r,wx=wx,wc=wc)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
! determine δτ/δρ using an approximate kinetic energy functional
    call xcifc(ktype,n,rho=rhoir_,grho2=grho2,dtdr=dtdr,dtdgr2=dtdgr2)
    call ggair_4(gvrho,dtdr,dtdgr2,wx,wc,vx,vc)
    wxcir_(:)=wx(:)+wc(:)
  end if
! hybrid functionals
  if (hybrid) then
    t1=1.d0-hybridc
! scale exchange part of energy
    exir_(:)=t1*exir_(:)
! scale exchange part of potential
    vxcir_(:)=t1*vx(:)+vc(:)
  else
    vxcir_(:)=vx(:)+vc(:)
  end if
end if
if (spinpol) then
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

