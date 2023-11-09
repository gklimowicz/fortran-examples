
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom. This file is distributed under the terms of the GNU
! General Public License. See the file COPYING for license details.

subroutine setupw90
use modmain
use modw90
implicit none
! local variables
integer is,ia,ias
real(8) real_lattice(3,3),recip_lattice(3,3)
! allocatable arrays
integer, allocatable :: proj_l(:),proj_m(:),proj_radial(:),proj_s(:)
integer, allocatable :: exclude_bands(:)
real(8), allocatable :: atoms_cart(:,:),proj_site(:,:)
real(8), allocatable :: proj_z(:,:),proj_x(:,:)
real(8), allocatable :: proj_zona(:),proj_s_qaxis(:,:)
character(256), allocatable :: atom_symbols(:)
! allocate global arrays
if (allocated(nnlist)) deallocate(nnlist)
allocate(nnlist(nkptnr,num_nnmax))
if (allocated(nncell)) deallocate(nncell)
allocate(nncell(3,nkptnr,num_nnmax))
! allocate local arrays
allocate(proj_l(num_bands),proj_m(num_bands))
allocate(proj_radial(num_bands),proj_s(num_bands))
allocate(exclude_bands(num_bands))
allocate(atoms_cart(3,natmtot),proj_site(3,num_bands))
allocate(proj_z(3,num_bands),proj_x(3,num_bands))
allocate(proj_zona(num_bands),proj_s_qaxis(3,num_bands))
allocate(atom_symbols(natmtot))
real_lattice=br_ang*transpose(avec)
recip_lattice=(1.d0/br_ang)*transpose(bvec)
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  atom_symbols(ias)=trim(spsymb(is))
  atoms_cart(:,ias)=br_ang*atposc(:,ia,is)
end do
call wannier_setup(seedname,ngridk,nkptnr,real_lattice,recip_lattice,vkl, &
 num_bands,natmtot,atom_symbols,atoms_cart,.false.,spinpol,nntot,nnlist, &
 nncell,num_bands,num_wann,proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x, &
 proj_zona,exclude_bands,proj_s,proj_s_qaxis)
deallocate(proj_l,proj_m)
deallocate(proj_radial,proj_s,exclude_bands)
deallocate(atoms_cart,proj_site,proj_z,proj_x)
deallocate(proj_zona,proj_s_qaxis,atom_symbols)
end subroutine

