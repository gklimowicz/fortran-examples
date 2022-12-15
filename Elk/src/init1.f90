
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init1
! !INTERFACE:
subroutine init1
! !USES:
use modmain
use moddftu
use modulr
use modtddft
use modtest
use modvars
! !DESCRIPTION:
!   Generates the $k$-point set and then allocates and initialises global
!   variables which depend on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
logical lsym(48)
integer is,ias,nppt
integer io,ilo,i1,i2,i3
integer ik,isym,jspn
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3,n
real(8) vl(3),vc(3),t1
real(8) boxl(3,0:3)
real(8) ts0,ts1
! external functions
complex(8), external :: gauntyry

call timesec(ts0)

!---------------------!
!     k-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) then
  ngridk(:)=1
  vkloff(:)=0.d0
  autokpt=.false.
end if
! store the point group symmetries for reducing the k-point set
if (reducek.eq.0) then
  nsymkpt=1
  symkpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    if (reducek.eq.2) then
! check symmetry is symmorphic
      if (.not.tv0symc(isym)) goto 10
! check also that the spin rotation is the same as the spatial rotation
      if (spinpol) then
        if (lspnsymc(isym).ne.lsplsymc(isym)) goto 10
      end if
    end if
    lsym(lsplsymc(isym))=.true.
10 continue
  end do
  nsymkpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymkpt=nsymkpt+1
      symkpt(:,:,nsymkpt)=symlat(:,:,isym)
    end if
  end do
end if
if (any(task.eq.[20,21,22,23])) then
! generate k-points along a path for band structure plots
  call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
  nkpt=npp1d
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkpt))
  do ik=1,nkpt
    vkl(:,ik)=vplp1d(:,ik)
    call r3mv(bvec,vkl(:,ik),vkc(:,ik))
  end do
  nkptnr=nkpt
else if (task.eq.25) then
! effective mass calculation
  nkpt=(2*ndspem+1)**3
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3,nkpt))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkpt))
! map vector to [0,1)
  call r3frac(epslat,vklem)
  ik=0
  do i3=-ndspem,ndspem
    do i2=-ndspem,ndspem
      do i1=-ndspem,ndspem
        ik=ik+1
        ivk(1,ik)=i1; ivk(2,ik)=i2; ivk(3,ik)=i3
        vc(1)=dble(i1); vc(2)=dble(i2); vc(3)=dble(i3)
        vc(:)=vc(:)*deltaem
        call r3mv(binv,vc,vl)
        vkl(:,ik)=vklem(:)+vl(:)
        call r3mv(bvec,vkl(:,ik),vkc(:,ik))
      end do
    end do
  end do
  nkptnr=nkpt
else
! determine the k-point grid automatically from radkpt if required
  if (autokpt) then
    t1=radkpt/twopi
    ngridk(:)=int(t1*sqrt(bvec(1,:)**2+bvec(2,:)**2+bvec(3,:)**2))+1
  end if
! set up the default k-point box
  boxl(:,0)=vkloff(:)/dble(ngridk(:))
  if (task.eq.102) boxl(:,0)=0.d0
  boxl(:,1)=boxl(:,0)
  boxl(:,2)=boxl(:,0)
  boxl(:,3)=boxl(:,0)
  boxl(1,1)=boxl(1,1)+1.d0
  boxl(2,2)=boxl(2,2)+1.d0
  boxl(3,3)=boxl(3,3)+1.d0
! k-point set and box for Fermi surface plots
  if (any(task.eq.[100,101,102])) then
    ngridk(:)=np3d(:)
    if (task.ne.102) boxl(:,:)=vclp3d(:,:)
  end if
! allocate the k-point set arrays
  if (allocated(ivkik)) deallocate(ivkik)
  allocate(ivkik(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
  if (allocated(ivkiknr)) deallocate(ivkiknr)
  allocate(ivkiknr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
  nkptnr=ngridk(1)*ngridk(2)*ngridk(3)
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3,nkptnr))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkptnr))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkptnr))
  if (allocated(wkpt)) deallocate(wkpt)
  allocate(wkpt(nkptnr))
! generate the k-point set
  call genppts(.false.,nsymkpt,symkpt,ngridk,nkptnr,epslat,bvec,boxl,nkpt, &
   ivkik,ivkiknr,ivk,vkl,vkc,wkpt,wkptnr)
! write to VARIABLES.OUT
  if (wrtvars) then
    call writevars('nsymkpt',iv=nsymkpt)
    call writevars('symkpt',nv=9*nsymkpt,iva=symkpt)
    call writevars('ngridk',nv=3,iva=ngridk)
    call writevars('vkloff',nv=3,rva=vkloff)
    call writevars('nkpt',iv=nkpt)
    call writevars('ivkik',nv=nkptnr,iva=ivkik)
    call writevars('ivk',nv=3*nkptnr,iva=ivk)
    call writevars('vkl',nv=3*nkptnr,rva=vkl)
    call writevars('wkpt',nv=nkpt,rva=wkpt)
  end if
end if
if (any(task.eq.[700,701,731,732,733,741,742,743,771,772,773])) then
! generate ultracell reciprocal lattice vectors if required
  call reciplat(avecu,bvecu,omegau,omegabzu)
! generate the kappa, k+kappa and Q-points if required
  call genkpakq
end if
! write the k-points to test file
call writetest(910,'k-points (Cartesian)',nv=3*nkpt,tol=1.d-8,rva=vkc)

!---------------------!
!     G+k-vectors     !
!---------------------!
if ((xctype(1).lt.0).or.tddos.or. &
 (any(task.eq.[5,10,205,300,600,620,630]))) then
  nppt=nkptnr
else
  nppt=nkpt
end if
! find the maximum number of G+k-vectors
call findngkmax(nkpt,vkc,nspnfv,vqcss,ngvc,vgc,gkmax,ngkmax)
! allocate the G+k-vector arrays
if (allocated(ngk)) deallocate(ngk)
allocate(ngk(nspnfv,nppt))
if (allocated(igkig)) deallocate(igkig)
allocate(igkig(ngkmax,nspnfv,nppt))
if (allocated(vgkl)) deallocate(vgkl)
allocate(vgkl(3,ngkmax,nspnfv,nppt))
if (allocated(vgkc)) deallocate(vgkc)
allocate(vgkc(3,ngkmax,nspnfv,nppt))
if (allocated(gkc)) deallocate(gkc)
allocate(gkc(ngkmax,nspnfv,nppt))
if (allocated(sfacgk)) deallocate(sfacgk)
allocate(sfacgk(ngkmax,natmtot,nspnfv,nppt))
do ik=1,nppt
  do jspn=1,nspnfv
    vl(:)=vkl(:,ik)
    vc(:)=vkc(:,ik)
! spin-spiral case
    if (spinsprl) then
      if (jspn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! generate the G+k-vectors
    call gengkvec(ngvc,ivg,vgc,vl,vc,gkmax,ngkmax,ngk(jspn,ik), &
     igkig(:,jspn,ik),vgkl(:,:,jspn,ik),vgkc(:,:,jspn,ik),gkc(:,jspn,ik))
! generate structure factors for G+k-vectors
    call gensfacgp(ngk(jspn,ik),vgkc(:,:,jspn,ik),ngkmax,sfacgk(:,:,jspn,ik))
  end do
end do
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('nspnfv',iv=nspnfv)
  call writevars('ngk',nv=nspnfv*nkpt,iva=ngk)
  do ik=1,nkpt
    do jspn=1,nspnfv
      call writevars('igkig',n1=jspn,n2=ik,nv=ngk(jspn,ik),iva=igkig(:,jspn,ik))
    end do
  end do
end if

!---------------------------------!
!     APWs and local-orbitals     !
!---------------------------------!
apwordmax=0
lorbordmax=0
nlomax=0
lolmax=0
do is=1,nspecies
  lmoapw(is)=0
  do l1=0,lmaxapw
! find the maximum APW order
    apwordmax=max(apwordmax,apword(l1,is))
! find total number of APW coefficients (l, m and order)
    lmoapw(is)=lmoapw(is)+(2*l1+1)*apword(l1,is)
    if (l1.eq.lmaxo) nlmwf(is)=lmoapw(is)
  end do
! find the maximum number of local-orbitals
  nlomax=max(nlomax,nlorb(is))
! find the maximum local-orbital order and angular momentum
  n=0
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    lolmax=max(lolmax,l1)
    lorbordmax=max(lorbordmax,lorbord(ilo,is))
    n=n+2*l1+1
  end do
! number of (l,m) components used for generating the muffin-tin wavefunctions
  nlmwf(is)=max(nlmwf(is),n)
end do
lolmmax=(lolmax+1)**2
! polynomial order used for APW and local-orbital radial derivatives
npapw=max(apwordmax+1,4)
nplorb=max(lorbordmax+1,4)
! set the APW and local-orbital linearisation energies to the default
if (allocated(apwe)) deallocate(apwe)
allocate(apwe(apwordmax,0:lmaxapw,natmtot))
if (allocated(lorbe)) deallocate(lorbe)
allocate(lorbe(lorbordmax,maxlorb,natmtot))
do ias=1,natmtot
  is=idxis(ias)
  do l1=0,lmaxapw
    do io=1,apword(l1,is)
      apwe(io,l1,ias)=apwe0(io,l1,is)
    end do
  end do
  do ilo=1,nlorb(is)
    do io=1,lorbord(ilo,is)
      lorbe(io,ilo,ias)=lorbe0(io,ilo,is)
    end do
  end do
end do
! generate the local-orbital index
call genidxlo
! allocate radial function arrays
if (allocated(apwfr)) deallocate(apwfr)
allocate(apwfr(nrmtmax,2,apwordmax,0:lmaxapw,natmtot))
if (allocated(apwdfr)) deallocate(apwdfr)
allocate(apwdfr(apwordmax,0:lmaxapw,natmtot))
if (allocated(lofr)) deallocate(lofr)
allocate(lofr(nrmtmax,2,nlomax,natmtot))

!-------------------------!
!     DFT+U variables     !
!-------------------------!
if (dftu.ne.0) then
! allocate energy arrays to calculate Slater integrals with Yukawa potential
  if (allocated(efdu)) deallocate(efdu)
  allocate(efdu(0:lmaxdm,natmtot))
! allocate radial functions to calculate Slater integrals with Yukawa potential
  if (allocated(fdufr)) deallocate(fdufr)
  allocate(fdufr(nrmtmax,0:lmaxdm,natmtot))
end if

!---------------------------------------!
!     eigenvalue equation variables     !
!---------------------------------------!
! total number of empty states (M. Meinert)
nempty=nint(nempty0*max(natmtot,1))
if (nempty.lt.1) nempty=1
! number of first-variational states
nstfv=int(chgval/2.d0)+nempty+1
! overlap and Hamiltonian matrix sizes
if (allocated(nmat)) deallocate(nmat)
allocate(nmat(nspnfv,nkpt))
nmatmax=0
do ik=1,nkpt
  do jspn=1,nspnfv
    n=ngk(jspn,ik)+nlotot
    if (nstfv.gt.n) then
      write(*,*)
      write(*,'("Error(init1): number of first-variational states larger than &
       &matrix size")')
      write(*,'("Increase rgkmax or decrease nempty")')
      write(*,*)
      stop
    end if
    nmat(jspn,ik)=n
    nmatmax=max(nmatmax,n)
  end do
end do
! number of second-variational states
nstsv=nstfv*nspinor
! allocate second-variational arrays
if (allocated(evalsv)) deallocate(evalsv)
allocate(evalsv(nstsv,nkpt))
if (allocated(occsv)) deallocate(occsv)
allocate(occsv(nstsv,nkpt))
occsv(:,:)=0.d0
! allocate overlap and Hamiltonian integral arrays
if (allocated(oalo)) deallocate(oalo)
allocate(oalo(apwordmax,nlomax,natmtot))
if (allocated(ololo)) deallocate(ololo)
allocate(ololo(nlomax,nlomax,natmtot))
if (allocated(haa)) deallocate(haa)
allocate(haa(lmmaxo,apwordmax,0:lmaxapw,apwordmax,0:lmaxapw,natmtot))
if (allocated(hloa)) deallocate(hloa)
allocate(hloa(lmmaxo,apwordmax,0:lmaxapw,nlomax,natmtot))
if (allocated(hlolo)) deallocate(hlolo)
allocate(hlolo(lmmaxo,nlomax,nlomax,natmtot))
! allocate and generate complex Gaunt coefficient array
if (allocated(gntyry)) deallocate(gntyry)
allocate(gntyry(lmmaxo,lmmaxapw,lmmaxapw))
do l1=0,lmaxapw
  do m1=-l1,l1
    lm1=l1*(l1+1)+m1+1
    do l3=0,lmaxapw
      do m3=-l3,l3
        lm3=l3*(l3+1)+m3+1
        do l2=0,lmaxo
          do m2=-l2,l2
            lm2=l2*(l2+1)+m2+1
            gntyry(lm2,lm3,lm1)=gauntyry(l1,l2,l3,m1,m2,m3)
          end do
        end do
      end do
    end do
  end do
end do
! write to VARIABLES.OUT
if (wrtvars) then
  call writevars('nempty',iv=nempty)
  call writevars('nstfv',iv=nstfv)
  call writevars('nlotot',iv=nlotot)
  call writevars('nstsv',iv=nstsv)
end if

call timesec(ts1)
timeinit=timeinit+ts1-ts0

end subroutine
!EOC

