
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: atom
! !INTERFACE:
subroutine atom(sol,ptnucl,zn,nst,n,l,k,occ,xctype,xcgrad,nr,r,eval,rho,vr,rwf)
! !USES:
use modxcifc
! !INPUT/OUTPUT PARAMETERS:
!   sol    : speed of light in atomic units (in,real)
!   ptnucl : .true. if the nucleus is a point particle (in,logical)
!   zn     : nuclear charge (in,real)
!   nst    : number of states to solve for (in,integer)
!   n      : priciple quantum number of each state (in,integer(nst))
!   l      : quantum number l of each state (in,integer(nst))
!   k      : quantum number k (l or l+1) of each state (in,integer(nst))
!   occ    : occupancy of each state (inout,real(nst))
!   xctype : exchange-correlation type (in,integer(3))
!   xcgrad : 1 for GGA functional, 0 otherwise (in,integer)
!   nr     : number of radial mesh points (in,integer)
!   r      : radial mesh (in,real(nr))
!   eval   : eigenvalue without rest-mass energy for each state (out,real(nst))
!   rho    : charge density (out,real(nr))
!   vr     : self-constistent potential (out,real(nr))
!   rwf    : major and minor components of radial wavefunctions for each state
!            (out,real(nr,2,nst))
! !DESCRIPTION:
!   Solves the Dirac-Kohn-Sham equations for an atom using the
!   exchange-correlation functional {\tt xctype} and returns the self-consistent
!   radial wavefunctions, eigenvalues, charge densities and potentials. Requires
!   the exchange-correlation interface routine {\tt xcifc}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Fixed s.c. convergence problem, October 2003 (JKD)
!   Added support for GGA functionals, June 2006 (JKD)
!
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: sol
logical, intent(in) :: ptnucl
real(8), intent(in) :: zn
integer, intent(in) :: nst
integer, intent(in) :: n(nst),l(nst),k(nst)
real(8), intent(inout) :: occ(nst)
integer, intent(in) :: xctype(3),xcgrad
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(out) :: eval(nst)
real(8), intent(out) :: rho(nr),vr(nr)
real(8), intent(out) :: rwf(nr,2,nst)
! local variables
integer, parameter :: maxscl=200
integer ir,ist,iscl
real(8), parameter :: fourpi=12.566370614359172954d0
! potential convergence tolerance
real(8), parameter :: eps=1.d-6
real(8) sm,dv,dvp,ze,beta,t1
! allocatable arrays
real(8), allocatable :: vn(:),vh(:),ex(:),ec(:),vx(:),vc(:),vrp(:)
real(8), allocatable :: ri(:),wpr(:,:),fr1(:),fr2(:),gr1(:),gr2(:)
real(8), allocatable :: grho(:),g2rho(:),g3rho(:)
if (nst.le.0) then
  write(*,*)
  write(*,'("Error(atom): invalid nst : ",I8)') nst
  write(*,*)
  stop
end if
! allocate local arrays
allocate(vn(nr),vh(nr),ex(nr),ec(nr),vx(nr),vc(nr),vrp(nr))
allocate(ri(nr),wpr(4,nr),fr1(nr),fr2(nr),gr1(nr),gr2(nr))
if (xcgrad.eq.1) then
  allocate(grho(nr),g2rho(nr),g3rho(nr))
end if
! find total electronic charge
ze=0.d0
do ist=1,nst
  ze=ze+occ(ist)
end do
! set up nuclear potential
call potnucl(ptnucl,nr,r,zn,vn)
do ir=1,nr
  ri(ir)=1.d0/r(ir)
! initialise the Kohn-Sham potential to the nuclear potential
  vr(ir)=vn(ir)
end do
! determine the weights for radial integration
call wsplintp(nr,r,wpr)
dvp=0.d0
vrp(:)=0.d0
! initialise mixing parameter
beta=0.5d0
! initialise eigenvalues to relativistic values (minus the rest mass energy)
do ist=1,nst
  t1=sqrt(dble(k(ist)**2)-(zn/sol)**2)
  t1=(dble(n(ist)-abs(k(ist)))+t1)**2
  t1=1.d0+((zn/sol)**2)/t1
  eval(ist)=sol**2/sqrt(t1)-sol**2
end do
! start of self-consistent loop
do iscl=1,maxscl
! solve the Dirac equation for each state
!$OMP PARALLEL DO DEFAULT(SHARED)
  do ist=1,nst
    call rdirac(sol,n(ist),l(ist),k(ist),nr,r,vr,eval(ist),rwf(:,1,ist), &
     rwf(:,2,ist))
  end do
!$OMP END PARALLEL DO
! compute the charge density
  do ir=1,nr
    sm=0.d0
    do ist=1,nst
      sm=sm+occ(ist)*(rwf(ir,1,ist)**2+rwf(ir,2,ist)**2)
    end do
    fr1(ir)=sm
    fr2(ir)=sm*ri(ir)
    rho(ir)=(1.d0/fourpi)*sm*ri(ir)**2
  end do
  call splintwp(nr,wpr,fr1,gr1)
  call splintwp(nr,wpr,fr2,gr2)
! find the Hartree potential
  t1=gr2(nr)
  do ir=1,nr
    vh(ir)=gr1(ir)*ri(ir)+t1-gr2(ir)
  end do
! normalise charge density and potential
  t1=ze/gr1(nr)
  rho(:)=t1*rho(:)
  vh(:)=t1*vh(:)
! compute the exchange-correlation energy and potential
  if (xcgrad.eq.1) then
! GGA functional
! |grad rho|
    call fderiv(1,nr,r,rho,grho)
! grad^2 rho
    call fderiv(2,nr,r,rho,g2rho)
    do ir=1,nr
      g2rho(ir)=g2rho(ir)+2.d0*ri(ir)*grho(ir)
    end do
! approximate (grad rho).(grad |grad rho|)
    do ir=1,nr
      g3rho(ir)=grho(ir)*g2rho(ir)
    end do
    call xcifc(xctype,nr,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex,ec=ec,&
     vx=vx,vc=vc)
  else
! LDA functional
    call xcifc(xctype,nr,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
  end if
! self-consistent potential
  vr(:)=vh(:)+vx(:)+vc(:)
! determine change in potential
  sm=0.d0
  do ir=1,nr
    sm=sm+(vr(ir)-vrp(ir))**2
  end do
  dv=sqrt(sm)/dble(nr)
  if (iscl.gt.2) then
! reduce beta if change in potential is diverging
    if (dv.gt.dvp) beta=beta*0.8d0
    beta=max(beta,0.01d0)
  end if
  dvp=dv
  do ir=1,nr
! mix old and new potentials
    vr(ir)=(1.d0-beta)*vrp(ir)+beta*vr(ir)
    vrp(ir)=vr(ir)
! add nuclear potential
    vr(ir)=vr(ir)+vn(ir)
  end do
! check for convergence
  if ((iscl.gt.2).and.(dv.lt.eps)) goto 10
! end self-consistent loop
end do
write(*,*)
write(*,'("Warning(atom): maximum iterations exceeded")')
10 continue
deallocate(vn,vh,ex,ec,vx,vc,vrp)
deallocate(ri,wpr,fr1,fr2,gr1,gr2)
if (xcgrad.eq.1) deallocate(grho,g2rho,g3rho)
return

contains

pure subroutine splintwp(n,wp,f,g)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wp(*),f(n)
real(8), intent(out) :: g(n)
! local variables
integer i,j
real(8) sm
g(1)=0.d0
sm=wp(5)*f(1)+wp(6)*f(2)+wp(7)*f(3)+wp(8)*f(4)
g(2)=sm
do i=2,n-2
  j=i*4+1
  sm=sm+wp(j)*f(i-1)+wp(j+1)*f(i)+wp(j+2)*f(i+1)+wp(j+3)*f(i+2)
  g(i+1)=sm
end do
j=(n-1)*4+1
g(n)=sm+wp(j)*f(n-3)+wp(j+1)*f(n-2)+wp(j+2)*f(n-1)+wp(j+3)*f(n)
end subroutine

end subroutine
!EOC

