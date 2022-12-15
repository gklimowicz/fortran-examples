
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotcoul
! !INTERFACE:
subroutine zpotcoul(nr,nri,np,ld1,rl,ngdg,igf,ngp,gpc,gclgp,ld2,jlgprmt,ylmgp, &
 sfacgp,zrhoir,ld3,zvclmt,zvclir)
! !USES:
use modmain
use modphonon
! !INPUT/OUTPUT PARAMETERS:
!   nr      : number of radial points for each species (in,integer(nspecies))
!   nri     : number of radial points on inner part (in,integer(nspecies))
!   np      : total number of points in muffin-tins (in,integer(nspecies))
!   ld1     : leading dimension (in,integer)
!   rl      : r^l on radial mesh for each species
!             (in,real(ld1,-lmaxo-1:lmaxo+2,nspecies))
!   ngdg    : G-vector grid sizes (in,integer(3))
!   igf     : map from G-vector index to FFT array (in,integer(*))
!   ngp     : number of G+p-vectors (in,integer)
!   gpc     : G+p-vector lengths (in,real(ngp))
!   gclgp   : Coulomb Green's function in G+p-space (in,real(ngp))
!   ld2     : leading dimension (in,integer)
!   jlgprmt : spherical Bessel functions for evergy G+p-vector and muffin-tin
!             radius (in,real(0:lnpsd,ld2,nspecies))
!   ylmgp   : spherical harmonics of the G+p-vectors (in,complex(lmmaxo,ngp))
!   sfacgp  : structure factors of the G+p-vectors (in,complex(ld2,natmtot))
!   zrhoir  : interstitial charge density (in,complex(*))
!   ld3     : leading dimension (in,integer)
!   zvclmt  : muffin-tin Coulomb potential, with the contribution from the
!             isolated muffin-tin density precalculated and passed in
!             (inout,complex(ld3,natmtot))
!   zvclir  : interstitial Coulomb potential (out,complex(*))
! !DESCRIPTION:
!   Calculates the Coulomb potential of a complex charge density by solving
!   Poisson's equation using the method of M. Weinert, {\it J. Math. Phys.}
!   {\bf 22}, 2433 (1981). First, the multipole moments of the muffin-tin charge
!   are determined for the $j$th atom of the $i$th species by
!   $$ q_{ij;lm}^{\rm MT}=\int_0^{R_i}r^{l+2}\rho_{ij;lm}(r)dr+z_{ij}Y_{00}
!   \,\delta_{l,0}\;, $$
!   where $R_i$ is the muffin-tin radius and $z_{ij}$ is a point charge located
!   at the atom center (usually the nuclear charge, which should be taken as
!   {\bf negative}). Next, the multipole moments of the continuation of the
!   interstitial density, $\rho^{\rm I}$, into the muffin-tin are found with
!   $$ q_{ij;lm}^{\rm I}=4\pi i^l R_i^{l+3}\sum_{\bf G}\frac{j_{l+1}(GR_i)}
!    {GR_i}\rho^{\rm I}({\bf G})\exp(i{\bf G}\cdot{\bf r}_{ij})Y_{lm}^*
!    (\hat{\bf G}), $$
!   remembering that
!   $$ \lim_{x\rightarrow 0}\frac{j_{l+n}(x)}{x^n}=\frac{1}{(2n+1)!!}
!    \delta_{l,0} $$
!   should be used for the case ${\bf G}=0$. A pseudocharge is now constructed
!   which is equal to the real density in the interstitial region and whose
!   multipoles are the difference between the real and interstitial muffin-tin
!   multipoles. This pseudocharge density is smooth in the sense that it can be
!   expanded in terms of the finite set of ${\bf G}$-vectors. In each muffin-tin
!   the pseudocharge has the form
!   $$ \rho_{ij}^{\rm P}({\bf r})=\rho^{\rm I}({\bf r}-{\bf r}_{ij})+\sum_{lm}
!    \rho_{ij;lm}^{\rm P}\frac{1}{R_i^{l+3}}\left(\frac{r}{R_i}\right)^l\left(1-
!    \frac{r^2}{R_i^2}\right)^{N_i}Y_{lm}(\hat{\bf r}) $$
!   where
!   $$ \rho_{ij;lm}^{\rm P}=\frac{(2l+2N_i+3)!!}{2^N_iN_i!(2l+1)!!}\left(
!    q_{ij;lm}^{\rm MT}-q_{ij;lm}^{\rm I}\right) $$
!   and $N_i\approx\frac{1}{4}R_iG_{\rm max}$ is generally a good choice.
!   The pseudocharge in reciprocal space is given by
!   $$ \rho^{\rm P}({\bf G})=\rho^{\rm I}({\bf G})+\sum_{ij;lm}2^{N_i}N_i!
!    \frac{4\pi(-i)^l}{\Omega R_i^l}\frac{j_{l+N_i+1}(GR_i)}{(GR_i)^{N_i+1}}
!    \rho_{ij;lm}^{\rm P}\exp(-i{\bf G}\cdot{\bf r}_{ij})Y_{lm}(\hat{\bf G}) $$
!   which may be used for solving Poisson's equation directly
!   $$ V^{\rm P}({\bf G})=\begin{cases}
!     4\pi\frac{\rho^{\rm P}({\bf G})}{G^2} & G>0 \\
!     0 & G=0 \end{cases}\;. $$
!   The usual Green's function approach is then employed to determine the
!   potential in the muffin-tin sphere due to charge in the sphere. In other
!   words
!   $$ V_{ij;lm}^{\rm MT}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r
!    \rho_{ij;lm}^{\rm MT}(r'){r'}^{l+2}dr'+r^l\int_r^{R_i}\frac{
!    \rho_{ij;lm}^{\rm MT}(r')}{{r'}^{l-1}}dr'\right)+\frac{1}{Y_{00}}
!    \frac{z_{ij}}{r}\delta_{l,0} $$
!   where the last term is the monopole arising from the point charge. All that
!   remains is to add the homogenous solution of Poisson's equation,
!   $$ V_{ij}^{\rm H}({\bf r})=\sum_{lm}V_{ij;lm}^{\rm H}\left(\frac{r}
!    {R_i}\right)^lY_{lm}(\hat{\bf r}), $$
!   to the muffin-tin potential so that it is continuous at the muffin-tin
!   boundary. Therefore the coefficients, $\rho_{ij;lm}^{\rm H}$, are given by
!   $$ V_{ij;lm}^{\rm H}=4\pi i^l\sum_{\bf G}j_{l}(Gr)V^{\rm P}({\bf G})
!    \exp(i{\bf G}\cdot{\bf r}_{ij})Y_{lm}^*(\hat{\bf G})-V_{ij;lm}^{\rm MT}
!    (R_i). $$
!   Finally note that the ${\bf G}$-vectors passed to the routine can represent
!   vectors with a non-zero offset, ${\bf G}+{\bf p}$ say, which is required for
!   calculating Coulomb matrix elements.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
integer, intent(in) :: ngdg(3),igf(*),ngp
real(8), intent(in) :: gpc(ngp),gclgp(ngp)
integer, intent(in) :: ld2
real(8), intent(in) :: jlgprmt(0:lnpsd,ld2,nspecies)
complex(8), intent(in) :: ylmgp(lmmaxo,ngp),sfacgp(ld2,natmtot)
complex(8), intent(in) :: zrhoir(*)
integer, intent(in) :: ld3
complex(8), intent(inout) :: zvclmt(ld3,natmtot)
complex(8), intent(out) :: zvclir(*)
! local variables
integer is,ia,ias,iro
integer l,lm,lma,lmb
integer ig,jg,i,i0,i1
real(8) t1,t2,t3
complex(8) z1,z2
! automatic arrays
real(8) rmtl(0:lmaxo+3,nspecies)
complex(8) qlm(lmmaxo,natmtot)
complex(8) zl(0:lmaxo),zlm(lmmaxo)
complex(8) zhmt(ld3)
! external functions
real(8), external :: factn2
! compute (R_mt)^l
do is=1,nspecies
  rmtl(0,is)=1.d0
  do l=1,lmaxo+3
    rmtl(l,is)=rmtl(l-1,is)*rmt(is)
  end do
end do
! compute the multipole moments from the muffin-tin potentials
t1=1.d0/fourpi
do ias=1,natmtot
  is=idxis(ias)
  i=np(is)-lmmaxo
  do l=0,lmaxo
    t2=t1*dble(2*l+1)*rmtl(l+1,is)
    lma=l**2+1; lmb=lma+2*l
    qlm(lma:lmb,ias)=t2*zvclmt(i+lma:i+lmb,ias)
  end do
end do
! Fourier transform density to G-space and store in zvclir
call zcopy(ngdg(1)*ngdg(2)*ngdg(3),zrhoir,1,zvclir,1)
call zfftifc(3,ngdg,-1,zvclir)
! subtract the multipole moments of the interstitial charge density
do is=1,nspecies
  do l=0,lmaxo
    zl(l)=fourpi*zil(mod(l,4))*rmtl(l+2,is)
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zlm(:)=0.d0
    do ig=1,ngp
      jg=igf(ig)
      if (gpc(ig).gt.epslat) then
        z1=zvclir(jg)*sfacgp(ig,ias)/gpc(ig)
        zlm(1)=zlm(1)+jlgprmt(1,ig,is)*z1*zl(0)*y00
        do l=1,lmaxo
          lma=l**2+1; lmb=lma+2*l
          z2=jlgprmt(l+1,ig,is)*z1*zl(l)
          zlm(lma:lmb)=zlm(lma:lmb)+z2*conjg(ylmgp(lma:lmb,ig))
        end do
      else
        t1=(fourpi/3.d0)*rmtl(3,is)*y00
        zlm(1)=zlm(1)+t1*zvclir(jg)
      end if
    end do
    qlm(:,ias)=qlm(:,ias)-zlm(:)
  end do
end do
! find the smooth pseudocharge within the muffin-tin whose multipoles are the
! difference between the real muffin-tin and interstitial multipoles
t1=(fourpi/omega)*factn2(2*lnpsd+1)
do ias=1,natmtot
  is=idxis(ias)
  do l=0,lmaxo
    t2=t1/(factn2(2*l+1)*rmtl(l,is))
    z1=t2*zilc(mod(l,4))
    lma=l**2+1; lmb=lma+2*l
    zlm(lma:lmb)=z1*qlm(lma:lmb,ias)
  end do
! add the pseudocharge and real interstitial densities in G-space
  do ig=1,ngp
    jg=igf(ig)
    if (gpc(ig).gt.epslat) then
      t2=gpc(ig)*rmt(is)
      t3=1.d0/t2**lnpsd
      z1=t3*zlm(1)*y00
      do l=1,lmaxo
        lma=l**2+1; lmb=lma+2*l
        t3=t3*t2
        z1=z1+t3*sum(zlm(lma:lmb)*ylmgp(lma:lmb,ig))
      end do
      z2=jlgprmt(lnpsd,ig,is)*conjg(sfacgp(ig,ias))
      zvclir(jg)=zvclir(jg)+z1*z2
    else
      t2=y00/factn2(2*lnpsd+1)
      zvclir(jg)=zvclir(jg)+t2*zlm(1)
    end if
  end do
end do
! solve Poisson's equation in G+p-space for the pseudocharge
do ig=1,ngp
  jg=igf(ig)
  zvclir(jg)=gclgp(ig)*zvclir(jg)
end do
! match potentials at muffin-tin boundary by adding homogeneous solution
do ias=1,natmtot
  is=idxis(ias)
  iro=nri(is)+1
! find the spherical harmonic expansion of the interstitial potential at the
! muffin-tin radius
  zlm(:)=0.d0
  do ig=1,ngp
    z1=fourpi*zvclir(igf(ig))*sfacgp(ig,ias)
    zlm(1)=zlm(1)+jlgprmt(0,ig,is)*z1*y00
    do l=1,lmaxo
      lma=l**2+1; lmb=lma+2*l
      z2=jlgprmt(l,ig,is)*z1*zil(mod(l,4))
      zlm(lma:lmb)=zlm(lma:lmb)+z2*conjg(ylmgp(lma:lmb,ig))
    end do
  end do
! calculate the homogenous solution
  i=np(is)-lmmaxo
  do l=0,lmaxi
    t1=1.d0/rmtl(l,is)
    do lm=l**2+1,(l+1)**2
      z1=t1*(zlm(lm)-zvclmt(i+lm,ias))
      i1=lmmaxi*(nri(is)-1)+lm
      zhmt(lm:i1:lmmaxi)=z1*rl(1:nri(is),l,is)
      i0=i1+lmmaxi
      i1=lmmaxo*(nr(is)-iro)+i0
      zhmt(i0:i1:lmmaxo)=z1*rl(iro:nr(is),l,is)
    end do
  end do
  do l=lmaxi+1,lmaxo
    t1=1.d0/rmtl(l,is)
    do lm=l**2+1,(l+1)**2
      z1=t1*(zlm(lm)-zvclmt(i+lm,ias))
      i0=lmmaxi*nri(is)+lm
      i1=lmmaxo*(nr(is)-iro)+i0
      zhmt(i0:i1:lmmaxo)=z1*rl(iro:nr(is),l,is)
    end do
  end do
  zvclmt(1:np(is),ias)=zvclmt(1:np(is),ias)+zhmt(1:np(is))
! store the nuclear potential without the self-term for the phonon dynamical
! matrix calculation
  if (tphdyn) then
    if (ias.eq.iasph) zvnmt(1:np(is))=zhmt(1:np(is))
  end if
end do
! Fourier transform interstitial potential to real-space
call zfftifc(3,ngdg,1,zvclir)
end subroutine
!EOC

