
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init4
use modmain
use modphonon
use modpw
use modvars
implicit none
! local variables
integer ik,jspn,n,i
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
real(8) vl(3),vc(3)
! external functions
real(8), external :: gaunt

!---------------------------!
!     H+k-vector arrays     !
!---------------------------!
if (any(task.eq.[135,170,171,172,173])) then
  if (task.eq.135) hkmax=0.5d0*gmaxvr-epslat
  call findngkmax(nkpt,vkc,nspnfv,vqcss,ngvec,vgc,hkmax,nhkmax)
! allocate the H+k-vector arrays
  if (allocated(nhk)) deallocate(nhk)
  allocate(nhk(nspnfv,nkpt))
  if (allocated(ihkig)) deallocate(ihkig)
  allocate(ihkig(nhkmax,nspnfv,nkpt))
  if (allocated(vhkl)) deallocate(vhkl)
  allocate(vhkl(3,nhkmax,nspnfv,nkpt))
  if (allocated(vhkc)) deallocate(vhkc)
  allocate(vhkc(3,nhkmax,nspnfv,nkpt))
  if (allocated(hkc)) deallocate(hkc)
  allocate(hkc(nhkmax,nspnfv,nkpt))
  if (allocated(sfachk)) deallocate(sfachk)
  allocate(sfachk(nhkmax,natmtot,nspnfv,nkpt))
! initialise H+k-vectors arrays
  do ik=1,nkpt
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
! generate the H+k-vectors
      call gengkvec(ngvec,ivg,vgc,vl,vc,hkmax,nhkmax,nhk(jspn,ik), &
       ihkig(:,jspn,ik),vhkl(:,:,jspn,ik),vhkc(:,:,jspn,ik),hkc(:,jspn,ik))
! generate structure factors for H+k-vectors
      call gensfacgp(nhk(jspn,ik),vhkc(:,:,jspn,ik),nhkmax,sfachk(:,:,jspn,ik))
    end do
  end do
! write to VARIABLES.OUT
  if (wrtvars) then
    call writevars('hkmax',rv=hkmax)
    call writevars('nhk',nv=nspnfv*nkpt,iva=nhk)
    do ik=1,nkpt
      do jspn=1,nspnfv
        call writevars('ihkig',n1=jspn,n2=ik,nv=nhk(jspn,ik), &
         iva=ihkig(:,jspn,ik))
      end do
    end do
  end if
end if

!-----------------------------!
!     G+k+q-vector arrays     !
!-----------------------------!
if (task.eq.205) then
  if (allocated(vkql)) deallocate(vkql)
  allocate(vkql(3,nkptnr))
  if (allocated(vkqc)) deallocate(vkqc)
  allocate(vkqc(3,nkptnr))
  if (allocated(ngkq)) deallocate(ngkq)
  allocate(ngkq(nspnfv,nkptnr))
  if (allocated(igkqig)) deallocate(igkqig)
  allocate(igkqig(ngkmax,nspnfv,nkptnr))
  if (allocated(vgkql)) deallocate(vgkql)
  allocate(vgkql(3,ngkmax,nspnfv,nkptnr))
  if (allocated(vgkqc)) deallocate(vgkqc)
  allocate(vgkqc(3,ngkmax,nspnfv,nkptnr))
  if (allocated(gkqc)) deallocate(gkqc)
  allocate(gkqc(ngkmax,nspnfv,nkptnr))
  if (allocated(sfacgkq)) deallocate(sfacgkq)
  allocate(sfacgkq(ngkmax,natmtot,nspnfv,nkptnr))
end if

!---------------------------!
!     G+q-vector arrays     !
!---------------------------!
if (task.eq.205) then
  if (allocated(vgqc)) deallocate(vgqc)
  allocate(vgqc(3,ngtot))
  if (allocated(gqc)) deallocate(gqc)
  allocate(gqc(ngtot))
  if (allocated(gclgq)) deallocate(gclgq)
  allocate(gclgq(ngvec))
  if (allocated(jlgqrmt)) deallocate(jlgqrmt)
  allocate(jlgqrmt(0:lnpsd,ngvec,nspecies))
  if (allocated(ylmgq)) deallocate(ylmgq)
  allocate(ylmgq(lmmaxo,ngvec))
  if (allocated(sfacgq)) deallocate(sfacgq)
  allocate(sfacgq(ngvec,natmtot))
  if (allocated(ffacgq)) deallocate(ffacgq)
  allocate(ffacgq(ngtot,nspecies))
  if (allocated(dcfunig)) deallocate(dcfunig)
  allocate(dcfunig(ngtot))
  if (allocated(dcfunir)) deallocate(dcfunir)
  allocate(dcfunir(ngtot))
end if

!-----------------------------------------------------------------!
!     phonon density functional perturbation theory variables     !
!-----------------------------------------------------------------!
if (task.eq.205) then
  if (allocated(drhomt)) deallocate(drhomt)
  allocate(drhomt(npmtmax,natmtot))
  if (allocated(drhoir)) deallocate(drhoir)
  allocate(drhoir(ngtot))
  if (allocated(dmagmt)) deallocate(dmagmt)
  if (allocated(dmagir)) deallocate(dmagir)
  if (spinpol) then
    allocate(dmagmt(npmtmax,natmtot,ndmag))
    allocate(dmagir(ngtot,ndmag))
  end if
  if (allocated(dvclmt)) deallocate(dvclmt)
  allocate(dvclmt(npmtmax,natmtot))
  if (allocated(dvclir)) deallocate(dvclir)
  allocate(dvclir(ngtot))
  if (allocated(zvnmt)) deallocate(zvnmt)
  allocate(zvnmt(npmtmax))
  if (allocated(gvsmt)) deallocate(gvsmt)
  allocate(gvsmt(npmtmax))
! combined target array for Kohn-Sham potential and magnetic field derivative
  if (allocated(dvsbs)) deallocate(dvsbs)
  n=npmtmax*natmtot+ngtot
  if (spinpol) n=n+(npcmtmax*natmtot+ngtot)*ndmag
  allocate(dvsbs(n))
! zero the array
  dvsbs(:)=0.d0
! associate pointer arrays with target
  dvsmt(1:npmtmax,1:natmtot)=>dvsbs(1:)
  i=npmtmax*natmtot+1
  dvsir(1:ngtot)=>dvsbs(i:)
  if (spinpol) then
    i=i+ngtot
    dbsmt(1:npcmtmax,1:natmtot,1:ndmag)=>dvsbs(i:)
    i=i+npcmtmax*natmtot*ndmag
    dbsir(1:ngtot,1:ndmag)=>dvsbs(i:)
  end if
  if (allocated(dvsig)) deallocate(dvsig)
  allocate(dvsig(ngtot))
  if (allocated(dsocfr)) deallocate(dsocfr)
  if (spinorb) then
    allocate(dsocfr(nrcmtmax,natmtot))
  end if
  if (allocated(dhaa)) deallocate(dhaa)
  allocate(dhaa(lmmaxo,apwordmax,0:lmaxapw,apwordmax,0:lmaxapw,natmtot))
  if (allocated(dhloa)) deallocate(dhloa)
  allocate(dhloa(lmmaxo,apwordmax,0:lmaxapw,nlomax,natmtot))
  if (allocated(dhlolo)) deallocate(dhlolo)
  allocate(dhlolo(lmmaxo,nlomax,nlomax,natmtot))
! allocate and generate real Gaunt coefficient array
  if (allocated(gntyyy)) deallocate(gntyyy)
  allocate(gntyyy(lmmaxo,lmmaxapw,lmmaxapw))
  do l1=0,lmaxapw
    do m1=-l1,l1
      lm1=l1*(l1+1)+m1+1
      do l3=0,lmaxapw
        do m3=-l3,l3
          lm3=l3*(l3+1)+m3+1
          do l2=0,lmaxo
            do m2=-l2,l2
              lm2=l2*(l2+1)+m2+1
              gntyyy(lm2,lm3,lm1)=gaunt(l1,l2,l3,m1,m2,m3)
            end do
          end do
        end do
      end do
    end do
  end do
  if (allocated(devalfv)) deallocate(devalfv)
  allocate(devalfv(nstfv,nspnfv,nkptnr))
  if (allocated(devalsv)) deallocate(devalsv)
  allocate(devalsv(nstsv,nkptnr))
  if (allocated(doccsv)) deallocate(doccsv)
  allocate(doccsv(nstsv,nkptnr))
end if

end subroutine

