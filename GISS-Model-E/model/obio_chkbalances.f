#include "rundeck_opts.h"

      subroutine obio_chkbalances(vrbos,nstep,i,j)
      !check balances globally

      USE obio_incom
      USE obio_com,  only: rhs_obio,chng_by

      implicit none

      integer nstep,i,j,nt
      real*8 balnc(15)
      real*8 tolrnc
      logical vrbos

      tolrnc =  1.d-10
!-----------------------------------------------------------------------
      balnc(:)=0.d0
!nitrates
      balnc(1) = rhs_obio(i,j,1,5)+rhs_obio(i,j,1,6)
     .                    + rhs_obio(i,j,1,7)+rhs_obio(i,j,1,8)
     .                    + rhs_obio(i,j,1,11)
     .                    +(rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .                     +rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13)
     .                     +rhs_obio(i,j,7,12))
     .                     *bn
     .                    + rhs_obio(i,j,2,5)+rhs_obio(i,j,2,6)
     .                    + rhs_obio(i,j,2,7)+rhs_obio(i,j,2,8)  !'Nchange by growth='
      balnc(2) = rhs_obio(i,j,10,16) /cnratio                    !'Nchange by detrital sinking='
      balnc(3) = rhs_obio(i,j,1,9)+rhs_obio(i,j,9,14)*bn !'Nchange by excretion of DON by zoo='
      balnc(4) = rhs_obio(i,j,1,15)
     .                   +(rhs_obio(i,j,5,14)+rhs_obio(i,j,5,15)
     .                    +rhs_obio(i,j,6,14)+rhs_obio(i,j,6,15)
     .                    +rhs_obio(i,j,7,14)+rhs_obio(i,j,7,15)
     .                    +rhs_obio(i,j,8,14)+rhs_obio(i,j,8,15))
     .                     *bn        !'Nchange by excretion of DON by phyto='
      balnc(5) = rhs_obio(i,j, 1,14)+rhs_obio(i,j,10,14)/cnratio !'Nchange by detrital breakdown='
      balnc(6) = rhs_obio(i,j, 1,10)+rhs_obio(i,j,10,10)/cnratio !'Nchange by detrital breakdown='
      balnc(7) = rhs_obio(i,j,10,5)/cnratio
     .                 +(rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .                 + rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8))
     .                 *bn                !'Nchange by phyto death='
      balnc(8) = rhs_obio(i,j,10,7)/cnratio
     .             +rhs_obio(i,j,9,11)*bn !'Nchange by zoo death='

      if (vrbos) then
      do nt=1,8
        if (abs(balnc(nt)).gt.tolrnc) then
           write(*,'(a,4i6,e18.8)') 
     .       'RHS imbalance, nitrates:',nstep,i,j,nt,balnc(nt)
          if (nt.eq.1) write(*,'(a,3i6,3e18.8)')'detail RHS imb, n:',
     .      nstep,i,j,balnc(nt),rhs_obio(i,j,1,11),rhs_obio(i,j,7,12)*bn
        endif
      enddo
      endif

!-----------------------------------------------------------------------
      balnc(:)=0.d0
!ammonium
      balnc(1) = rhs_obio(i,j,1,5)+rhs_obio(i,j,1,6)
     .                    + rhs_obio(i,j,1,7)+rhs_obio(i,j,1,8)
     .                    + rhs_obio(i,j,1,11)
     .                    +(rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .                     +rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13)
     .                     +rhs_obio(i,j,7,12))
     .                     *bn
     .                    + rhs_obio(i,j,2,5)+rhs_obio(i,j,2,6)
     .                    + rhs_obio(i,j,2,7)+rhs_obio(i,j,2,8)    !'Achange by growth='
      balnc(2) = rhs_obio(i,j,2,9)/bn+rhs_obio(i,j,9,10)           !'Achange by regener from zoopl='
      balnc(3) = rhs_obio(i,j,2,10)/bn+rhs_obio(i,j,10,9)/(cnratio*bn)  !'Achange by regener from detrt='
      if (vrbos) then
      do nt=1,3
        if (abs(balnc(nt)).gt.tolrnc) write(*,'(a,4i6,e18.8)') 
     .       'RHS imbalance, ammonium:',nstep,i,j,nt,balnc(nt)
      enddo
      endif


!-----------------------------------------------------------------------
      balnc(:)=0.d0
!silicate 
      balnc(1) = rhs_obio(i,j,3,5) +rhs_obio(i,j,5,13)*bs ! 'Schange by growth=', 
      balnc(2) = rhs_obio(i,j,11,16)/csratio              ! 'Schange by detrital sinking=',
      balnc(3) = rhs_obio(i,j, 3,11)+rhs_obio(i,j,11,11)  ! 'Schange by reminr=', 
      balnc(4) = rhs_obio(i,j, 3,15)   
     .                   +(rhs_obio(i,j,5,14)+rhs_obio(i,j,5,15))
     .                     *bs                            ! 'Schange by excretion&respir=', 
      balnc(5) = rhs_obio(i,j,11, 5)/bs+rhs_obio(i,j,5,5) ! 'Schange by phyt death=', 
      balnc(6) = rhs_obio(i,j,11, 9)/bs+rhs_obio(i,j,5,9) ! 'Schange by zoo  death=', 
      if (vrbos) then
      do nt=1,6
        if (abs(balnc(nt)).gt.tolrnc) write(*,'(a,4i6,e18.8)') 
     .       'RHS imbalance, silicate:',nstep,i,j,nt,balnc(nt)
      enddo
      endif

!-----------------------------------------------------------------------
      balnc(:)=0.d0
!iron conservation
      balnc(1) = rhs_obio(i,j,4,5)+rhs_obio(i,j,4,6)  
     .                    +rhs_obio(i,j,4,7)+rhs_obio(i,j,4,8)
     .                    +(rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .                     +rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13)
     .                     +rhs_obio(i,j,7,12))
     .                     *bf        ! 'Ichange by growth=', 
      balnc(2) = rhs_obio(i,j,12,16)/cfratio                 ! 'Ichange by detrital sinking=',
      balnc(3) = rhs_obio(i,j,4,9)/bf+rhs_obio(i,j,9,10)  ! 'Ichange by grazing eff=', 
      balnc(4) = rhs_obio(i,j,4,11)+rhs_obio(i,j,12,7)    ! 'Ichange by regeneration=', 
      balnc(5) = rhs_obio(i,j,12,9)/bf+rhs_obio(i,j,9,11) ! 'Ichange by regeneration=', 
      balnc(6) = rhs_obio(i,j,4,10)/bf+rhs_obio(i,j,9,14) ! 'Ichange by excrtn zoo=',
      balnc(7) = rhs_obio(i,j,4,12)+rhs_obio(i,j,12,12)   ! 'Ichange by reminerlzn=',
      balnc(8) = rhs_obio(i,j,12,6)/bf+rhs_obio(i,j,9,12) ! 'Ichange by zoo death =',
      balnc(9) = rhs_obio(i,j, 4,13)+rhs_obio(i,j,12,4)   ! 'Ichange by scavenging',
      balnc(10) = rhs_obio(i,j, 4,15)
     .                   +(rhs_obio(i,j,5,14)+rhs_obio(i,j,5,15)
     .                    +rhs_obio(i,j,6,14)+rhs_obio(i,j,6,15)
     .                    +rhs_obio(i,j,7,14)+rhs_obio(i,j,7,15)
     .                    +rhs_obio(i,j,8,14)+rhs_obio(i,j,8,15))
     .                     *bf                            ! 'Ichange by excrtnresp=',
      balnc(11) = rhs_obio(i,j,12, 5)/bf
     .                 +rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .                 + rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8) !'Ichange by phyt death=',
      balnc(12) = rhs_obio(i,j, 4, 4)+rhs_obio(i,j,4,14)     !'Ichange by Ideposition and bottom sink'

      if (vrbos) then
      do nt=1,12
        if (abs(balnc(nt)).gt.tolrnc) write(*,'(a,4i6,e18.8)') 
     .       'RHS imbalance, iron:',nstep,i,j,nt,balnc(nt)
        if (nt.eq.12) write(*,'(a,4i6,3e18.8)') 
     .       'RHS imbalance, iron, detail:',nstep,i,j,nt,
     .        rhs_obio(i,j, 4, 4),rhs_obio(i,j,4,14),balnc(nt)
      enddo
      endif

!-----------------------------------------------------------------------
      balnc(:)=0.d0
!carbon 
      balnc(1) = (rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .               +  rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13)
     .                       +rhs_obio(i,j,7,12))
     .               *mgchltouMC
     .               + rhs_obio(i,j,14,6)               !'chng by growth  ='
      balnc(2) = ((rhs_obio(i,j,5,9)+rhs_obio(i,j,6,9)
     .            +rhs_obio(i,j,7,9)+rhs_obio(i,j,8,9))
     .            +  rhs_obio(i,j,9,9)) * mgchltouMC    !'chng by grazing ='
      balnc(3) = (rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .               +  rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8))
     .               *mgchltouMC
     .               +  rhs_obio(i,j,10,5)/uMtomgm3     !'chng by phyt dea='
      balnc(4) = (rhs_obio(i,j,5,14)+rhs_obio(i,j,6,14)
     .               +  rhs_obio(i,j,7,14)+rhs_obio(i,j,8,14))
     .               *mgchltouMC
     .               +  rhs_obio(i,j,13,5)              !'chng by DOC prod='
      balnc(5) = (rhs_obio(i,j,5,15)+rhs_obio(i,j,6,15)
     .               +  rhs_obio(i,j,7,15)+rhs_obio(i,j,8,15))
     .               *mgchltouMC
     .               +  rhs_obio(i,j,14,5)              !'chng by DIC prod='
      balnc(6) = rhs_obio(i,j,9,10)*mgchltouMC 
     .                            + rhs_obio(i,j,13,9)  !'chng by zoo graz='
      balnc(7) = rhs_obio(i,j,9,11)*mgchltouMC
     .                            + rhs_obio(i,j,10,7)/uMtomgm3  !'chng by zoo dea ='
      balnc(8) = rhs_obio(i,j,9,12)*mgchltouMC
     .                            + rhs_obio(i,j,10,8)/uMtomgm3  !'chng by zoo dea ='
      balnc(9) = rhs_obio(i,j,9,14)*mgchltouMC
     .                            + rhs_obio(i,j,13,12)     !'chng by excz    ='
      balnc(10) = rhs_obio(i,j,9,15)*mgchltouMC
     .                            + rhs_obio(i,j,14,15)    !'chng by zoo prod='
      balnc(11) = rhs_obio(i,j,10,14)/uMtomgm3
     .                            + rhs_obio(i,j,13,13)    !'chng by docdet  ='
      balnc(12) = rhs_obio(i,j,10,9)/uMtomgm3
     .                            + rhs_obio(i,j,13,10)    !'chng by regen   ='
      balnc(13) = rhs_obio(i,j,10,10)/uMtomgm3
     .                            + rhs_obio(i,j,14,10)    !'chng by reminer ='
      balnc(14) = rhs_obio(i,j,13,14)
     .                            + rhs_obio(i,j,14,14)    !'chng by docbac  ='
 
      if (vrbos) then
      do nt=1,14
        if (abs(balnc(nt)).gt.tolrnc) write(*,'(a,4i6,e18.8)') 
     .       'RHS imbalance, carbon:',nstep,i,j,nt,balnc(nt)
      enddo
      endif

!-----------------------------------------------------------------------
#ifdef TRACERS_Alkalinity
      balnc(:)=0.d0
!alkalinity
      balnc(1) =  rhs_obio(i,j,15,1) + rhs_obio(i,j,15,5)      

      do nt=1,1
        if (abs(balnc(nt)).gt.tolrnc) write(*,'(a,4i6,3e18.8)') 
     .       'RHS imbalance, alkalinity:',nstep,i,j,nt,
     .       rhs_obio(i,j,15,1),rhs_obio(i,j,15,5),balnc(nt)
      enddo
#endif

!-----------------------------------------------------------------------
      if (vrbos) then
      write(*,*) '_____________________________________________'
      write(*,*) '---------------------------------------------'
      write(*,'(a,3i10)') 'N Conserv diagn, at ',nstep,i,j
      write(*,'(a,3i10)') 'all units in mili-mol,N/m2/s'   
!note here: N/C-detritus is in micro,grams/lt so use cnratio to get to uM,N
!all other terms which are in uM,C divide by (cnratio*12)
      print*, 'Nchange by excretion of DON by zoo=',
     .               rhs_obio(i,j,1,9)+rhs_obio(i,j,9,14)*bn
      print*, 'Nchange by excretion of DON by phyto=', 
     .                     rhs_obio(i,j,1,15)
     .                   +(rhs_obio(i,j,5,14)+rhs_obio(i,j,5,15)
     .                    +rhs_obio(i,j,6,14)+rhs_obio(i,j,6,15)
     .                    +rhs_obio(i,j,7,14)+rhs_obio(i,j,7,15)
     .                    +rhs_obio(i,j,8,14)+rhs_obio(i,j,8,15))
     .                     *bn
      print*, 'Nchange by detrital breakdown=',
     .                 rhs_obio(i,j, 1,14)+rhs_obio(i,j,10,14)/cnratio
      print*, 'Nchange by remineralization=',
     .                 rhs_obio(i,j, 1,10)+rhs_obio(i,j,10,10)/cnratio
      print*, 'Nchange by detrital sinking=',
     .                 rhs_obio(i,j,10,16) /cnratio
      print*, 'Nchange by nfixatn=',
     .                 rhs_obio(i,j,1,11)
     .                +rhs_obio(i,j,7,12)*bn
      print*, 'Nchange by nfixatn2=',
     .                 rhs_obio(i,j,1,11),rhs_obio(i,j,7,12)*bn
      print*, 'Nchange by nutrient uptake=',
     .                     (rhs_obio(i,j,1,5)+rhs_obio(i,j,1,6)
     .                    + rhs_obio(i,j,1,7)+rhs_obio(i,j,1,8))/bn
     .               + rhs_obio(i,j,14,6)/mgchltouMC
     .                    +(rhs_obio(i,j,2,5)+rhs_obio(i,j,2,6)
     .                    + rhs_obio(i,j,2,7)+rhs_obio(i,j,2,8))/bn
      print*, 'Nchange by growth=', 
     .                      rhs_obio(i,j,1,5)+rhs_obio(i,j,1,6)
     .                    + rhs_obio(i,j,1,7)+rhs_obio(i,j,1,8)
     .                    +(rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .                     +rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13))
     .                     *bn
     .                    + rhs_obio(i,j,2,5)+rhs_obio(i,j,2,6)
     .                    + rhs_obio(i,j,2,7)+rhs_obio(i,j,2,8)
      write(*,'(a,12e12.4)')'detail Nchange by growth',
     . rhs_obio(i,j,1,5:8),rhs_obio(i,j,2,5:8),rhs_obio(i,j,5:8,13)*bn
      print*, 'Nchange by phyto death=', 
     .                   rhs_obio(i,j,10,5)/cnratio
     .                 +(rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .                 + rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8))
     .                 *bn
      print*, 'Nchange by zoo death=', 
     .              rhs_obio(i,j,10,7)/cnratio
     .             +rhs_obio(i,j,9,11)*bn

 
!ammonium conservation
      print*, '   '
      print*, 'ammo conservation'
      print*, 'Achange by growth=', 
     .                      rhs_obio(i,j,1,5)+rhs_obio(i,j,1,6)
     .                    + rhs_obio(i,j,1,7)+rhs_obio(i,j,1,8)
     .                    +(rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .                     +rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13))
     .                     *bn
     .                    + rhs_obio(i,j,2,5)+rhs_obio(i,j,2,6)
     .                    + rhs_obio(i,j,2,7)+rhs_obio(i,j,2,8)
      print*, 'Achange by regener from zoopl=',
     .              rhs_obio(i,j,2,9)/bn+rhs_obio(i,j,9,10)
      print*, 'Achange by regener from detrt=',
     .       rhs_obio(i,j,2,10)/bn+rhs_obio(i,j,10,9)/(cnratio*bn)  


!silicate conservation
      print*, '   '
      print*, 'sili conservation'
      print*, 'Schange by growth=', 
     .                     rhs_obio(i,j,3,5)
     .                   +rhs_obio(i,j,5,13)*bs
      print*, 'Schange by reminr=', 
     .                     rhs_obio(i,j, 3,11)+rhs_obio(i,j,11,11) 
      print*, 'Schange by excretion&respir=', 
     .                     rhs_obio(i,j, 3,15)   
     .                   +(rhs_obio(i,j,5,14)+rhs_obio(i,j,5,15))
     .                     *bs
      print*, 'Schange by phyt death=', 
     .                   rhs_obio(i,j,11, 5)/bs+rhs_obio(i,j,5,5)
      print*, 'Schange by zoo  death=', 
     .                rhs_obio(i,j,11, 9)/bs+rhs_obio(i,j,5,9)   
      print*, 'Schange by detrital sinking=',
     .                  rhs_obio(i,j,11,16)/csratio

!iron conservation
      print*, '   '
      print*, 'iron conservation'
      print*, 'Ichange by deposition=', 
     .                     rhs_obio(i,j, 4, 4)   
      print*, 'Ichange by growth=', 
     .                     rhs_obio(i,j,4,5)+rhs_obio(i,j,4,6)  
     .                    +rhs_obio(i,j,4,7)+rhs_obio(i,j,4,8)
     .                    +(rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .                     +rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13))
     .                     *bf
      print*, 'Ichange by grazing eff=', 
     .          rhs_obio(i,j,4,9)/bf+rhs_obio(i,j,9,10)
      print*, 'Ichange by regeneration=', 
     .          rhs_obio(i,j,4,11)+rhs_obio(i,j,12,7)
      print*, 'Ichange by regeneration=', 
     .          rhs_obio(i,j,12,9)/bf+rhs_obio(i,j,9,11)
      print*, 'Ichange by excrtn zoo=',
     .          rhs_obio(i,j,4,10)/bf+rhs_obio(i,j,9,14)
      print*, 'Ichange by reminerlzn=',
     .          rhs_obio(i,j,4,12)+rhs_obio(i,j,12,12)
      print*, 'Ichange by zoo death =',
     .          rhs_obio(i,j,12,6)/bf+rhs_obio(i,j,9,12)
      print*, 'Ichange by scavenging',
     .          rhs_obio(i,j, 4,13)+rhs_obio(i,j,12,4)
      print*, 'Ichange by excrtnresp=',
     .          rhs_obio(i,j, 4,15)
     .                   +(rhs_obio(i,j,5,14)+rhs_obio(i,j,5,15)
     .                    +rhs_obio(i,j,6,14)+rhs_obio(i,j,6,15)
     .                    +rhs_obio(i,j,7,14)+rhs_obio(i,j,7,15)
     .                    +rhs_obio(i,j,8,14)+rhs_obio(i,j,8,15))
     .                     *bf
      print*, 'Ichange by phyt death=',
     .        rhs_obio(i,j,12, 5)/bf
     .                 +rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .                 + rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8)
      print*, 'Ichange by detrital sinking=',
     .                  rhs_obio(i,j,12,16)/cfratio
      print*, 'Ichange by bottom sink=',
     .                  rhs_obio(i,j,4,14)

!carbon conservation
      !convert all to mili-mol,C/m2
!     rhs_obio(i,j,5,:) = rhs_obio(i,j,5,:) * mgchltouMC
!     rhs_obio(i,j,6,:) = rhs_obio(i,j,6,:) * mgchltouMC
!     rhs_obio(i,j,7,:) = rhs_obio(i,j,7,:) * mgchltouMC
!     rhs_obio(i,j,8,:) = rhs_obio(i,j,8,:) * mgchltouMC
!     rhs_obio(i,j,9,:) = rhs_obio(i,j,9,:) * mgchltouMC
!     rhs_obio(i,j,10,:) = rhs_obio(i,j,10,:) / uMtomgm3

      chng_by(i,j,1) =((rhs_obio(i,j,5,9)+rhs_obio(i,j,6,9)
     .                  +rhs_obio(i,j,7,9)+rhs_obio(i,j,8,9))
     .               +  rhs_obio(i,j,9,9)) * mgchltouMC
      chng_by(i,j,2) = (rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .               +  rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8))
     .               *mgchltouMC
     .               +  rhs_obio(i,j,10,5)/uMtomgm3
      chng_by(i,j,3) = (rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .               +  rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13)
     .                       +rhs_obio(i,j,7,12))
     .               *mgchltouMC
     .               + rhs_obio(i,j,14,6)
      chng_by(i,j,4) = (rhs_obio(i,j,5,14)+rhs_obio(i,j,6,14)
     .               +  rhs_obio(i,j,7,14)+rhs_obio(i,j,8,14))
     .               *mgchltouMC
     .               +  rhs_obio(i,j,13,5)
      chng_by(i,j,5) = (rhs_obio(i,j,5,15)+rhs_obio(i,j,6,15)
     .               +  rhs_obio(i,j,7,15)+rhs_obio(i,j,8,15))
     .               *mgchltouMC
     .               +  rhs_obio(i,j,14,5)
      chng_by(i,j,6) =  rhs_obio(i,j,9,10)*mgchltouMC 
     .                            + rhs_obio(i,j,13,9)
      chng_by(i,j,7) = rhs_obio(i,j,9,11)*mgchltouMC
     .                            + rhs_obio(i,j,10,7)/uMtomgm3 
      chng_by(i,j,8) = rhs_obio(i,j,9,12)*mgchltouMC
     .                            + rhs_obio(i,j,10,8)/uMtomgm3
      chng_by(i,j,9) = rhs_obio(i,j,9,14)*mgchltouMC
     .                            + rhs_obio(i,j,13,12) 
      chng_by(i,j,10)= rhs_obio(i,j,9,15)*mgchltouMC
     .                            + rhs_obio(i,j,14,15)
      chng_by(i,j,11)= rhs_obio(i,j,10,14)/uMtomgm3
     .                            + rhs_obio(i,j,13,13)
      chng_by(i,j,12)= rhs_obio(i,j,10,9)/uMtomgm3
     .                            + rhs_obio(i,j,13,10)
      chng_by(i,j,13)= rhs_obio(i,j,10,10)/uMtomgm3
     .                            + rhs_obio(i,j,14,10)
      chng_by(i,j,14)= rhs_obio(i,j,13,14)
     .                            + rhs_obio(i,j,14,14)
 

      print*, '   '
      write(*,'(a,3i10)') 'Conserv diagn, at ',nstep,i,j
      write(*,'(a,3i10)') 'all units in mili-mol,C/m2/s'   
      write(*,*)'chng by grazing =', chng_by(i,j,1)
      write(*,*)'chng by phyt dea=', chng_by(i,j,2)
      write(*,*)'chng by growth  =', chng_by(i,j,3)
      write(*,*)'chng by DOC prod=', chng_by(i,j,4)
      write(*,*)'chng by DIC prod=', chng_by(i,j,5)
      write(*,*)'chng by zoo graz=', chng_by(i,j,6)
      write(*,*)'chng by zoo dea =', chng_by(i,j,7)
      write(*,*)'chng by zoo dea =', chng_by(i,j,8)
      write(*,*)'chng by excz    =', chng_by(i,j,9)
      write(*,*)'chng by zoo prod=', chng_by(i,j,10)
      write(*,*)'chng by docdet  =', chng_by(i,j,11)
      write(*,*)'chng by regen   =', chng_by(i,j,12)
      write(*,*)'chng by reminer =', chng_by(i,j,13)
      write(*,*)'chng by docbac  =', chng_by(i,j,14)
      write(*,*) '____________________________________________'
      write(*,*) '--------------------------------------------'

!phytoplankton terms
      write(*,'(a,2i7,6e18.8)')'Diat tendencies:',i,j,
     .  rhs_obio(i,j,5,5),rhs_obio(i,j,5,9),rhs_obio(i,j,5,13),
     .  rhs_obio(i,j,5,14),rhs_obio(i,j,5,15),rhs_obio(i,j,5,16)
      write(*,'(a,2i7,6e18.8)')'Chlo tendencies:',i,j,
     .  rhs_obio(i,j,6,6),rhs_obio(i,j,6,9),rhs_obio(i,j,6,13),
     .  rhs_obio(i,j,6,14),rhs_obio(i,j,6,15),rhs_obio(i,j,6,16)
      write(*,'(a,2i7,7e18.8)')'Cyan tendencies:',i,j,
     .  rhs_obio(i,j,7,7),rhs_obio(i,j,7,9),rhs_obio(i,j,7,13),
     .  rhs_obio(i,j,7,14),rhs_obio(i,j,7,15),rhs_obio(i,j,7,16),
     .  rhs_obio(i,j,7,12)
      write(*,'(a,2i7,6e18.8)')'Cocc tendencies:',i,j,
     .  rhs_obio(i,j,8,8),rhs_obio(i,j,8,9),rhs_obio(i,j,8,13),
     .  rhs_obio(i,j,8,14),rhs_obio(i,j,8,15),rhs_obio(i,j,8,16)

      endif    ! vrbos
      end subroutine obio_chkbalances
