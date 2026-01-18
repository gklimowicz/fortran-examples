#include "rundeck_opts.h"
c Family chemistry calculations: Equilibrium values are production/loss
c from reactions *within* family only:


      SUBROUTINE Oxfam(Lmax,I,J)
!@sum Oxfam Find O,O1D & Ox initial conc assuming equilibrium with O3
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only : LM
      USE TRACER_COM, only : n_CH4, n_Ox, nn_Ox, nn_CH4
      use photolysis, only: rj
      USE TRCHEM_Shindell_COM, only:ss,rr,y,nO2,nM,nH2O,nO,nO1D,nO3,pOx
     &                             ,rrbi,rrtri
      IMPLICIT NONE

C**** Local parameters and variables and arguments:

!@var az,bz,P1 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var Lmax maximum altitude for chemistry
      integer, intent(IN)   :: Lmax,I,J
      integer               :: L
      real*8                :: az, bz, P1

      do L=1,Lmax
c       P1 = [O3]/[Ox] 
c       P2=[O]/[Ox] 
c       P3 = [O1D]/[Ox] 
c       bz = [O1D]/[O3] = P3/P1
c       az = [O]/[O3] = P2/P1
c       P1+P2+P3=1 (total Ox = sum of parts) = P1+P1az+P1bz; so P1 = 1/(1+az+bz)
c
c       for concentration of O(1D):
        bz=ss(rj%O3__O1D_O2,L,I,J)
     &    /(rr(rrbi%O1D_O2__O_O2,L)*y(nO2,L)
     &      +rr(rrbi%O1D_M__O_M,L)*y(nM,L)
     &      +rr(rrbi%O1D_H2O__OH_OH,L)*y(nH2O,L)
     &      +rr(rrbi%O1D_CH4__OH_CH3O2,L)*y(nn_CH4,L)
     &      +rr(rrbi%O1D_CH4__HCHO_H2,L)*y(nn_CH4,L))
        ! here we USED TO tune bz with a pressure criterion
c       for concentration of O:
        if (y(nO2,L) > 0.d0) then
          az=(ss(rj%O3__O1D_O2,L,I,J)
     &        +ss(rj%O3__O_O2,L,I,J))
     &      /(rr(rrtri%O_O2__O3_M,L)*y(nO2,L))
          P1=1.d0/(1.d0+az+bz)
        else
          az=0.d0 ! actually az=infinite but multiplied by P1(=0) below:
          P1=0.d0
        end if
        y(nO,L)=P1*az*y(nn_Ox,L)
        y(nO1D,L)=P1*bz*y(nn_Ox,L)
        y(nO3,L)=y(nn_Ox,L)-y(nO,L)-y(nO1D,L)
        if(y(nO,L) < 0.)  y(nO,L)  =0.d0
        if(y(nO1D,L) < 0.)y(nO1D,L)=0.d0
        if(y(nO3,L) < 1.) y(nO3,L) =1.d0
        if(y(nn_Ox,L) < 1.)y(nn_Ox,L)=1.d0
        pOx(I,J,L)=y(nO3,L)/y(nn_Ox,L)
      end do
c
      return
      END SUBROUTINE Oxfam



      SUBROUTINE NOxfam(Lmax,I,J)
!@sum NOxfam Find NOx family (NO,NO2,NO3,HONO) partitioning assuming
!@+   equilibrium at a given concentration of NOx. Only called during
!@+   daylight, then assume NO3=HONO=1.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:

      USE RESOLUTION, only         : LS1=>LS1_NOMINAL
      USE ATM_COM, only            : LTROPO
      USE TRACER_COM, only         : n_NOx,nn_NOx,nn_Alkenes
      use photolysis, only: rj
      USE TRCHEM_Shindell_COM, only:rr,y,yNO3,nO3,nHO2,nO,nC2O3,nCH3O2,
     & pNO3,ta,nXO2,ss,nNO,nNO2,pNOx,nNO3,nHONO,which_trop,nClO,nOClO,
     & nBrO,rrbi,rrtri

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var b,c,p1,p2 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var maxT LTROPO(I,J) or LS1-1, depending upon what_trop variable
!@+ Or the top layer of chemistry in the unlikely event that is lower.
!@+ Note in that case, loops like L=maxT+1,Lmax will do nothing.
      integer             :: L,maxT
      integer, intent(IN) :: Lmax,I,J
      real*8              :: b,c,p1,p2,d
      
      select case(which_trop)
      case(0); maxT=min(Ltropo(I,J),Lmax)
      case(1); maxT=min(ls1-1,Lmax)
      case default; call stop_model('which_trop problem 8',255)
      end select

      do L=1,Lmax
        ! here we USED TO set NO3 back to zero at dawn.
c       B is for NO->NO2 reactions :
        B=rr(rrbi%O3_NO__NO2_O2,L)*y(nO3,L)
     &    +rr(rrbi%HO2_NO__OH_NO2,L)*y(nHO2,L)
     &    +rr(rrtri%NO_O__NO2_M,L)*y(nO,L)

        if(L <= maxT)then  ! Troposphere:
          B=B
     &      +rr(rrbi%CH3O2_NO__HCHO_NO2,L)*y(nCH3O2,L)
     &      +rr(rrbi%C2O3_NO__HCHO_NO2,L)*y(nC2O3,L)
     &      +rr(rrbi%XO2_NO__NO2_M,L)*y(nXO2,L)
        else               ! Stratosphere:
          B=B
     &      +rr(rrbi%ClO_NO__NO2_Cl,L)*y(nClO,L)
     &      +rr(rrbi%NO_OClO__NO2_ClO,L)*y(nOClO,L)
     &      +rr(rrbi%BrO_NO__Br_NO2,L)*y(nBrO,L)
        end if

C       C is for NO2->NO reactions :
        C=ss(rj%NO2__NO_O,L,I,J)
     &    +rr(rrbi%O_NO2__NO_O2,L)*y(nO,L)
        ! below forms NO3, assume some goes to NO:
        C=C
     &    +rr(rrbi%NO2_O3__NO3_O2,L)*y(nO3,L)
     &      *ss(rj%NO3__NO_O2,L,I,J)
     &      /(ss(rj%NO3__NO_O2,L,I,J)
     &        +ss(rj%NO3__NO2_O,L,I,J))
        p2=B/(B+C)
        p1=1-p2

C       Set NO3: D is loss rxns NO3->NO2 or NO
        D=ss(rj%NO3__NO_O2,L,I,J)
     &    +ss(rj%NO3__NO2_O,L,I,J)
     &    +rr(rrbi%NO3_NO__NO2_NO2,L)*p1*y(nn_NOx,L)
     &    +rr(rrbi%NO2_NO3__NO_NO2,L)*p2*y(nn_NOx,L)
     &    +rr(rrbi%NO3_NO3__NO2_NO2,L)*yNO3(I,J,L)
     &    +rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L)

        yNO3(I,J,L)=(rr(rrbi%NO2_O3__NO3_O2,L)*y(nO3,L)*p2
     &    *y(nn_NOx,L))/D
        if(yNO3(I,J,L).ge.1.d-1*y(nn_NOx,L))
     &    yNO3(I,J,L)=1.d-1*y(nn_NOx,L)
        y(nNO,L)= p1*(y(nn_NOx,L)-yNO3(I,J,L))
        y(nNO2,L)=p2*(y(nn_NOx,L)-yNO3(I,J,L))

C       Set limits on NO, NO2, NOx:
        if(y(nNO,L)   < 1.)   y(nNO,L) = 1.d0
        if(y(nNO2,L)  < 1.)  y(nNO2,L) = 1.d0
        if(y(nn_NOx,L) < 1.) y(nn_NOx,L) = 1.d0
        pNOx(I,J,L)=y(nNO2,L)/y(nn_NOx,L)
        pNO3(I,J,L)=yNO3(I,J,L)/y(nn_NOx,L)
        y(nNO3,L)=yNO3(I,J,L)
        y(nHONO,L)=1.d0
      end do

      return
      END SUBROUTINE NOxfam


      SUBROUTINE HOxfam(Lmax,I,J)
!@sum HOxfam Find HOx family (OH,HO2) partitioning assuming equilibrium
!@+   concentration of HOx.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only : LS1=>LS1_NOMINAL
      USE RESOLUTION, only : LM
      USE GEOM, only : LAT2D_DG
      USE ATM_COM, only: LTROPO,PMIDL00

      USE TRACER_COM, only : rsulf1,rsulf2,rsulf4
      USE TRACER_COM, only : nn_CH4,nn_HNO3,nn_CH3OOH,nn_H2O2,nn_HCHO,
     &                       nn_CO,nn_Paraffin,nn_Alkenes,nn_Isoprene,
     &                       nn_AlkylNit,nn_Terpenes,
     &                       nn_HBr,nn_HOCl,nn_HCl

      use photolysis, only: rj
      USE TRCHEM_Shindell_COM, only:pHOx,rr,y,nNO2,nNO,nH2O,nO3,nCH3O2,
     &                        nO2,nM,nHO2,nOH,nH2,nAldehyde,nXO2,nXO2N,
     &                        ta,ss,nC2O3,nROR,yso2,ydms,which_trop,nO1D
     &         ,OxlossbyH,dt2,nBrO,nClO,nOClO,nBr,nCl,SF3,nO
     &         ,rrbi,rrtri,yNO3

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var aqqz,bqqz,cqqz,cz,dz,sqroot,temp_yHOx,ratio dummy vars
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var Lmax maximum altitude for chemistry
!@var rHprod,rHspecloss,rkzero,rktot temporary var during OH->H rxns
!@var maxT LTROPO(I,J) or LS1-1, depending upon what_trop variable
!@+ Or the top layer of chemistry in the unlikely event that is lower.
!@+ Note in that case, loops like L=maxT+1,Lmax will do nothing.

      integer             :: L, maxT 
      integer, intent(IN) :: Lmax,I,J
      real*8              :: aqqz, bqqz, cqqz, cz, dz, sqroot, 
     &   temp_yHOx,ratio,rHprod,rHspecloss,rkzero,rktot,
     &   yAtomicH

      select case(which_trop)
      case(0); maxT=min(ltropo(I,J),Lmax)
      case(1); maxT=min(ls1-1,Lmax)
      case default; call stop_model('which_trop problem 9',255)
      end select

C ------------ Troposphere ------------
      do L=1,maxT   ! >> beginning of altitude loop <<
c First calculate equilibrium amount of HOx:
c A: loss rxns with HOx**2
c B: loss rxns linear in HOx
c C: prod equations
c all: in terms of HO2 (so *pHOx when OH is reactant)
        aqqz=2.d0*(pHOx(I,J,L)*rr(rrbi%OH_HO2__H2O_O2,L)
     &    +pHOx(I,J,L)*pHOx(I,J,L)
     &      *(rr(rrbi%OH_OH__H2O_O,L)
     &        +rr(rrtri%OH_OH__H2O2_M,L))
     &      +rr(rrbi%HO2_HO2__H2O2_O2,L))

        bqqz=pHOx(I,J,L)
     &    *(rr(rrbi%CH4_OH__H2O_CH3O2,L)*y(nn_CH4,L)
     &      +rr(rrbi%OH_HNO3__H2O_NO3,L)*y(nn_HNO3,L)
     &      +rr(rrtri%OH_NO2__HNO3_M,L)*y(nNO2,L)
     &      +rr(rrtri%OH_NO__HONO_M,L)*y(nNO,L)
     &      +rr(rrbi%CH3OOH_OH__CH3O2_H2O,L)*y(nn_CH3OOH,L))
     &    +rr(rrbi%CH3O2_HO2__CH3OOH_O2,L)*y(nCH3O2,L)
     &    +pHOx(I,J,L)
     &    *(rr(rrbi%Aldehyde_OH__C2O3_M,L)*y(nAldehyde,L)
     &      +rr(rrbi%Paraffin_OH__HO2_M,L)*y(nn_Paraffin,L)*0.89d0
!!!!!&      +rr(rrbi%Alkenes_OH__HCHO_HO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_OH__HCHO_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.15d0
     &      +rr(rrbi%AlkylNit_OH__NO2_XO2,L)*y(nn_AlkylNit,L)
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_OH__HCHO_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.15d0
#endif  /* TRACERS_TERP */
     &    )
     &    +rr(rrbi%XO2_HO2__CH3OOH_O2,L)*y(nXO2,L)
     &    +rr(rrbi%XO2N_HO2__CH3OOH_O2,L)*y(nXO2N,L)
     &    +pHOx(I,J,L)*(rsulf1(i,j,l)*ydms(i,j,l) ! oxidation of DMS
     &      +rsulf2(i,j,l)*ydms(i,j,l)) ! oxidation of SO2

        cqqz=2.d0*ss(rj%H2O2__OH_OH,L,I,J)*y(nn_H2O2,L)
     &    +ss(rj%HNO3__OH_NO2,L,I,J)*y(nn_HNO3,L)
     &    +2.d0*ss(rj%HCHO__CO_HO2,L,I,J)*y(nn_HCHO,L)
     &    +2.d0*ss(rj%CH3OOH__HCHO_HO2,L,I,J)*y(nn_CH3OOH,L)
     &    +(rr(rrbi%CH3O2_NO__HCHO_NO2,L)*y(nNO,L)
     &      +2.d0*rr(rrbi%CH3O2_CH3O2__HCHO_HCHO,L)*y(nCH3O2,L)
     &    +rr(rrbi%ClO_CH3O2__Cl_HCHO,L)*y(nClO,L)
     &    )*y(nCH3O2,L)

        ! 1.66/1.31 accounts for HOx production via O(1D)+CH4-->CH3O2 path:
        cqqz=cqqz
     &    +((2.d0*rr(rrbi%O1D_H2O__OH_OH,L)*y(nH2O,L)
     &      +(1.66d0/1.31d0)*rr(rrbi%O1D_CH4__OH_CH3O2,L)*y(nn_CH4,L)
     &    )*y(nO1D,L))
     &    +ss(rj%Aldehyde__HCHO_CO,L,I,J)*y(nAldehyde,L)*2.d0
     &    +(rr(rrbi%C2O3_NO__HCHO_NO2,L)*y(nNO,L)
     &    +rr(rrbi%C2O3_C2O3__HCHO_HCHO,L)*y(nC2O3,L)*2.d0)*y(nC2O3,L)
     &    +(rr(rrbi%ROR_M__Aldehyde_HO2,L)*y(nM,L)*0.94d0
     &      +rr(rrbi%ROR_M__HO2_M,L))*y(nROR,L)
     &    +rr(rrbi%Alkenes_O3__HCHO_CO,L)*y(nn_Alkenes,L)*y(nO3,L)
     &      *0.65d0
     &    +rr(rrbi%Isoprene_O3__HCHO_Alkenes,L)*y(nn_Isoprene,L)
     &      *y(nO3,L)*0.58d0
     &    +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &      *yNO3(I,J,L)*0.9d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__HCHO_Alkenes,L)*y(nn_Terpenes,L)
     &      *y(nO3,L)*0.58d0
#endif  /* TRACERS_TERP */

        sqroot=sqrt(bqqz*bqqz+4.d0*aqqz*cqqz)
        y(nHO2,L)=(sqroot-bqqz)/(2.d0*aqqz)
        y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
        temp_yHOx=y(nOH,L)+y(nHO2,L)

c Now partition HOx into OH and HO2:
        ! CZ: OH->HO2 reactions :
        cz=rr(rrbi%OH_O3__HO2_O2,L)*y(nO3,L)
     &    +rr(rrbi%CO_OH__HO2_O2,L)*y(nn_CO,L)
     &    +rr(rrbi%OH_H2O2__H2O_HO2,L)*y(nn_H2O2,L)
     &    +rr(rrbi%H2_OH__HO2_H2O,L)*y(nH2,L)
     &    +rr(rrbi%HCHO_OH__HO2_CO,L)*y(nn_HCHO,L)
     &    +rr(rrbi%Paraffin_OH__HO2_M,L)*y(nn_Paraffin,L)
     &      *0.11d0
     &    +rr(rrbi%Isoprene_OH__HCHO_Alkenes,L)*y(nn_Isoprene,L)*0.85d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_OH__HCHO_Alkenes,L)*y(nn_Terpenes,L)*0.85d0
#endif  /* TRACERS_TERP */
     &    +rr(rrbi%Alkenes_OH__HCHO_HO2,L)*y(nn_Alkenes,L)
     &    +rsulf4(i,j,l)*yso2(i,j,l) ! SO2 oxidation: 

        ! DZ: HO2->OH reactions :
        dz=rr(rrbi%HO2_O3__OH_O2,L)*y(nO3,L)
     &    +rr(rrbi%HO2_NO__OH_NO2,L)*y(nNO,L)
     &    +rr(rrbi%C2O3_HO2__HCHO_HO2,L)*(0.79d0*y(nC2O3,L))
c Previous few lines represent additional OH production via reaction 41
c which also produces HO2 and R15 then S4/(S4+S14) fraction.

        if(cz+dz > 0.)then
          y(nOH,L)=(dz/(cz+dz))*temp_yHOx
          if(y(nOH,L) > temp_yHOx) y(nOH,L)=temp_yHOx-1.d0
          ! here we USED TO cap OH as a function of latitude
        else
          y(nOH,L)=1.d0
        end if
        y(nHO2,L)=(temp_yHOx-y(nOH,L))
        ! Some limits on OH, HO2:
        if(y(nOH,L) < 1.d0)y(nOH,L)=1.d0
        if(y(nHO2,L) < 1.d0)y(nHO2,L)=1.d0
        if(y(nHO2,L) > 1.d9)y(nHO2,L)=1.d9
        pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)

      end do  ! >> end of troposphere loop <<

C ------------ Stratosphere ------------
      do L=maxT+1,Lmax
c First calculate equilibrium amount of HOx:
c A: loss rxns with HOx**2
c B: loss rxns linear in HOx
c C: prod equations
c all: in terms of HO2 (so *pHOx when OH is reactant)
        aqqz=2.d0*(pHOx(I,J,L)*rr(rrbi%OH_HO2__H2O_O2,L)
     &    +pHOx(I,J,L)*pHOx(I,J,L)
     &      *(rr(rrbi%OH_OH__H2O_O,L)
     &        +rr(rrtri%OH_OH__H2O2_M,L))
     &      +rr(rrbi%HO2_HO2__H2O2_O2,L))

        bqqz=pHOx(I,J,L)
     &    *(rr(rrbi%CH4_OH__H2O_CH3O2,L)*y(nn_CH4,L)
     &      +rr(rrbi%OH_HNO3__H2O_NO3,L)*y(nn_HNO3,L)
     &      +rr(rrtri%OH_NO2__HNO3_M,L)*y(nNO2,L)
     &      +rr(rrtri%OH_NO__HONO_M,L)*y(nNO,L))
     &    +rr(rrbi%OH_HCl__H2O_Cl,L)*y(nn_HCl,L)*pHOx(I,J,L)
     &    +rr(rrbi%OH_HOCl__H2O_ClO,L)*y(nn_HOCl,L)*pHOx(I,J,L)
     &    +rr(rrbi%OClO_OH__HOCl_O2,L)*y(nOClO,L)*pHOx(I,J,L)
     &    +rr(rrbi%Cl_HO2__HCl_O2,L)*y(nCl,L)
     &    +rr(rrbi%ClO_OH__HCl_O2,L)*y(nClO,L)*pHOx(I,J,L)
     &    +rr(rrbi%ClO_HO2__HOCl_O2,L)*y(nClO,L)
     &    +rr(rrbi%HBr_OH__H2O_Br,L)*y(nn_HBr,L)*pHOx(I,J,L)
     &    +rr(rrbi%Br_HO2__HBr_O2,L)*y(nBr,L)
     &    +rr(rrbi%BrO_HO2__HOBr_O2,L)*y(nBrO,L)
     &    +rr(rrbi%BrO_OH__HBr_O2,L)*y(nBrO,L)*pHOx(I,J,L)
     &    +pHOx(I,J,L)*(rsulf1(i,j,l)*ydms(i,j,l) ! oxidation of SO2
     &    +rsulf2(i,j,l)*ydms(i,j,l)) ! oxidation of DMS

        ! Use OH production without O1D explicitly:
        cqqz=2.d0*ss(rj%H2O2__OH_OH,L,i,j)*y(nn_H2O2,L)
     &    +ss(rj%HNO3__OH_NO2,L,i,j)*y(nn_HNO3,L)
     &    +ss(rj%HOCl__OH_Cl,L,i,j)*y(nn_HOCl,L) 
     &    +rr(rrbi%O_HCl__OH_Cl,L)*y(nn_HCl,L)*y(nO,L)
     &    +rr(rrbi%O_HOCl__OH_ClO,L)*y(nn_HOCl,L)*y(nO,L)
     &    +rr(rrbi%Cl_HOCl__Cl2_OH,L)*y(nn_HOCl,L)*y(nCl,L)
     &    +rr(rrbi%Cl_H2O2__HCl_HO2,L)*y(nCl,L)*y(nn_H2O2,L)
     &    +rr(rrbi%Cl_H2__HCl_HO2,L)*y(nCl,L)*y(nH2,L)
     &    +rr(rrbi%Br_H2O2__HBr_HO2,L)*y(nBr,L)*y(nn_H2O2,L)
     &    +rr(rrbi%O_HBr__OH_Br,L)*y(nn_HBr,L)*y(nO,L)
     &    +rr(rrbi%ClO_CH3O2__Cl_HCHO,L)*y(nClO,L)*y(nCH3O2,L)
     
        ! water vapor photolysis in SRBs:
        if(PMIDL00(L) < 10.d0) cqqz = cqqz + 0.5d0*SF3(I,J,L)*y(nH2O,L)

        ! production from O1D NO LONGER limited to O1D amount or
        ! O1D fraction via r10,11,s2:
        cqqz=cqqz
     &    +(2.d0*rr(rrbi%O1D_H2O__OH_OH,L)*y(nH2O,L)
     &      +(1.66d0/1.31d0)*rr(rrbi%O1D_CH4__OH_CH3O2,L)*y(nn_CH4,L)
     &    )*y(nO1D,L)

        sqroot=sqrt(bqqz*bqqz+4.d0*aqqz*cqqz)
        y(nHO2,L)=(sqroot-bqqz)/(2.d0*aqqz)
        y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
        temp_yHOx=y(nOH,L)+y(nHO2,L)

c Include loss of OH into atomic H using production
c via OH + O -> O2 + H, loss via H + O3 -> OH + O2 and
c H + O2 + M -> HO2 + M , and affects on OH/HO2 and Ox
        rHprod=rr(rrbi%O_OH__O2_H,L)*y(nOH,L)*y(nO,L)
        rHspecloss=y(nO3,L)*1.4d-10*exp(-470./ta(L))
        rkzero=y(nM,L)*4.4d-32*((ta(L)/300.d0)**(-1.3))
        rktot=(rkzero/(1+(rkzero/(7.5d-11*(ta(L)/300.d0)**0.2))))
        rktot=y(nO2,L)*rktot
        rHspecloss=rHspecloss+rktot
        if(rHspecloss==0. .or. rHspecloss+rktot==0.)
     &  call stop_model('rHspecloss or rHspecloss+rktot=0.',255)
        yAtomicH=rHprod/rHspecloss
        OxlossbyH(L)= yAtomicH*rHspecloss*dt2

c Now partition HOx into OH and HO2:
c CZ: OH->HO2 reactions :
        cz=rr(rrbi%OH_O3__HO2_O2,L)*y(nO3,L)
     &    +rr(rrbi%CO_OH__HO2_O2,L)*y(nn_CO,L)
     &    +rr(rrbi%OH_H2O2__H2O_HO2,L)*y(nn_H2O2,L)
     &    +rr(rrbi%H2_OH__HO2_H2O,L)*y(nH2,L)
     &    +rr(rrbi%HCHO_OH__HO2_CO,L)*y(nn_HCHO,L)
     &    +rr(rrbi%ClO_OH__HO2_Cl,L)*y(nClO,L)
     &    +rr(rrbi%BrO_OH__Br_HO2,L)*y(nBrO,L)
     &    +rr(rrbi%O_OH__O2_H,L)*y(nO,L)*rktot/(rHspecloss+rktot)
     &    +rsulf4(i,j,l)*yso2(i,j,l) ! SO2 oxidation: 

        dz=rr(rrbi%HO2_O3__OH_O2,L)*y(nO3,L)
     &    +rr(rrbi%HO2_NO__OH_NO2,L)*y(nNO,L)
     &    +rr(rrbi%Cl_HO2__OH_ClO,L)*y(nCl,L)
     &    +rr(rrbi%O_HO2__OH_O2,L)*y(nO,L)

        if(cz+dz > 0)then
          y(nOH,L)=(dz/(cz+dz))*temp_yHOx-yAtomicH
          if(y(nOH,L) > temp_yHOx)y(nOH,L)=temp_yHOx-1.d0
        else
          y(nOH,L)=1.d0
        end if
        y(nHO2,L)=(temp_yHOx-y(nOH,L))

        if(y(nOH,L)  < 1)   y(nOH,L)  = 1.d0
        if(y(nHO2,L) < 1)   y(nHO2,L) = 1.d0
        if(y(nHO2,L) > 1.d9)y(nHO2,L) = 1.d9
        pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)

      end do  ! end of stratosphere loop

      return
      END SUBROUTINE HOxfam



      SUBROUTINE ClOxfam(Lmax,I,J)
!@sum ClOxfam Find ClOx family (Cl,ClO,OClO,Cl2,Cl2O2) partitioning 
!@+   assuming equilibrium within ClOx.
!@auth Drew Shindell

C**** GLOBAL parameters and variables:

      USE ATM_COM, only   : LTROPO
      USE RESOLUTION, only : LS1=>LS1_NOMINAL
      USE TRACER_COM, only : n_ClOx,n_HOCl,n_ClONO2,n_HCl,n_H2O2,n_CH4
      USE TRACER_COM, only : nn_ClOx,nn_HOCl,nn_ClONO2,nn_HCl,nn_H2O2,
     &    nn_CH4
      use photolysis, only : sza,rj
      USE TRCHEM_Shindell_COM, only:pClOx,rr,y,nClO,nOClO,nCl,nCl2O2,
     &    ta,ss,nO3,nHO2,nNO3,nO,nNO,nBr,nOH,nBrO,nCH3O2,nM,nCl2,nH2,
     &    dt2,pClx,pOClOx,nNO2,which_trop,yCl2,yCl2O2,ClOx_old,
     &    rrmono,rrbi,rrtri

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var A,B,C,D,F,G,Q,V,X,YY,dClOx,ww,p1,p2,p3,dOClO,ratioc dummy vars
!@var ClOx_old total ClOx at start of chemical timestep
!@var rnormnum is temporary var for conservation of Cl
!@var destCl, prodCl temporary vars for simple steady state Cl calc
!@var maxT LTROPO(I,J) or LS1-1, depending upon what_trop variable
!@+ Or the top layer of chemistry in the unlikely event that is lower.
!@+ Note in that case, loops like L=maxT+1,Lmax will do nothing.
      integer               :: L,maxT
      integer, intent(IN)   :: Lmax,I,J
      real*8                :: A,B,C,D,F,G,Q,V,X,YY,dClOx,ww,p1,p2,p3,
     &                         dOClO,ratioc,rnormnum,destCl,prodCl

      select case(which_trop)
      case(0); maxT=min(ltropo(I,J),Lmax)
      case(1); maxT=min(ls1-1,Lmax)
      case default; call stop_model('which_trop problem 10',255)
      end select

      do L=maxT+1,Lmax ! stratosphere
c Set Cl2 and default Cl2O2:
        y(nCl2,L)=yCl2(I,J,L)       ! set in chemstep
        yCl2O2(I,J,L)=0.d0          ! non-zero at low temp, see below
        y(nCl2O2,L)=yCl2O2(I,J,L)   ! non-zero at low temp, see below

c Full ClOxfam code from offline photochemistry:
        y(nClO,L)=y(nn_ClOx,L)*pClOx(I,J,L)
        y(nOClO,L)=y(nn_ClOx,L)*pOClOx(I,J,L)
        y(nCl,L)=y(nn_ClOx,L)*pClx(I,J,L)
        if(y(nn_ClOx,L) == 0) CYCLE

c Low temperature stabilizes ClO dimer, use [Cl2O2] only for
c calculating Cl amount, otherwise ignore:
        if(ta(L) <= 220.)then
          y(nCl2O2,L)=
     &       (rr(rrtri%ClO_ClO__Cl2O2_M,L)*y(nClO,L)*y(nClO,L))
     &      /(rr(rrmono%Cl2O2_M__ClO_ClO,L)
     &       +ss(rj%Cl2O2__Cl_Cl,L,i,j))
          if(y(nCl2O2,L) > y(nClO,L)) y(nCl2O2,L)=y(nClO,L)
          y(nClO,L)=y(nClO,L)-y(nCl2O2,L)
        end if

! Below is some documentation based on some very old notes from Drew. A word
! of caution, they are from the mid-1990s and the current code might be a bit
! (or a lot) different, but it is the best piece of documentation we can have.
!
! ClO equilibrium reactions
! A: Cl to ClO
! B: OClO to ClO
! C: ClO to Cl or OClO
!
! OClO equilibrium reactions
! D: OClO production
! E: OClO loss
!
! Note that some reactions are going to appear in both the ClO and OClO
! equilibrium sections
!
! F: non-family (i.e. non-ClOx) to ClO
! G: ClO to non-family
!
! Cl equilibrium reactions
! V: Cl production within ClOx family
! W: Cl loss within ClOx family
! X: Cl production from outside ClOx family
! YY: Cl loss from outside ClOx family

        A=y(nO3,L)*rr(rrbi%Cl_O3__ClO_O2,L)
     &    +y(nOClO,L)*rr(rrbi%Cl_OClO__ClO_ClO,L)
     &    +y(nHO2,L)*rr(rrbi%Cl_HO2__OH_ClO,L)
        B=y(nCl,L)*rr(rrbi%Cl_OClO__ClO_ClO,L)
     &    +y(nO,L)*rr(rrbi%O_OClO__ClO_O2,L)
     &    +y(nNO,L)*rr(rrbi%NO_OClO__NO2_ClO,L)
     &    +y(nBr,L)*rr(rrbi%Br_OClO__BrO_ClO,L)
     &    +ss(rj%OClO__O_ClO,L,i,j)
        C=y(nO,L)*rr(rrbi%ClO_O__Cl_O2,L)
     &    +y(nO3,L)*(rr(rrbi%ClO_O3__Cl_O2,L)
     &      +rr(rrbi%ClO_O3__OClO_O2,L))
     &    +y(nOH,L)*rr(rrbi%ClO_OH__HO2_Cl,L)
     &    +y(nNO,L)*rr(rrbi%ClO_NO__NO2_Cl,L)
     &    +y(nBrO,L)*(rr(rrbi%BrO_ClO__OClO_Br,L)
     &      +rr(rrbi%BrO_ClO__Br_Cl,L))
     &    +ss(rj%ClO__Cl_O,L,i,j)
     &    +rr(rrbi%ClO_CH3O2__Cl_HCHO,L)*y(nCH3O2,L)
     &    +y(nClO,L)*(1.d-12*exp(-1590./TA(L))
     &      +3.d-11*exp(-2450./TA(L))
     &      +3.5d-13*exp(-1370./TA(L))) 
        D=y(nO3,L)*rr(rrbi%ClO_O3__OClO_O2,L)
     &    +y(nBrO,L)*rr(rrbi%BrO_ClO__OClO_Br,L)
        F=(rr(rrbi%OH_HOCl__H2O_ClO,L)*y(nOH,L)*y(nn_HOCl,L)
     &      +rr(rrbi%O_HOCl__OH_ClO,L)*y(nO,L)*y(nn_HOCl,L)
     &      +rr(rrbi%ClONO2_O__ClO_NO3,L)*y(nn_ClONO2,L)*y(nO,L)
     &    )/y(nn_ClOx,L)
        G=rr(rrbi%ClO_OH__HCl_O2,L)*y(nOH,L)
     &    +rr(rrbi%ClO_HO2__HOCl_O2,L)*y(nHO2,L)
     &    +2.d0*rr(rrtri%ClO_ClO__Cl2O2_M,L)*y(nClO,L)
     &    +rr(rrtri%ClO_NO2__ClONO2_M,L)*y(nNO2,L)
        Q=rr(rrbi%OClO_OH__HOCl_O2,L)*y(nOH,L)
        V=C-D
        X=rr(rrbi%OH_Cl2__HOCl_Cl,L)*y(nOH,L)*y(nCl2,L)
     &    +rr(rrbi%O_HCl__OH_Cl,L)*y(nO,L)*y(nn_HCl,L)
     &    +2.d0*ss(rj%Cl2__Cl_Cl,L,i,j)*y(nCl2,L)
     &    +2.d0*ss(rj%Cl2O2__Cl_Cl,L,i,j)*y(nCl2O2,L)
     &    +ss(rj%HOCl__OH_Cl,L,i,j)*y(nn_HOCl,L)
     &    +ss(rj%ClONO2__Cl_NO3,L,i,j)*y(nn_ClONO2,L)
        X=X/y(nn_ClOx,L)
        YY=rr(rrbi%Cl_HOCl__Cl2_OH,L)*y(nn_HOCl,L)
     &    +rr(rrbi%Cl_H2O2__HCl_HO2,L)*y(nn_H2O2,L)
     &    +rr(rrbi%Cl_HO2__HCl_O2,L)*y(nHO2,L)
     &    +rr(rrbi%Cl_CH4__HCl_CH3O2,L)*y(nn_CH4,L)
     &    +rr(rrbi%Cl_H2__HCl_HO2,L)*y(nHO2,L)
        if((dt2*y(nn_ClOx,L)) /= 0)then
          dClOx=(y(nn_ClOx,L)-ClOx_old(L))/(dt2*y(nn_ClOx,L))
        else
          dClOx=0.d0
        end if
        ww=A+yy+dClOx

        if(v /= 0)then 
          p1=(-(ww/v)*(B+C+G+dClOx))+A-B
          if(p1 /= 0)then
            p1=abs(-((x/v)*(B+C+G+dClOx)+F+B)/p1)
            if((B+D-Q-dClOx) /= 0)then
              p3=abs((-p1*D+D)/(B+D-Q-dClOx))
            else 
              p3=0.d0
            end if   
            p2=1-p1-p3
          else 
            p2=1.d0
            p1=0.d0
            p3=0.d0  
          end if  
        else 
          p2=1.d0 
          p1=0.d0
          p3=0.d0  
        end if          
       
        if(SZA > 90.)then
          dOClO=(y(nClO,L)*D-y(nOClO,L)*(B+Q))*dt2
          y(nOClO,L)=y(nOClO,L)+dOClO
          if(y(nOClO,L) < 0.) y(nOClO,L)=0.d0
        end if 
       
        y(nCl,L)=p1*y(nn_ClOx,L)
        y(nClO,L)=p2*y(nn_ClOx,L)
  
        if(SZA < 90.) y(nOClO,L)=p3*y(nn_ClOx,L)  
        if(y(nCl,L) < 0)   y(nCl,L)   =0.d0
        if(y(nClO,L) < 1)  y(nClO,L)  =0.d0
        if(y(nOClO,L) < 1) y(nOClO,L) =0.d0
        if(y(nn_ClOx,L) < 1)y(nn_ClOx,L)=1.0
           
c Normalize so that amount of ClOx doesn't change:
        if((y(nCl,L)+y(nClO,L)+y(nOClO,L)) > 0.)then
          rnormnum=y(nn_ClOx,L)/(y(nCl,L)+y(nClO,L)+y(nOClO,L))
          y(nCl,L)=y(nCl,L)*rnormnum
          y(nClO,L)=y(nClO,L)*rnormnum
          y(nOClO,L)=y(nOClO,L)*rnormnum
        end if

        pClOx(I,J,L)=y(nClO,L)/y(nn_ClOx,L)
        pClx(I,J,L)=y(nCl,L)/y(nn_ClOx,L)
        pOClOx(I,J,L)=y(nOClO,L)/y(nn_ClOx,L)

      end do  ! end of stratosphere loop

      yCl2O2(I,J,:)=y(nCl2O2,:)

      return
      END SUBROUTINE ClOxfam



      SUBROUTINE BrOxfam(Lmax,I,J)
!@sum BrOxfam Find BrOx family (Br,BrO) partitioning 
!@+   assuming equilibrium within BrOx.
!@auth Drew Shindell

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only : LS1=>LS1_NOMINAL
      USE ATM_COM, only    : LTROPO
      USE TRACER_COM, only : n_BrOx,n_H2O2,n_HBr,n_HOBr,n_BrONO2
      USE TRACER_COM, only : nn_BrOx,nn_H2O2,nn_HBr,nn_HOBr,nn_BrONO2
      use photolysis, only: rj
      USE TRCHEM_Shindell_COM, only:rr,y,nO3,nClO,nOClO,nNO,nO,nBr,nOH,
     &    nBrO,ss,nHO2,nNO2,pBrOx,which_trop,rrbi,rrtri

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var a,b,c,d,eq,f,p2,p1 dummy vars
!@var maxT LTROPO(I,J) or LS1-1, depending upon what_trop variable
!@+ Or the top layer of chemistry in the unlikely event that is lower.
!@+ Note in that case, loops like L=maxT+1,Lmax will do nothing.
      integer             :: L,maxT
      integer, intent(IN) :: Lmax,I,J
      real*8              :: a,b,c,d,eq,f,p2,p1
      
      select case(which_trop)
      case(0); maxT=min(ltropo(I,J),Lmax)
      case(1); maxT=min(ls1-1,Lmax)
      case default; call stop_model('which_trop problem 11',255)
      end select

      do L=maxT+1,Lmax  ! stratosphere
        a=y(nO3,L)*rr(rrbi%Br_O3__BrO_O2,L)
     &    +y(nOClO,L)*rr(rrbi%Br_OClO__BrO_ClO,L)
        b=y(nO,L)*rr(rrbi%BrO_O__Br_O2,L)
     &    +y(nNO,L)*rr(rrbi%BrO_NO__Br_NO2,L)
     &    +y(nClO,L)
     &      *(rr(rrbi%BrO_ClO__OClO_Br,L)
     &      +rr(rrbi%BrO_ClO__Br_Cl,L))
     &    +2*y(nBrO,L)*rr(rrbi%BrO_BrO__Br_Br,L)
     &    +y(nOH,L)*rr(rrbi%BrO_OH__Br_HO2,L)
     &    +ss(rj%BrO__Br_O,L,i,j)
        c=rr(rrbi%BrO_HO2__HOBr_O2,L)*y(nHO2,L)
     &    +rr(rrbi%BrO_OH__HBr_O2,L)*y(nOH,L)
     &    +rr(rrtri%BrO_NO2__BrONO2_M,L)*y(nNO2,L)    
        d=rr(rrbi%Br_HO2__HBr_O2,L)*y(nHO2,L)
     &    +rr(rrbi%Br_H2O2__HBr_HO2,L)*y(nn_H2O2,L)
        eq=rr(rrbi%HBr_OH__H2O_Br,L)*y(nn_HBr,L)*y(nOH,L)
     &    +rr(rrbi%O_HBr__OH_Br,L)*y(nn_HBr,L)*y(nO,L)
     &    +ss(rj%HOBr__Br_OH,L,i,j)*y(nn_HOBr,L)
        f=ss(rj%BrONO2__BrO_NO2,L,i,j)*y(nn_BrONO2,L)
        if(a+b /= 0)then
          p2=a/(a+b)
        else
          p2=0.d0
        end if
        if(p2 < 0)p2=0.d0
        if(p2 > 1)p2=1.d0
        p1=1-p2    
        y(nBr,L)=p1*y(nn_BrOx,L)
        y(nBrO,L)=p2*y(nn_BrOx,L)
        if(y(nBr,L) < 1)   y(nBr,L)   =0.d0
        if(y(nBrO,L) < 1)  y(nBrO,L)  =0.d0
        if(y(nBr,L) > 1d9) y(nBr,L)   =0.d0
        if(y(nBrO,L) > 1d9)y(nBrO,L)  =0.d0
        if(y(nn_BrOx,L) < 1)y(nn_BrOx,L)=1.d0
        pBrOx(I,J,L)=y(nBrO,L)/y(nn_BrOx,L)
      end do ! end of stratosphere loop

      return
      END SUBROUTINE BrOxfam


