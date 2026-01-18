#include "rundeck_opts.h"

      SUBROUTINE OADVT3 (MA,RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     DT,QLIMIT, OIJL)
!@sum  OADVT advects tracers using the quadratic upstream scheme.
C****
C**** Input:  MB (kg) = mass before advection
C****          DT (s) = time step
C****       MU (kg/s) = west to east mass flux
C****       MV (kg/s) = south to north mass flux
C****       MW (kg/s) = downward vertical mass flux
C****          QLIMIT = whether tracer is non-negative
C**** Input/Output: RM (kg) = tracer mass
C****   RX,RY,RZ (kg) = first moments of tracer mass
C****   RXX,RYY,RZZ,RXY,RYZ,RZX (kg) = second moments of tracer mass
C****       OIJL (kg) = diagnostic accumulation of tracer mass flux
C****
      USE OCEAN, only : im,jm,lmo
      USE OCEAN_DYN, only : mb=>mmi,smu,smv,smw

      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),     DIMENSION
     &     (IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     MA,RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO,3) :: OIJL
      INTEGER I,J,L, J_0H
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT

      logical :: HAVE_NORTH_POLE  ! ,HAVE_SOUTH_POLE
         call getDomainBounds(grid, J_STRT_HALO=J_0H)
         call getDomainBounds(grid, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C        call getDomainBounds(grid, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

C****
C**** Load mass after advection from mass before advection
C****
      MA = MB
C****
C**** Advect the tracer using the Slopes Scheme
C****
      if(qlimit) then
      CALL OADVTX4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMU,5d-1*DT,OIJL(1,J_0H,1,1))
      CALL OADVTY4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMV,     DT,OIJL(1,J_0H,1,2))
      CALL OADVTZ4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMW,     DT,OIJL(1,J_0H,1,3))
      CALL OADVTX4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMU,5d-1*DT,OIJL(1,J_0H,1,1))
      else
      CALL OADVTX3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMU,5d-1*DT,OIJL(1,J_0H,1,1))
      CALL OADVTY3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMV,     DT,OIJL(1,J_0H,1,2))
      CALL OADVTZ3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMW,     DT,OIJL(1,J_0H,1,3))
      CALL OADVTX3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MA,SMU,5d-1*DT,OIJL(1,J_0H,1,1))
      endif

C**** Fill in values at the poles
      if (HAVE_NORTH_POLE) then
      DO 20 L=1,LMO
      DO 20 I=1,IM
      RM(I,JM,L) = RM(1,JM,L)
   20 RZ(I,JM,L) = RZ(1,JM,L)
      end if
C     if (HAVE_SOUTH_POLE) then
C     DO 30 L=1,LMO
C     DO 30 I=1,IM
C     RM(I, 1,L) = RM(IM,1,L)
C  30 RZ(I, 1,L) = RZ(IM,1,L)
C     end if

      RETURN
      END SUBROUTINE OADVT3

      Subroutine OADVTX4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MM,MU,DT,OIJL)
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOU=>LMU
      USE OCEAN, only : nbyzu,i1yzu,i2yzu, nbyzm,i1yzm,i2yzm
      USE OCEANR_DIM, only : ogrid
      Implicit None
C**** Interface variables
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  MU
      Real*8,Intent(InOut),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX, MM,OIJL
C**** Local variables
      Integer :: I,J,L, J1P,JNP, N,NC,NCOURANT
      real*8, dimension(im) :: mudt
      Real*8 :: AM,AM2,A,FM,FX,FY,FZ,FXX,FYY,FZZ,FXY,FYZ,FZX,
     &     MMnew,bymnew, zCOURANT,mcheck,courmax
      Real*8 :: AMim1,FMim1,FMim10,AM_im,FM_im,FM_im0,
     &     FXim1,FYim1,FZim1,FXXim1,FYYim1,FZZim1,FXYim1,FYZim1,FZXim1,
     &     FX_im,FY_im,FZ_im,FXX_im,FYY_im,FZZ_im,FXY_im,FYZ_im,FZX_im

      REAL*8 :: RM0,FM0,FM_PASS,FX_PASS,FXX_PASS

C**** Extract domain decomposition band parameters
      J1P = Max (oGRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (oGRID%J_STOP, JM-1)  !  Exclude north pole

C****
C**** Loop over layers and latitudes
C****
      Do L=1,LMO
      Do J=J1P,JNP
        if(nbyzu(j,l) == 0) cycle


c
c adjust mass fluxes so that courant numbers are never greater than 1
c
        courmax = 0.
        mudt(1:2) = mu(1:2,j,l)*dt
        mudt(im) = mu(im,j,l)*dt
        i=1
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(im)-mudt(i))
          else
            mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        do n=1,nbyzu(j,l)
          i = i1yzu(n,j,l)
          if(i.gt.1) mudt(i-1) = 0. ! zero mudt at west boundary
          mudt(i) = mu(i,j,l)*dt
          do i=max(2,i1yzu(n,j,l)),min(i2yzu(n,j,l),im-1)
            mudt(i+1) = mu(i+1,j,l)*dt
            if(mudt(i).ge.0.) then
              mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
            else
              mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
            endif
            courmax = max(courmax,mudt(i)/mcheck)
          enddo
        enddo
        i = im
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
          else
            mcheck = -mm(1,j,l)-min(0d0,mudt(i)-mudt(1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        if(courmax.gt.1.) then
          ncourant = 1+int(courmax)
          write(6,*) 'ncourant ',j,l,ncourant
          zcourant = 1d0/ncourant
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              mudt(i) = mudt(i)*zcourant
            enddo
          enddo
        else
          ncourant = 1
        endif

        do nc=1,ncourant

c when flow out both sides would cause negative tracer mass, modify moments
          i=1
          if(l <= lmom(i,j)) then
            if(mudt(i).gt.0. .and. mudt(im).lt.0.) then
              call checkfluxo(mudt(im),mudt(i),mm(i,j,l),
     &             rm(i,j,l),rx(i,j,l),rxx(i,j,l))
            endif
          endif
          do n=1,nbyzm(j,l)
          if(i1yzm(n,j,l)>1 .and. i1yzm(n,j,l)==i2yzm(n,j,l)) cycle
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            if(mudt(i).gt.0. .and. mudt(i-1).lt.0.) then
              call checkfluxo(mudt(i-1),mudt(i),mm(i,j,l),
     &             rm(i,j,l),rx(i,j,l),rxx(i,j,l))
            endif
          enddo
          enddo

c
c calculate fluxes at the dateline
c
          i = im
          if(l <= lmou(i,j)) then
            AM = mudt(I)
            if(am.ge.0.) then   ! mass flux is positive or zero
              A  = AM / MM(I,J,L)
              FM = A*(RM(I,J,L)+(1-A)*(RX(I,J,L)+(1-2*A)*RXX(I,J,L)))
              FX  = AM*(A*A*(RX(I,J,L) + 3*(1-A)*RXX(I,J,L)) - 3*FM)
              FXX = AM*(AM*A**3*RXX(I,J,L) - 5*(AM*FM+FX))
              FY = A * (RY(I,J,L) + (1-A)*RXY(I,J,L))
              FXY = AM * (A*A*RXY(I,J,L) - 3*FY)
              FZ = A * (RZ(I,J,L) + (1-A)*RZX(I,J,L))
              FZX = AM * (A*A*RZX(I,J,L) - 3*FZ)
              FYY = A*RYY(I,J,L)
              FZZ = A*RZZ(I,J,L)
              FYZ = A*RYZ(I,J,L)

              fm0 = fm
              if(fm.lt.0.) then
                fm=0.
                fx_pass=0
                fxx_pass=0.
              elseif(fm.gt.rm(i,j,l)) then
                fm=rm(i,j,l)
                fx_pass=am*(-3.*fm)
                fxx_pass=am*(-5.*(am*fm+fx_pass))
              else
                fx_pass = fx
                fxx_pass = fxx
              endif
              fm_pass = fm

            else                ! mass flux is negative
              A  = AM / MM(1,J,L)
              FM = A*(RM(1,J,L) -(1+A)*(
     &             RX(1,J,L) -(1+2*A)*RXX(1,J,L)))
              FX  = AM*(A*A*(RX(1,J,L) -3*(1+A)*RXX(1,J,L)) -3*FM)
              FXX = AM*(AM*A**3*RXX(1,J,L) - 5*(AM*FM+FX))
              FY = A * (RY(1,J,L) - (1+A)*RXY(1,J,L))
              FXY = AM * (A*A*RXY(1,J,L) - 3*FY)
              FZ = A * (RZ(1,J,L) - (1+A)*RZX(1,J,L))
              FZX = AM * (A*A*RZX(1,J,L) - 3*FZ)
              FYY = A*RYY(1,J,L)
              FZZ = A*RZZ(1,J,L)
              FYZ = A*RYZ(1,J,L)

              fm_pass = fm
              fx_pass = fx
              fxx_pass = fxx
              if(fm.gt.0.) then
                fm=0.
                fx=0.
                fxx=0.
              elseif(fm.lt.-rm(1,j,l)) then
                fm=-rm(1,j,l)
                fx=am*(-3.*fm)
                fxx=am*(-5.*(am*fm+fx))
              endif
              fm0 = fm

            endif
          else
            am = 0.
            fm = 0.
            fx = 0.
            fy = 0.
            fz = 0.
            fxx = 0.
            fyy = 0.
            fzz = 0.
            fxy = 0.
            fyz = 0.
            fzx = 0.
            fm0 = 0.
            fm_pass = 0.
            fx_pass = 0.
            fxx_pass = 0.
          endif
          amim1 = am
          fmim1 = fm
          fmim10 = fm_pass
          fxim1 = fx_pass
          fyim1 = fy
          fzim1 = fz
          fxxim1 = fxx_pass
          fyyim1 = fyy
          fzzim1 = fzz
          fxyim1 = fxy
          fyzim1 = fyz
          fzxim1 = fzx
          am_im = am
          fm_im = fm
          fm_im0 = fm0
          fx_im = fx
          fy_im = fy
          fz_im = fz
          fxx_im = fxx
          fyy_im = fyy
          fzz_im = fzz
          fxy_im = fxy
          fyz_im = fyz
          fzx_im = fzx

c
c loop over basins
c
          do n=1,nbyzm(j,l)
            if(i1yzm(n,j,l)>1 .and. i1yzm(n,j,l)==i2yzm(n,j,l)) cycle
            do i=i1yzm(n,j,l),min(i2yzm(n,j,l),im-1)
              AM = mudt(I)
              if(am.ge.0.) then ! mass flux is positive or zero
                A  = AM / MM(I,J,L)
                FM = A*(RM(I,J,L)+(1-A)*(RX(I,J,L)+(1-2*A)*RXX(I,J,L)))
                FX  = AM*(A*A*(RX(I,J,L) + 3*(1-A)*RXX(I,J,L)) - 3*FM)
                FXX = AM*(AM*A**3*RXX(I,J,L) - 5*(AM*FM+FX))
                FY = A * (RY(I,J,L) + (1-A)*RXY(I,J,L))
                FXY = AM * (A*A*RXY(I,J,L) - 3*FY)
                FZ = A * (RZ(I,J,L) + (1-A)*RZX(I,J,L))
                FZX = AM * (A*A*RZX(I,J,L) - 3*FZ)
                FYY = A*RYY(I,J,L)
                FZZ = A*RZZ(I,J,L)
                FYZ = A*RYZ(I,J,L)

                fm0 = fm
                if(fm.lt.0.) then
                  fm=0.
                  fx_pass=0
                  fxx_pass=0.
                elseif(fm.gt.rm(i,j,l)) then
                  fm=rm(i,j,l)
                  fx_pass=am*(-3.*fm)
                  fxx_pass=am*(-5.*(am*fm+fx_pass))
                else
                  fx_pass = fx
                  fxx_pass = fxx
                endif
                fm_pass = fm

              else              ! mass flux is negative
                A  = AM / MM(I+1,J,L)
                FM = A*(RM(I+1,J,L) -(1+A)*(
     &               RX(I+1,J,L) -(1+2*A)*RXX(I+1,J,L)))
                FX  = AM*(A*A*(RX(I+1,J,L) -3*(1+A)*RXX(I+1,J,L)) -3*FM)
                FXX = AM*(AM*A**3*RXX(I+1,J,L) - 5*(AM*FM+FX))
                FY = A * (RY(I+1,J,L) - (1+A)*RXY(I+1,J,L))
                FXY = AM * (A*A*RXY(I+1,J,L) - 3*FY)
                FZ = A * (RZ(I+1,J,L) - (1+A)*RZX(I+1,J,L))
                FZX = AM * (A*A*RZX(I+1,J,L) - 3*FZ)
                FYY = A*RYY(I+1,J,L)
                FZZ = A*RZZ(I+1,J,L)
                FYZ = A*RYZ(I+1,J,L)

                fm_pass = fm
                fx_pass = fx
                fxx_pass = fxx
                if(fm.gt.0.) then
                  fm=0.
                  fx=0.
                  fxx=0.
                elseif(fm.lt.-rm(i+1,j,l)) then
                  fm=-rm(i+1,j,l)
                  fx=am*(-3.*fm)
                  fxx=am*(-5.*(am*fm+fx))
                endif
                fm0 = fm

              endif
              MMnew = MM(I,J,L) + (AMim1-AM)
              bymnew = 1d0/mmnew
              AM2 = AMim1+AM
              RM0 = RM(I,J,L) +  (FMim10-FM0)
              RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
c              RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
c     &             3d0*(AM2*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))*bymnew
c              RXX(I,J,L) = (RXX(I,J,L)*MM(I,J,L)**2
c     &             +2.5*RM(I,J,L)*(MM(I,J,L)**2 -MMnew**2 -3.*AM2*AM2)
c     &             +5.*(MM(I,J,L)*(MM(I,J,L)*(FMim1-FM)-FXim1-FX)
c     &             +AM2*RX(I,J,L)*MMnew)
c     &             +(FXXim1-FXX) ) * (bymnew*bymnew)
              RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &             3d0*(AM2*RM0-MM(I,J,L)*(FMim10+FM0)))*bymnew
              RXX(I,J,L) = (RXX(I,J,L)*MM(I,J,L)**2
     &             +2.5*RM0*(MM(I,J,L)**2 -MMnew**2 -3.*AM2*AM2)
     &             +5.*(MM(I,J,L)*(MM(I,J,L)*(FMim10-FM0)-FXim1-FX)
     &             +AM2*RX(I,J,L)*MMnew)
     &             +(FXXim1-FXX) ) * (bymnew*bymnew)
              RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
              RXY(I,J,L) = (RXY(I,J,L)*MM(I,J,L) + (FXYim1-FXY) +
     &             3d0*(AM2*RY(I,J,L)-MM(I,J,L)*(FYim1+FY)))*bymnew
              RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
              RZX(I,J,L) = (RZX(I,J,L)*MM(I,J,L) + (FZXim1-FZX) +
     &             3d0*(AM2*RZ(I,J,L)-MM(I,J,L)*(FZim1+FZ)))*bymnew
              RYY(I,J,L) = RYY(I,J,L) + (FYYim1-FYY)
              RZZ(I,J,L) = RZZ(I,J,L) + (FZZim1-FZZ)
              RYZ(I,J,L) = RYZ(I,J,L) + (FYZim1-FYZ)
              MM(I,J,L) = MMnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
! clean up roundoff errors
              if(rm(i,j,l).le.0d0) then
                rm(i,j,l)=0d0
                rx(i,j,l)=0d0
                ry(i,j,l)=0d0
                rz(i,j,l)=0d0
                rxx(i,j,l)=0d0
                ryy(i,j,l)=0d0
                rzz(i,j,l)=0d0
                rxy(i,j,l)=0d0
                ryz(i,j,l)=0d0
                rzx(i,j,l)=0d0
              endif
              amim1 = am
              fmim1 = fm
              fmim10 = fm_pass
              fxim1 = fx_pass
              fyim1 = fy
              fzim1 = fz
              fxxim1 = fxx_pass
              fyyim1 = fyy
              fzzim1 = fzz
              fxyim1 = fxy
              fyzim1 = fyz
              fzxim1 = fzx
            enddo
          enddo

c
c update at dateline
c
          i = im
          if(l <= lmom(i,j)) then
            am = am_im
            fm = fm_im
            fm0 = fm_im0
            fx = fx_im
            fy = fy_im
            fz = fz_im
            fxx = fxx_im
            fyy = fyy_im
            fzz = fzz_im
            fxy = fxy_im
            fyz = fyz_im
            fzx = fzx_im
            MMnew = MM(I,J,L) + (AMim1-AM)
            bymnew = 1d0/mmnew
            AM2 = AMim1+AM
            RM0 = RM(I,J,L) +  (FMim10-FM0)
            RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
c            RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
c     &           3d0*(AM2*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))*bymnew
c            RXX(I,J,L) = (RXX(I,J,L)*MM(I,J,L)**2
c     &           +2.5*RM(I,J,L)*(MM(I,J,L)**2 -MMnew**2 -3.*AM2*AM2)
c     &           +5.*(MM(I,J,L)*(MM(I,J,L)*(FMim1-FM)-FXim1-FX)
c     &           +AM2*RX(I,J,L)*MMnew)
c     &           +(FXXim1-FXX) ) * (bymnew*bymnew)
            RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &           3d0*(AM2*RM0-MM(I,J,L)*(FMim10+FM0)))*bymnew
            RXX(I,J,L) = (RXX(I,J,L)*MM(I,J,L)**2
     &           +2.5*RM0*(MM(I,J,L)**2 -MMnew**2 -3.*AM2*AM2)
     &           +5.*(MM(I,J,L)*(MM(I,J,L)*(FMim10-FM0)-FXim1-FX)
     &           +AM2*RX(I,J,L)*MMnew)
     &           +(FXXim1-FXX) ) * (bymnew*bymnew)
            RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
            RXY(I,J,L) = (RXY(I,J,L)*MM(I,J,L) + (FXYim1-FXY) +
     &           3d0*(AM2*RY(I,J,L)-MM(I,J,L)*(FYim1+FY)))*bymnew
            RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
            RZX(I,J,L) = (RZX(I,J,L)*MM(I,J,L) + (FZXim1-FZX) +
     &           3d0*(AM2*RZ(I,J,L)-MM(I,J,L)*(FZim1+FZ)))*bymnew
            RYY(I,J,L) = RYY(I,J,L) + (FYYim1-FYY)
            RZZ(I,J,L) = RZZ(I,J,L) + (FZZim1-FZZ)
            RYZ(I,J,L) = RYZ(I,J,L) + (FYZim1-FYZ)
            MM(I,J,L) = MMnew
            OIJL(I,J,L) = OIJL(I,J,L) + FM
! clean up roundoff errors
            if(rm(i,j,l).le.0d0) then
              rm(i,j,l)=0d0
              rx(i,j,l)=0d0
              ry(i,j,l)=0d0
              rz(i,j,l)=0d0
              rxx(i,j,l)=0d0
              ryy(i,j,l)=0d0
              rzz(i,j,l)=0d0
              rxy(i,j,l)=0d0
              ryz(i,j,l)=0d0
              rzx(i,j,l)=0d0
            endif
          endif
        enddo ! nc
      enddo ! j
      enddo ! l

      Return
      End Subroutine OADVTX4

      Subroutine OADVTY4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MO,MV,DT,OIJL)
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOV=>LMV
      USE OCEAN, only : nbyzm,i1yzm,i2yzm, nbyzv,i1yzv,i2yzv
      USE OCEANR_DIM, only : grid=>ogrid
      use DOMAIN_DECOMP_1D, only : halo_update, 
     &     hasSouthPole, hasNorthPole
      Implicit None
C**** Interface variables
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  MV
      Real*8,Intent(InOut),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX, MO,OIJL
C**** Local variables
      logical :: qnp
      Integer :: I,J,L, J1H,J1P,JNP,JN, N, imin,imax
      Real*8 :: BM,BM2,B,FM,FX,FY,FZ,FXX,FYY,FZZ,FXY,FYZ,FZX,Mnew,bymnew
      REAL*8 :: RM0,FM0,FM_PASS,FY_PASS,FYY_PASS
      real*8, dimension(im) :: BMjm1,FMjm1,FMjm10,
     &     FXjm1,FYjm1,FZjm1,FXXjm1,FYYjm1,FZZjm1,FXYjm1,FYZjm1,FZXjm1

C**** Extract domain decomposition band parameters
      J1H = Max (GRID%J_STRT-1, 1)   !  Exclude south pole
      J1P = Max (GRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (GRID%J_STOP, JM-1)  !  Exclude north pole
      JN = GRID%J_STOP
      QNP = hasNorthPole(grid)

      do l=1,lmo
      do j=j1p,jnp
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
        if(mv(i,j,l).gt.0. .and. mv(i,j-1,l).lt.0.) then
          call checkfluxo(mv(i,j-1,l)*dt,mv(i,j,l)*dt,mo(i,j,l),
     &         rm(i,j,l),ry(i,j,l),ryy(i,j,l))
        endif
      enddo
      enddo
      enddo
      enddo

      call halo_update(grid,mo)
      call halo_update(grid,rm)
      call halo_update(grid,rx)
      call halo_update(grid,ry)
      call halo_update(grid,rz)
      call halo_update(grid,rxx)
      call halo_update(grid,ryy)
      call halo_update(grid,rzz)
      call halo_update(grid,rxy)
      call halo_update(grid,ryz)
      call halo_update(grid,rzx)

      do l=1,lmo

        do i=1,im
          bmjm1(i) = 0.
          fmjm1(i) = 0.
          fmjm10(i) = 0.
          fxjm1(i) = 0.
          fyjm1(i) = 0.
          fzjm1(i) = 0.
          fxxjm1(i) = 0.
          fyyjm1(i) = 0.
          fzzjm1(i) = 0.
          fxyjm1(i) = 0.
          fyzjm1(i) = 0.
          fzxjm1(i) = 0.
        enddo

        if(qnp) then
          j = jm
          do i=2,im
            mo(i,j,l) = mo(1,j,l)
            rm(i,j,l) = rm(1,j,l)
            rz(i,j,l) = rz(1,j,l)
            rzz(i,j,l) = rzz(1,j,l)
          enddo
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
            rxx(i,j,l) = 0.
            ryy(i,j,l) = 0.
            rxy(i,j,l) = 0.
            ryz(i,j,l) = 0.
            rzx(i,j,l) = 0.
          enddo
        endif

        j = j1h
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            bm = mv(i,j,l)*dt
            if(bm.ge.0.) then   ! mass flux is positive or zero
              B  = BM / MO(I,J,L)
              FM = B*(RM(I,J,L)+(1-B)*(RY(I,J,L)+(1-2*B)*RYY(I,J,L)))
              FY  = BM*(B*B*(RY(I,J,L) + 3*(1-B)*RYY(I,J,L)) - 3*FM)
              FYY = BM*(BM*B**3*RYY(I,J,L) - 5*(BM*FM+FY))
              FX = B * (RX(i,j,L) + (1-B)*RXY(i,j,L))
              FXY = BM * (B*B*RXY(i,j,L) - 3*FX)
              FZ = B * (RZ(i,j,L) + (1-B)*RYZ(i,j,L))
              FYZ = BM * (B*B*RYZ(i,j,L) - 3*FZ)
              FXX = B*RXX(I,J,L)
              FZZ = B*RZZ(I,J,L)
              FZX = B*RZX(I,J,L)

              fm0 = fm
              if(fm.lt.0.) then
                fm=0.
                fy_pass=0
                fyy_pass=0.
              elseif(fm.gt.rm(i,j,l)) then
                fm=rm(i,j,l)
                fy_pass=bm*(-3.*fm)
                fyy_pass=bm*(-5.*(bm*fm+fy_pass))
              else
                fy_pass = fy
                fyy_pass = fyy
              endif
              fm_pass = fm

            else                ! mass flux is negative
              B  = BM / MO(i,j+1,L)
              FM = B*(RM(I,J+1,L) -(1+B)*(
     &             RY(I,J+1,L) -(1+2*B)*RYY(I,J+1,L)))
              FY  = BM*(B*B*(RY(I,J+1,L) -3*(1+B)*RYY(I,J+1,L)) -3*FM)
              FYY = BM*(BM*B**3*RYY(I,J+1,L) - 5*(BM*FM+FY))
              FX = B * (RX(i,j+1,L) - (1+B)*RXY(i,j+1,L))
              FXY = BM * (B*B*RXY(i,j+1,L) - 3*FX)
              FZ = B * (RZ(i,j+1,L) - (1+B)*RYZ(i,j+1,L))
              FYZ = BM * (B*B*RYZ(i,j+1,L) - 3*FZ)
              FXX = B*RXX(i,j+1,L)
              FZZ = B*RZZ(i,j+1,L)
              FZX = B*RZX(i,j+1,L)

              fm_pass = fm
              fy_pass = fy
              fyy_pass = fyy
              if(fm.gt.0.) then
                fm=0.
                fy=0.
                fyy=0.
              elseif(fm.lt.-rm(i,j+1,l)) then
                fm=-rm(i,j+1,l)
                fy=bm*(-3.*fm)
                fyy=bm*(-5.*(bm*fm+fy))
              endif
              fm0 = fm

            endif
            bmjm1(i) = bm
            fmjm1(i) = fm
            fmjm10(i) = fm_pass
            fxjm1(i) = fx
            fyjm1(i) = fy_pass
            fzjm1(i) = fz
            fxxjm1(i) = fxx
            fyyjm1(i) = fyy_pass
            fzzjm1(i) = fzz
            fxyjm1(i) = fxy
            fyzjm1(i) = fyz
            fzxjm1(i) = fzx
          enddo
        enddo

        do j=j1p,jn
          do n=1,nbyzm(j,l)
            imin=i1yzm(n,j,l)
            imax=i2yzm(n,j,l)
            if(j == jm) imax=im
            do i=imin,imax
              bm = mv(i,j,l)*dt
              if(bm.ge.0.) then ! mass flux is positive or zero
                B  = BM / MO(I,J,L)
                FM = B*(RM(I,J,L)+(1-B)*(RY(I,J,L)+(1-2*B)*RYY(I,J,L)))
                FY  = BM*(B*B*(RY(I,J,L) + 3*(1-B)*RYY(I,J,L)) - 3*FM)
                FYY = BM*(BM*B**3*RYY(I,J,L) - 5*(BM*FM+FY))
                FX = B * (RX(i,j,L) + (1-B)*RXY(i,j,L))
                FXY = BM * (B*B*RXY(i,j,L) - 3*FX)
                FZ = B * (RZ(i,j,L) + (1-B)*RYZ(i,j,L))
                FYZ = BM * (B*B*RYZ(i,j,L) - 3*FZ)
                FXX = B*RXX(I,J,L)
                FZZ = B*RZZ(I,J,L)
                FZX = B*RZX(I,J,L)

                fm0 = fm
                if(fm.lt.0.) then
                  fm=0.
                  fy_pass=0
                  fyy_pass=0.
                elseif(fm.gt.rm(i,j,l)) then
                  fm=rm(i,j,l)
                  fy_pass=bm*(-3.*fm)
                  fyy_pass=bm*(-5.*(bm*fm+fy_pass))
                else
                  fy_pass = fy
                  fyy_pass = fyy
                endif
                fm_pass = fm

              else              ! mass flux is negative
                B  = BM / MO(i,j+1,L)
                FM = B*(RM(I,J+1,L) -(1+B)*(
     &               RY(I,J+1,L) -(1+2*B)*RYY(I,J+1,L)))
                FY  = BM*(B*B*(RY(I,J+1,L) -3*(1+B)*RYY(I,J+1,L)) -3*FM)
                FYY = BM*(BM*B**3*RYY(I,J+1,L) - 5*(BM*FM+FY))
                FX = B * (RX(i,j+1,L) - (1+B)*RXY(i,j+1,L))
                FXY = BM * (B*B*RXY(i,j+1,L) - 3*FX)
                FZ = B * (RZ(i,j+1,L) - (1+B)*RYZ(i,j+1,L))
                FYZ = BM * (B*B*RYZ(i,j+1,L) - 3*FZ)
                FXX = B*RXX(i,j+1,L)
                FZZ = B*RZZ(i,j+1,L)
                FZX = B*RZX(i,j+1,L)

                fm_pass = fm
                fy_pass = fy
                fyy_pass = fyy
                if(fm.gt.0.) then
                  fm=0.
                  fy=0.
                  fyy=0.
                elseif(fm.lt.-rm(i,j+1,l)) then
                  fm=-rm(i,j+1,l)
                  fy=bm*(-3.*fm)
                  fyy=bm*(-5.*(bm*fm+fy))
                endif
                fm0 = fm

              endif
              Mnew = MO(I,J,L) + (BMjm1(i)-BM)
              bymnew = 1d0/mnew
              BM2 = BMjm1(i)+BM
              RM0 = RM(I,J,L) +  (FMjm10(i)-FM0)
              RM(I,J,L) = RM(I,J,L) +  (FMjm1(i)-FM)
c              RY(I,J,L) = (RY(I,J,L)*MO(I,J,L) + (FYjm1(i)-FY) + 3d0*
c     &             (BM2*RM(I,J,L)-MO(I,J,L)*(FMjm1(i)+FM)))*bymnew
c              RYY(I,J,L) = (RYY(I,J,L)*MO(I,J,L)**2
c     &             +2.5*RM(I,J,L)*(MO(I,J,L)**2 -Mnew**2 -3.*BM2*BM2)
c     &             +5.*(MO(I,J,L)*(MO(I,J,L)*(FMjm1(i)-FM)-FYjm1(i)-FY)
c     &             +BM2*RY(I,J,L)*Mnew)
c     &             +(FYYjm1(i)-FYY) ) * (bymnew*bymnew)
              RY(I,J,L) = (RY(I,J,L)*MO(I,J,L) + (FYjm1(i)-FY) + 3d0*
     &             (BM2*RM0-MO(I,J,L)*(FMjm10(i)+FM0)))*bymnew
              RYY(I,J,L) = (RYY(I,J,L)*MO(I,J,L)**2
     &             +2.5*RM0*(MO(I,J,L)**2 -Mnew**2 -3.*BM2*BM2)
     &            +5.*(MO(I,J,L)*(MO(I,J,L)*(FMjm10(i)-FM0)-FYjm1(i)-FY)
     &             +BM2*RY(I,J,L)*Mnew)
     &             +(FYYjm1(i)-FYY) ) * (bymnew*bymnew)
              RX(I,J,L) = RX(I,J,L) + (FXjm1(i)-FX)
              RXY(I,J,L) = (RXY(I,J,L)*MO(I,J,L) +(FXYjm1(i)-FXY) + 3d0*
     &             (BM2*RX(I,J,L)-MO(I,J,L)*(FXjm1(i)+FX)))*bymnew
              RZ(I,J,L) = RZ(I,J,L) + (FZjm1(i)-FZ)
              RYZ(I,J,L) = (RYZ(I,J,L)*MO(I,J,L) +(FYZjm1(i)-FYZ) + 3d0*
     &             (BM2*RZ(I,J,L)-MO(I,J,L)*(FZjm1(i)+FZ)))*bymnew
              RXX(I,J,L) = RXX(I,J,L) + (FXXjm1(i)-FXX)
              RZZ(I,J,L) = RZZ(I,J,L) + (FZZjm1(i)-FZZ)
              RZX(I,J,L) = RZX(I,J,L) + (FZXjm1(i)-FZX)
              MO(I,J,L) = Mnew

! clean up roundoff errors
              if(rm(i,j,l).le.0d0) then
                rm(i,j,l)=0d0
                rx(i,j,l)=0d0
                ry(i,j,l)=0d0
                rz(i,j,l)=0d0
                rxx(i,j,l)=0d0
                ryy(i,j,l)=0d0
                rzz(i,j,l)=0d0
                rxy(i,j,l)=0d0
                ryz(i,j,l)=0d0
                rzx(i,j,l)=0d0
              endif

              OIJL(I,J,L) = OIJL(I,J,L) + FM
              bmjm1(i) = bm
              fmjm1(i) = fm
              fmjm10(i) = fm_pass
              fxjm1(i) = fx
              fyjm1(i) = fy_pass
              fzjm1(i) = fz
              fxxjm1(i) = fxx
              fyyjm1(i) = fyy_pass
              fzzjm1(i) = fzz
              fxyjm1(i) = fxy
              fyzjm1(i) = fyz
              fzxjm1(i) = fzx
            enddo
          enddo
        enddo

c
c average the pole
c
        if(qnp) then
        if(l<=lmom(1,jm)) then
          j = jm
          mo(:,j,l) = sum(mo(:,j,l))/im
          rm(:,j,l) = sum(rm(:,j,l))/im
          rz(:,j,l) = sum(rz(:,j,l))/im
          rzz(:,j,l) = sum(rzz(:,j,l))/im
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
            rxx(i,j,l) = 0.
            ryy(i,j,l) = 0.
            rxy(i,j,l) = 0.
            ryz(i,j,l) = 0.
            rzx(i,j,l) = 0.
          enddo
        endif
        endif

      enddo ! l

      return
      end subroutine oadvty4

      SUBROUTINE OADVTZ4(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MO,MW,DT,OIJL)
C****
      USE CONSTANT, only : by3,by12
      USE OCEAN, only : im,jm,lmo,lmm
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      use domain_decomp_1d, only : getDomainBounds,halo_update
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX, OIJL, MO
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     MW
      REAL*8, INTENT(IN) :: DT
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     CMUP,FMUP,FXUP,FYUP,FZUP,FXXUP,FYYUP,FZZUP,FXYUP,FYZUP,FZXUP,
     &     FMUP0
      REAL*8 :: CM,C,CM2,FM,FX,FY,FZ,FXX,FYY,FZZ,FXY,FYZ,FZX,MNEW,BYMNEW
      REAL*8 :: RM0,FM0,FM_PASS,FZ_PASS,FZZ_PASS
      INTEGER :: I,J,L,N

      INTEGER :: J_0,J_1,J_0S,J_1S

      logical :: done,mw_sv_alloc
      integer :: iter
      real*8, dimension(:,:,:), allocatable :: mw_sv

! Notes on the treatment of courant numbers c exceeding unity
! (usually due to a long timestep):
! 1. The strategy assumes that this occurs rarely, and thus
!    avoids any pre-checks and other operations that would
!    slow execution for representative circumstances.
! 2. Within the transport loop, c is checked. Where it exceeds unity,
!    a mass flux corresponding to c=1 is used, and the remainder
!    is applied the following substep (via the do-while loop).
! 3. This strategy has not yet been verified to work well where c>=2.
!    More work will be done to check the overall realism of such situations.

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)

      iter = 0
      done = .false.
      mw_sv_alloc = .false.

      do while(.not.done)

      iter = iter + 1
      done = .true.

      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            cmup(i,j) = 0.
            fmup(i,j) = 0.
            fmup0(i,j) = 0.
            fxup(i,j) = 0.
            fyup(i,j) = 0.
            fzup(i,j) = 0.
            fxxup(i,j) = 0.
            fyyup(i,j) = 0.
            fzzup(i,j) = 0.
            fxyup(i,j) = 0.
            fyzup(i,j) = 0.
            fzxup(i,j) = 0.
          enddo
        enddo
      enddo

      do l=2,lmo-1
      do j=j_0,j_1
      do n=1,nbyzm(j,l+1)
      do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
        if(mw(i,j,l).gt.0. .and. mw(i,j,l-1).lt.0.) then
          call checkfluxo(mw(i,j,l-1)*dt,mw(i,j,l)*dt,mo(i,j,l),
     &         rm(i,j,l),rz(i,j,l),rzz(i,j,l))
        endif
      enddo
      enddo
      enddo
      enddo

      do l=1,lmo
        do j=j_0,j_1
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              CM = DT*MW(I,J,L)
              if(cm.ge.0.) then ! mass flux is downward or zero
                C  = CM/MO(I,J,L)
                if(c.gt.1d0) then ! see notes above for this caes
                  write(6,*) 'c>1:',i,j,l,c,mo(i,j,l)
                  done = .false.
                  if(.not.mw_sv_alloc) then
                    allocate(mw_sv(im,j_0:j_1,lmo))
                    mw_sv(:,j_0:j_1,:) = mw(:,j_0:j_1,:)
                    mw_sv_alloc = .true.
                  endif
                  cm = mo(i,j,l)
                  c  = 1d0
                  mw(i,j,l) = mw(i,j,l) - cm/dt ! remainder for next substep
                endif
                FM = C*(RM(I,J,L)+(1-C)*(RZ(I,J,L)+(1-2*C)*RZZ(I,J,L)))
                FZ  = CM*(C*C*(RZ(I,J,L) + 3*(1-C)*RZZ(I,J,L)) - 3*FM)
                FZZ = CM*(CM*C**3*RZZ(I,J,L) - 5*(CM*FM+FZ))
                FX = C*(RX(I,J,L)+(1d0-C)*RZX(I,J,L))
                FZX = CM*(C*C*RZX(I,J,L)-3d0*FX)
                FY = C*(RY(I,J,L)+(1d0-C)*RYZ(I,J,L))
                FYZ = CM*(C*C*RYZ(I,J,L)-3d0*FY)
                FXX = C*RXX(I,J,L)
                FYY = C*RYY(I,J,L)
                FXY = C*RXY(I,J,L)

                fm0 = fm
                if(fm.lt.0.) then
                  fm=0.
                  fz_pass=0
                  fzz_pass=0.
                elseif(fm.gt.rm(i,j,l)) then
                  fm=rm(i,j,l)
                  fz_pass=cm*(-3.*fm)
                  fzz_pass=cm*(-5.*(cm*fm+fz_pass))
                else
                  fz_pass = fz
                  fzz_pass = fzz
                endif
                fm_pass = fm

              else              ! mass flux is upward
                C  = CM/MO(I,J,L+1)
                if(c.lt.-1d0) then ! see notes above for this caes
                  write(6,*) 'c<-1:',i,j,l,c,mo(i,j,l+1)
                  done = .false.
                  if(.not.mw_sv_alloc) then
                    allocate(mw_sv(im,j_0:j_1,lmo))
                    mw_sv(:,j_0:j_1,:) = mw(:,j_0:j_1,:)
                    mw_sv_alloc = .true.
                  endif
                  cm = -mo(i,j,l+1)
                  c  = -1d0
                  mw(i,j,l) = mw(i,j,l) - cm/dt ! remainder for next substep
                endif
                FM = C*(RM(I,J,L+1) -(1+C)*(
     &               RZ(I,J,L+1) -(1+2*C)*RZZ(I,J,L+1)))
                FZ  = CM*(C*C*(RZ(I,J,L+1) -3*(1+C)*RZZ(I,J,L+1)) -3*FM)
                FZZ = CM*(CM*C**3*RZZ(I,J,L+1) - 5*(CM*FM+FZ))
                FX = C*(RX(I,J,L+1)-(1d0+C)*RZX(I,J,L+1))
                FZX = CM*(C*C*RZX(I,J,L+1)-3d0*FX)
                FY = C*(RY(I,J,L+1)-(1d0+C)*RYZ(I,J,L+1))
                FYZ = CM*(C*C*RYZ(I,J,L+1)-3d0*FY)
                FXX = C*RXX(I,J,L+1)
                FYY = C*RYY(I,J,L+1)
                FXY = C*RXY(I,J,L+1)

                fm_pass = fm
                fz_pass = fz
                fzz_pass = fzz
                if(fm.gt.0.) then
                  fm = 0.
                  fz = 0.
                  fzz = 0.
                elseif(fm.lt.-rm(i,j,l+1)) then
                  fm=-rm(i,j,l+1)
                  fz=cm*(-3.*fm)
                  fzz=cm*(-5.*(cm*fm+fz))
                endif
                fm0 = fm

              endif
              mnew = MO(I,J,L) + (CMUP(I,J)-CM)
              bymnew = 1d0/mnew
              CM2 = CMUP(I,J)+CM
              RM0 = RM(I,J,L) + (FMUP0(I,J)-FM0)
              RM(I,J,L) = RM(I,J,L) + (FMUP(I,J)-FM)
c              RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZUP(I,J)-FZ) +3d0*
c     &             (CM2*RM(I,J,L)-MO(I,J,L)*(FMUP(I,J)+FM)))*bymnew
c              RZZ(I,J,L) = (RZZ(I,J,L)*MO(I,J,L)**2
c     &            +2.5*RM(I,J,L)*(MO(I,J,L)**2 -Mnew**2 -3.*CM2*CM2)
c     &            +5.*(MO(I,J,L)*(MO(I,J,L)*(FMUP(i,j)-FM)-FZUP(i,j)-FZ)
c     &            +CM2*RZ(I,J,L)*Mnew)
c     &            +(FZZUP(i,j)-FZZ) ) * (bymnew*bymnew)
              RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZUP(I,J)-FZ) +3d0*
     &             (CM2*RM0-MO(I,J,L)*(FMUP0(I,J)+FM0)))*bymnew
              RZZ(I,J,L) = (RZZ(I,J,L)*MO(I,J,L)**2
     &            +2.5*RM0*(MO(I,J,L)**2 -Mnew**2 -3.*CM2*CM2)
     &          +5.*(MO(I,J,L)*(MO(I,J,L)*(FMUP0(i,j)-FM0)-FZUP(i,j)-FZ)
     &            +CM2*RZ(I,J,L)*Mnew)
     &            +(FZZUP(i,j)-FZZ) ) * (bymnew*bymnew)
              RX(I,J,L) = RX(I,J,L) + (FXUP(I,J)-FX)
              RZX(I,J,L) = (RZX(I,J,L)*MO(I,J,L) +(FZXUP(I,J)-FZX) +3d0*
     &             (CM2*RX(I,J,L)-MO(I,J,L)*(FXUP(I,J)+FX)))*bymnew
              RY(I,J,L) = RY(I,J,L) + (FYUP(I,J)-FY)
              RYZ(I,J,L) = (RYZ(I,J,L)*MO(I,J,L) +(FYZUP(I,J)-FYZ) +3d0*
     &             (CM2*RY(I,J,L)-MO(I,J,L)*(FYUP(I,J)+FY)))*bymnew
              RXX(I,J,L) = RXX(I,J,L) + (FXXUP(I,J)-FXX)
              RYY(I,J,L) = RYY(I,J,L) + (FYYUP(I,J)-FYY)
              RXY(I,J,L) = RXY(I,J,L) + (FXYUP(I,J)-FXY)
c
              MO(I,J,L) = mnew

! clean up roundoff errors
              if(rm(i,j,l).le.0d0) then
                rm(i,j,l)=0d0
                rx(i,j,l)=0d0
                ry(i,j,l)=0d0
                rz(i,j,l)=0d0
                rxx(i,j,l)=0d0
                ryy(i,j,l)=0d0
                rzz(i,j,l)=0d0
                rxy(i,j,l)=0d0
                ryz(i,j,l)=0d0
                rzx(i,j,l)=0d0
              endif

              cmup(i,j) = cm
              fmup(i,j) = fm
              fmup0(i,j) = fm_pass
              fxup(i,j) = fx
              fyup(i,j) = fy
              fzup(i,j) = fz_pass
              fxxup(i,j) = fxx
              fyyup(i,j) = fyy
              fzzup(i,j) = fzz_pass
              fxyup(i,j) = fxy
              fyzup(i,j) = fyz
              fzxup(i,j) = fzx
              OIJL(I,J,L) = OIJL(I,J,L) + FM
            enddo
          enddo
        enddo
      enddo

      if(mw_sv_alloc .and. iter.eq.1) then ! zero out the |c|<=1 mw
        do l=1,lmo
          do j=j_0,j_1
            do n=1,nbyzm(j,l)
              do i=i1yzm(n,j,l),i2yzm(n,j,l)
                if(mw(i,j,l).ne.mw_sv(i,j,l)) cycle ! |c|>1 point
                mw(i,j,l) = 0.
              enddo
            enddo
          enddo
        enddo
      endif

      enddo ! iterations over substeps

      if(mw_sv_alloc) then ! restore input mw
        mw(:,j_0:j_1,:) = mw_sv(:,j_0:j_1,:)
      endif

      return
      END SUBROUTINE OADVTZ4

      subroutine checkfluxo(aml,amr,m,rm,rxm,rxxm)
      implicit none
      real*8 :: aml,amr,m,rm,rxm,rxxm
      real*8 :: a,fl,fr
c flux out the right side
      A = AMR / M
      FR = A*(RM + (1.-A)*(RXM + (1.-2.*A)*RXXM))
c flux out the left side
      A = AML / M
      FL = A*(RM - (1.+A)*(RXM - (1.+2.*A)*RXXM))
c
      if(rm+fl-fr.le.0.) then
        rxm = 0.
        rxxm = 0.
      endif
      return
      end subroutine checkfluxo

      Subroutine OADVTX3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MM,MU,DT,OIJL)
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOU=>LMU
      USE OCEAN, only : nbyzu,i1yzu,i2yzu, nbyzm,i1yzm,i2yzm
      USE OCEANR_DIM, only : ogrid
      Implicit None
C**** Interface variables
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  MU
      Real*8,Intent(InOut),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX, MM,OIJL
C**** Local variables
      Integer :: I,J,L, J1P,JNP, N,NC,NCOURANT
      real*8, dimension(im) :: mudt
      Real*8 :: AM,AM2,A,FM,FX,FY,FZ,FXX,FYY,FZZ,FXY,FYZ,FZX,
     &     MMnew,bymnew, zCOURANT,mcheck,courmax
      Real*8 :: AMim1,FMim1,AM_im,FM_im,
     &     FXim1,FYim1,FZim1,FXXim1,FYYim1,FZZim1,FXYim1,FYZim1,FZXim1,
     &     FX_im,FY_im,FZ_im,FXX_im,FYY_im,FZZ_im,FXY_im,FYZ_im,FZX_im

C**** Extract domain decomposition band parameters
      J1P = Max (oGRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (oGRID%J_STOP, JM-1)  !  Exclude north pole

C****
C**** Loop over layers and latitudes
C****
      Do L=1,LMO
      Do J=J1P,JNP
        if(nbyzu(j,l) == 0) cycle


c
c adjust mass fluxes so that courant numbers are never greater than 1
c
        courmax = 0.
        mudt(1:2) = mu(1:2,j,l)*dt
        mudt(im) = mu(im,j,l)*dt
        i=1
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(im)-mudt(i))
          else
            mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        do n=1,nbyzu(j,l)
          i = i1yzu(n,j,l)
          if(i.gt.1) mudt(i-1) = 0. ! zero mudt at west boundary
          mudt(i) = mu(i,j,l)*dt
          do i=max(2,i1yzu(n,j,l)),min(i2yzu(n,j,l),im-1)
            mudt(i+1) = mu(i+1,j,l)*dt
            if(mudt(i).ge.0.) then
              mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
            else
              mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
            endif
            courmax = max(courmax,mudt(i)/mcheck)
          enddo
        enddo
        i = im
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
          else
            mcheck = -mm(1,j,l)-min(0d0,mudt(i)-mudt(1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        if(courmax.gt.1.) then
          ncourant = 1+int(courmax)
          write(6,*) 'ncourant ',j,l,ncourant
          zcourant = 1d0/ncourant
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              mudt(i) = mudt(i)*zcourant
            enddo
          enddo
        else
          ncourant = 1
        endif

        do nc=1,ncourant

c
c calculate fluxes at the dateline
c
          i = im
          if(l <= lmou(i,j)) then
            AM = mudt(I)
            if(am.ge.0.) then   ! mass flux is positive or zero
              A  = AM / MM(I,J,L)
              FM = A*(RM(I,J,L)+(1-A)*(RX(I,J,L)+(1-2*A)*RXX(I,J,L)))
              FX  = AM*(A*A*(RX(I,J,L) + 3*(1-A)*RXX(I,J,L)) - 3*FM)
              FXX = AM*(AM*A**3*RXX(I,J,L) - 5*(AM*FM+FX))
              FY = A * (RY(I,J,L) + (1-A)*RXY(I,J,L))
              FXY = AM * (A*A*RXY(I,J,L) - 3*FY)
              FZ = A * (RZ(I,J,L) + (1-A)*RZX(I,J,L))
              FZX = AM * (A*A*RZX(I,J,L) - 3*FZ)
              FYY = A*RYY(I,J,L)
              FZZ = A*RZZ(I,J,L)
              FYZ = A*RYZ(I,J,L)

            else                ! mass flux is negative
              A  = AM / MM(1,J,L)
              FM = A*(RM(1,J,L) -(1+A)*(
     &             RX(1,J,L) -(1+2*A)*RXX(1,J,L)))
              FX  = AM*(A*A*(RX(1,J,L) -3*(1+A)*RXX(1,J,L)) -3*FM)
              FXX = AM*(AM*A**3*RXX(1,J,L) - 5*(AM*FM+FX))
              FY = A * (RY(1,J,L) - (1+A)*RXY(1,J,L))
              FXY = AM * (A*A*RXY(1,J,L) - 3*FY)
              FZ = A * (RZ(1,J,L) - (1+A)*RZX(1,J,L))
              FZX = AM * (A*A*RZX(1,J,L) - 3*FZ)
              FYY = A*RYY(1,J,L)
              FZZ = A*RZZ(1,J,L)
              FYZ = A*RYZ(1,J,L)

            endif
          else
            am = 0.
            fm = 0.
            fx = 0.
            fy = 0.
            fz = 0.
            fxx = 0.
            fyy = 0.
            fzz = 0.
            fxy = 0.
            fyz = 0.
            fzx = 0.
          endif
          amim1 = am
          fmim1 = fm
          fxim1 = fx
          fyim1 = fy
          fzim1 = fz
          fxxim1 = fxx
          fyyim1 = fyy
          fzzim1 = fzz
          fxyim1 = fxy
          fyzim1 = fyz
          fzxim1 = fzx
          am_im = am
          fm_im = fm
          fx_im = fx
          fy_im = fy
          fz_im = fz
          fxx_im = fxx
          fyy_im = fyy
          fzz_im = fzz
          fxy_im = fxy
          fyz_im = fyz
          fzx_im = fzx

c
c loop over basins
c
          do n=1,nbyzm(j,l)
            if(i1yzm(n,j,l)>1 .and. i1yzm(n,j,l)==i2yzm(n,j,l)) cycle
            do i=i1yzm(n,j,l),min(i2yzm(n,j,l),im-1)
              AM = mudt(I)
              if(am.ge.0.) then ! mass flux is positive or zero
                A  = AM / MM(I,J,L)
                FM = A*(RM(I,J,L)+(1-A)*(RX(I,J,L)+(1-2*A)*RXX(I,J,L)))
                FX  = AM*(A*A*(RX(I,J,L) + 3*(1-A)*RXX(I,J,L)) - 3*FM)
                FXX = AM*(AM*A**3*RXX(I,J,L) - 5*(AM*FM+FX))
                FY = A * (RY(I,J,L) + (1-A)*RXY(I,J,L))
                FXY = AM * (A*A*RXY(I,J,L) - 3*FY)
                FZ = A * (RZ(I,J,L) + (1-A)*RZX(I,J,L))
                FZX = AM * (A*A*RZX(I,J,L) - 3*FZ)
                FYY = A*RYY(I,J,L)
                FZZ = A*RZZ(I,J,L)
                FYZ = A*RYZ(I,J,L)

              else              ! mass flux is negative
                A  = AM / MM(I+1,J,L)
                FM = A*(RM(I+1,J,L) -(1+A)*(
     &               RX(I+1,J,L) -(1+2*A)*RXX(I+1,J,L)))
                FX  = AM*(A*A*(RX(I+1,J,L) -3*(1+A)*RXX(I+1,J,L)) -3*FM)
                FXX = AM*(AM*A**3*RXX(I+1,J,L) - 5*(AM*FM+FX))
                FY = A * (RY(I+1,J,L) - (1+A)*RXY(I+1,J,L))
                FXY = AM * (A*A*RXY(I+1,J,L) - 3*FY)
                FZ = A * (RZ(I+1,J,L) - (1+A)*RZX(I+1,J,L))
                FZX = AM * (A*A*RZX(I+1,J,L) - 3*FZ)
                FYY = A*RYY(I+1,J,L)
                FZZ = A*RZZ(I+1,J,L)
                FYZ = A*RYZ(I+1,J,L)

              endif
              MMnew = MM(I,J,L) + (AMim1-AM)
              bymnew = 1d0/mmnew
              AM2 = AMim1+AM
              RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
              RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &             3d0*(AM2*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))*bymnew
              RXX(I,J,L) = (RXX(I,J,L)*MM(I,J,L)**2
     &             +2.5*RM(I,J,L)*(MM(I,J,L)**2 -MMnew**2 -3.*AM2*AM2)
     &             +5.*(MM(I,J,L)*(MM(I,J,L)*(FMim1-FM)-FXim1-FX)
     &             +AM2*RX(I,J,L)*MMnew)
     &             +(FXXim1-FXX) ) * (bymnew*bymnew)
              RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
              RXY(I,J,L) = (RXY(I,J,L)*MM(I,J,L) + (FXYim1-FXY) +
     &             3d0*(AM2*RY(I,J,L)-MM(I,J,L)*(FYim1+FY)))*bymnew
              RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
              RZX(I,J,L) = (RZX(I,J,L)*MM(I,J,L) + (FZXim1-FZX) +
     &             3d0*(AM2*RZ(I,J,L)-MM(I,J,L)*(FZim1+FZ)))*bymnew
              RYY(I,J,L) = RYY(I,J,L) + (FYYim1-FYY)
              RZZ(I,J,L) = RZZ(I,J,L) + (FZZim1-FZZ)
              RYZ(I,J,L) = RYZ(I,J,L) + (FYZim1-FYZ)
              MM(I,J,L) = MMnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
              amim1 = am
              fmim1 = fm
              fxim1 = fx
              fyim1 = fy
              fzim1 = fz
              fxxim1 = fxx
              fyyim1 = fyy
              fzzim1 = fzz
              fxyim1 = fxy
              fyzim1 = fyz
              fzxim1 = fzx
            enddo
          enddo

c
c update at dateline
c
          i = im
          if(l <= lmom(i,j)) then
            am = am_im
            fm = fm_im
            fx = fx_im
            fy = fy_im
            fz = fz_im
            fxx = fxx_im
            fyy = fyy_im
            fzz = fzz_im
            fxy = fxy_im
            fyz = fyz_im
            fzx = fzx_im
            MMnew = MM(I,J,L) + (AMim1-AM)
            bymnew = 1d0/mmnew
            AM2 = AMim1+AM
            RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
            RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &           3d0*(AM2*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))*bymnew
            RXX(I,J,L) = (RXX(I,J,L)*MM(I,J,L)**2
     &           +2.5*RM(I,J,L)*(MM(I,J,L)**2 -MMnew**2 -3.*AM2*AM2)
     &           +5.*(MM(I,J,L)*(MM(I,J,L)*(FMim1-FM)-FXim1-FX)
     &           +AM2*RX(I,J,L)*MMnew)
     &           +(FXXim1-FXX) ) * (bymnew*bymnew)
            RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
            RXY(I,J,L) = (RXY(I,J,L)*MM(I,J,L) + (FXYim1-FXY) +
     &           3d0*(AM2*RY(I,J,L)-MM(I,J,L)*(FYim1+FY)))*bymnew
            RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
            RZX(I,J,L) = (RZX(I,J,L)*MM(I,J,L) + (FZXim1-FZX) +
     &           3d0*(AM2*RZ(I,J,L)-MM(I,J,L)*(FZim1+FZ)))*bymnew
            RYY(I,J,L) = RYY(I,J,L) + (FYYim1-FYY)
            RZZ(I,J,L) = RZZ(I,J,L) + (FZZim1-FZZ)
            RYZ(I,J,L) = RYZ(I,J,L) + (FYZim1-FYZ)
            MM(I,J,L) = MMnew
            OIJL(I,J,L) = OIJL(I,J,L) + FM
          endif
        enddo ! nc
      enddo ! j
      enddo ! l

      Return
      End Subroutine OADVTX3

      Subroutine OADVTY3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MO,MV,DT,OIJL)
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOV=>LMV
      USE OCEAN, only : nbyzm,i1yzm,i2yzm, nbyzv,i1yzv,i2yzv
      USE OCEANR_DIM, only : grid=>ogrid
      use DOMAIN_DECOMP_1D, only : halo_update, 
     &     hasSouthPole, hasNorthPole
      Implicit None
C**** Interface variables
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  MV
      Real*8,Intent(InOut),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX, MO,OIJL
C**** Local variables
      logical :: qnp
      Integer :: I,J,L, J1H,J1P,JNP,JN, N, imin,imax
      Real*8 :: BM,BM2,B,FM,FX,FY,FZ,FXX,FYY,FZZ,FXY,FYZ,FZX,Mnew,bymnew
      real*8, dimension(im) :: BMjm1,FMjm1,
     &     FXjm1,FYjm1,FZjm1,FXXjm1,FYYjm1,FZZjm1,FXYjm1,FYZjm1,FZXjm1

C**** Extract domain decomposition band parameters
      J1H = Max (GRID%J_STRT-1, 1)   !  Exclude south pole
      J1P = Max (GRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (GRID%J_STOP, JM-1)  !  Exclude north pole
      JN = GRID%J_STOP
      QNP = hasNorthPole(grid)

      call halo_update(grid,mo)
      call halo_update(grid,rm)
      call halo_update(grid,rx)
      call halo_update(grid,ry)
      call halo_update(grid,rz)
      call halo_update(grid,rxx)
      call halo_update(grid,ryy)
      call halo_update(grid,rzz)
      call halo_update(grid,rxy)
      call halo_update(grid,ryz)
      call halo_update(grid,rzx)

      do l=1,lmo

        do i=1,im
          bmjm1(i) = 0.
          fmjm1(i) = 0.
          fxjm1(i) = 0.
          fyjm1(i) = 0.
          fzjm1(i) = 0.
          fxxjm1(i) = 0.
          fyyjm1(i) = 0.
          fzzjm1(i) = 0.
          fxyjm1(i) = 0.
          fyzjm1(i) = 0.
          fzxjm1(i) = 0.
        enddo

        if(qnp) then
          j = jm
          do i=2,im
            mo(i,j,l) = mo(1,j,l)
            rm(i,j,l) = rm(1,j,l)
            rz(i,j,l) = rz(1,j,l)
            rzz(i,j,l) = rzz(1,j,l)
          enddo
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
            rxx(i,j,l) = 0.
            ryy(i,j,l) = 0.
            rxy(i,j,l) = 0.
            ryz(i,j,l) = 0.
            rzx(i,j,l) = 0.
          enddo
        endif

        j = j1h
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            bm = mv(i,j,l)*dt
            if(bm.ge.0.) then   ! mass flux is positive or zero
              B  = BM / MO(I,J,L)
              FM = B*(RM(I,J,L)+(1-B)*(RY(I,J,L)+(1-2*B)*RYY(I,J,L)))
              FY  = BM*(B*B*(RY(I,J,L) + 3*(1-B)*RYY(I,J,L)) - 3*FM)
              FYY = BM*(BM*B**3*RYY(I,J,L) - 5*(BM*FM+FY))
              FX = B * (RX(i,j,L) + (1-B)*RXY(i,j,L))
              FXY = BM * (B*B*RXY(i,j,L) - 3*FX)
              FZ = B * (RZ(i,j,L) + (1-B)*RYZ(i,j,L))
              FYZ = BM * (B*B*RYZ(i,j,L) - 3*FZ)
              FXX = B*RXX(I,J,L)
              FZZ = B*RZZ(I,J,L)
              FZX = B*RZX(I,J,L)
            else                ! mass flux is negative
              B  = BM / MO(i,j+1,L)
              FM = B*(RM(I,J+1,L) -(1+B)*(
     &             RY(I,J+1,L) -(1+2*B)*RYY(I,J+1,L)))
              FY  = BM*(B*B*(RY(I,J+1,L) -3*(1+B)*RYY(I,J+1,L)) -3*FM)
              FYY = BM*(BM*B**3*RYY(I,J+1,L) - 5*(BM*FM+FY))
              FX = B * (RX(i,j+1,L) - (1+B)*RXY(i,j+1,L))
              FXY = BM * (B*B*RXY(i,j+1,L) - 3*FX)
              FZ = B * (RZ(i,j+1,L) - (1+B)*RYZ(i,j+1,L))
              FYZ = BM * (B*B*RYZ(i,j+1,L) - 3*FZ)
              FXX = B*RXX(i,j+1,L)
              FZZ = B*RZZ(i,j+1,L)
              FZX = B*RZX(i,j+1,L)
            endif
            bmjm1(i) = bm
            fmjm1(i) = fm
            fxjm1(i) = fx
            fyjm1(i) = fy
            fzjm1(i) = fz
            fxxjm1(i) = fxx
            fyyjm1(i) = fyy
            fzzjm1(i) = fzz
            fxyjm1(i) = fxy
            fyzjm1(i) = fyz
            fzxjm1(i) = fzx
          enddo
        enddo

        do j=j1p,jn
          do n=1,nbyzm(j,l)
            imin=i1yzm(n,j,l)
            imax=i2yzm(n,j,l)
            if(j == jm) imax=im
            do i=imin,imax
              bm = mv(i,j,l)*dt
              if(bm.ge.0.) then ! mass flux is positive or zero
                B  = BM / MO(I,J,L)
                FM = B*(RM(I,J,L)+(1-B)*(RY(I,J,L)+(1-2*B)*RYY(I,J,L)))
                FY  = BM*(B*B*(RY(I,J,L) + 3*(1-B)*RYY(I,J,L)) - 3*FM)
                FYY = BM*(BM*B**3*RYY(I,J,L) - 5*(BM*FM+FY))
                FX = B * (RX(i,j,L) + (1-B)*RXY(i,j,L))
                FXY = BM * (B*B*RXY(i,j,L) - 3*FX)
                FZ = B * (RZ(i,j,L) + (1-B)*RYZ(i,j,L))
                FYZ = BM * (B*B*RYZ(i,j,L) - 3*FZ)
                FXX = B*RXX(I,J,L)
                FZZ = B*RZZ(I,J,L)
                FZX = B*RZX(I,J,L)
              else              ! mass flux is negative
                B  = BM / MO(i,j+1,L)
                FM = B*(RM(I,J+1,L) -(1+B)*(
     &               RY(I,J+1,L) -(1+2*B)*RYY(I,J+1,L)))
                FY  = BM*(B*B*(RY(I,J+1,L) -3*(1+B)*RYY(I,J+1,L)) -3*FM)
                FYY = BM*(BM*B**3*RYY(I,J+1,L) - 5*(BM*FM+FY))
                FX = B * (RX(i,j+1,L) - (1+B)*RXY(i,j+1,L))
                FXY = BM * (B*B*RXY(i,j+1,L) - 3*FX)
                FZ = B * (RZ(i,j+1,L) - (1+B)*RYZ(i,j+1,L))
                FYZ = BM * (B*B*RYZ(i,j+1,L) - 3*FZ)
                FXX = B*RXX(i,j+1,L)
                FZZ = B*RZZ(i,j+1,L)
                FZX = B*RZX(i,j+1,L)
              endif
              Mnew = MO(I,J,L) + (BMjm1(i)-BM)
              bymnew = 1d0/mnew
              BM2 = BMjm1(i)+BM
              RM(I,J,L) = RM(I,J,L) +  (FMjm1(i)-FM)
              RY(I,J,L) = (RY(I,J,L)*MO(I,J,L) + (FYjm1(i)-FY) + 3d0*
     &             (BM2*RM(I,J,L)-MO(I,J,L)*(FMjm1(i)+FM)))*bymnew
              RYY(I,J,L) = (RYY(I,J,L)*MO(I,J,L)**2
     &             +2.5*RM(I,J,L)*(MO(I,J,L)**2 -Mnew**2 -3.*BM2*BM2)
     &             +5.*(MO(I,J,L)*(MO(I,J,L)*(FMjm1(i)-FM)-FYjm1(i)-FY)
     &             +BM2*RY(I,J,L)*Mnew)
     &             +(FYYjm1(i)-FYY) ) * (bymnew*bymnew)
              RX(I,J,L) = RX(I,J,L) + (FXjm1(i)-FX)
              RXY(I,J,L) = (RXY(I,J,L)*MO(I,J,L) +(FXYjm1(i)-FXY) + 3d0*
     &             (BM2*RX(I,J,L)-MO(I,J,L)*(FXjm1(i)+FX)))*bymnew
              RZ(I,J,L) = RZ(I,J,L) + (FZjm1(i)-FZ)
              RYZ(I,J,L) = (RYZ(I,J,L)*MO(I,J,L) +(FYZjm1(i)-FYZ) + 3d0*
     &             (BM2*RZ(I,J,L)-MO(I,J,L)*(FZjm1(i)+FZ)))*bymnew
              RXX(I,J,L) = RXX(I,J,L) + (FXXjm1(i)-FXX)
              RZZ(I,J,L) = RZZ(I,J,L) + (FZZjm1(i)-FZZ)
              RZX(I,J,L) = RZX(I,J,L) + (FZXjm1(i)-FZX)
              MO(I,J,L) = Mnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
              bmjm1(i) = bm
              fmjm1(i) = fm
              fxjm1(i) = fx
              fyjm1(i) = fy
              fzjm1(i) = fz
              fxxjm1(i) = fxx
              fyyjm1(i) = fyy
              fzzjm1(i) = fzz
              fxyjm1(i) = fxy
              fyzjm1(i) = fyz
              fzxjm1(i) = fzx
            enddo
          enddo
        enddo

c
c average the pole
c
        if(qnp) then
        if(l<=lmom(1,jm)) then
          j = jm
          mo(:,j,l) = sum(mo(:,j,l))/im
          rm(:,j,l) = sum(rm(:,j,l))/im
          rz(:,j,l) = sum(rz(:,j,l))/im
          rzz(:,j,l) = sum(rzz(:,j,l))/im
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
            rxx(i,j,l) = 0.
            ryy(i,j,l) = 0.
            rxy(i,j,l) = 0.
            ryz(i,j,l) = 0.
            rzx(i,j,l) = 0.
          enddo
        endif
        endif

      enddo ! l

      return
      end subroutine oadvty3

      SUBROUTINE OADVTZ3(RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX,
     &     MO,MW,DT,OIJL)
C****
      USE CONSTANT, only : by3,by12
      USE OCEAN, only : im,jm,lmo,lmm
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      use domain_decomp_1d, only : getDomainBounds,halo_update
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ,RXX,RYY,RZZ,RXY,RYZ,RZX, OIJL, MO
      REAL*8, INTENT(IN),
     &  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     MW
      REAL*8, INTENT(IN) :: DT
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     CMUP,FMUP,FXUP,FYUP,FZUP,FXXUP,FYYUP,FZZUP,FXYUP,FYZUP,FZXUP
      REAL*8 :: CM,C,CM2,FM,FX,FY,FZ,FXX,FYY,FZZ,FXY,FYZ,FZX,MNEW,BYMNEW
      INTEGER :: I,J,L,N

      INTEGER :: J_0,J_1,J_0S,J_1S

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)

      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            cmup(i,j) = 0.
            fmup(i,j) = 0.
            fxup(i,j) = 0.
            fyup(i,j) = 0.
            fzup(i,j) = 0.
            fxxup(i,j) = 0.
            fyyup(i,j) = 0.
            fzzup(i,j) = 0.
            fxyup(i,j) = 0.
            fyzup(i,j) = 0.
            fzxup(i,j) = 0.
          enddo
        enddo
      enddo

      do l=1,lmo
        do j=j_0,j_1
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              CM = DT*MW(I,J,L)
              if(cm.ge.0.) then ! mass flux is downward or zero
                C  = CM/MO(I,J,L)
                IF(C.GT.1d0)  WRITE (6,*) 'C>1:',I,J,L,C,MO(I,J,L)
                FM = C*(RM(I,J,L)+(1-C)*(RZ(I,J,L)+(1-2*C)*RZZ(I,J,L)))
                FZ  = CM*(C*C*(RZ(I,J,L) + 3*(1-C)*RZZ(I,J,L)) - 3*FM)
                FZZ = CM*(CM*C**3*RZZ(I,J,L) - 5*(CM*FM+FZ))
                FX = C*(RX(I,J,L)+(1d0-C)*RZX(I,J,L))
                FZX = CM*(C*C*RZX(I,J,L)-3d0*FX)
                FY = C*(RY(I,J,L)+(1d0-C)*RYZ(I,J,L))
                FYZ = CM*(C*C*RYZ(I,J,L)-3d0*FY)
                FXX = C*RXX(I,J,L)
                FYY = C*RYY(I,J,L)
                FXY = C*RXY(I,J,L)
              else              ! mass flux is upward
                C  = CM/MO(I,J,L+1)
                IF(C.LT.-1d0)  WRITE (6,*) 'C<-1:',I,J,L,C,MO(I,J,L+1)
                FM = C*(RM(I,J,L+1) -(1+C)*(
     &               RZ(I,J,L+1) -(1+2*C)*RZZ(I,J,L+1)))
                FZ  = CM*(C*C*(RZ(I,J,L+1) -3*(1+C)*RZZ(I,J,L+1)) -3*FM)
                FZZ = CM*(CM*C**3*RZZ(I,J,L+1) - 5*(CM*FM+FZ))
                FX = C*(RX(I,J,L+1)-(1d0+C)*RZX(I,J,L+1))
                FZX = CM*(C*C*RZX(I,J,L+1)-3d0*FX)
                FY = C*(RY(I,J,L+1)-(1d0+C)*RYZ(I,J,L+1))
                FYZ = CM*(C*C*RYZ(I,J,L+1)-3d0*FY)
                FXX = C*RXX(I,J,L+1)
                FYY = C*RYY(I,J,L+1)
                FXY = C*RXY(I,J,L+1)
              endif
              mnew = MO(I,J,L) + (CMUP(I,J)-CM)
              bymnew = 1d0/mnew
              CM2 = CMUP(I,J)+CM
              RM(I,J,L) = RM(I,J,L) + (FMUP(I,J)-FM)
              RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZUP(I,J)-FZ) +3d0*
     &             (CM2*RM(I,J,L)-MO(I,J,L)*(FMUP(I,J)+FM)))*bymnew
              RZZ(I,J,L) = (RZZ(I,J,L)*MO(I,J,L)**2
     &            +2.5*RM(I,J,L)*(MO(I,J,L)**2 -Mnew**2 -3.*CM2*CM2)
     &            +5.*(MO(I,J,L)*(MO(I,J,L)*(FMUP(i,j)-FM)-FZUP(i,j)-FZ)
     &            +CM2*RZ(I,J,L)*Mnew)
     &            +(FZZUP(i,j)-FZZ) ) * (bymnew*bymnew)
              RX(I,J,L) = RX(I,J,L) + (FXUP(I,J)-FX)
              RZX(I,J,L) = (RZX(I,J,L)*MO(I,J,L) +(FZXUP(I,J)-FZX) +3d0*
     &             (CM2*RX(I,J,L)-MO(I,J,L)*(FXUP(I,J)+FX)))*bymnew
              RY(I,J,L) = RY(I,J,L) + (FYUP(I,J)-FY)
              RYZ(I,J,L) = (RYZ(I,J,L)*MO(I,J,L) +(FYZUP(I,J)-FYZ) +3d0*
     &             (CM2*RY(I,J,L)-MO(I,J,L)*(FYUP(I,J)+FY)))*bymnew
              RXX(I,J,L) = RXX(I,J,L) + (FXXUP(I,J)-FXX)
              RYY(I,J,L) = RYY(I,J,L) + (FYYUP(I,J)-FYY)
              RXY(I,J,L) = RXY(I,J,L) + (FXYUP(I,J)-FXY)
c
              MO(I,J,L) = mnew
              cmup(i,j) = cm
              fmup(i,j) = fm
              fxup(i,j) = fx
              fyup(i,j) = fy
              fzup(i,j) = fz
              fxxup(i,j) = fxx
              fyyup(i,j) = fyy
              fzzup(i,j) = fzz
              fxyup(i,j) = fxy
              fyzup(i,j) = fyz
              fzxup(i,j) = fzx
              OIJL(I,J,L) = OIJL(I,J,L) + FM
            enddo
          enddo
        enddo
      enddo

      return
      END SUBROUTINE OADVTZ3

      subroutine relax_zmoms(M,RM,KWT,RZ,RZZ)
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, DZO_orig=>DZO, BYDZOE
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      USE OCEANR_DIM, only : ogrid
      Implicit None
      Real*8,Intent(In),
     &     Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     &     M,RM,KWT
      Real*8,Intent(InOut),
     &     Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     &     RZ,RZZ
C**** Local variables
      Integer :: I,J,L,N, J1P,JNP
      real*8, dimension(IM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,LMO)
     &     :: r
      real*8, dimension(IM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,0:LMO)
     &     :: rzgrad

      real*8, dimension(0:lmo+1) :: dzo
      real*8, dimension(lmo) :: by6arr,dzby12

      dzo(1:lmo) = dzo_orig(1:lmo)
      dzo(0) = dzo(1)
      dzo(lmo+1) = dzo(lmo)
      do l=1,lmo
        by6arr(l) = .5d0*dzo(l)/sum(dzo(l-1:l+1))
        dzby12(l) = .5d0*by6arr(l)*dzo(l)
      enddo


C**** Extract domain decomposition band parameters
      J1P = oGRID%J_STRT
      JNP = oGRID%J_STOP

C****
C**** Loop over layers and latitudes
C****
      do l=1,lmo
      do j=j1p,jnp
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            r(i,j,l) = rm(i,j,l)/m(i,j,l)
          enddo
        enddo
      enddo
      enddo

      do l=1,lmo-1
      do j=j1p,jnp
        do n=1,nbyzm(j,l+1)
          do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
            rzgrad(i,j,l) = (r(i,j,l+1)-r(i,j,l))*bydzoe(l)
          enddo
        enddo
      enddo
      enddo
      do j=j1p,jnp
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            rzgrad(i,j,0) = rzgrad(i,j,1)
            l = lmom(i,j)
            rzgrad(i,j,l) = rzgrad(i,j,l-1)
          enddo
        enddo
      enddo

      do l=1,lmo
      do j=j1p,jnp
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            rz(i,j,l) = kwt(i,j,l)*rz(i,j,l) + (1.-kwt(i,j,l))*
     &           m(i,j,l)*by6arr(l)*(
     &           dzo(l+1)*rzgrad(i,j,l-1)+dzo(l-1)*rzgrad(i,j,l)
     &           + .5d0*dzo(l)*(rzgrad(i,j,l-1)+rzgrad(i,j,l))  )
            rzz(i,j,l) = kwt(i,j,l)*rzz(i,j,l) + (1.-kwt(i,j,l))*
     &           m(i,j,l)*(dzby12(l)*(rzgrad(i,j,l)-rzgrad(i,j,l-1)))
          enddo
        enddo
      enddo
      enddo
      return
      end subroutine relax_zmoms

      subroutine diffuse_moms(dz,k,q,qz,qzz,n,dt)
      use tridiag_mod, only : tridiag
      implicit none
      integer :: n
      real*8, dimension(n) :: dz,k,q,qz,qzz
      real*8 :: dt
c
      integer :: i,i3,n3
      real*8, dimension(:), allocatable :: q3,a,b,c,r
      n3 = 3*n
      allocate(q3(n3),a(n3),b(n3),c(n3),r(n3))
c store dt/dz in a, k/dze in c
      i3 = 1
      do i=1,n
        c(i3:i3+1) = 3./dz(i)
        a(i3:i3+2) = dt*c(i3)
        i3 = i3+3
      enddo
      i3 = 3
      do i=1,n-1
        c(i3) = k(i)*3./(.5*(dz(i)+dz(i+1)))
        i3 = i3+3
      enddo
c      c(1) = (k(1)/3.)*c(1)
c      c(2) = (k(1)*2./3.)*c(2)
      c(1) = k(1)*c(1)
      c(2) = k(1)*c(2)
      i3 = 2
      do i=2,n
        i3 = i3 + 2
        c(i3) = (k(i-1)*2./3. + k(i)/3.)*c(i3)
        i3 = i3 + 1
        c(i3) = (k(i-1)/3. + k(i)*2./3.)*c(i3)
      enddo
      i3 = 0
      do i=1,n
        i3 = i3 + 1
        r(i3) = q(i) - 2.*qz(i)/3. + 2.*qzz(i)/9.
        i3 = i3 + 1
        r(i3) = q(i)               - 4.*qzz(i)/9.
        i3 = i3 + 1
        r(i3) = q(i) + 2.*qz(i)/3. + 2.*qzz(i)/9.
      enddo

      i3=n3
        c(i3) = 0.
        a(i3) = -c(i3-1)*a(i3)
        b(i3) = 1d0 - a(i3) - c(i3)
      do i3=n3-1,2,-1
        c(i3) = -c(i3)*a(i3)
        a(i3) = -c(i3-1)*a(i3)
        b(i3) = 1d0 - a(i3) - c(i3)
      enddo
      i3=1
        c(i3) = -c(i3)*a(i3)
        a(i3) = 0.
        b(i3) = 1d0 - a(i3) - c(i3)
      call tridiag(a,b,c,r,q3,n3)
      i3 = 1
      do i=1,n
        q(i) = sum(q3(i3:i3+2))/3.
        qz(i) = (q3(i3+2)-q3(i3))*.75
        qzz(i) = (q(i)-q3(i3+1))*2.25
        i3 = i3 + 3
      enddo
      deallocate(q3,a,b,c,r)
      return
      end subroutine diffuse_moms
