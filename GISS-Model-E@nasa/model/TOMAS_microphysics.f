#include "rundeck_opts.h"


!@sum multicoag  :  performs coagulation on the aerosol size distribution
!@+  defined by Nk and Mk (number and mass).  See "An Efficient
!@+   Numerical Solution to the Stochastic Collection Equation", S.
!@+   Tzivion, G. Feingold, and Z. Levin, J Atmos Sci, 44, no 21, 3139-
!@+   3149, 1987.  Unless otherwise noted, all equation references refer
!@+   to this paper.  Some equations are taken from "Atmospheric Chemistry
!@+   and Physics: From Air Pollution to Climate Change" by Seinfeld
!@+   and Pandis (S&P).  

!@+   This routine uses a "moving sectional" approach in which the
!@+   aerosol size bins are defined in terms of dry aerosol mass.
!@+   Addition or loss of water, therefore, does not affect which bin
!@+   a particle falls into.  As a result, this routine does not
!@+   change Mk(water), although water masses are needed to compute
!@+   particle sizes and, therefore, coagulation coefficients.  Aerosol
!@+   water masses in each size bin will need to be updated later
!@+   (in another routine) to reflect changes that result from
!@+   coagulation.
!@+   The user must supply the mass and number distributions, Mk and Nk,
!@+   as well as the time step, dt.

!@auth Yunha Lee

#if (defined TOMAS_12_10NM) || (defined TOMAS_12_3NM)

      SUBROUTINE multicoag(dt)


      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk
      USE CONSTANT,   only:  pi,gasc   
      IMPLICIT NONE

      real dt !time step [s]
      integer n,c,bh,bl,dts_old,count ! for 15 size bins lumping
      integer k,j,i,jj,kk    !counters
      real*8 dNdt(ibins), dMdt(ibins,icomp-idiag)
!@var   dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
!@+   the mass boundaries of the size bins.  xbar is the average mass
!@+   of a given size bin (it varies with time in this algorithm).  phi
!@+   and eff are defined in the reference, equations 13a and b.

      real*8 xbar(ibins), phi(ibins), eff(ibins)
      real*8 Nki(ibins), Mki(ibins, icomp) 
!@var kij represents the coagulation coefficient (cm3/s) normalized by the
!@+   volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
      real kij(ibins,ibins)
      real Dpk(ibins)             !diameter (m) of particles in bin k
      real Dk(ibins)              !Diffusivity (m2/s) of bin k particles
      real ck(ibins)              !Mean velocity (m/2) of bin k particles
      real olddiff                !used to iterate to find diffusivity
      real density                !density (kg/m3) of particles
      real mu                     !viscosity of air (kg/m s)
      real mfp                    !mean free path of air molecule (m)
      real Kn                     !Knudsen number of particle
      real*8 mp         !particle mass (kg)
      real beta                   !correction for coagulation coeff.
      real*8 aerodens
      external aerodens
      real*8 Mktot      !total mass of aerosol

      !temporary summation variables
      real*8 k1m(icomp-idiag),k1mx(icomp-idiag)
      real*8 k1mx2(icomp-idiag)
      real*8 k1mtot,k1mxtot
      real*8 sk2mtot, sk2mxtot
      real*8 sk2m(icomp-idiag), sk2mx(icomp-idiag)
      real*8 sk2mx2(icomp-idiag)
      real*8 High_in,dNtot,dMtot
      real*8 mtotal

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl   

      real zeta                      !see reference, eqn 6
      real tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
      real dts  !internal time step (<dt for stability)
      real tsum !time so far
      real*8 Neps !minimum value for Nk
      real*8 mi, mf   !initial and final masses
      logical is_nan
      external is_nan
      parameter(zeta=1.28125, dtlimit=0.25, itlimit=10)
      real kB  !kB is Boltzmann constant (J/K)
      parameter (kB=1.38e-23, Neps=1.0e-3)

C-----CODE--------------------------------------------------------------

      tsum = 0.0
      dts_old=10.
C If any Nk are zero, then set them to a small value to avoid division by zero
      do k=1,ibins
         if (Nk(k) .lt. Neps) then
            Nk(k)=Neps
            Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
            do j=1,icomp
               if (j.ne.srtso4)then
                  Mk(k,j)=0.d0
               endif
            enddo
         endif
cyhl To check whether mass is conserved during coagulation. 
c         print*,'Nk',k,Nk(k),Mk(k,2)
      enddo

      Nki(:)=Nk(:)
      Mki(:,:) =Mk(:,:)

C Calculate air viscosity and mean free path

      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6

      do k=1,ibins

         mso4=Mk(k,srtso4) 
         mnacl=Mk(k,srtna)
         mno3=0.e0
         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
         mnh4=Mk(k,srtnh4) !0.1875*mso4  !assume ammonium bisulfate
         mecob=Mk(k,srtecob)
         mecil=Mk(k,srtecil)
         mocil=Mk(k,srtocil)
         mocob=Mk(k,srtocob)
         mdust=Mk(k,srtdust)          
         mh2o=Mk(k,srth2o)   

         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 

         Mktot=0.d0
         do j=1,icomp
            Mktot=Mktot+Mk(k,j)
         enddo
         mp=Mktot/Nk(k)
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
          ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
      enddo

C Calculate coagulation coefficients

      do i=1,ibins
         do j=1,ibins
            Kn=4.0*(Dk(i)+Dk(j))          
     &        /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
            beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
            kij(j,i)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
            kij(j,i)=kij(j,i)*1.0e6/boxvol  !normalize by grid cell volume
         enddo
      enddo
 10   continue     !repeat process here if multiple time steps are needed

C Calculate xbar, phi and eff

      do k=1,ibins

         xbar(k)=0.0
         do j=1,icomp-idiag
            xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
         enddo

       if(k.lt.ibins-1)then !from 1 to 10 bins

         eff(k)=2./9.*Nk(k)/xk(k)
     *        *(4.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
         phi(k)=2./9.*Nk(k)/xk(k)
     *        *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2./3.*Nk(k)/xk(k)
            phi(k)=0.0

         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2./3.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
      else                      ! from 11 bins to 12 bins
            eff(k)=2./31./31.*Nk(k)/xk(k)
     *           *(32.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
            phi(k)=2./31./31.*Nk(k)/xk(k)
     *           *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2./31.*Nk(k)/xk(k)
            phi(k)=0.0

         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2./31.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
      endif

      enddo

C Necessary initializations
         sk2mtot=0.0
         sk2mxtot=0.0
         do j=1,icomp-idiag
            sk2m(j)=0.0
            sk2mx(j)=0.0
            sk2mx2(j)=0.0
         enddo

C Calculate rates of change for Nk and Mk

      do k=1,ibins

         !Initialize to zero
         do j=1,icomp-idiag
            k1m(j)=0.0
            k1mx(j)=0.0
            k1mx2(j)=0.0
         enddo
         High_in=0.0
         k1mtot=0.0
         k1mxtot=0.0

         !Calculate sums
         do j=1,icomp-idiag
            if (k .gt. 1.and.k.lt.ibins) then
               do i=1,k-1
                  k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                  k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)*zeta
                  k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2.
     *                 *zeta**3.
               enddo
            elseif(k.eq.ibins)then
                  k1m(j)= sk2m(j)+kij(k,k-1)*Mk(k-1,j)
                  k1mx(j)=sk2mx(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)*4.754
                  k1mx2(j)=sk2mx2(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)**2.
     *                 *107.4365
            endif
            k1mtot=k1mtot+k1m(j)
            k1mxtot=k1mxtot+k1mx(j)
         enddo
         if(k.lt.ibins)then
            do i=k+1, ibins
               High_in=High_in+Nk(i)*kij(k,i)
            enddo
         endif

         if(k.lt.ibins-1)then

         dNdt(k)= -Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.125
     *    -(phi(k)*k1mtot+(eff(k)-phi(k))/6./xk(k)*k1mxtot)
     *    -kij(k,k)*(phi(k)/3.*xbar(k)*Nk(k)+(eff(k)-phi(k))/18.
     *        /xk(k)*zeta*xbar(k)*xbar(k)*Nk(k))

         if (k .gt. 1) then      
cyhl Nk*low_in changes to -0.5*Kij*Nk**2.
         dNdt(k)=dNdt(k)+0.625*kij(k-1,k-1)*Nk(k-1)**2 
     * +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1)
     *        *sk2mxtot)
     * +kij(k-1,k-1)*(phi(k-1)/3.*xbar(k-1)*Nk(k-1)+(eff(k-1)
     *  -phi(k-1))/18./xk(k-1)*zeta*xbar(k-1)*xbar(k-1)*Nk(k-1))


cyhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
cyhl It results in much closer result as 30 bins. Apr.27.08
 
      call nanstop(dNdt(k),246,k,0)

         endif

         do j=1,icomp-idiag
                
            dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     *  -(phi(k)*xk(k+1)*k1m(j)+(eff(k)+2.*phi(k))/6.*k1mx(j)
     *           +(eff(k)-phi(k))/6./xk(k)*k1mx2(j)) ! term3
     *  - kij(k,k)*Nk(k)*Mk(k,j)/3. ! I assume 1/2Nk and 2/3Mk for half bin
     *  - kij(k,k)*(phi(k)*xk(k+1)*Mk(k,j)/3.
     *           +(eff(k)+2.*phi(k))/6.*zeta*xbar(k)*Mk(k,j)/3.
     * +(eff(k)-phi(k))/6./xk(k)*zeta**3.*xbar(k)**2.*Mk(k,j)/3.)
      
cyhl  Term9(-kij(k,k)*Nk(k)*Mk(k,j)) is cancled out by term6 (k)  
          
            if (k .gt. 1) then
               dMdt(k,j)=dMdt(k,j)
     *  +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)+2.*phi(k-1))/6.*sk2mx(j)
     *        +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) !term1
     *  +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3.
     *  +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3.
     *    +(eff(k-1)+2.*phi(k-1))/6.*zeta
     *    *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6.
     *    /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)

               call nanstop(dMdt(k,J),265,k,J)    
            endif
         enddo
      else if (k.eq.ibins-1)then

         dNdt(k)=0.625*kij(k-1,k-1)*Nk(k-1)**2 
     * +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1)
     *        *sk2mxtot)
     * +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)/3.*(phi(k-1)+(eff(k-1)
     *  -phi(k-1))/6./xk(k-1)*zeta*xbar(k-1))

cyhl updated the following
         dNdt(k)=dNdt(k)-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.02
     *    -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot)
     *    -kij(k,k)*xbar(k)*Nk(k)*0.484*(phi(k)+(eff(k)-phi(k))/62.
     *        /xk(k)*4.754*xbar(k))

cyhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
cyhl It results in much closer result as 30 bins. Apr.27.08

         do j=1,icomp-idiag
            dMdt(k,j)=
     *  +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)+2.*phi(k-1))/6.*sk2mx(j)
     *        +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) !term1
     *  +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3.
     *  +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3.
     *    +(eff(k-1)+2.*phi(k-1))/6.*zeta
     *    *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6.
     *    /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)


cyhl updated the following
            dMdt(k,j)= dMdt(k,j)+Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     * -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k))
     *    *k1mx(j)+(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) ! term3
     *  - kij(k,k)*Nk(k)*Mk(k,j)*0.103226 ! I assume 1/2Nk and 2/3Mk for half bin
     *  - kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62.
     *   +0.484*phi(k))*4.754*xbar(k)
     * +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)
         enddo

      else if (k.eq.ibins)then
         dNdt(k)=-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.103226
     *    -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot)
     *    -kij(k,k)*0.484*xbar(k)*Nk(k)*(phi(k)+(eff(k)-phi(k))
     *       /62./xk(k)*4.754*xbar(k))
     *    +0.52*kij(k-1,k-1)*Nk(k-1)**2 
     * +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mxtot)
     * +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)*0.484*(phi(k-1)+(eff(k-1)
     *  -phi(k-1))/62./xk(k-1)*4.754*xbar(k-1))

         do j=1,icomp-idiag
            dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     * -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k))*k1mx(j)
     *   +(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) ! term3
     *  - kij(k,k)*Nk(k)*Mk(k,j)*0.103226 ! I assume 1/2Nk and 2/3Mk for half bin
     *  - kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62.
     *   +0.484*phi(k))*4.754*xbar(k)
     * +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)

     *  +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)/62.+0.484*phi(k-1))
     *    *sk2mx(j)+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mx2(j)) !term1
     *  +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)*0.103226
     *  +kij(k-1,k-1)*Mk(k-1,j)*0.484*(phi(k-1)*xk(k)
     * +(eff(k-1)/62.+0.484*phi(k-1))*4.754*xbar(k-1)
     * +(eff(k-1)-phi(k-1))/62./xk(k-1)*107.4365*xbar(k-1)**2.)
         enddo
      endif

         !Save the summations that are needed for the next size bin
         sk2mtot=k1mtot
         sk2mxtot=k1mxtot
         do j=1,icomp-idiag
            sk2m(j)=k1m(j)
            sk2mx(j)=k1mx(j)
            sk2mx2(j)=k1mx2(j)
         enddo
c         print*,'dNdt',Nk(k),dNdt(k)
      enddo  !end of main k loop


C Update Nk and Mk according to rates of change and time step

      !If any Mkj are zero, add a small amount to achieve finite
      !time steps
      do k=1,ibins
         do j=1,icomp-idiag
            if (Mk(k,j) .eq. 0.d0) then
               !add a small amount of mass
               mtotal=0.d0
               do jj=1,icomp-idiag
                  mtotal=mtotal+Mk(k,jj)
               enddo
               Mk(k,j)=1.d-10*mtotal
            endif
         enddo
      enddo

      !Choose time step
      dts=dt-tsum      !try to take entire remaining time step
cdbg      limit='comp'
      do k=1,ibins
         if (Nk(k) .gt. Neps) then
            !limit rates of change for this bin
            if (dNdt(k) .lt. 0.0) tlimit=dtlimit
            if (dNdt(k) .gt. 0.0) tlimit=itlimit
            if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then 
               dts=Nk(k)*tlimit/abs(dNdt(k))
!               if(dts.eq.0) print*,'dts Nk',k,dts,dNdt(k),Nk(k)
            endif
            do j=1,icomp-idiag
               if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
               if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
               if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then 
                  mtotal=0.d0
                  do jj=1,icomp-idiag
                     mtotal=mtotal+Mk(k,jj)
                  enddo         !only use this criteria if this species is significant

                  if ((Mk(k,j)/mtotal) .gt. 5.d-4) then
                     dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
 
                     if(dts.eq.0)print*,'dts MK',k,j,dts,dt,tsum,Mk(k,j)
     &                    ,nk(k),mtotal,dMdt(k,j),xk(k),xk(k+1)
                     if(dts.eq.0.and.dts_old.eq.0.) then
                       open (1044,file='debug_coag.dat',access='append',
     &                      status='unknown')
                       write(1044,*)'start',k,j,dts,dt,tsum,Mk(k,j)
     &                    ,nk(k),mtotal,dMdt(k,j),xk(k),xk(k+1)
                       write(1044,*)'dts MK, T,P',temp,pres
                       write(1044,*)'dts Mk, BOX',boxvol,boxmass,rh
                       do kk=1,ibins
                       write(1044,*)'dts MK, Nk=',Nk(kk),Nki(kk)
                       enddo
                       do jj=1,icomp
                         write(1044,*)'comp',jj
                       do kk=1,ibins
            write(1044,*)'dts MK, Mk=',Mk(kk,jj),Mki(kk,jj)
                       enddo
                       enddo
                     endif

                  else
                     if (dMdt(k,j) .lt. 0.0) then !set dmdt to 0 to avoid very small mk going negative
                        dMdt(k,j)=0.0
                     endif
                  endif
c                  print*,'mk',k,j,Mk(k,j)/mtotal,dMdt(k,j)
               endif
            enddo
c            print*,'dNdt',k,Nk(k),dNdt(k),Mk(k,1),dMdt(k,1) 
         else
            !nothing in this bin - don't let it affect time step
            Nk(k)=Neps
            Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
            if (dNdt(k) .lt. 0.0) dNdt(k)=0.0  !make sure mass/number don't go negative
            do j=1,icomp-idiag
               if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
            enddo
         endif
      enddo

c      if (dts .lt. 20.) write(*,*) 'dts<20. in multicoag',dts,tsum
!       if (dts .eq. 0.) then
!          write(*,*) 'time step is 0. dts=',dts
!          call stop_model('dts=0 in multicoag',255)
!       endif

!YUNHA (Sep, 2012) this is newly added to prevent an occasional crash. 
      if(dts.eq.0) then
! When dts is small and this is due to the mass/number is out of size range, 
! it calls mnfix. Doing this, it will help to avoid dts=0 case. 
        Nki(:)=Nk(:)
        Mki(:,:)=Mk(:,:)
        call mnfix(Nki,Mki)       !YUNHA LEE (08/28/2012) 
        count=0
        do k=1,ibins
          if(Nki(k).eq.Nk(k)) then
            count=count+1
          endif
        enddo

        if(count.lt.ibins) then
          print*,'dts=0 but mnfix worked',dts,tsum,count ! count<5 is random choice
          Nk(:)=Nki(:)
          Mk(:,:)=Mki(:,:)
        endif
      endif
     
      if(dts.eq.0.and.count.eq.ibins)then
!This case, next dts will be zero and model should be stopped. 
        write(*,*) 'time step is 0',count
        call stop_model('dts=0 in multicoag',255)
      endif

      do k=1,ibins
         Nk(k)=Nk(k)+dNdt(k)*dts
         do j=1,icomp-idiag
            Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
         enddo
      enddo

      if(dts.lt.1e-10.and.dts_old.lt.1e-10) then
! When dts is small in two sequently, then mnfix is called. It might help to get a larger dts next time.  However, I am not sure if it is helpful   !YUNHA LEE (09/05/2012) 
        Nki(:)=Nk(:)
        Mki(:,:)=Mk(:,:)
        call mnfix(Nki,Mki)     
        Nk(:)=Nki(:)
        Mk(:,:)=Mki(:,:)
      endif

      tsum=tsum+dts
      dts_old=dts
      if (tsum .lt. dt) goto 10

      RETURN
      END SUBROUTINE multicoag
#endif 


c$$$  This is multicoag for TOMAS-30 model
c$$$
c$$$!@+   **************************************************
c$$$!@+   *  multicoag                                     *
c$$$!@+   **************************************************
c$$$
c$$$!@auth   Peter Adams, June 1999
c$$$!@+   Modified to allow for multicomponent aerosols, February 2000
c$$$
c$$$!@sum   :  performs coagulation on the aerosol size distribution
c$$$!@+   defined by Nk and Mk (number and mass).  See "An Efficient
c$$$!@+   Numerical Solution to the Stochastic Collection Equation", S.
c$$$!@+   Tzivion, G. Feingold, and Z. Levin, J Atmos Sci, 44, no 21, 3139-
c$$$!@+   3149, 1987.  Unless otherwise noted, all equation references refer
c$$$!@+   to this paper.  Some equations are taken from "Atmospheric Chemistry
c$$$!@+   and Physics: From Air Pollution to Climate Change" by Seinfeld
c$$$!@+   and Pandis (S&P).  
c$$$
c$$$!@+   This routine uses a "moving sectional" approach in which the
c$$$!@+   aerosol size bins are defined in terms of dry aerosol mass.
c$$$!@+   Addition or loss of water, therefore, does not affect which bin
c$$$!@+   a particle falls into.  As a result, this routine does not
c$$$!@+   change Mk(water), although water masses are needed to compute
c$$$!@+   particle sizes and, therefore, coagulation coefficients.  Aerosol
c$$$!@+   water masses in each size bin will need to be updated later
c$$$!@+   (in another routine) to reflect changes that result from
c$$$!@+   coagulation.
c$$$
c$$$C-----INPUTS------------------------------------------------------------
c$$$
c$$$!@+   The user must supply the mass and number distributions, Mk and Nk,
c$$$!@+   as well as the time step, dt.
c$$$
c$$$C-----OUTPUTS-----------------------------------------------------------
c$$$
c$$$!@+   The program updates Nk and Mk.
c$$$
c$$$      SUBROUTINE multicoag(dt)
c$$$
c$$$      USE TOMAS_AEROSOL
c$$$      USE TRACER_COM, only : xk
c$$$      USE CONSTANT,   only:  pi,gasc   
c$$$      IMPLICIT NONE
c$$$
c$$$C-----ARGUMENT DECLARATIONS---------------------------------------------
c$$$
c$$$c      real dt         !time step (s)
c$$$      real dt
c$$$
c$$$C-----VARIABLE DECLARATIONS---------------------------------------------
c$$$
c$$$      integer k,j,i,jj,kk    !counters
c$$$      real*8 dNdt(ibins), dMdt(ibins,icomp-idiag)
c$$$      real*8 xbar(ibins), phi(ibins), eff(ibins)
c$$$
c$$$C kij represents the coagulation coefficient (cm3/s) normalized by the
c$$$C volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
c$$$      real kij(ibins,ibins)
c$$$      real Dpk(ibins)             !diameter (m) of particles in bin k
c$$$      real Dk(ibins)              !Diffusivity (m2/s) of bin k particles
c$$$      real ck(ibins)              !Mean velocity (m/2) of bin k particles
c$$$      real olddiff                !used to iterate to find diffusivity
c$$$      real density                !density (kg/m3) of particles
c$$$      real mu                     !viscosity of air (kg/m s)
c$$$      real mfp                    !mean free path of air molecule (m)
c$$$      real Kn                     !Knudsen number of particle
c$$$      real*8 mp         !particle mass (kg)
c$$$      real beta                   !correction for coagulation coeff.
c$$$      real aerodens
c$$$      external aerodens
c$$$      real*8 Mktot      !total mass of aerosol
c$$$
c$$$      !temporary summation variables
c$$$      real*8 k1m(icomp-idiag),k1mx(icomp-idiag)
c$$$      real*8 k1mx2(icomp-idiag)
c$$$      real*8 k1mtot,k1mxtot
c$$$      real*8 sk2mtot, sk2mxtot
c$$$      real*8 sk2m(icomp-idiag), sk2mx(icomp-idiag)
c$$$      real*8 sk2mx2(icomp-idiag)
c$$$      real*8 in
c$$$      real*8 mtotal
c$$$
c$$$      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
c$$$      real mecil,mecob,mocil,mocob
c$$$      real mdust,mnacl   
c$$$
c$$$      real zeta                      !see reference, eqn 6
c$$$      real tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
c$$$      real dts  !internal time step (<dt for stability)
c$$$      real tsum !time so far
c$$$      real*8  Neps              !minimum value for Nk
c$$$cdbg      character*12 limit        !description of what limits time step
c$$$
c$$$      real*8 mi, mf             !initial and final masses
c$$$      logical is_nan
c$$$      external is_nan
c$$$      
c$$$!@+   VARIABLE COMMENTS...
c$$$
c$$$!@+   dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
c$$$!@+   the mass boundaries of the size bins.  xbar is the average mass
c$$$!@+   of a given size bin (it varies with time in this algorithm).  phi
c$$$!@+   and eff are defined in the reference, equations 13a and b.
c$$$
c$$$C-----ADJUSTABLE PARAMETERS---------------------------------------------
c$$$
c$$$      parameter(zeta=1.0625, dtlimit=0.25, itlimit=10.)
c$$$      real kB  !kB is Boltzmann constant (J/K)
c$$$      parameter (kB=1.38e-23, Neps=1.0e-3)
c$$$
c$$$ 1    format(16E15.3)
c$$$
c$$$C-----CODE--------------------------------------------------------------
c$$$
c$$$      tsum = 0.0
c$$$
c$$$C If any Nk are zero, then set them to a small value to avoid division by zero
c$$$      do k=1,ibins
c$$$         if (Nk(k) .lt. Neps) then
c$$$            Nk(k)=Neps
c$$$            Mk(k,srtso4)=Neps*sqrt(xk(k+1)*xk(k)) !make the added particles SO4
c$$$            do j=1,icomp
c$$$               if (j.ne.srtso4)then
c$$$                  Mk(k,j)=0.d0
c$$$               endif
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$
c$$$C Calculate air viscosity and mean free path
c$$$
c$$$      mu=2.5277e-7*temp**0.75302
c$$$      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6
c$$$!      call nanstop(mfp,124, 0,0)
c$$$C Calculate particle sizes and diffusivities
c$$$      do k=1,ibins
c$$$
c$$$         mso4=Mk(k,srtso4) 
c$$$         mnacl=Mk(k,srtna)
c$$$         mno3=0.e0
c$$$         if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
c$$$         mnh4=0.1875*mso4  !assume ammonium bisulfate
c$$$         mecob=Mk(k,srtecob)
c$$$         mecil=Mk(k,srtecil)
c$$$         mocil=Mk(k,srtocil)
c$$$         mocob=Mk(k,srtocob)
c$$$         mdust=Mk(k,srtdust)          
c$$$         mh2o=Mk(k,srth2o)   
c$$$
c$$$         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
c$$$     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
c$$$
c$$$c$$$         density=aerodens(Mk(k,srtso4),0.d0,Mk(k,srtnh4),
c$$$c$$$     &        Mk(k,srtna),Mk(k,srtecil),Mk(k,srtecob),
c$$$c$$$     &        Mk(k,srtocil),Mk(k,srtocob),Mk(k,srtdust),
c$$$c$$$     &        Mk(k,srth2o))    !assume bisulfate
c$$$!         call nanstop(density,130,k,0)
c$$$
c$$$Ckpc  Add 0.2x first for ammonium, and then add 1.0x in the loop
c$$$         Mktot=0.d0
c$$$         do j=1,icomp
c$$$            Mktot=Mktot+Mk(k,j)
c$$$         enddo
c$$$         mp=Mktot/Nk(k)
c$$$
c$$$         Dpk(k)=((mp/density)*(6./pi))**(0.333)
c$$$      call nanstop(dpk(k),139,k,0)
c$$$         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
c$$$
c$$$         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
c$$$     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
c$$$
c$$$         ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
c$$$
c$$$      enddo
c$$$
c$$$C Calculate coagulation coefficients
c$$$
c$$$      do i=1,ibins
c$$$         do j=1,ibins
c$$$            Kn=4.0*(Dk(i)+Dk(j))          
c$$$     &        /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
c$$$
c$$$            beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
c$$$
c$$$            !This is S&P eqn 12.46 with non-continuum correction, beta
c$$$            kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
c$$$
c$$$            kij(i,j)=kij(i,j)*1.0e6/boxvol  !normalize by grid cell volume
c$$$
c$$$         enddo
c$$$      enddo
c$$$
c$$$ 10   continue     !repeat process here if multiple time steps are needed
c$$$
c$$$C Calculate xbar, phi and eff
c$$$
c$$$      do k=1,ibins
c$$$
c$$$         xbar(k)=0.0
c$$$         do j=1,icomp-idiag
c$$$            xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
c$$$         enddo
c$$$
c$$$         eff(k)=2.*Nk(k)/xk(k)*(2.-xbar(k)/xk(k))    !eqn 13a
c$$$      call nanstop(eff(k),180,k,0)
c$$$         phi(k)=2.*Nk(k)/xk(k)*(xbar(k)/xk(k)-1.)    !eqn 13b
c$$$      call nanstop(phi(k),182,k,0)   
c$$$         !Constraints in equation 15
c$$$         if (xbar(k) .lt. xk(k)) then
c$$$            eff(k)=2.*Nk(k)/xk(k)
c$$$            phi(k)=0.0
c$$$         else if (xbar(k) .gt. xk(k+1)) then
c$$$            phi(k)=2.*Nk(k)/xk(k)
c$$$            eff(k)=0.0
c$$$         endif
c$$$      enddo
c$$$
c$$$C Necessary initializations
c$$$         sk2mtot=0.0
c$$$         sk2mxtot=0.0
c$$$         do j=1,icomp-idiag
c$$$            sk2m(j)=0.0
c$$$            sk2mx(j)=0.0
c$$$            sk2mx2(j)=0.0
c$$$         enddo
c$$$
c$$$C Calculate rates of change for Nk and Mk
c$$$
c$$$      do k=1,ibins
c$$$
c$$$         !Initialize to zero
c$$$         do j=1,icomp-idiag
c$$$            k1m(j)=0.0
c$$$            k1mx(j)=0.0
c$$$            k1mx2(j)=0.0
c$$$         enddo
c$$$         in=0.0
c$$$         k1mtot=0.0
c$$$         k1mxtot=0.0
c$$$
c$$$         !Calculate sums
c$$$         do j=1,icomp-idiag
c$$$            if (k .gt. 1) then
c$$$               do i=1,k-1
c$$$
c$$$                  k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
c$$$                  k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)
c$$$                  k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2
c$$$               enddo
c$$$            endif
c$$$            k1mtot=k1mtot+k1m(j)
c$$$            k1mxtot=k1mxtot+k1mx(j)
c$$$         enddo
c$$$         if (k .lt. ibins) then
c$$$            do i=k+1,ibins
c$$$               in=in+Nk(i)*kij(k,i)
c$$$            enddo
c$$$         endif
c$$$
c$$$         !Calculate rates of change
c$$$         dNdt(k)= 
c$$$     &           -kij(k,k)*Nk(k)**2
c$$$     &           -phi(k)*k1mtot
c$$$     &           -zeta*(eff(k)-phi(k))/(2*xk(k))*k1mxtot
c$$$     &           -Nk(k)*in
c$$$      call nanstop(dNdt(k),240,k,0)
c$$$         if (k .gt. 1) then
c$$$         dNdt(k)=dNdt(k)+
c$$$     &           0.5*kij(k-1,k-1)*Nk(k-1)**2
c$$$     &           +phi(k-1)*sk2mtot
c$$$     &           +zeta*(eff(k-1)-phi(k-1))/(2*xk(k-1))*sk2mxtot
c$$$      call nanstop(dNdt(k),246,k,0)
c$$$         endif
c$$$
c$$$         do j=1,icomp-idiag
c$$$            dMdt(k,j)= 
c$$$     &           +Nk(k)*k1m(j)
c$$$     &           -kij(k,k)*Nk(k)*Mk(k,j)
c$$$     &           -Mk(k,j)*in
c$$$     &           -phi(k)*xk(k+1)*k1m(j)
c$$$     &           -0.5*zeta*eff(k)*k1mx(j)
c$$$     &           +zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
c$$$      call nanstop(dMdt(k,j),257,k,j)
c$$$C      write(79,*) '243,k,j,dMdt',k,j,dMdt(k,j)
c$$$            if (k .gt. 1) then
c$$$               dMdt(k,j)=dMdt(k,j)+
c$$$     &           kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)
c$$$     &           +phi(k-1)*xk(k)*sk2m(j)
c$$$     &           +0.5*zeta*eff(k-1)*sk2mx(j)
c$$$     &           -zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
c$$$      call nanstop(dMdt(k,j),265,k,j)
c$$$C      write(79,*) '243,k,j,dMdt',k,j,dMdt(k,j)
c$$$            endif
c$$$cdbg            if (j. eq. srtso4) then
c$$$cdbg               if (k. gt. 1) then
c$$$cdbg                  write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j),
c$$$cdbg     &               Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j),
c$$$cdbg     &               0.5*zeta*eff(k)*k1mx(j),
c$$$cdbg     &               zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j),
c$$$cdbg     &               kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j),
c$$$cdbg     &               phi(k-1)*xk(k)*sk2m(j),
c$$$cdbg     &               0.5*zeta*eff(k-1)*sk2mx(j),
c$$$cdbg     &               zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
c$$$cdbg               else
c$$$cdbg                  write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j),
c$$$cdbg     &               Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j),
c$$$cdbg     &               0.5*zeta*eff(k)*k1mx(j),
c$$$cdbg     &               zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
c$$$cdbg               endif
c$$$cdbg            endif
c$$$         enddo
c$$$
c$$$cdbg         write(*,*) 'k,dNdt,dMdt: ', k, dNdt(k), dMdt(k,srtso4)
c$$$
c$$$         !Save the summations that are needed for the next size bin
c$$$         sk2mtot=k1mtot
c$$$         sk2mxtot=k1mxtot
c$$$         do j=1,icomp-idiag
c$$$            sk2m(j)=k1m(j)
c$$$            sk2mx(j)=k1mx(j)
c$$$            sk2mx2(j)=k1mx2(j)
c$$$         enddo
c$$$
c$$$      enddo  !end of main k loop
c$$$
c$$$C Update Nk and Mk according to rates of change and time step
c$$$
c$$$      !If any Mkj are zero, add a small amount to achieve finite
c$$$      !time steps
c$$$      do k=1,ibins
c$$$         do j=1,icomp-idiag
c$$$            if (Mk(k,j) .eq. 0.d0) then
c$$$               !add a small amount of mass
c$$$               mtotal=0.d0
c$$$               do jj=1,icomp-idiag
c$$$                  mtotal=mtotal+Mk(k,jj)
c$$$               enddo
c$$$               Mk(k,j)=1.d-10*mtotal
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$
c$$$      !Choose time step
c$$$      dts=dt-tsum      !try to take entire remaining time step
c$$$cdbg      limit='comp'
c$$$      do k=1,ibins
c$$$         if (Nk(k) .gt. Neps) then
c$$$            !limit rates of change for this bin
c$$$            if (dNdt(k) .lt. 0.0) tlimit=dtlimit
c$$$            if (dNdt(k) .gt. 0.0) tlimit=itlimit
c$$$            if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then 
c$$$               dts=Nk(k)*tlimit/abs(dNdt(k))
c$$$
c$$$C      write(79,*) 'k,tlimit,dts',k,tlimit,dts
c$$$cdbg               limit='number'
c$$$cdbg               write(limit(8:9),'(I2)') k
c$$$cdbg               write(*,*) Nk(k), dNdt(k)
c$$$            endif
c$$$            do j=1,icomp-idiag
c$$$               if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
c$$$               if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
c$$$               if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then 
c$$$               mtotal=0.d0
c$$$               do jj=1,icomp-idiag
c$$$                  mtotal=mtotal+Mk(k,jj)
c$$$               enddo
c$$$               !only use this criteria if this species is significant
c$$$               if ((Mk(k,j)/mtotal) .gt. 1.d-5) then
c$$$                  dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
c$$$
c$$$C      write(79,*) 'k,j,tlimit,dts',k,j,tlimit,dts
c$$$               else
c$$$                  if (dMdt(k,j) .lt. 0.0) then
c$$$                     !set dmdt to 0 to avoid very small mk going negative
c$$$                     dMdt(k,j)=0.0
c$$$                  endif
c$$$               endif
c$$$cdbg                  limit='mass'
c$$$cdbg                  write(limit(6:7),'(I2)') k
c$$$cdbg                  write(limit(9:9),'(I1)') j
c$$$cdbg                  write(*,*) Mk(k,j), dMdt(k,j)
c$$$               endif
c$$$            enddo
c$$$         else
c$$$            !nothing in this bin - don't let it affect time step
c$$$            Nk(k)=Neps
c$$$            Mk(k,srtso4)=Neps*1.4*xk(k) !make the added particles SO4
c$$$            !make sure mass/number don't go negative
c$$$            if (dNdt(k) .lt. 0.0) dNdt(k)=0.0
c$$$            do j=1,icomp-idiag
c$$$               if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
c$$$            enddo
c$$$         endif
c$$$      enddo
c$$$c      if (dts .lt. 20.) write(*,*) 'dts<20. in multicoag'
c$$$       if (dts .eq. 0.) then
c$$$          write(*,*) 'time step is 0'
c$$$C          pause
c$$$          stop
c$$$C       go to 20
c$$$       endif
c$$$
c$$$      !Change Nk and Mk
c$$$cdbg      write(*,*) 't=',tsum+dts,' ',limit
c$$$      do k=1,ibins
c$$$         Nk(k)=Nk(k)+dNdt(k)*dts
c$$$         do j=1,icomp-idiag
c$$$            Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
c$$$         enddo
c$$$      enddo
c$$$c      print *, 'tsum=', tsum, 'dts=', dts
c$$$      !Update time and repeat process if necessary
c$$$      tsum=tsum+dts
c$$$c      print *, 'tsum=', tsum
c$$$      if (tsum .lt. dt) goto 10
c$$$
c$$$      RETURN
c$$$      END subroutine multicoag



!@sum cond_nuc   :  calculates the change in the aerosol size distribution
!@+   due to so4 condensation and binary/ternary nucleation during the
!@+   overal microphysics timestep.

!@auth   Jeff Pierce, May 2007

!@+   Initial values of
!@+   =================

!@var   Nki(ibins) - number of particles per size bin in grid cell
!@var   Nnuci - number of nucleation size particles per size bin in grid cell
!@var   Mnuci - mass of given species in nucleation pseudo-bin (kg/grid cell)
!@var   Mki(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!@+                  gas phase except water
!@var   H2SO4rate - rate of H2SO4 chemical production [kg s^-1]
!@var   dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

!@var   Nkf, Mkf, Gcf - same as above, but final values
!@var   Nknuc, Mknuc - same as above, final values from just nucleation
!@var   Nkcond, Mkcond - same as above, but final values from just condensation
!@var   fn, fn1

      SUBROUTINE cond_nuc(Nki,Mki,Gci,Nkf,Mkf,Gcf,fnavg,fn1avg,
     &     H2SO4rate,dti,num_iter,Nknuc,Mknuc,Nkcond,Mkcond,lev)            

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk
      USE DOMAIN_DECOMP_ATM, only : am_i_root
      IMPLICIT NONE

      real*8 Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      real*8  Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      real*8 Nknuc(ibins), Mknuc(ibins, icomp)
      real*8 Nkcond(ibins),Mkcond(ibins,icomp)
      real*8 H2SO4rate
      real dti
      real*8 fnavg    ! nucleation rate of clusters cm-3 s-1
      real*8 fn1avg   ! formation rate of particles to first size bin cm-3 s-1
      integer, intent(in) :: lev
      real*8 dt
      integer i,j,k,c,jc,l,n,flag  !dmw         ! counters
      real*8 fn       ! nucleation rate of clusters cm-3 s-1
      real*8 fn1      ! formation rate of particles to first size bin cm-3 s-1
      real*8 CSi,CSa   ! intial and average condensation sinks
      real*8 CS1,CS2       ! guesses for condensation sink [s^-1]
      real*8 CStest   !guess for condensation sink
      real*8 Nk1(ibins), Mk1(ibins, icomp), Gc1(icomp-1)
      real*8 Nk2(ibins), Mk2(ibins, icomp), Gc2(icomp-1)
      real*8 Nk3(ibins), Mk3(ibins, icomp), Gc3(icomp-1)
      logical nflg ! returned from nucleation, says whether nucleation occurred or not
      real*8 mcond,mcond1    !mass to condense [kg]
      real*8 tol      !tolerance
      real*8 eps      !small number
      real*8 sinkfrac(ibins) !fraction of condensation sink coming from bin k
      real*8 totmass  !the total mass of H2SO4 generated during the timestep
      real*8 tmass
      real*8 CSch     !fractional change in condensation sink
      real*8 CSch_tol !tolerance in change in condensation sink
      real*8 addt     !adaptive timestep time
      real*8 time_rem !time remaining
      integer num_iter !number of iteration
      real*8 sumH2SO4 !used for finding average H2SO4 conc over timestep
      integer iter ! number of iteration
      real*8 rnuc !critical radius [nm]
      real*8 gasConc  !gas concentration [kg]
      real*8 mass_change !change in mass during nucleation.f
      real*8 total_nh4_1,total_nh4_2
      real*8 min_tstep !minimum timestep [s]
      integer nuc_bin           ! the nucleation bin
      real*8 sumfn, sumfn1 ! used for getting average nucleation rates
      real*8 mcond_soa
      parameter(eps=1E-30)
      parameter(CSch_tol=0.01)
      parameter(min_tstep=1.0d0)


      dt = dble(dti)


C Initialize values of Nkf, Mkf, Gcf, and time
      do j=1,icomp-1
         Gc1(j)=Gci(j)
      enddo
      do k=1,ibins
         Nk1(k)=Nki(k)
         Nknuc(k)=Nki(k)
         Nkcond(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
            Mknuc(k,j)=Mki(k,j)
            Mkcond(k,j)=Mki(k,j)
         enddo
      enddo


C     Get initial condensation sink
      CS1 = 0.d0
      call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac)

C     Get initial H2SO4 concentration guess (assuming no nucleation)
C     Make sure that H2SO4 concentration doesn't exceed the amount generated
C     during that timestep (this will happen when the condensation sink is very low)

C     get the steady state H2SO4 concentration
      call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,lev)
      Gc1(srtso4) = gasConc
      addt = min_tstep

      totmass = H2SO4rate*addt*96.d0/98.d0

C     Get change size distribution due to nucleation with initial guess  
      call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,nuc_bin,
     &     addt,lev)          

! for so4
      mass_change = 0.d0

      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
      enddo
      mcond = totmass-mass_change ! mass of h2so4 to condense

      if (mcond.lt.0.d0)then
         tmass = 0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tmass = tmass + Mk2(k,j)
            enddo
         enddo

         if (abs(mcond).gt.totmass*1.0d-8) then
            if (-mcond.lt.Mk2(nuc_bin,srtso4)) then

               tmass = 0.d0
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(nuc_bin,j)
               enddo
               Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
               Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
               mcond = 0.d0
            else
               print*,'budget fudge 2 in cond_nuc.f'
               do k=2,ibins
                  Nk2(k) = Nk1(k)
                  Mk2(k,srtso4) = Mk1(k,srtso4)
               enddo
               Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
               Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
               mcond = 0.d0        

            endif
         else
            mcond = 0.d0
         endif
      endif
      
      mass_change = 0.d0

      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtocil)-Mk1(k,srtocil))
      enddo

      mcond_soa = SOArate*addt-mass_change ! mass of soa to condense
      tmass = 0.d0
      do k=1,ibins-1
         do j=1,icomp-idiag
            tmass = tmass + Mk2(k,j)
         enddo
      enddo
      if (mcond_soa.gt.tmass)then ! limit soa
         mcond_soa=tmass
      endif


C     Get guess for condensation
      call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3)

      if(mcond_soa.eq.0) goto 17
      do k=1,ibins
         Nk2(k)=Nk3(k)
         do j=1,icomp
            Mk2(k,j)=Mk3(k,j)
         enddo
      enddo
      call ezcond(Nk2,Mk2,mcond_soa,srtocil,Nk3,Mk3)

 17   continue

      Gc3(srtnh4) = Gc1(srtnh4)   
      if (Gc1(srtnh4) .lt. 0.0) then
         print*, 'less than zero Gc', Gc1(srtnh4)
      endif



      call eznh3eqm(Gc3,Mk3)
      call ezwatereqm(Mk3)


! check to see how much condensation sink changed
      call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac) !problem here. 


      if(CS2.eq.CS1.or.CS1.le.0.) then
!MODELE-TOMAS : Somehow CSch becomes NaN, which mean no condensation and nucleation occurs.
! I take whole timestep in this case. Need to look this further in future. 
        time_rem=dt 
        addt=dt  !This will let it go only one. 

      else
         
         CSch = abs(CS2 - CS1)/CS1    
       
c      if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
         ! get starting adaptive timestep to not allow condensationk sink
         ! to change that much
         addt = addt*CSch_tol/CSch/2
         addt = min(addt,dt)
         addt = max(addt,min_tstep)
         time_rem = dt ! time remaining
         
      endif
         num_iter = 0
         sumH2SO4=0.d0
         sumfn = 0.d0
         sumfn1 = 0.d0
         ! do adaptive timesteps
         do while (time_rem .gt. 0.d0)
            num_iter = num_iter + 1
C     get the steady state H2SO4 concentration
            if (num_iter.gt.1)then ! no need to recalculate for first step
               call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,lev)
               Gc1(srtso4) = gasConc
            endif

            sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
            totmass = H2SO4rate*addt*96.d0/98.d0

            call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,
     &           nuc_bin,addt,lev) 
            sumfn = sumfn + fn*addt
            sumfn1 = sumfn1 + fn1*addt

            mass_change = 0.d0
            do k=1,ibins
               mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
            enddo
            mcond = totmass-mass_change ! mass of h2so4 to condense

            if (mcond.lt.0.d0)then
               tmass = 0.d0
               do k=1,ibins
                  do j=1,icomp-idiag
                     tmass = tmass + Mk2(k,j)
                  enddo
               enddo
               if (abs(mcond).gt.totmass*1.0D-8) then
                  if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
c                     if (CS1.gt.1.0D-5)then
c                        print*,'budget fudge 1 in cond_nuc.f'
c                     endif
                     tmass = 0.d0
                     do j=1,icomp-idiag
                        tmass = tmass + Mk2(nuc_bin,j)
                     enddo
                     Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                     Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                     mcond = 0.d0
                  else
                     print*,'budget fudge 2 in cond_nuc.f'
                     do k=2,ibins
                        Nk2(k) = Nk1(k)
                        Mk2(k,srtso4) = Mk1(k,srtso4)
                     enddo
                     Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
                     Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
                     mcond = 0.d0 
                  endif
               else
                  mcond = 0.d0
               endif
            endif

            mass_change = 0.d0

            do k=1,ibins
               mass_change=mass_change+(Mk2(k,srtocil)-Mk1(k,srtocil))
            enddo

            mcond_soa = SOArate*addt-mass_change ! mass of soa to condense
            tmass = 0.d0
            do k=1,ibins-1
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(k,j)
               enddo
            enddo
            if (mcond_soa.gt.tmass)then  ! limit soa
               mcond_soa=tmass
            endif
            
            do k=1,ibins
               Nknuc(k) = Nknuc(k)+Nk2(k)-Nk1(k)
               do j=1,icomp-idiag
                  Mknuc(k,j)=Mknuc(k,j)+Mk2(k,j)-Mk1(k,j)
               enddo
            enddo          
           
            call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3)
            do k=1,ibins
               Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
               do j=1,icomp-idiag
                  Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
               enddo
            enddo

            if(mcond_soa.eq.0) goto 19
            do k=1,ibins
               Nk2(k)=Nk3(k)
               do j=1,icomp
                  Mk2(k,j)=Mk3(k,j)
               enddo
            enddo

            call ezcond(Nk2,Mk2,mcond_soa,srtocil,Nk3,Mk3)
            do k=1,ibins
               Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
               do j=1,icomp-idiag
                  Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
               enddo
            enddO
 19         continue

            Gc3(srtnh4) = Gc1(srtnh4)

            call eznh3eqm(Gc3,Mk3)
            call ezwatereqm(Mk3)
            
! check to see how much condensation sink changed
            call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac)  

            time_rem = time_rem - addt
            if (time_rem .gt. 0.d0) then
               
!MODELE-TOMAS : newly added to prevent NaN for CSch. 
!Somehow, mcond >0 but not nucleation and no condensation occurs. 
!same size distribution causes problems. 

               if(CS2.eq.CS1.or.CS1.le.0.) then
                  addt = time_rem
                  addt = max(addt,min_tstep)
               else
                  CSch = abs(CS2 - CS1)/CS1 !TOMAS - CS1 = ZERO?
                  
                  
                  addt = min(addt*CSch_tol/CSch,addt*1.5d0) ! allow adaptive timestep to change
                  addt = min(addt,time_rem) ! allow adaptive timestep to change
                  addt = max(addt,min_tstep)
               endif
Cjrp               endif
               CS1 = CS2
               Gc1(srtnh4)=Gc3(srtnh4)
               do k=1,ibins
                  Nk1(k)=Nk3(k)
                  do j=1,icomp
                     Mk1(k,j)=Mk3(k,j)
                  enddo
               enddo         
            endif
         enddo
         Gcf(srtso4)=sumH2SO4/dt
         fnavg = sumfn/dt
         fn1avg = sumfn1/dt

      do k=1,ibins
         Nkf(k)=Nk3(k)
         do j=1,icomp
            Mkf(k,j)=Mk3(k,j)
         enddo
      enddo      

      Gcf(srtnh4)=Gc3(srtnh4)

      return
      end


!@sum getCondSink_kerm :  calculates the condensation sink (first order loss
!@+    rate of condensing gases) from the aerosol size distribution.
!@+    This is the cond sink in kerminen et al 2004 Parameterization for 
!@+    new particle formation AS&T Eqn 6.
!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   Nk(ibins) - number of particles per size bin in grid cell
!@var   Nnuc - number of particles per size bin in grid cell
!@var   Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
!@var   Mk(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   spec - number of the species we are finding the condensation sink for

!   Output 
!   =================
!@var   CS - condensation sink [s^-1]
!@var   sinkfrac(ibins) - fraction of condensation sink from a bin

      SUBROUTINE getCondSink_kerm(Nko,Mko,CS,Dpmean,Dp1,dens1)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi,gasc  
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 Nko(ibins), Mko(ibins, icomp)
      real*8 CS
      real*8 Dpmean             ! the number mean diameter [m]
      real*8 Dp1                ! the size of the first size bin [m]
      real*8 dens1              ! the density of the first size bin [kg/m3]

      integer i,j,k,c           ! counters
      real*8 mu                  !viscosity of air (kg/m s)
      real*8 mfp                 !mean free path of air molecule (m)

      real Di                   !diffusivity of gas in air (m2/s)
      real*8 Neps               !tolerance for number
      real density              !density [kg m^-3]
      real*8 mp                 !mass per particle [kg]
      real*8 Dpk(ibins)         !diameter of particle [m]
      real*8 Kn                 !Knudson number
      real*8 beta(ibins)        !non-continuum correction factor
      real*8 Mktot              !total mass in bin [kg]
      real*8 Dtot,Ntot          ! used on getting the number mean diameter
      
      real mso4, mh2o, mno3, mnh4 !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      real gasdiff
      external aerodens         !!, gasdiff
      
      parameter(Neps=1.0d10)
      real*8 alpha(icomp)       ! accomodation coef  
      data alpha/0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65/
      real Sv(icomp)            !parameter used for estimating diffusivity
      data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88,42.88,
     &     42.88/
      
      
C     get some parameters  
!!      mu=2.5277e-7*temp**0.75302
!!      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6
!for old debugging      mfp=2.0*mu/(pres*sqrt(8.0*0.6589/(pi*gasc*temp)))  !S&P eqn 8.6


cyhl the following should be used instead of above two lines!
      Di=gasdiff(temp,pres,98.0,Sv(srtso4)) ! Di is diffusivity of condensing gas in air [m2/s]
      mfp=2.d0*Di/sqrt(8.0*gasc*temp/(pi*(molwt(srtso4)+2.)/1000.))   !0.098 for H2SO4 m.w. [kg/mol]
!     the denominator is mean speed. sqrt(8.0*R*temp/(pi*molwt(srtso4)/1000.)) !S&P2 eqn 9.2 ms=mean speed [m/s]
         
cyhl but this needs some tuning before it actually uses! 

c      Di=gasdiff(temp,pres,98.0,Sv(srtso4))
c      print*,'Di',Di

C     get size dependent values
      CS = 0.d0
      Ntot = 0.d0
      Dtot = 0.d0
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
               Mktot=Mktot+Mko(k,j)
            enddo

            mso4=Mko(k,srtso4) 
            mnacl=Mko(k,srtna)
            mno3=0.e0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mko(k,srtnh4) !0.1875*mso4    !assume ammonium bisulfate
            mecob=Mko(k,srtecob)
            mecil=Mko(k,srtecil)
            mocil=Mko(k,srtocil)
            mocob=Mko(k,srtocob)
            mdust=Mko(k,srtdust)          
            mh2o=Mko(k,srth2o)   
            
            density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *           ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate                  
            mp=Mktot/Nko(k)
          else
!nothing in this bin - set to "typical value"
            density=1500.
            mp=sqrt(xk(k+1)*xk(k))
          endif
          Dpk(k)=((mp/density)*(6./pi))**(0.333)
          Kn=2.0*mfp/Dpk(k)     !S&P eqn 11.35 (text)
          CS=CS+0.5d0*(Dpk(k)*Nko(k)/(boxvol*1.0D-6)*(1+Kn))/
     &         (1.d0+0.377d0*Kn+1.33d0*Kn*(1+Kn))
          Ntot = Ntot + Nko(k)
          Dtot = Dtot + Nko(k)*Dpk(k)
          if (k.eq.1)then
            Dp1=Dpk(k)
            dens1 = density
          endif
        enddo      
        
        if (Ntot.gt.1D15)then
          Dpmean = Dtot/Ntot
        else
          Dpmean = 150.d0
        endif
        
        return
        end
      
      
!@sum getCondSink :  calculates the condensation sink (first order loss
!@+   rate of condensing gases) from the aerosol size distribution.
!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   Nk(ibins) - number of particles per size bin in grid cell
!@var   Nnuc - number of particles per size bin in grid cell
!@var   Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
!@var   Mk(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   spec - number of the species we are finding the condensation sink for

!   Output
!   =================
!@var   CS - condensation sink [s^-1]
!@var   sinkfrac(ibins) - fraction of condensation sink from a bin

      SUBROUTINE getCondSink(Nko,Mko,spec,CS,sinkfrac)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi,gasc  
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 Nko(ibins), Mko(ibins, icomp)
      real*8 CS, sinkfrac(ibins)
      integer spec

      integer i,j,k,c           ! counters
      real*8 mu                 !viscosity of air (kg/m s)
      real*8 mfp                !mean free path of air molecule (m)
      real Di                   !diffusivity of gas in air (m2/s), and molecular weight (kg/mol)
      real*8 Neps               !tolerance for number
      real*8 density            !density [kg m^-3]
      real mw
      real*8 mp                 !mass per particle [kg]
      real*8 Dpk(ibins)         !diameter of particle [m]
      real*8 Kn                 !Knudson number
      real*8 beta(ibins)        !non-continuum correction factor
      real*8 Mktot              !total mass in bin [kg]

      real mso4, mh2o, mno3, mnh4 !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      real gasdiff
      external aerodens         !, gasdiff

      parameter(Neps=1.0d10)
      real*8 alpha(icomp)       ! accomodation coef  
      data alpha/0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65/
      real Sv(icomp)            !parameter used for estimating diffusivity
      data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88,42.88,
     &         42.88/


cyhl the following two lines are commented out (not accurate mfp)
C     get some parameters
c$$$      mu=2.5277e-7*temp**0.75302 
c$$$      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6 !bug?


cyhl mfp is now for the condensing gas in the air   10/17/2010
      Di=gasdiff(temp,pres,98.0,Sv(spec))  ! YHL(10/17/2010) - this is not accurate for SOA, but leave this for now. 
      mfp=2.d0*Di/sqrt(8.0*gasc*temp/(pi*(molwt(srtso4)+2.)/1000.))   
!     molwt(srtso4) is 96, so adding 2 will make 98. 
!     the denominator is mean speed. sqrt(8.0*R*temp/(pi*molwt(srtso4)/1000.)) !S&P2 eqn 9.2 ms=mean speed [m/s]

C     get size dependent values
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
                  Mktot=Mktot+Mko(k,j)
            enddo

            mso4=Mko(k,srtso4) 
            mnacl=Mko(k,srtna)
            mno3=0.e0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mko(k,srtnh4)!0.1875*mso4    !assume ammonium bisulfate
            mecob=Mko(k,srtecob)
            mecil=Mko(k,srtecil)
            mocil=Mko(k,srtocil)
            mocob=Mko(k,srtocob)
            mdust=Mko(k,srtdust)          
            mh2o=Mko(k,srth2o)   
            
            density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *           ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate                  
            mp=Mktot/Nko(k)
          else
!nothing in this bin - set to "typical value"
            density=1500.d0
            mp=sqrt(xk(k+1)*xk(k))
          endif
          Dpk(k)=((mp/density)*(6.d0/pi))**(1.d0/3.d0)
          Kn=2.0*mfp/Dpk(k)     !S&P eqn 11.35 (text)
          beta(k)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(spec)) !S&P eqn 11.35
        enddo      
C     get condensation sink
        CS = 0.d0
        surf_area = 0.d0
        do k=1,ibins
          CS = CS + Dpk(k)*Nko(k)*beta(k)
          surf_area = surf_area+Nko(k)*pi*(Dpk(k)*1.0D6)**2
        enddo
        do k=1,ibins
          sinkfrac(k) = Dpk(k)*Nko(k)*beta(k)/CS
        enddo
c     CS = 2.d0*pi*Di*CS/(boxvol*1D-6)      
        CS = 2.d0*pi*dble(Di)*CS/(boxvol*1D-6)
        surf_area = surf_area/boxvol
        return
        end
      

!@sum getCoagLoss  :  calculates the first order loss rate of particles
!@+   of a given size with respect to coagulation.
!@auth   Jeff Pierce, April 2007

C-----INPUTS------------------------------------------------------------

!@var   d1: diameter in [m] of particle

C-----OUTPUTS-----------------------------------------------------------

!@var   ltc: first order loss rate with respect to coagulation [s-1]
!@var   sinkfrac: fraction of sink from each bin

      SUBROUTINE getCoagLoss(d1,ltc,Nko,Mko,sinkfrac)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi,gasc  
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 d1       ! diameter of the particle [m]
      real*8 ltc    ! first order loss rate [s-1]
      real*8 Nko(ibins), Mko(ibins, icomp)
      real*8 sinkfrac(ibins)

      integer i,k,j
      real*8 density  ! density of particles [kg/m3]
      real*8 density1 ! density of particles in first bin
      real*8 MW, kB
      real*8 mu       !viscosity of air (kg/m s)
      real*8 mfp      !mean free path of air molecule (m)
      real*8 Kn       !Knudsen number of particle
      real*8 Mktot    !total mass of aerosol
      real*8 kij(ibins)
      real*8 Dpk(ibins) !diameter (m) of particles in bin k
      real*8 Dk(ibins),Dk1 !Diffusivity (m2/s) of bin k particles
      real*8 ck(ibins),ck1 !Mean velocity (m/2) of bin k particles
      real*8 neps
      real*8 meps
      real*8 mp       ! mass of the particle (kg)
      real*8 beta     !correction for coagulation coeff.
      real*8 kij_self !coagulation coefficient for self coagulation

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      external aerodens

      parameter (kB= 1.38E-23) !pi and gas constant (J/mol K)
      parameter (neps=1E8, meps=1E-8)

      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gasc*temp)))  !S&P eqn 8.6

C Calculate particle sizes and diffusivities
      do k=1,ibins
         Mktot = 0.d0
         do j=1, icomp
            Mktot=Mktot+Mko(k,j)
         enddo
         Mktot=Mktot+2.d0*Mko(k,srtso4)/96.d0-Mko(k,srtnh4)/18.d0 ! account for h+

         if (Mktot.gt.meps)then
            mso4=Mko(k,srtso4) 
            mnacl=Mko(k,srtna)
            mno3=0.e0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mko(k,srtnh4)!0.1875*mso4    !assume ammonium bisulfate
            mecob=Mko(k,srtecob)
            mecil=Mko(k,srtecil)
            mocil=Mko(k,srtocil)
            mocob=Mko(k,srtocob)
            mdust=Mko(k,srtdust)          
            mh2o=Mko(k,srth2o)   
            
         density=aerodens(mso4,mno3,mnh4 !mno3 taken off!
     *        ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
         else
            density = 1400.
         endif
         if(Nko(k).gt.neps .and. Mktot.gt.meps)then
            mp=Mktot/Nko(k)
         else
            mp=sqrt(xk(k)*xk(k+1))
         endif
         if (k.eq.1) density1 = density
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
         ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
      enddo
      
      Kn = 2.0*mfp/d1
      Dk1 = kB*temp/(3.0*pi*mu*d1)             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
      mp = 4.d0/3.d0*pi*(d1/2.d0)**3.d0*density1
      ck1=sqrt(8.0*kB*temp/(pi*mp)) !S&P Table 12.1      

      do i=1,ibins
         Kn=4.0*(Dk(i)+Dk1)          
     &        /(sqrt(ck(i)**2+ck1**2)*(Dpk(i)+d1)) !S&P eqn 12.51
         beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
                                !This is S&P eqn 12.46 with non-continuum correction, beta
         kij(i)=2.0*pi*(Dpk(i)+d1)*(Dk(i)+Dk1)*beta
         kij(i)=kij(i)*1.0e6/boxvol !normalize by grid cell volume
c         kij(i)=kij(i)*1.0e6 !cm3/s
      enddo
Cjrp
Cjrp      !self coagulation
Cjrp      Kn=4.0*(Dk1+Dk1)          
Cjrp     &     /(sqrt(ck1**2+ck1**2)*(d1+d1)) !S&P eqn 12.51
Cjrp      beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn)) !S&P eqn 12.50
Cjrp                                !This is S&P eqn 12.46 with non-continuum correction, beta
Cjrp      kij_self=2.0*pi*(d1+d1)*(Dk1+Dk1)*beta
Cjrp      print*,'kij_self',kij_self*1E6
Cjrp      kij_self=kij_self*1.0e6/boxvol !normalize by grid cell volume      
         
      ltc = 0.d0

      do i=1,ibins
         ltc = ltc + kij(i)*Nko(i)
      enddo
      do i=1,ibins
         sinkfrac(i) = kij(i)*Nko(i)/ltc
      enddo

      RETURN
      END


!@sum NH3_GISStoTOMAS   :  puts ammonia to the particle phase until 
!@+   there is 2 moles of ammonium per mole of sulfate and the remainder
!@+   of ammonia is left in the gas phase.
!@auth   Jeff Pierce, April 2007

      SUBROUTINE NH3_GISStoTOMAS(giss_nh3g,giss_nh4a,Gce,Mke)

      USE TOMAS_AEROSOL

      IMPLICIT NONE

      real*8 Gce(icomp-1)
      real*8 Mke(ibins,icomp)
      real*8 giss_nh3g, giss_nh4a

      integer k
      real*8 tot_nh3  !total kmoles of ammonia
      real*8 tot_so4  !total kmoles of so4
      real*8 sfrac    !fraction of sulfate that is in that bin

      ! get the total number of kmol nh3
      tot_nh3 = giss_nh3g/17.d0 + giss_nh4a/18.d0

      ! get the total number of kmol so4
      tot_so4=0.d0
      do k=1,ibins
         tot_so4 = tot_so4 + Mke(k,srtso4)/96.d0
      enddo

      ! see if there is free ammonia
      if (tot_nh3/2.d0.lt.tot_so4)then ! no free ammonia
        Gce(srtnh4) = 0.d0      ! no gas phase ammonia
        do k=1,ibins
          sfrac = Mke(k,srtso4)/96.d0/tot_so4
          Mke(k,srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
        enddo
      else                      ! free ammonia
c     Mnuce(srtnh4) = Mnuce(srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
        do k=1,ibins
          Mke(k,srtnh4) = Mke(k,srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
        enddo
        Gce(srtnh4) = (tot_nh3 - tot_so4*2.d0)*17.d0 ! put whats left over in the gas phase
      endif
      

      RETURN
      END


!@sum getNucRate  :  calls the Vehkamaki 2002 and Napari 2002 nucleation
!@+   parameterizations and gets the binary and ternary nucleation rates.
!@auth   Jeff Pierce, April 2007

!   Initial values of
!   =================

!@var   Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!@+                  gas phase except water

!-----OUTPUTS-----------------------------------------------------------

!@var   fn - nucleation rate [# cm-3 s-1]
!@var   rnuc - radius of nuclei [nm]
!@var   nflg - says if nucleation happend

      SUBROUTINE getNucRate(Gci,fn,mnuc,nflg,l)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi  
      USE TRACER_COM, only : xk

      IMPLICIT NONE

      integer j,i,k,l
      real*8 Gci(icomp-1)
      real*8 fn       ! nucleation rate to first bin cm-3 s-1
      real*8 mnuc     !mass of nucleating particle [kg]
      logical nflg

      real*8 nh3ppt   ! gas phase ammonia in pptv
      real*8 h2so4    ! gas phase h2so4 in molec cc-1
      real*8 gtime    ! time to grow to first size bin [s]
      real*8 ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      real*8 Mktot    ! total mass in bin
      real*8 neps
      real*8 meps
      real*8 density  ! density of particle [kg/m3]
      real*8 frac     ! fraction of particles growing into first size bin
      real*8 d1,d2    ! diameters of particles [m]
      real*8 mp       ! mass of particle [kg]
      real*8 mold     ! saved mass in first bin
      real*8 rnuc     ! critical nucleation radius [nm]
      real*8 sinkfrac(ibins) ! fraction of loss to different size bins
      real*8 nadd     ! number to add
      real*8 CS       ! kerminan condensation sink [m-2]
      real*8 Dpmean   ! the number wet mean diameter of the existing aerosol
      real*8 Dp1      ! the wet diameter of bin 1
      real*8 dens1    ! density in bin 1 [kg m-3]
      real*8 GR       ! growth rate [nm hr-1]
      real*8 gamma,eta ! used in kerminen 2004 parameterzation
      real*8 drymass,wetmass,WR
      real*8 fn_c     ! barrierless nucleation rate
      real*8 h1,h2,h3,h4,h5,h6
      parameter (neps=1E8, meps=1E-8)


      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             pres/101325.*273./temp ! corrected for pressure (because this should be concentration)

      fn = 0.d0
      rnuc = 0.d0

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then

            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc 
            if (ion_nuc.eq.1.and.ionrate.ge.1.d0) then
               call ion_nucl(h2so4,surf_area,temp,ionrate,rh,h1,h2,h3,
     &              h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
            nflg=.true.
          elseif (bin_nuc.eq.1) then

            if((actv_nuc.eq.1).and.(l.le.3))then
              call bl_nucl(h2so4,fn,rnuc)
            else
              call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
            endif
              if ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
                call ion_nucl(h2so4,surf_area,temp,ionrate,rh,h1,h2,h3,
     &               h4,h5,h6)
              else
                h1=0.d0
              endif
              
              if (h1.gt.fn)then
                fn=h1
                rnuc=h5
              endif
              if (fn.gt.1.0d-6)then
                nflg=.true.
              else
                fn = 0.d0
                nflg=.false.
              endif
         elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
            call ion_nucl(h2so4,surf_area,temp,ionrate,rh,h1,h2,h3,
     &           h4,h5,h6)
            fn=h1
            rnuc=h5
         else
            nflg=.false.
         endif
         call cf_nucl(temp,rh,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max for ternary
         fn = min(fn,fn_c)    
      else
         nflg=.false.
      endif

      if (fn.gt.0.d0) then
         call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1)
         d1 = rnuc*2.d0*1D-9
         drymass = 0.d0
         do j=1,icomp-idiag
            drymass = drymass + Mk(1,j)
         enddo
         wetmass = 0.d0
         do j=1,icomp
            wetmass = wetmass + Mk(1,j)
         enddo
         WR = wetmass/drymass
         
         call getGrowthTime(d1,Dp1,Gci(srtso4)*(1.d0+soa_amp)*WR,temp,
     &        boxvol,dens1,gtime)
         GR = (Dp1-d1)*1D9/gtime*3600.d0 ! growth rate, nm hr-1
         
         gamma = 0.23d0*(d1*1.0d9)**(0.2d0)*(Dp1*1.0d9/3.d0)**0.075d0*
     &        (Dpmean*1.0d9/150.d0)**0.048d0*(dens1*1.0d-3)**
     &        (-0.33d0)*(temp/293.d0) ! equation 5 in kerminen
         eta = gamma*CS/GR

         fn = fn*exp(eta/(Dp1*1.0D9)-eta/(d1*1.0D9))

         mnuc = sqrt(xk(1)*xk(2))
      endif

      return
      end
      

!@sum getH2SO4conc :  uses newtons method to solve for the steady state 
!@+   H2SO4 concentration when nucleation is occuring.

!@+   It solves for H2SO4 in 0 = P - CS*H2SO4 - M(H2SO4)

!@+   where P is the production rate of H2SO4, CS is the condensation sink
!@+   and M(H2SO4) is the loss of mass towards making new particles.

!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   H2SO4rate - H2SO4 generation rate [kg box-1 s-1]
!@var   CS - condensation sink [s-1]
!@var   NH3conc - ammonium in box [kg box-1]
!@var   prev - logical flag saying if a previous guess should be used or not
!@var   gasConc_prev - the previous guess [kg/box] (not used if prev is false)

!-----OUTPUTS-----------------------------------------------------------

!@var   gasConc - gas H2SO4 [kg/box]

      SUBROUTINE getH2SO4conc(H2SO4rate,CS,NH3conc,gasConc,level)

      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 H2SO4rate
      real*8 CS
      real*8 NH3conc
      real*8 gasConc
      logical prev
      real*8 gasConc_prev

      integer, intent(in) :: level  ! vertical layer 
      integer i,j,k,c           ! counters
      real*8 fn, rnuc           ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      real*8 mnuc, mnuc1        ! mass of nucleated particle [kg]
      real*8 fn1, rnuc1         ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      real*8 res                ! d[H2SO4]/dt, need to find the solution where res = 0
      real*8 massnuc            ! mass being removed by nucleation [kg s-1 box-1]
      real*8 gasConc1           ! perturbed gasConc
      real*8 gasConc_hi, gasConc_lo
      real*8 res1               ! perturbed res
      real*8 res_new            ! new guess for res
      real*8 dresdgasConc       ! derivative for newtons method
      real*8 Gci(icomp-1)       !array to carry gas concentrations
      logical nflg              !says if nucleation occured

      real*8 H2SO4min           !minimum H2SO4 concentration in parameterizations (molec/cm3)
!     real*8 pi
      integer iter,iter1
      real*8 CSeps              ! low limit for CS
      real*8 max_H2SO4conc      !maximum H2SO4 concentration in parameterizations (kg/box)
      real*8 nh3ppt             !ammonia concentration in ppt
     
      parameter(H2SO4min=1.D4) !molecules cm-3
      parameter(CSeps=1.0d-20)


      do i=1,icomp-1
         Gci(i)=0.d0
      enddo
      Gci(srtnh4)=NH3conc
! make sure CS doesn't equal zero
! some specific stuff for napari vs. vehk
      if ((bin_nuc.eq.1).or.(tern_nuc.eq.1))then
         nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             pres/101325.*273./temp ! corrected for pressure (because this should be concentration)
         if ((nh3ppt.gt.1.0d0).and.(tern_nuc.eq.1))then
            max_H2SO4conc=1.0D9*boxvol/1000.d0*98.d0/6.022d23
         elseif (bin_nuc.eq.1)then
            max_H2SO4conc=1.0D11*boxvol/1000.d0*98.d0/6.022d23
         else
            max_H2SO4conc = 1.0D100
         endif
      else
         max_H2SO4conc = 1.0D100
      endif
      
C     Checks for when condensation sink is very small
      if (CS.gt.CSeps) then
         gasConc = H2SO4rate/CS
      else
         if ((bin_nuc.eq.1).or.(tern_nuc.eq.1)) then
            gasConc = max_H2SO4conc
         else
            print*,'condesation sink too small in getH2SO4conc.f'
            call stop_model('CS too small getH2SO4',255)
         endif
      endif
      
      gasConc = min(gasConc,max_H2SO4conc)
      Gci(srtso4) = gasConc
      call getNucRate(Gci,fn,mnuc,nflg,level)
      
      if (fn.gt.0.d0) then      ! nucleation occured
         gasConc_lo = H2SO4min*boxvol/(1000.d0/98.d0*6.022d23) !convert to kg/box
         
C     Test to see if gasConc_lo gives a res < 0 (this means ANY nucleation is too high)
         Gci(srtso4) = gasConc_lo*1.000001d0
         call getNucRate(Gci,fn1,mnuc1,nflg,level)
         if (fn1.gt.0.d0) then
            massnuc = mnuc1*fn1*(1.d0/(1.d0+soa_amp))*boxvol*98.d0/96.d0
c     massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1350.*fn1*boxvol*
c     massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1800.*fn1*boxvol*
c     &           98.d0/96.d0

            res = H2SO4rate - CS*gasConc_lo - massnuc
            if (res.lt.0.d0) then ! any nucleation too high

               gasConc = gasConc_lo*1.000001 ! have nucleation occur and fix mass balance after
               return
            endif
         endif
         
         gasConc_hi = gasConc   ! we know this must be the upper limit (since no nucleation)
                                !take density of nucleated particle to be 1350 kg/m3
         massnuc = mnuc*fn*(1.d0/(1.d0+soa_amp))*boxvol*98.d0/96.d0
         res = H2SO4rate - CS*gasConc - massnuc
         
                                ! check to make sure that we can get solution
         if (res.gt.0.d0) then
            return
         endif
         
         iter = 0

         do while ((abs(res/H2SO4rate).gt.1.D-4).and.(iter.lt.40))
            iter = iter+1
            if (res .lt. 0.d0) then ! H2SO4 concentration too high, must reduce
               gasConc_hi = gasConc ! old guess is new upper bound
            elseif (res .gt. 0.d0) then ! H2SO4 concentration too low, must increase
               gasConc_lo = gasConc ! old guess is new lower bound
            endif

            gasConc = sqrt(gasConc_hi*gasConc_lo) ! take new guess as logmean
            Gci(srtso4) = gasConc
            call getNucRate(Gci,fn,mnuc,nflg,level)
            massnuc = mnuc*fn*(1.d0/(1.d0+soa_amp))*boxvol*98.d0/96.d0
            res = H2SO4rate - CS*gasConc - massnuc

            if (iter.eq.30.and.CS.gt.1.0D-5)then
               print*,'getH2SO4conc iter break'
               print*,'H2SO4rate',H2SO4rate,'CS',CS
               print*,'gasConc',gasConc,'massnuc',massnuc
               print*,'res/H2SO4rate',res/H2SO4rate
            endif
         enddo
         
      else                      ! nucleation didn't occur
      endif
      
      return
      end
      



!@sum getGrowthTime  :  calculates the time it takes for a particle to grow
!@+   from one size to the next by condensation of sulfuric acid (and
!@+   associated NH3 and water) onto particles.

!@+   This subroutine  assumes that the growth happens entirely in the kinetic
!@+   regine such that the dDp/dt is not size dependent.  The time for growth 
!@+   to the first size bin may then be approximated by the time for growth via
!@+   sulfuric acid (not including nh4 and water) to the size of the first size bin
!@+   (not including nh4 and water).

!@auth   Jeff Pierce, April 2007

C-----INPUTS------------------------------------------------------------

!@var   d1: intial diameter [m]
!@var   d2: final diameter [m]
!@var     h2so4: h2so4 ammount [kg]
!@var     temp: temperature [K]
!@var     boxvol: box volume [cm3]

C-----OUTPUTS-----------------------------------------------------------

!@var  gtime: the time it takes the particle to grow to first size bin [s]

      SUBROUTINE getGrowthTime(d1,d2,h2so4,temp,boxvol,density,gtime)

      USE CONSTANT,   only:  pi,gasc  
      IMPLICIT NONE

      real*8 d1,d2    ! initial and final diameters [m]
      real*8 h2so4    ! h2so4 ammount [kg]
      real*8 temp     ! temperature [K]
      real*8 boxvol   ! box volume [cm3]
      real*8 gtime    ! the time it will take the particle to grow 
                                ! to first size bin [s]
      real*8 density  ! density of particles in first bin [kg/m3]

      real*8 MW
      real*8 csulf    ! concentration of sulf acid [kmol/m3]
      real*8 mspeed   ! mean speed of molecules [m/s]
      real*8 alpha    ! accomidation coef

      parameter(MW=98.d0) ! density [kg/m3], mol wgt sulf [kg/kmol]
      parameter(alpha=0.65)


      csulf = h2so4/MW/(boxvol*1d-6) ! SA conc. [kmol/m3]
      mspeed = sqrt(8.d0*gasc*temp*1000.d0/(pi*MW))

C   Kinetic regime expression (S&P 11.25) solved for T
      gtime = (d2-d1)/(4.d0*MW/density*mspeed*alpha*csulf)

      RETURN
      END




!@sum nucleation :  calls the Vehkamaki 2002 and Napari 2002 nucleation
!@+   parameterizations and gets the binary and ternary nucleation rates.
!@+   The number of particles added to the first size bin is calculated
!@+   by comparing the growth rate of the new particles to the coagulation
!@+   sink.

!@auth   Jeff Pierce, April 2007

C-----INPUTS------------------------------------------------------------

!   Initial values of
!   =================

!@var   Nki(ibins) - number of particles per size bin in grid cell
!@var   Mki(ibins, icomp) - mass of a given species per size bin/grid cell
!@var   Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!@+                  gas phase except water
!@var   dt - total model time step to be taken (s)

!-----OUTPUTS-----------------------------------------------------------

!@var   Nkf, Mkf, Gcf - same as above, but final values
!@var   fn, fn1

      SUBROUTINE nucleation(Nki,Mki,Gci,Nkf,Mkf,Gcf,fn,fn1,totsulf,
     &     nuc_bin,dt,l)


      USE TOMAS_AEROSOL
      USE CONSTANT,   only:  pi 
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      integer j,i,k,l
      real*8 Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      real*8 Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      real*8 totsulf
      integer nuc_bin
      real*8 dt
      real*8 fn       ! nucleation rate of clusters cm-3 s-1
      real*8 fn1      ! formation rate of particles to first size bin cm-3 s-1

      real*8 nh3ppt   ! gas phase ammonia in pptv
      real*8 h2so4    ! gas phase h2so4 in molec cc-1
      real*8 rnuc     ! critical nucleation radius [nm]
      real*8 gtime    ! time to grow to first size bin [s]
      real*8 ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      real*8 Mktot    ! total mass in bin
      real*8 neps
      real*8 meps
      real*8 density  ! density of particle [kg/m3]
      real*8 frac     ! fraction of particles growing into first size bin
      real*8 d1,d2    ! diameters of particles [m]
      real*8 mp       ! mass of particle [kg]
      real*8 Dpk(ibins) !diameter of particle [m]
      real*8 mold     ! saved mass in first bin
      real*8 mnuc     !mass of nucleation
      real*8 sinkfrac(ibins) ! fraction of loss to different size bins
      real*8 nadd     ! number to add
      real*8 CS       ! kerminan condensation sink [m-2]
      real*8 Dpmean   ! the number wet mean diameter of the existing aerosol
      real*8 Dp1      ! the wet diameter of bin 1
      real*8 dens1    ! density in bin 1 [kg m-3]
      real*8 GR       ! growth rate [nm hr-1]
      real*8 gamma,eta ! used in kerminen 2004 parameterzation
      real*8 drymass,wetmass,WR

      real mso4, mh2o, mno3, mnh4  !mass of each component (kg/grid box)
      real mecil,mecob,mocil,mocob
      real mdust,mnacl  
      real*8 aerodens
      external aerodens
      real*8 fn_c     ! barrierless nucleation rate
      real*8 h1,h2,h3,h4,h5,h6
      parameter (neps=1E8, meps=1E-8)


      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             pres/101325.*273./temp ! corrected for pressure (because this should be concentration)

      fn = 0.d0
      fn1 = 0.d0
      rnuc = 0.d0
      gtime = 0.d0

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then

            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
            if (ion_nuc.eq.1.and.ionrate.ge.1.d0) then
               call ion_nucl(h2so4,surf_area,temp,ionrate,rh,h1,h2,h3,
     &              h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
         elseif (bin_nuc.eq.1) then

            if((actv_nuc.eq.1).and.(l.le.3))then
              call bl_nucl(h2so4,fn,rnuc)
            else
              call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
            endif
            if ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
               call ion_nucl(h2so4,surf_area,temp,ionrate,rh,h1,h2,h3,
     &              h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
            if (fn.lt.1.0d-6)then
               fn = 0.d0
            endif
         elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
            call ion_nucl(h2so4,surf_area,temp,ionrate,rh,h1,h2,h3,
     &           h4,h5,h6)
            fn=h1
            rnuc=h5
         endif    
         call cf_nucl(temp,rh,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max
         fn = min(fn,fn_c) 
      endif

      d1 = rnuc*2.d0*1D-9
      do k=1,ibins
         if (Nki(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
               Mktot=Mktot+Mki(k,j)
            enddo

            mso4=Mki(k,srtso4) 
            mnacl=Mki(k,srtna)
            mno3=0.0
            if ((mso4+mno3) .lt. 1.e-8) mso4=1.e-8
            mnh4=Mki(k,srtnh4)!0.1875*mso4    !assume ammonium bisulfate
            mecob=Mki(k,srtecob)
            mecil=Mki(k,srtecil)
            mocil=Mki(k,srtocil)
            mocob=Mki(k,srtocob)
            mdust=Mki(k,srtdust)          
            mh2o=Mki(k,srth2o)   
            
            density=aerodens(mso4,mno3,mnh4 
     *           ,mnacl,mecil,mecob,mocil,mocob,mdust,mh2o) !assume bisulfate 
            
            mp=Mktot/Nki(k)
         else
                                !nothing in this bin - set to "typical value"
            density=1500.d0
            mp=sqrt(xk(k+1)*xk(k))
         endif
         Dpk(k)=((mp/density)*(6.d0/pi))**(1.d0/3.d0)
      enddo

C     if nucleation occured, see how many particles grow to join the first size
C     section
      if (fn.gt.0.d0.and.Dpk(1).gt.d1) then

         call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1)
         d1 = rnuc*2.d0*1D-9
         drymass = 0.d0
         do j=1,icomp-idiag
            drymass = drymass + Mk(1,j)
         enddo
         wetmass = 0.d0
         do j=1,icomp
            wetmass = wetmass + Mk(1,j)
         enddo
         WR = wetmass/drymass
         
         call getGrowthTime(d1,Dp1,Gci(srtso4)*(1.d0+soa_amp)*WR,temp,
     &        boxvol,dens1,gtime)
         GR = (Dp1-d1)*1D9/gtime*3600.d0 ! growth rate, nm hr-1
         
         gamma = 0.23d0*(d1*1.0d9)**(0.2d0)*(Dp1*1.0d9/3.d0)**0.075d0*
     &        (Dpmean*1.0d9/150.d0)**0.048d0*(dens1*1.0d-3)**
     &        (-0.33d0)*(temp/293.d0) ! equation 5 in kerminen
         eta = gamma*CS/GR

         fn1 = fn*exp(eta/(Dp1*1.0D9)-eta/(d1*1.0D9))

         mnuc = sqrt(xk(1)*xk(2))
 
         nadd = fn1
         
         nuc_bin = 1
         
         mold = Mki(nuc_bin,srtso4)
         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+nadd*mnuc*
     &        boxvol*dt/(1.d0+soa_amp)
         Mkf(nuc_bin,srtocil) = Mki(nuc_bin,srtocil)+nadd*mnuc*
     &        boxvol*dt*(1.d0-1.d0/(1.d0+soa_amp))
         Nkf(nuc_bin) = Nki(nuc_bin)+nadd*boxvol*dt
         Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
         Gcf(srtnh4) = Gci(srtnh4)
         
         do k=1,ibins
            if (k .ne. nuc_bin)then
               Nkf(k) = Nki(k)
               do i=1,icomp
                  Mkf(k,i) = Mki(k,i)
               enddo
            else
               do i=1,icomp
                  if (i.ne.srtso4) then
                     Mkf(k,i) = Mki(k,i)
                  endif
               enddo
            endif
         enddo
         
         do k=1,ibins
            if (Nkf(k).lt.1.d0) then
               Nkf(k) = 0.d0
               do j=1,icomp
                  Mkf(k,j) = 0.d0
               enddo
            endif
         enddo
         call mnfix(Nkf,Mkf)
         
C     there is a chance that Gcf will go less than zero because we are artificially growing
C     particles into the first size bin.  don't let it go less than zero.

         
      else if (fn.gt.0.d0.and.d1.gt.Dpk(1)) then
         fn1=fn
         k=1
         do while (d1.ge.Dpk(k+1))
            k=k+1
         enddo

         nuc_bin=k
         mnuc=sqrt(xk(nuc_bin)*xk(nuc_bin+1))

         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+fn1*mnuc*
     &        boxvol*dt/(1.d0+soa_amp)
         Mkf(nuc_bin,srtocil) = Mki(nuc_bin,srtocil)+fn1*mnuc*
     &        boxvol*dt*(1.d0-1.d0/(1.d0+soa_amp))
         Nkf(nuc_bin) = Nki(nuc_bin)+fn1*boxvol*dt
         Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
         Gcf(srtnh4) = Gci(srtnh4)

         if(k.gt.5)then
            print*,'big rnuc',Dpk(k),d1,Dpk(k+1)
         endif

         do k=1,ibins
            if (k .ne. nuc_bin)then
               Nkf(k) = Nki(k)
               do i=1,icomp
                  Mkf(k,i) = Mki(k,i)
               enddo
            else
               do i=1,icomp
                  if (i.ne.srtso4) then
                     Mkf(k,i) = Mki(k,i)
                  endif
               enddo
            endif
         enddo
         
         do k=1,ibins
            if (Nkf(k).lt.1.d0) then
               Nkf(k) = 0.d0
               do j=1,icomp
                  Mkf(k,j) = 0.d0
               enddo
            endif
         enddo
         call mnfix(Nkf,Mkf)
         
      else
         
         do k=1,ibins
            Nkf(k) = Nki(k)
            do i=1,icomp
               Mkf(k,i) = Mki(k,i)
            enddo
         enddo
         
      endif
      
      RETURN
      END


!@sum ezcond  :  takes a given amount of mass and condenses it
!@+   across the bins accordingly.  
!@auth   Jeff Pierce, May 2007

!   Initial values of
!   =================

!@var   Nki(ibins) - number of particles per size bin in grid cell
!@var   Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
!@var   mcond - mass of species to condense [kg/grid cell]
!@var   spec - the number of the species to condense

!-----OUTPUTS-----------------------------------------------------------

!     Nkf, Mkf - same as above, but final values

      SUBROUTINE ezcond(Nki,Mki,mcondi,spec,Nkf,Mkf)

      USE TOMAS_AEROSOL
      USE DOMAIN_DECOMP_ATM, only : am_i_root 
      USE TRACER_COM, only : xk
      IMPLICIT NONE

      real*8 Nki(ibins), Mki(ibins, icomp)
      real*8 Nkf(ibins), Mkf(ibins, icomp)
      real*8 mcondi
      integer spec

      integer i,j,k,c,n,jc,flag,l           ! counters
      real*8 mcond
      real*8 CS                 ! condensation sink [s^-1]
      real*8 sinkfrac(ibins+1)  ! fraction of CS in size bin
      real*8 Nk1(ibins), Mk1(ibins, icomp)
      real*8 Nk2(ibins), Mk2(ibins, icomp)
      real*8 madd               ! mass to add to each bin [kg]
      real*8 maddp(ibins)       ! mass to add per particle [kg]
      real*8 mconds             ! mass to add per step [kg]
      integer nsteps            ! number of condensation steps necessary
      integer floor, ceil       ! the floor and ceiling (temporary)
      real*8 eps                ! small number
      real*8 tdt                !the value 2/3
      real*8 mpo,mpw            !dry and "wet" mass of particle
      real*8 WR                 !wet ratio
      real*8 tau(ibins)         !driving force for condensation
      real*8 totsinkfrac        ! total sink fraction not including nuc bin
      real*8 CSeps              ! lower limit for condensation sink
      real*8 tot_m,tot_s        !total mass, total sulfate mass
      real*8 ratio              ! used in mass correction
      real*8 fracch(ibins,icomp)
      real*8 totch
      real*8 zeros(icomp)
      real*8 tot_i,tot_f,tot_fa ! used for conservation of mass check
      
      parameter(eps=1.d-40)
      parameter(CSeps=1.d-20)


      tdt=2.d0/3.d0

      mcond=mcondi

! initialize variables
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
      enddo




      call mnfix(Nk1,Mk1)


! get the sink fractions
      call getCondSink(Nk1,Mk1,spec,CS,sinkfrac) ! set Nnuc to zero for this calc

! make sure that condensation sink isn't too small
      if (CS.lt.CSeps) then     ! just make particles in first bin
         Mkf(1,spec) = Mk1(1,spec) + mcond
         Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
         do j=1,icomp
            if (icomp.ne.spec) then
               Mkf(1,j) = Mk1(1,j)
            endif
         enddo
         do k=2,ibins
            Nkf(k) = Nk1(k)
            do j=1,icomp
               Mkf(k,j) = Mk1(k,j)
            enddo
         enddo
         return
      endif
      
! determine how much mass to add to each size bin
! also determine how many condensation steps we need
      totsinkfrac = 0.d0
      do k=1,ibins
        totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
      enddo
      nsteps = 1
      do k=1,ibins
         if (sinkfrac(k).lt.1.0D-20)then
            madd = 0.d0
         else
            madd = mcond*sinkfrac(k)/totsinkfrac
         endif
         mpo=0.0
         do j=1,icomp-idiag
            mpo=mpo + Mk1(k,j)
         enddo
         floor = int(madd*0.00001/mpo)
         ceil = floor + 1
         nsteps = max(nsteps,ceil) ! don't let the mass increase by more than 10%
      enddo

! mass to condense each step
      mconds = mcond/nsteps
      
! do steps of condensation
      do i=1,nsteps
         if (i.ne.1) then
            call getCondSink(Nk1,Mk1,spec,
     &        CS,sinkfrac)      ! set Nnuc to zero for this calculation
            totsinkfrac = 0.d0
            do k=1,ibins
              totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
            enddo
         endif      
         
         tot_m=0.d0
         tot_s=0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tot_m = tot_m + Mk1(k,j)
               if (j.eq.spec) then
                  tot_s = tot_s + Mk1(k,j)
               endif
            enddo
         enddo

         if (mcond.gt.tot_m*1.0D-3) then

            do k=1,ibins
               mpo=0.0
               mpw=0.0
!WIN'S CODE MODIFICATION 6/19/06
!THIS MUST CHANGED WITH THE NEW dmdt_int.f
               do j=1,icomp-idiag
                  mpo = mpo+Mk1(k,j) !accumulate dry mass
               enddo
               do j=1,icomp
                  mpw = mpw+Mk1(k,j) ! have wet mass include amso4
               enddo
               WR = mpw/mpo     !WR = wet ratio = total mass/dry mass
               if (Nk1(k) .gt. 0.d0) then
                  maddp(k) = mconds*sinkfrac(k)/totsinkfrac/Nk1(k)
                  mpw=mpw/Nk1(k)
                  tau(k)=1.5d0*((mpw+maddp(k)*WR)**tdt-mpw**tdt) !added WR to moxid term (win, 5/15/06)
c     tau(k)=0.d0
c     maddp(k)=0.d0
               else
                                !nothing in this bin - set tau to zero
                  tau(k)=0.d0
                  maddp(k) = 0.d0
               endif
            enddo

            call mnfix(Nk1,Mk1)
                                ! do condensation

            call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,maddp)
c     call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec)
C     jrp         totch=0.0
C     jrp         do k=1,ibins
C     jrp            do j=1,icomp
C     jrp               fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
C     jrp               totch = totch + (Mk2(k,j)-Mk1(k,j))
C     jrp            enddo
C     jrp         enddo

         elseif (mcond.gt.tot_s*1.0D-12) then
            do k=1,ibins
               if (Nk1(k) .gt. 0.d0) then
                  maddp(k) = mconds*sinkfrac(k)/totsinkfrac
               else
                  maddp(k) = 0.d0
               endif
               Mk2(k,spec)=Mk1(k,spec)+maddp(k)
               do j=1,icomp
                  if (j.ne.spec) then
                     Mk2(k,j)=Mk1(k,j)
                  endif
               enddo
               Nk2(k)=Nk1(k)
            enddo
            call mnfix(Nk2,Mk2)
         else ! do nothing
            mcond = 0.d0
            do k=1,ibins
               Nk2(k)=Nk1(k)
               do j=1,icomp
                  Mk2(k,j)=Mk1(k,j)
               enddo
            enddo
         endif
         if (i.ne.nsteps)then
            do k=1,ibins
               Nk1(k)=Nk2(k)
               do j=1,icomp
                  Mk1(k,j)=Mk2(k,j)
               enddo
            enddo            
         endif

      enddo

      do k=1,ibins
         Nkf(k)=Nk2(k)
         do j=1,icomp
            Mkf(k,j)=Mk2(k,j)
         enddo
      enddo

! check for conservation of mass
      tot_i = 0.d0
      tot_fa = mcond
      tot_f = 0.d0
      do k=1,ibins
         tot_i=tot_i+Mki(k,spec)
         tot_f=tot_f+Mkf(k,spec)
         tot_fa=tot_fa+Mki(k,spec)
      enddo

      if (mcond.gt.0.d0.and.
     &    abs((mcond-(tot_f-tot_i))/mcond).gt.0.d0) then
         if (abs((mcond-(tot_f-tot_i))/mcond).lt.1.d0) then
            ! do correction of mass
            ratio = (tot_f-tot_i)/mcond
            do k=1,ibins
               Mkf(k,spec)=Mki(k,spec)+
     &              (Mkf(k,spec)-Mki(k,spec))/ratio
            enddo
            call mnfix(Nkf,Mkf)
         else
            if(am_i_root())then
            print*,'ERROR in ezcond',spec
            print*,'Condensation error',(mcond-(tot_f-tot_i))/mcond
            print*,'mcond',mcond,'change',tot_f-tot_i
            print*,'tot_i',tot_i,'tot_fa',tot_fa,'tot_f',tot_f
            print*,'Nki',Nki
            print*,'Nkf',Nkf
            print*,'Mki',Mki
            print*,'Mkf',Mkf
            call STOP_model('ERROR in ezcond',255)
            endif
         endif
      endif

! check for conservation of mass
      tot_i = 0.d0
      tot_f = 0.d0
      do k=1,ibins
         tot_i=tot_i+Mki(k,srtnh4)
         tot_f=tot_f+Mkf(k,srtnh4)
      enddo
      if(tot_i.gt.0)then !YUNHA LEE  - I added this to avoid floating invalid in MODELE-TOMAS.
      if (abs(tot_f-tot_i)/tot_i.gt.1.0D-8)then
         print*,'No N conservation in ezcond.f'
         print*,'initial',tot_i
         print*,'final',tot_f
      endif
      endif

      return
      end


!@auth  Peter Adams/Modified by Yunha Lee Mar 2008 

!@sum TMCOND :  CONDENSATION Based on Tzivion, Feingold, Levin, JAS 1989 and 
!@+             Stevens, Feingold, Cotton, JAS 1996

!@+  The supersaturation is calculated outside of the routine and assumed
!@+  to be constant at its average value over the timestep.
!@+  
!@+  The method has three basic components:
!@+  (1) first a top hat representation of the distribution is construced
!@+      in each bin and these are translated according to the analytic
!@+      solutions
!@+  (2) The translated tophats are then remapped to bins.  Here if a 
!@+      top hat entirely or in part lies below the lowest bin it is 
!@+      not counted.
!@+  

!@+   Additional notes (Peter Adams)

!@+   I have changed the routine to handle multicomponent aerosols.  The
!@+   arrays of mass moments are now two dimensional (size and species).
!@+   Only a single component (CSPECIES) is allowed to condense during
!@+   a given call to this routine.  Multicomponent condensation/evaporation
!@+   is accomplished via multiple calls.  Variables YLC and YUC are
!@+   similar to YL and YU except that they refer to the mass of the 
!@+   condensing species, rather than total aerosol mass.

!@+   I have removed ventilation variables (VSW/VNTF) from the subroutine
!@+   call.  They still exist internally within this subroutine, but
!@+   are initialized such that they do nothing.

!@+   I have created a new variable, AMKDRY, which is the total mass in
!@+   a size bin (sum of all chemical components excluding water).  I
!@+   have also created WR, which is the ratio of total wet mass to 
!@+   total dry mass in a size bin.

!@+   AMKC(k,j) is the total amount of mass after condensation of species
!@+   j in particles that BEGAN in bin k.  It is used as a diagnostic
!@+   for tracking down numerical errors.

!@+   End of my additional notes

!@var  TAU(k) ......... Forcing for diffusion = (2/3)*CPT*ETA_BAR*DELTA_T
!@var  X(K) ........ Array of bin limits in mass space
!@var  AMKD(K,J) ... Input array of mass moments
!@var  ANKD(K) ..... Input array of number moments
!@var  AMK(K,J) .... Output array of mass moments
!@var  ANK(K) ...... Output array of number moments
!@var  CSPECIES .... Index of chemical species that is condensing

      SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,CSPECIES,moxd)

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk
      USE DOMAIN_DECOMP_ATM, only : am_i_root

      IMPLICIT NONE

      INTEGER L,I,J,K,IMN,CSPECIES,kk
      real*8 DN,DM,DYI,TAU(ibins),XL,XU,YL,YLC,YU,YUC
      real*8 TEPS,NEPS,EX2,ZERO
      real*8 XI,XX,XP,YM,WTH,W1,W2,WW,AVG
      real*8 VSW,VNTF(ibins)
      real*8 TAU_L, maxtau
      real*8 X(ibins+1),AMKD(ibins,icomp),ANKD(ibins)
      real*8 AMK(ibins,icomp),ANK(ibins)
      real*8 AMKDRY(ibins), AMKWET(ibins), WR(ibins)
      real*8 DMDT_INT
      real*8 AMKD_tot
      PARAMETER (TEPS=1.0d-40,NEPs=1.0d-20)
      PARAMETER (EX2=2.d0/3.d0,ZERO=0.0d0)
      real*8 moxd(ibins)        !moxid/Nact (win, 5/25/06)
      real*8 c1, c2             !correction factor (win, 5/25/06)
      real*8 xk_hi,tmpvar,xk_lo,frac_lo_n,frac_lo_m
!      external DMDT_INT

 3    format(I4,200E20.11)



C If any ANKD are zero, set them to a small value to avoid division by zero
      do k=1,ibins
         if (ANKD(k) .lt. NEPS) then
            ANKD(k)=NEPS
            AMKD(k,srtso4)=NEPS*sqrt(X(k)*X(k+1)) !make the added particles SO4
            do j=1,icomp
               if (j.ne.srtso4)then
                  AMKD(k,j)=0.d0
               endif
            enddo
         endif
      enddo

Cpja Sometimes, after repeated condensation calls, the average bin mass
Cpja can be just above the bin boundary - in that case, transfer a some
Cpja to the next highest bin

!TOMAS - If AMKD_tot/ANKD(k) >X (k+1), it will still great after this part. 
! mp is same by moving 10% mass and number.  That's the point here??
      do k=1,ibins-1
        AMKD_tot=0.0d0
        do kk=1,icomp-idiag
        AMKD_tot=AMKD_tot+AMKD(k,kk)
        enddo
         if (AMKD_tot/ANKD(k).gt.X(k+1)) then
            do j=1,icomp
               AMKD(k+1,j)=AMKD(k+1,j)+0.1d0*AMKD(k,j)
               AMKD(k,j)=AMKD(k,j)*0.9d0
            enddo
            ANKD(k+1)=ANKD(k+1)+0.1d0*ANKD(k)
            ANKD(k)=ANKD(k)*0.9d0
         endif
      enddo

Cpja Initialize ventilation variables so they don't do anything
      VSW=0.0d0
      DO L=1,ibins
         VNTF(L)=0.0d0
      ENDDO

Cpja Initialize AMKDRY and WR
      DO L=1,ibins
         AMKDRY(L)=0.0d0
         AMKWET(L)=0.0d0
         DO J=1,icomp-idiag
            AMKDRY(L)=AMKDRY(L)+AMKD(L,J)
         ENDDO
         DO J=1,icomp
            AMKWET(L)=AMKWET(L)+AMKD(L,J)
         ENDDO
         WR(L)=AMKWET(L)/AMKDRY(L)
      ENDDO

Cpja Initialize X() array of particle masses based on xk()
      DO L=1,ibins
         X(L)=xk(L)
      ENDDO

c
c Only solve when significant forcing is available
c
      maxtau=0.0d0
      do l=1,ibins
         maxtau=max(maxtau,abs(TAU(l)))
      enddo
      IF(ABS(maxtau).LT.TEPS)THEN
         DO L=1,ibins
            DO J=1,icomp
               AMK(L,J)=AMKD(L,J)
            ENDDO
            ANK(L)=ANKD(L)
         ENDDO
      ELSE
         DO L=1,ibins
            DO J=1,icomp
               AMK(L,J)=0.d0
            ENDDO
            ANK(L)=0.d0
         ENDDO
         WW=0.5d0
c        IF(TAU.LT.0.)WW=.5d0
c
c identify tophats and do lagrangian growth
c
         DO L=1,ibins
            IF(ANKD(L).EQ.0.)GOTO 200

            !if tau is zero, leave everything in same bin
            IF (TAU(L) .EQ. 0.) THEN
               ANK(L)=ANK(L)+ANKD(L)
               DO J=1,icomp
                  AMK(L,J)=AMK(L,J)+AMKD(L,J)
               ENDDO
            ENDIF
            IF (TAU(L) .EQ. 0.) GOTO 200

Cpja Limiting AVG, the average particle size to lie within the size
Cpja bounds causes particles to grow or shrink arbitrarily and is
Cpja wreacking havoc with choosing condensational timesteps and
Cpja conserving mass.  I have turned them off.
c            AVG=MAX(X(L),MIN(X(L+1),AMKDRY(L)/(NEPS+ANKD(L))))
            AVG=AMKDRY(L)/ANKD(L)
            XX=X(L)/AVG
cyhl Eq below needs to be changed for 15 size bins. 
cyhl double check this equation!!!
            if(l.lt.ibins-1)then
               XI=.5d0 + XX*(2.5d0 - 2.0d0*XX)

cyhl this is the old equation for 30 bins: XI=.5d0 + XX*(1.5d0 - XX)
            if (XI .LT. 1.d0) then
               !W1 will have sqrt of negative number
               write(*,*)'ERROR: tmcond - XI<1 for bin: ',L,XK(L)
               write(*,*)'lower limit is',X(L),TAU(L)
               write(*,*)'AVG is ',AVG
               write(*,*)'Nk is ', ANKD(L)
               write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
               write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
               call stop_model('ERROR in tmcond',255)
            endif
cyhl            W1 =SQRT(12.d0*(XI-1.d0))*AVG
cyhl            W2 =MIN(X(L+1)-AVG,AVG-X(L))
cyhl            WTH=W1*WW+W2*(1.d0-WW)
            W1 =SQRT(12.d0*(XI-1.d0))*AVG/4.0d0 ! cyhl 4.0=xk(k+1)/xk(k)
            W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0d0
            
            else                   ! 32 = xk(k+1)/xk(k)
            XI=.5d0 + XX*(16.5d0 - 16.0d0*XX)
            if (XI .LT. 1.d0) then
               !W1 will have sqrt of negative number
               write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
               write(*,*)'lower limit is',X(L)
               write(*,*)'AVG is ',AVG
               write(*,*)'Nk is ', ANKD(L)
               write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
               write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
               call stop_model('ERROR in tmcond',255)
            endif
               W1 =SQRT(12.d0*(XI-1.d0))*AVG/32.0d0 ! cyhl 32.0=xk(k+1)/xk(k)
               W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0d0
            endif

            WTH=W1*WW+W2*(1.d0-WW)  ! average of w1 and w2
            IF(WTH.GT.1.) then
               write(*,*)'WTH>1 in cond.f, bin #',L,W1,W2
               call stop_model('ERROR in tmcond',255)
            ENDIF
            XU=AVG+WTH*.5d0
            XL=AVG-WTH*.5d0
c Ventilation added bin-by-bin
            TAU_L=TAU(l)*MAX(1.d0,VNTF(L)*VSW)
            IF(TAU_L/TAU(l).GT. 6.) THEN
               PRINT *,'TAU..>6.',TAU(l),TAU_L,VSW,L
            ENDIF
            IF(TAU_L.GT.TAU(l)) THEN 
               PRINT *,'TAU...',TAU(l),TAU_L,VSW,L
            ENDIF
! prior to 5/25/06 (win)
!            YU=DMDT_INT(XU,TAU_L,WR(L))
!            YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
!            IF (YU .GT. X(ibins+1) ) THEN
!               YUC=YUC*X(ibins+1)/YU
!               YU=X(ibins+1)
!            ENDIF
!            YL=DMDT_INT(XL,TAU_L,WR(L)) 
!            YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
!add new correction factor to YU and YL (win, 5/25/06)
            YU=DMDT_INT(XU,TAU_L,WR(L))
            YL=DMDT_INT(XL,TAU_L,WR(L)) 
            if(moxd(L).eq.0d0) then
               c1=1.d0          !for so4cond call, without correction factor.
            else
               c1 = moxd(L)*2.d0/(YU+YL-XU-XL)
            endif
            c2 = c1 - (c1-1.d0)*(XU+XL)/(YU+YL)
            YU = YU*c2
            YL = YL*c2
!end part for fudging to get higher AVG 

            YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
            IF (YU .GT. X(ibins+1) ) THEN
               YUC=YUC*X(ibins+1)/YU
               YU=X(ibins+1)
            ENDIF
            YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
            DYI=1.d0/(YU-YL)

c            print*,'XL',XL,'YL',YL,'XU',XU,'YU',YU
c
c deal with portion of distribution that lies below lowest gridpoint
c
            IF(YL.LT.X(1))THEN
cpja Instead of the following, I will just add all new condensed
cpja mass to the same size bin
c               if ((YL/XL-1.d0) .LT. 1.d-3) then
c                  !insignificant growth - leave alone
c                  ANK(L)=ANK(L)+ANKD(L)
c                  DO J=1,icomp-1
c                     AMK(L,J)=AMK(L,J)+AMKD(L,J)
c                  ENDDO
c                  GOTO 200
c               else
c                  !subtract out lower portion
c                  write(*,*)'ERROR in cond - low portion subtracted'
c                  write(*,*) 'Nk,Mk: ',ANKD(L),AMKD(L,1),AMKD(L,2)
c                  write(*,*) 'TAU: ', TAU_L
c                  write(*,*) 'XL, YL, YLC: ',XL,YL,YLC
c                  write(*,*) 'XU, YU, YUC: ',XU,YU,YUC
c                  ANKD(L)=ANKD(L)*MAX(ZERO,(YU-X(1)))*DYI
c                  YL=X(1)
c                  YLC=X(1)*AMKD(1,CSPECIES)/AMKDRY(1)
c                  DYI=1.d0/(YU-YL)
c               endif
               ANK(L)=ANK(L)+ANKD(L)
               do j=1,icomp
                  if (J.EQ.CSPECIES) then
                     AMK(L,J)=AMK(L,J)+(YUC+YLC)*.5d0*ANKD(L)
                  else
                     AMK(L,J)=AMK(L,J)+AMKD(L,J)
                  endif
               enddo
               GOTO 200
            ENDIF
            IF(YU.LT.X(1))GOTO 200
c
c Begin remapping (start search at present location if condensation)
c
            IMN=1
            IF(TAU(l).GT.0.)IMN=L
            DO I=IMN,ibins
               IF(YL.LT.X(I+1))THEN
                  IF(YU.LE.X(I+1))THEN
                     DN=ANKD(L)
                     do j=1,icomp
                        DM=AMKD(L,J)
                        IF (J.EQ.CSPECIES) THEN
                           AMK(I,J)=(YUC+YLC)*.5d0*DN+AMK(I,J)
                        ELSE
                           AMK(I,J)=AMK(I,J)+DM
                        ENDIF
                     enddo
                     ANK(I)=ANK(I)+DN
                  ELSE
                     DN=ANKD(L)*(X(I+1)-YL)*DYI
                     do j=1,icomp
                        DM=AMKD(L,J)*(X(I+1)-YL)*DYI
                        IF (J.EQ.CSPECIES) THEN
                           XP=DMDT_INT(X(I+1),-1.0d0*TAU_L,WR(L))
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(I+1)-XP
                           AMK(I,J)=DN*(YM+YLC)*0.5d0+AMK(I,J)
                        ELSE
                           AMK(I,J)=AMK(I,J)+DM
                        ENDIF
                     enddo
                     ANK(I)=ANK(I)+DN
                     DO K=I+1,ibins
                        IF(YU.LE.X(K+1))GOTO 100
                        DN=ANKD(L)*(X(K+1)-X(K))*DYI
                        do j=1,icomp
                           DM=AMKD(L,J)*(X(K+1)-X(K))*DYI
                           IF (J.EQ.CSPECIES) THEN
                              XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))
                              YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP
                              AMK(K,J)=DN*1.5d0*YM+AMK(K,J)
                           ELSE
                              AMK(K,J)=AMK(K,J)+DM
                           ENDIF
                        enddo
                        ANK(K)=ANK(K)+DN
                     ENDDO
                     print*,'Trying to put stuff in bin ibins+1'
                     call stop_model('ERROR in tmcond,',255)
 100                 CONTINUE
                     DN=ANKD(L)*(YU-X(K))*DYI
                     do j=1,icomp
                        DM=AMKD(L,J)*(YU-X(K))*DYI
                        IF (J.EQ.CSPECIES) THEN
                           XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP
                           AMK(K,J)=DN*(YUC+YM)*0.5d0+AMK(K,J)
                        ELSE
                           AMK(K,J)=AMK(K,J)+DM
                        ENDIF
                     enddo
                     ANK(K)=ANK(K)+DN
                  ENDIF  !YU.LE.X(I+1)
                  GOTO200
               ENDIF   !YL.LT.X(I+1)
            ENDDO !I loop
 200        CONTINUE
         ENDDO    !L loop
      ENDIF

      RETURN
      END


c----------------------------------------------------------------------
!@sum Function DMDT_INT:  Here we apply the analytic solution to the
!@+  droplet growth equation in mass space for a given scale length which
!@+  mimics the inclusion of gas kinetic effects
!@+ Reference: Stevens et al. 1996, Elements of the Microphysical Structure
!@+           of Numerically Simulated Nonprecipitating Stratocumulus,
!@+           J. Atmos. Sci., 53(7),980-1006. 

!@+  This calculates a solution for m(t+dt) using eqn.(A3) from the reference

!@+  Comments by Peter Adams :
!@+  I have changed the length scale.  Non-continuum effects are
!@+  assumed to be taken into account in choice of tau (in so4cond
!@+  subroutine).

!@+  I have also added another argument to the function call, WR.  This
!@+  is the ratio of wet mass to dry mass of the particle.  I use this
!@+  information to calculate the amount of growth of the wet particle,
!@+  but then return the resulting dry mass.  This is the appropriate
!@+  way to implement the condensation algorithm in a moving sectional
!@+  framework.

!@+  End of comments

!@var  M0 ......... initial mass
!@var  L0 ......... length scale
!@var  Tau ........ forcing from vapor field

      real*8 FUNCTION DMDT_INT(M0,TAU,WR)
 
      IMPLICIT NONE
      real*8 M0,TAU,X,L0,C,ZERO,WR,MH2O
      PARAMETER (C=2.d0/3.d0,L0=0.0d0,ZERO=0.0d0)
 
      MH2O=(WR-1.d0)*M0
      X=((M0+MH2O)**C+L0)
      X=MAX(ZERO,SQRT(MAX(ZERO,C*TAU+X))-L0)
!win,5/14/06      DMDT_INT=X*X*X-MH2O
      DMDT_INT = X*X*X/WR  !<step5.2> change calculation to keep WR 
                           !constant after condensation/evap (win, 5/14/06)

Cpja Perform some numerical checks on dmdt_int
      if ((tau .gt. 0.0) .and. (dmdt_int .lt. m0)) dmdt_int=m0
      if ((tau .lt. 0.0) .and. (dmdt_int .gt. m0)) dmdt_int=m0

      RETURN
      END


!@sum cf_nucl  :  calculates the barrierless nucleation rate and radius of the 
!@+   critical nucleation cluster using the parameterization of...
!@+   Reference : Clement and Ford (1999) Atmos. Environ. 33:489-499
!@auth   Jeff Pierce, April 2007

      SUBROUTINE cf_nucl(temp,rh,cna,nh3ppt,fn)

      IMPLICIT NONE

      real*8 cna      ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 nh3ppt   ! mixing ratio of ammonia in ppt
      real*8 fn                   ! nucleation rate [cm-3 s-1]
      real*8 rnuc                 ! critical cluster radius [nm]

      real*8 temp                 ! temperature of air [K]
      real*8 rh                   ! relative humidity of air as a fraction
      real*8 alpha1


      if (nh3ppt .lt. 0.1) then
         alpha1=4.276e-10*sqrt(temp/293.15) ! For sulfuric acid
      else
         alpha1=3.684e-10*sqrt(temp/293.15) ! For ammonium sulfate
      endif
      fn = alpha1*cna**2*3600.
c sensitivity       fn = 1.e-3 * fn ! 10^-3 tuner
      if (fn.gt.1.0e9) fn=1.0e9 ! For numerical conversion

 10   return
      end




!@sum vehk_nucl    :  calculates the binary nucleation rate and radius of the 
!@+   critical nucleation cluster using the parameterization of 

!@+   Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck, 
!@+   M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric 
!@+   Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions." 
!@+   Journal of Geophysical Research-Atmospheres 107, no. D22 (2002).

!@auth   Jeff Pierce, April 2007

      SUBROUTINE vehk_nucl(temp,rh,cnai,fn,rnuc)

      IMPLICIT NONE

      real*8 cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 fn                   ! nucleation rate [cm-3 s-1]
      real*8 rnuc                 ! critical cluster radius [nm]

      real*8 fb0(10),fb1(10),fb2(10),fb3(10),fb4(10),fb(10)
      real*8 gb0(10),gb1(10),gb2(10),gb3(10),gb4(10),gb(10) ! set parameters
      real*8 temp                 ! temperature of air [K]
      real*8 rh                   ! relative humidity of air as a fraction
      real*8 cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 xstar                ! mole fraction sulfuric acid in cluster
      real*8 ntot                 ! total number of molecules in cluster
      integer i                 ! counter

c     Nucleation Rate Coefficients
c
      data fb0 /0.14309, 0.117489, -0.215554, -3.58856, 1.14598,
     $          2.15855, 1.6241, 9.71682, -1.05611, -0.148712        /
      data fb1 /2.21956, 0.462532, -0.0810269, 0.049508, -0.600796,
     $       0.0808121, -0.0160106, -0.115048, 0.00903378, 0.00283508/
      data fb2 /-0.0273911, -0.0118059, 0.00143581, -0.00021382, 
     $       0.00864245, -0.000407382, 0.0000377124, 0.000157098,
     $       -0.0000198417, -9.24619d-6 /
      data fb3 /0.0000722811, 0.0000404196, -4.7758d-6, 3.10801d-7,
     $       -0.0000228947, -4.01957d-7, 3.21794d-8, 4.00914d-7,
     $       2.46048d-8, 5.00427d-9 /
      data fb4 /5.91822, 15.7963, -2.91297, -0.0293333, -8.44985,
     $       0.721326, -0.0113255, 0.71186, -0.0579087, -0.0127081  / 

c     Coefficients of total number of molecules in cluster     
      data gb0 /-0.00295413, -0.00205064, 0.00322308, 0.0474323,
     $         -0.0125211, -0.038546, -0.0183749, -0.0619974,
     $         0.0121827, 0.000320184 /
      data gb1 /-0.0976834, -0.00758504, 0.000852637, -0.000625104,
     $         0.00580655, -0.000672316, 0.000172072, 0.000906958,
     $         -0.00010665, -0.0000174762 /      
      data gb2 /0.00102485, 0.000192654, -0.0000154757, 2.65066d-6,
     $         -0.000101674, 2.60288d-6, -3.71766d-7, -9.11728d-7,
     $         2.5346d-7, 6.06504d-8 /
      data gb3 /-2.18646d-6, -6.7043d-7, 5.66661d-8, -3.67471d-9,
     $         2.88195d-7, 1.19416d-8, -5.14875d-10, -5.36796d-9,
     $         -3.63519d-10, -1.42177d-11 /
      data gb4 /-0.101717, -0.255774, 0.0338444, -0.000267251,
     $         0.0942243, -0.00851515, 0.00026866, -0.00774234,
     $         0.000610065, 0.000135751 /


      cna=cnai/5.

c     Respect the limits of the parameterization
      if (cna .lt. 1.d4) then ! limit sulf acid conc
         fn = 0.
         rnuc = 1.
c         print*,'cna < 1D4', cna
         goto 10
      endif
      if (cna .gt. 1.0d11) cna=1.0e11 ! limit sulfuric acid conc  
      if (temp .lt. 230.15) temp=230.15 ! limit temp
      if (temp .gt. 305.15) temp=305.15 ! limit temp
      if (rh .lt. 1d-4) rh=1d-4 ! limit rh
      if (rh .gt. 1.) rh=1. ! limit rh
c
c     Mole fraction of sulfuric acid
      xstar=0.740997-0.00266379*temp-0.00349998*log(cna)
     &   +0.0000504022*temp*log(cna)+0.00201048*log(rh)
     &   -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2.
     &   -0.0000179059*temp*(log(rh))**2.
     &   +0.000184403*(log(rh))**3.
     &   -1.50345d-6*temp*(log(rh))**3.
c 
c     Nucleation rate coefficients 
      do i=1, 10
         fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2.
     &        +fb3(i)*temp**3.+fb4(i)/xstar
      enddo
c
c     Nucleation rate (1/cm3-s)
      fn = exp(fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2.
     &    +fb(4)*(log(rh))**3.+fb(5)*log(cna)
     &    +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna)
     &    +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2.
     &    +fb(10)*(log(cna))**3.)

c
c   Cap at 10^6 particles/s, limit for parameterization
      if (fn.gt.1.0d6) then
         fn=1.0d6
      endif
c
c     Coefficients of total number of molecules in cluster 
      do i=1, 10
         gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2.
     &        +gb3(i)*temp**3.+gb4(i)/xstar
      enddo
c     Total number of molecules in cluster
      ntot=exp(gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2.
     &    +gb(4)*(log(rh))**3.+gb(5)*log(cna)
     &    +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna)
     &    +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2.
     &    +gb(10)*(log(cna))**3.)

c     cluster radius
      rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]

 10   return
      end




!@sum  bl_nucl :  calculates a simple binary nucleation rate of 1 nm
!@+   particles.
!@+       j_1nm = A * [H2SO4]
!@auth   Jeff Pierce, April 2007

      SUBROUTINE bl_nucl(cnai,fn,rnuc)

      IMPLICIT NONE


      real*8,intent(in) :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]

      real*8 fn                 ! nucleation rate [cm-3 s-1]
      real*8 rnuc               ! critical cluster radius [nm]

      real*8 cna                ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 A                  ! prefactor... empirical
      parameter(A=2.0D-6)


      cna=cnai

      fn=A*cna
      rnuc=0.5d0                ! particle diameter of 1 nm

      return
      end SUBROUTINE bl_nucl
            

!@sum  napa_nucl :  calculates the ternary nucleation rate and radius of the 
!@+   critical nucleation cluster using the parameterization of 
!@+     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of 
!@+    Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical 
!@+     Research-Atmospheres 107, no. D19 (2002).

!@auth   Jeff Pierce, April 2007

      SUBROUTINE napa_nucl(temp,rh,cnai,nh3ppti,fn,rnuc)

      IMPLICIT NONE

      real*8 cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 nh3ppti              ! concentration of gas phase ammonia

      real*8 fn                   ! nucleation rate [cm-3 s-1]
      real*8 rnuc                 ! critical cluster radius [nm]

      real*8 aa0(20),a1(20),a2(20),a3(20),fa(20) ! set parameters
      real*8 fnl                  ! natural log of nucleation rate
      real*8 temp                 ! temperature of air [K]
      real*8 rh                   ! relative humidity of air as a fraction
      real*8 cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8 nh3ppt               ! concentration of gas phase ammonia
      integer i                 ! counter

      data aa0 /-0.355297, 3.13735, 19.0359, 1.07605, 6.0916,
     $         0.31176, -0.0200738, 0.165536,
     $         6.52645, 3.68024, -0.066514, 0.65874,
     $         0.0599321, -0.732731, 0.728429, 41.3016,
     $         -0.160336, 8.57868, 0.0530167, -2.32736        /

      data a1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378,
     $         1.64009, -0.752115, 3.26623, -0.258002, -0.204098,
     $         -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736,
     $         -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/
     
      data a2 /0.34536, 0.00561204, 0.000479808, -0.00796052,
     $         0.00939836, -0.00343852, 0.00525813, -0.0489703,
     $         0.00143456, 0.00106259, 0.0122938, -0.00165718,
     $         -0.0362432, 0.000147186, -0.027422, 0.000904383,
     $         -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/
     
      data a3 /-0.000824007, -9.74576d-06, -4.14699d-07, 7.61229d-06,
     $         -1.74927d-05, -1.09753d-05, -8.98038d-06, 0.000146967,
     $         -2.02036d-06, -1.2656d-06, 6.18554d-05, 3.41744d-06,
     $         4.93337d-05, -2.37711d-07, 4.93478d-05, -5.73788d-07, 
     $         8.39522d-08, -6.48365d-07, -2.93564d-05, 8.0459d-08   /


      cna=cnai
      nh3ppt=nh3ppti

c     Napari's parameterization is only valid within limited area
      if ((cna .lt. 1.d4).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
         fn = 0.
         rnuc = 1
         goto 10
      endif  
      if (cna .gt. 1.0d9) cna=1.0d9 ! limit sulfuric acid conc
      if (nh3ppt .gt. 100.) nh3ppt=100. ! limit temp  
      if (temp .lt. 240.) temp=240. ! limit temp
      if (temp .gt. 300.) temp=300. ! limit temp
      if (rh .lt. 0.05) rh=0.05 ! limit rh 
      if (rh .gt. 0.95) rh=0.95 ! limit rh

      do i=1,20
         fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
      enddo

      fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2.
     &  +fa(4)*log(nh3ppt)+fa(5)*(log(nh3ppt))**2.+fa(6)*rh
     &  +fa(7)*log(rh)+fa(8)*log(nh3ppt)/log(cna)+fa(9)*log(nh3ppt)
     &  *log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna)
     &  +fa(12)*rh
     &  *log(nh3ppt)+fa(13)*log(rh)/log(cna)+fa(14)*log(rh)
     &  *log(nh3ppt)+fa(15)*(log(nh3ppt))**2./log(cna)+fa(16)*log(cna)
     &  *(log(nh3ppt))**2.+fa(17)*(log(cna))**2.*log(nh3ppt)
     &  +fa(18)*rh
     &  *(log(nh3ppt))**2.+fa(19)*rh*log(nh3ppt)/log(cna)+fa(20)
     &  *(log(cna))**2.*(log(nh3ppt))**2.
c
c
      fn=exp(fnl)
c   Cap at 10^6 particles/cm3-s, limit for parameterization
      if (fn.gt.1.0d6) then
        fn=1.0d6
        fnl=log(fn)
      endif

      rnuc=0.141027-0.00122625*fnl-7.82211d-6*fnl**2.
     &     -0.00156727*temp-0.00003076*temp*fnl
     &     +0.0000108375*temp**2.

 10   return
      end



!@sum gasdiff: This function returns the diffusion constant of a species in
!@+   air (m2/s).  It uses the method of Fuller, Schettler, and
!@+   Giddings as described in Perry's Handbook for Chemical
!@+   Engineers.
!@auth   Peter Adams, May 2000

      real FUNCTION gasdiff(temp,pres,mw,Sv)

      IMPLICIT NONE

      real*8 temp, pres  !temperature (K) and pressure (Pa) of air
      real mw          !molecular weight (g/mol) of diffusing species
      real Sv          !sum of atomic diffusion volumes of diffusing species

      real mwair, Svair   !same as above, but for air
      real mwf, Svf

      parameter(mwair=28.9, Svair=20.1)

      mwf=sqrt((mw+mwair)/(mw*mwair))
      Svf=(Sv**(1./3.)+Svair**(1./3.))**2.
      gasdiff=1.0e-7*temp**1.75*mwf/pres*1.0e5/Svf

      RETURN
      END

!@sum mnfix  :  examines the mass and number distributions and
!@+   determines if any bins have an average mass outside their normal
!@+   range.  I have seen this happen because the GCM advection seems
!@+   to treat the mass and number tracers inconsistently.  If any bins
!@+   are out of range, I shift some mass and number to a new bin in
!@+   a way that conserves both.

!@auth   Peter Adams, September 2000/Jeff Pierce, August 2007

!@var   Nkx and Mkx are the number and mass distributions


      SUBROUTINE mnfix(Nkx,Mkx)

      USE TOMAS_AEROSOL
      USE TRACER_COM, only : xk

      IMPLICIT NONE

      real*8 Nkx(ibins), Mkx(ibins,icomp)

      integer k,j,L,jj          !counters
      real*8 tot_mass           ! total dry mass
      real*8 xk_hi, xk_lo       ! geometric mean mass of bins that mass is moving to
      real*8 xk_hi1, xk_hi2     ! used as tests to see if mass is way to low or high
      real*8 Neps,Meps
      real*8 tmpvar             ! temporary variable
      real*8 frac_lo_n, frac_lo_m ! fraction of the number and mass going to the lower bin
      real*8 totmass
      parameter(Neps=1.d-20,Meps=1.d0-50)


      ! first check for negative tracers
      do k=1,ibins
         if (Nkx(k).lt.Neps)then
            if (Nkx(k).gt.-1.0d10)then
               Nkx(k)=Neps
               do j=1,icomp
                  if (j.eq.srtso4)then
                     Mkx(k,j)=Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j)=0.d0
                  endif
               enddo
            else
               print*,'Negative tracer in mnfix'
               print*,'Nk',Nkx
               print*,'Mkx',Mkx
               call stop_model('neg tracer in mnfix - bin loop',255)
            endif
         endif
         totmass=0.d0
         do j=1,icomp-idiag
            totmass=totmass+Mkx(k,j)
         enddo
         if (totmass.lt.Meps) then
            Nkx(k)=Neps
            do jj=1,icomp
               if (jj.eq.srtso4)then
                  Mkx(k,jj)=Neps*sqrt(xk(k+1)*xk(k))
               else
                  Mkx(k,jj)=0.d0
               endif
            enddo
         endif            
         do j=1,icomp-idiag
            if (Mkx(k,j).lt.0.d0)then
               if (Mkx(k,j).gt.-25.0d0)then  !dmw trying to relax this # a bit
                  Nkx(k)=Neps
                  do jj=1,icomp
                     if (jj.eq.srtso4)then
                        Mkx(k,jj)=Neps*sqrt(xk(k+1)*xk(k))
                     else
                        Mkx(k,jj)=0.d0
                     endif
                  enddo
               else
                  print*,'Negative tracer in mnfix'
                  print*,'Nk',Nkx
                  print*,'Mkx',Mkx
                  call stop_model('neg tracer in mnfix - comp loop',255)
               endif
            endif
         enddo
      enddo

      do k=1,ibins
         tot_mass=0.0d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
         if (tot_mass/Nkx(k).gt.xk(k+1).or.
     &        tot_mass/Nkx(k).lt.xk(k)) then
c            print*,'out of bounts in mnfix, fixing'
c            print*,k,'AVG',tot_mass/Nkx(k),'lo',xk(k),'hi',xk(k+1)
            ! figure out which bins to redistribute to
            xk_hi1 = sqrt(xk(2)*xk(1))
            xk_hi2 = sqrt(xk(ibins+1)*xk(ibins))
            if (xk_hi1.gt.tot_mass/Nkx(k)) then
               !mass per particle very low
               !conserve mass at expense of number
               tmpvar = Nkx(k)
               Nkx(k)=Neps
               Nkx(1)=Nkx(1)+tot_mass/sqrt(xk(2)*xk(1))
               do j=1,icomp
                  tmpvar=Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j) = Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j)=0.d0
                  endif
                  Mkx(1,j)=Mkx(1,j)+tmpvar
               enddo
            elseif (xk_hi2.lt.tot_mass/Nkx(k)) then
               !mass per particle very high
               !conserve mass at expernse of number
               Nkx(k)=Neps
              Nkx(ibins)=Nkx(ibins)+tot_mass/sqrt(xk(ibins+1)*xk(ibins))
               do j=1,icomp
                  tmpvar=Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j) = Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j)=0.d0
                  endif
                  Mkx(ibins,j)=Mkx(ibins,j)+tmpvar
               enddo               
            else ! mass of particle is somewhere within the bins
               L = 2
               xk_hi = sqrt(xk(L+1)*xk(L))
               do while (xk_hi .lt. tot_mass/Nkx(k))
                  L=L+1
                  xk_hi = sqrt(xk(L+1)*xk(L))
               enddo
               xk_lo = sqrt(xk(L)*xk(L-1))
                                ! figure out how much of the number to distribute to the lower bin
               frac_lo_n = (tot_mass - Nkx(k)*xk_hi)/
     &              (Nkx(k)*(xk_lo-xk_hi))
               frac_lo_m = frac_lo_n*Nkx(k)*xk_lo/tot_mass

               tmpvar = Nkx(k)
               Nkx(k) = Neps
               Nkx(L-1) = Nkx(L-1) + frac_lo_n*tmpvar
               Nkx(L) = Nkx(L) + (1-frac_lo_n)*tmpvar
               do j=1,icomp
                  tmpvar = Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j) = Neps*sqrt(xk(k+1)*xk(k))
                  else
                     Mkx(k,j) = 0.d0
                  endif
                  Mkx(L-1,j) = Mkx(L-1,j) + frac_lo_m*tmpvar
                  Mkx(L,j) = Mkx(L,j) + (1-frac_lo_m)*tmpvar
               enddo
               tot_mass=0.0d0
               do j=1,icomp-idiag
                  tot_mass=tot_mass+Mkx(k,j)
               enddo
            endif 
         endif
      enddo

      do k=1,ibins
         tot_mass=0.0d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
         if (tot_mass/Nkx(k).gt.xk(k+1).or.
     &        tot_mass/Nkx(k).lt.xk(k)) then      
            print*,'ERROR in mnfix'
            print*,'bin',k,'lo',xk(k),'hi',xk(k+1),'avg',tot_mass/Nkx(k)
            print*,'tot_mass',tot_mass
            print*,'Nkx(k)',Nkx(k)
            print*,'Mkx'
            do j=1,icomp
               print*,Mkx(k,j)
            enddo
            call stop_model('out of range in mnfix',255)
         endif
      enddo

      return
      end



!@sum ion_nucl : Ion nucleation Rate calculation
!@+   from Modgil et al.(2005), JGR, vol. 110, D19205 
c
c--------------------------------------------------------------------
!@+  Mathematical Expressions for particle nucleation rate (h1,cm-3 s-1), 
!@+  nucleating H2SO4 flux (h2, cm-3 s-1), number of H2SO4 in average 
!@+  nucleating cluster (h3), number of H2O in average nucleating cluster (h4),
!@+  radius of average nucleating cluster (h5, nm), and first order loss of 
!@+  H2SO4 to particles (h6) are given below.

      subroutine ion_nucl(h2so4i,sai,t,qi,rh,h1,h2,h3,h4,h5,h6)

      implicit none
!@var  h2so4 = h2so4 concentration [molec cm^-3]
!@var  sa = aerosol surface area [um^2 cm^-3]
!@var  t = temperature [K]
!@var  q = ion formation rate [ion pairs cm^-3 s-1]
!@var  rh = relative humidity as a fraction

      real*8 h2so4i,qi,sai
      real*8 t,rh,h2so4,q,sa
      real*8 h1,h2,h3,h4,h6,h5
      
      h2so4=h2so4i
      q=qi
      sa=sai
      
      if (h2so4 .lt. 1.d5) then ! changed to 2E5 because function was diverging
         h1=0.d0
         h2=0.d0
         h3=0.d0
         h4=0.d0
         h5=0.d0
         h6=0.d0
         return
      endif
      if (h2so4 .gt. 1d8) h2so4 = 1d8
      if (sa .lt. 2.d0) sa = 2.d0
      if (sa .gt. 100.d0) sa = 100.d0
      if (t .lt. 190.d0) t = 190.d0
      if (t .gt. 300.d0) t = 300.d0
      if (rh .lt. 0.05d0) rh = 0.05d0
      if (rh .gt. 0.95d0) rh = 0.95d0
      if (q .lt. 1.d0) q = 1.d0
      if (q .gt. 50.d0) q = 50.d0      

      h6=0.000026859579119003205*SA + 1.7477354270484002d-8*q*SA + 
     &     1.5718068902491457d-8*SA**2 + 8.060796806911441d-8*SA*T + 
     &     3.904048293417882d-7*SA*Log(H2SO4) + 
     &     2.727259306977938d-7*SA*Log(RH)
      

      h3=-198.8039518313554 + 3357.132963009284*h6 - 130510.31325149858*
     &  h6**2 - 0.7093715033997716*q - 10.713505046150196*h6*q + 
     &  1103.4737682776713*h6**2*q + 0.0052565148186649*q**2 - 
     &  0.20195850414426988*h6*q**2 + 10.961027676935213*h6**2*q**2 - 
     &  26.553841269634976*RH + 2913.499196899548*h6*RH - 
     &  7558.996305824136*h6**2*RH + 0.050092880471591994*q*RH + 
     &  0.39840936335061017*h6*q*RH + 16.140386509938388*h6**2*q*RH - 
     &  0.0008159572217147427*q**2*RH - 0.02492462618304389*h6*q**2*RH + 
     &  3.2372842210428825*RH**2 + 1709.7485838150235*h6*RH**2 - 
     &  4016.182638678486*h6**2*RH**2 - 0.022142010235491123*q*RH**2 - 
     &  1.620063009925805*h6*q*RH**2-0.00028477984814528825*q**2*RH**2 + 
     &  22.136724153656015*RH**3 + 170.8982375938333*h6*RH**3 - 
     &  0.01881686723215867*q*RH**3 + 2.6974144456100957*T - 
     &  96.60591604496483*h6*T + 1772.137264721083*h6**2*T - 
     &  0.0009432251807207652*q*T + 0.06072064184950673*h6*q*T - 
     &  2.5196932894429502*h6**2*q*T - 0.000013848768113392552*q**2*T - 
     &  0.0001948394841164792*h6*q**2*T + 0.1828636512279507*RH*T - 
     &  55.135341874839185*h6*RH*T + 164.02631709083576*h6**2*RH*T + 
     &  0.001745921048607296*q*RH*T + 0.035017713828742754*h6*q*RH*T + 
     &  4.057082638293583d-6*q**2*RH*T - 0.3900441693913758*RH**2*T - 
     &  8.955078982582657*h6*RH**2*T+0.00021434974336412074*q*RH**2*T - 
     &  0.14947568974964962*RH**3*T - 0.022748394377623382*T**2 + 
     &  0.7227721843282815*h6*T**2 - 5.386480671871808*h6**2*T**2 - 
     &  0.000035250836279611095*q*T**2-0.0003363405774846326*h6*q*T**2 + 
     &  2.9254973516794257d-8*q**2*T**2 + 0.003994529829421164*RH*T**2 + 
     &  0.2074067980035454*h6*RH*T**2 - 5.136172472264946d-6*q*RH*T**2 + 
     &  0.0020603328018819816*RH**2*T**2+0.000042019279193164354*T**3 - 
     &  0.002100661388749787*h6*T**3 + 7.309966632740304d-8*q*T**3 - 
     &  0.000016323969052607556*RH*T**3+ 12.330627568298462*Log(H2SO4) + 
     &  768.3961789008589*h6*Log(H2SO4)-11568.47324943553*h6**2*
     &  Log(H2SO4)+0.14349043416922366*q*Log(H2SO4)+0.8946851157223353*
     &  h6*q*Log(H2SO4) - 68.46004143191098*h6**2*q*Log(H2SO4) - 
     &  0.0006241121793370407*q**2*Log(H2SO4)+0.011897674833907721*h6*
     &  q**2*Log(H2SO4) + 1.7860574934328677*RH*Log(H2SO4) + 
     &  316.5406316191978*h6*RH*Log(H2SO4) - 2036.825340216443*h6**2*RH*
     &  Log(H2SO4) - 0.026323914507434605*q*RH*Log(H2SO4) - 
     &  0.37505804775954393*h6*q*RH*Log(H2SO4) + 0.00003454867680790666*
     &  q**2*RH*Log(H2SO4) + 2.844877302874606*RH**2*Log(H2SO4) - 
     &  1.5845895178086176*h6*RH**2*Log(H2SO4) + 0.001732608008714275*q*
     &  RH**2*Log(H2SO4) + 0.5611003862827533*RH**3*Log(H2SO4) + 
     &  0.18033151281768975*T*Log(H2SO4) - 7.807090214680351*h6*T*
     &  Log(H2SO4) + 52.76241348342321*h6**2*T*Log(H2SO4) + 
     &  0.0011535134888316242*q*T*Log(H2SO4) + 0.0068466874844708295*h6*
     &  q*T*Log(H2SO4) - 4.231224168766194d-7*q**2*T*Log(H2SO4) - 
     &  0.10505349775719895*RH*T*Log(H2SO4) - 1.9241106452950727*h6*RH*
     &  T*Log(H2SO4) + 2.0451815715440337d-6*q*RH*T*Log(H2SO4) - 
     &  0.015299483302534183*RH**2*T*Log(H2SO4) + 0.000775633115370002*
     &  T**2*Log(H2SO4) + 0.04228608723566267*h6*T**2*Log(H2SO4) - 
     &  9.572221299803945d-7*q*T**2*Log(H2SO4) + 0.0002669785990812474*
     &  RH*T**2*Log(H2SO4) - 1.7595742533055222d-6*T**3*Log(H2SO4) - 
     &  2.7165200489046812*Log(H2SO4)**2 + 1.963036665672805*h6*
     &  Log(H2SO4)**2 + 33.88545004559797*h6**2*Log(H2SO4)**2 - 
     &  0.01647722099703982*q*Log(H2SO4)**2 - 0.08498050322218324*h6*q*
     &  Log(H2SO4)**2 + 0.00003320475358154802*q**2*Log(H2SO4)**2 + 
     &  0.50626436977659*RH*Log(H2SO4)**2 + 3.914586690682404*h6*RH*
     &  Log(H2SO4)**2 + 0.0006654515705980484*q*RH*Log(H2SO4)**2 - 
     &  0.025122425208873058*RH**2*Log(H2SO4)**2 - 0.021380868797539664*
     &  T*Log(H2SO4)**2 - 0.35405514523597137*h6*T*Log(H2SO4)**2 - 
     &  0.000020639290758942666*q*T*Log(H2SO4)**2+0.0003587951811915662*
     &  RH*T*Log(H2SO4)**2 + 7.620708111644729d-6*T**2*Log(H2SO4)**2 + 
     &  0.2196641696573127*Log(H2SO4)**3 + 1.7291708055805226*h6*
     &  Log(H2SO4)**3 + 0.00038146321602426414*q*Log(H2SO4)**3 - 
     &  0.011997306640487447*RH*Log(H2SO4)**3 + 0.0003857955558500776*T*
     &  Log(H2SO4)**3 - 0.004779827937902779*Log(H2SO4)**4

      h3=EXP(h3)



      h1=456229.3726785317 - 696754.0061755505/h3 - 
     &  8.954389043957226d7*h6 + (1.4677717736521986d8*h6)/h3 + 
     &  1867.5296995211318*q - (2798.172491398116*q)/h3 + 
     &  1500.05530404756*h6*q - (171625.68387665015*h6*q)/h3 - 
     &  5924.898937400813*T + (7657.054762118453*T)/h3 + 
     &  1.1074613167489376d6*h6*T - (1.556115054313954d6*h6*T)/h3 - 
     &  28.49469750360138*q*T + (50.46283217861981*q*T)/h3 + 
     &  1753.7251886642075*h6*q*T - (975.4050494746682*h6*q*T)/h3 + 
     &  25.45577410331681*T**2 - (26.190051353137182*T**2)/h3 - 
     &  4516.4884641323815*h6*T**2 + (5030.318058537955*h6*T**2)/h3 + 
     &  0.14403307295844858*q*T**2 - (0.2897062226196713*q*T**2)/h3 - 
     &  15.561717210345154*h6*q*T**2+(18.061492698402766*h6*q*T**2)/h3- 
     &  0.03621852394282133*T**3 + (0.02659567895141871*T**3)/h3 + 
     &  6.0696427040893886*h6*T**3 - (4.585112604562606*h6*T**3)/h3 - 
     &  0.00024142548506534756*q*T**3 + 
     &  (0.0005400351895947918*q*T**3)/h3+0.0341853829283136*h6*q*T**3- 
     &  (0.04459905758423761*h6*q*T**3)/h3-81941.3895777097*Log(H2SO4)+ 
     &  (90462.54571516594*Log(H2SO4))/h3 + 
     &  1.6483362741354546d7*h6*Log(H2SO4) - 
     &  (2.3542460879418086d7*h6*Log(H2SO4))/h3 - 
     &  319.8790300996385*q*Log(H2SO4) + 
     &  (226.84336756604796*q*Log(H2SO4))/h3 - 
     &  26455.51148372377*h6*q*Log(H2SO4) + 
     &  (60326.79915791394*h6*q*Log(H2SO4))/h3 + 
     &  1063.8621634534027*T*Log(H2SO4) - 
     &  (891.6122304868614*T*Log(H2SO4))/h3 - 
     &  203718.41653995324*h6*T*Log(H2SO4) + 
     &  (243001.1509313233*h6*T*Log(H2SO4))/h3 + 
     &  4.963508394566835*q*T*Log(H2SO4) - 
     &  (5.401955915120605*q*T*Log(H2SO4))/h3 + 
     &  25.811543727271804*h6*q*T*Log(H2SO4) - 
     &  (355.292262810142*h6*q*T*Log(H2SO4))/h3 - 
     &  4.5691538346459195*T**2*Log(H2SO4) + 
     &  (2.468423948121569*T**2*Log(H2SO4))/h3 + 
     &  830.195816015575*h6*T**2*Log(H2SO4) - 
     &  (745.2009340266254*h6*T**2*Log(H2SO4))/h3 - 
     &  0.025450884793884514*q*T**2*Log(H2SO4) + 
     &  (0.03593327468289806*q*T**2*Log(H2SO4))/h3 + 
     &  1.3302486952619124*h6*q*T**2*Log(H2SO4) - 
     &  (0.320851018492291*h6*q*T**2*Log(H2SO4))/h3 + 
     &  0.006498242231697865*T**3*Log(H2SO4) - 
     &  (0.0013540108127227987*T**3*Log(H2SO4))/h3 - 
     &  1.1150374268629073*h6*T**3*Log(H2SO4) + 
     &  (0.5953715882656012*h6*T**3*Log(H2SO4))/h3 + 
     &  0.000043146778270078744*q*T**3*Log(H2SO4) - 
     &  (0.0000735656618351653*q*T**3*Log(H2SO4))/h3 - 
     &  0.004059155580949069*h6*q*T**3*Log(H2SO4) + 
     &  (0.0028886014749212705*h6*q*T**3*Log(H2SO4))/h3 + 
     &  4912.2066512397205*Log(H2SO4)**2 - 
     &  (3022.0044642322964*Log(H2SO4)**2)/h3 - 
     &  1.0095170363345986d6*h6*Log(H2SO4)**2 + 
     &  (1.1998313113698033d6*h6*Log(H2SO4)**2)/h3 + 
     &  18.018115617541998*q*Log(H2SO4)**2 + 
     &  (7.561352267717425*q*Log(H2SO4)**2)/h3 + 
     &  3133.188395383639*h6*q*Log(H2SO4)**2 - 
     &  (3338.848456178725*h6*q*Log(H2SO4)**2)/h3 - 
     &  63.77612570995598*T*Log(H2SO4)**2 + 
     &  (20.100257708705428*T*Log(H2SO4)**2)/h3 + 
     &  12469.209965937996*h6*T*Log(H2SO4)**2 - 
     &  (11883.666207080118*h6*T*Log(H2SO4)**2)/h3 - 
     &  0.2853187730266578*q*T*Log(H2SO4)**2 + 
     &  (0.03880844916086884*q*T*Log(H2SO4)**2)/h3 - 
     &  21.708027380830988*h6*q*T*Log(H2SO4)**2 + 
     &  (28.111166854183463*h6*q*T*Log(H2SO4)**2)/h3 + 
     &  0.27388018029546984*T**2*Log(H2SO4)**2 + 
     &  (0.005568448436161079*T**2*Log(H2SO4)**2)/h3 - 
     &  50.78398945757657*h6*T**2*Log(H2SO4)**2 + 
     &  (33.25484929904544*h6*T**2*Log(H2SO4)**2)/h3 + 
     &  0.0014875097131454307*q*T**2*Log(H2SO4)**2 - 
     &  (0.0008810798671059988*q*T**2*Log(H2SO4)**2)/h3 + 
     &  0.007045532214236005*h6*q*T**2*Log(H2SO4)**2 - 
     &  (0.05521643331038819*h6*q*T**2*Log(H2SO4)**2)/h3 - 
     &  0.00038943745406964085*T**3*Log(H2SO4)**2 - 
     &  (0.00015308439633013156*T**3*Log(H2SO4)**2)/h3 + 
     &  0.06817761559886615*h6*T**3*Log(H2SO4)**2 - 
     &  (0.019520711865759842*h6*T**3*Log(H2SO4)**2)/h3 - 
     &  2.554240267270117d-6*q*T**3*Log(H2SO4)**2 + 
     &  (2.526957991768314d-6*q*T**3*Log(H2SO4)**2)/h3 + 
     &  0.00011978108924107446*h6*q*T**3*Log(H2SO4)**2 - 
     &  98.09290363960072*Log(H2SO4)**3 + 
     &  (3.8697214915504805*Log(H2SO4)**3)/h3 + 
     &  20559.743351544734*h6*Log(H2SO4)**3 - 
     &  (18687.867717989076*h6*Log(H2SO4)**3)/h3 - 
     &  0.33118998078526996*q*Log(H2SO4)**3 - 
     &  (0.6955616584928909*q*Log(H2SO4)**3)/h3 - 
     &  92.45613590327787*h6*q*Log(H2SO4)**3 + 
     &  (8.948488896745264*h6*q*Log(H2SO4)**3)/h3 + 
     &  1.273863980203035*T*Log(H2SO4)**3 + 
     &  (0.38051947439108685*T*Log(H2SO4)**3)/h3 - 
     &  253.8228380167286*h6*T*Log(H2SO4)**3 + 
     &  (171.19178070061787*h6*T*Log(H2SO4)**3)/h3 + 
     &  0.005376622635136611*q*T*Log(H2SO4)**3 + 
     &  (0.006658042079999949*q*T*Log(H2SO4)**3)/h3 + 
     &  0.8230892248689082*h6*q*T*Log(H2SO4)**3 - 
     &  (0.058021973134188776*h6*q*T*Log(H2SO4)**3)/h3 - 
     &  0.005471105747328322*T**2*Log(H2SO4)**3 - 
     &  (0.003699167378319841*T**2*Log(H2SO4)**3)/h3 + 
     &  1.0332570198673408*h6*T**2*Log(H2SO4)**3 - 
     &  (0.38600322062471826*h6*T**2*Log(H2SO4)**3)/h3 - 
     &  0.000028584962254144285*q*T**2*Log(H2SO4)**3 - 
     &  (0.000016062832568931284*q*T**2*Log(H2SO4)**3)/h3 - 
     &  0.001819022861961261*h6*q*T**2*Log(H2SO4)**3 + 
     &  7.779696142330516d-6*T**3*Log(H2SO4)**3 + 
     &  (8.513832613126146d-6*T**3*Log(H2SO4)**3)/h3 - 
     &  0.0013866861532058*h6*T**3*Log(H2SO4)**3 + 
     &  4.981099027173453d-8*q*T**3*Log(H2SO4)**3 + 
     &  178821.264938151*Log(RH) - (634125.5419575685*Log(RH))/h3 - 
     &  2.310391119522488d7*h6*Log(RH) + 
     &  (3.564836811852359d7*h6*Log(RH))/h3 - 
     &  744.9440395010726*q*Log(RH)+(1019.4951846825716*q*Log(RH))/h3- 
     &  28008.93484108728*h6*q*Log(RH) - 
     &  (7136.405275749431*h6*q*Log(RH))/h3 - 
     &  2275.393426900043*T*Log(RH) + (9231.8244493922*T*Log(RH))/h3 + 
     &  237571.37589421016*h6*T*Log(RH) - 
     &  (295411.0553218543*h6*T*Log(RH))/h3 + 
     &  4.451622613453289*q*T*Log(RH) - 
     &  (4.150562571935793*q*T*Log(RH))/h3 + 
     &  595.8398944235868*h6*q*T*Log(RH) - 
     &  (210.09033695512161*h6*q*T*Log(RH))/h3 + 
     &  9.678343416541734*T**2*Log(RH) - 
     &  (45.69882758602488*T**2*Log(RH))/h3 - 
     &  747.9878102876322*h6*T**2*Log(RH) + 
     &  (649.4424346377774*h6*T**2*Log(RH))/h3 + 
     &  0.005153649541788936*q*T**2*Log(RH) - 
     &  (0.0023655719118541533*q*T**2*Log(RH))/h3 - 
     &  2.848932474900831*h6*q*T**2*Log(RH) + 
     &  (0.2889680076326506*h6*q*T**2*Log(RH))/h3 - 
     &  0.013824587420637588*T**3*Log(RH) + 
     &  (0.07673955581386811*T**3*Log(RH))/h3 + 
     &  0.6813253055163031*h6*T**3*Log(RH) - 
     &  (0.2663343569482817*h6*T**3*Log(RH))/h3 - 
     &  0.00004364745106737275*q*T**3*Log(RH) - 
     &  (9.434880576275503d-6*q*T**3*Log(RH))/h3 + 
     &  0.0030430214338531205*h6*q*T**3*Log(RH) + 
     &  (0.00371157673745342*h6*q*T**3*Log(RH))/h3 - 
     &  31842.277946766168*Log(H2SO4)*Log(RH) + 
     &  (68134.57926920953*Log(H2SO4)*Log(RH))/h3 + 
     &  3.493985612508899d6*h6*Log(H2SO4)*Log(RH) - 
     &  (4.587173661378969d6*h6*Log(H2SO4)*Log(RH))/h3 + 
     &  190.19310993267845*q*Log(H2SO4)*Log(RH) - 
     &  (233.68205749312318*q*Log(H2SO4)*Log(RH))/h3 - 
     &  1547.3248018686509*h6*q*Log(H2SO4)*Log(RH) + 
     &  (4045.133810971437*h6*q*Log(H2SO4)*Log(RH))/h3 + 
     &  405.31725032582267*T*Log(H2SO4)*Log(RH) - 
     &  (1057.655396923827*T*Log(H2SO4)*Log(RH))/h3 - 
     &  33903.90722480216*h6*T*Log(H2SO4)*Log(RH) + 
     &  (32965.76708786528*h6*T*Log(H2SO4)*Log(RH))/h3 - 
     &  1.5034587825111836*q*T*Log(H2SO4)*Log(RH) + 
     &  (1.106270806420077*q*T*Log(H2SO4)*Log(RH))/h3 - 
     &  25.549725610666396*h6*q*T*Log(H2SO4)*Log(RH) + 
     &  (20.09072468984517*h6*q*T*Log(H2SO4)*Log(RH))/h3 - 
     &  1.7253372988643991*T**2*Log(H2SO4)*Log(RH) + 
     &  (5.612615323997249*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  95.62445916373963*h6*T**2*Log(H2SO4)*Log(RH) - 
     &  (48.96291394056088*h6*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  0.0019300142376745193*q*T**2*Log(H2SO4)*Log(RH) + 
     &  (0.00015584582354154648*q*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  0.19177967246162606*h6*q*T**2*Log(H2SO4)*Log(RH) - 
     &  (0.16871358892301552*h6*q*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  0.00246781658651071*T**3*Log(H2SO4)*Log(RH) - 
     &  (0.01009444649935961*T**3*Log(H2SO4)*Log(RH))/h3 - 
     &  0.06589616928988701*h6*T**3*Log(H2SO4)*Log(RH) - 
     &  (0.01596274185896012*h6*T**3*Log(H2SO4)*Log(RH))/h3 + 
     &  4.109041091272658d-6*q*T**3*Log(H2SO4)*Log(RH) + 
     &  (2.041743340493583d-7*q*T**3*Log(H2SO4)*Log(RH))/h3 - 
     &  0.0001572125066332154*h6*q*T**3*Log(H2SO4)*Log(RH) + 
     &  1932.4011188760162*Log(H2SO4)**2*Log(RH) - 
     &  (781.4228432455719*Log(H2SO4)**2*Log(RH))/h3 - 
     &  181291.39888342103*h6*Log(H2SO4)**2*Log(RH) + 
     &  (210184.24050416853*h6*Log(H2SO4)**2*Log(RH))/h3 - 
     &  13.666906354559005*q*Log(H2SO4)**2*Log(RH) + 
     &  (15.930151624229865*q*Log(H2SO4)**2*Log(RH))/h3 + 
     &  305.7703887043172*h6*q*Log(H2SO4)**2*Log(RH) - 
     &  (391.71164595591085*h6*q*Log(H2SO4)**2*Log(RH))/h3 - 
     &  24.621398404913865*T*Log(H2SO4)**2*Log(RH) + 
     &  (20.079316822912197*T*Log(H2SO4)**2*Log(RH))/h3 + 
     &  1658.5845815289838*h6*T*Log(H2SO4)**2*Log(RH) - 
     &  (1404.5530139788664*h6*T*Log(H2SO4)**2*Log(RH))/h3 + 
     &  0.11823318365930788*q*T*Log(H2SO4)**2*Log(RH) - 
     &  (0.08056635668752316*q*T*Log(H2SO4)**2*Log(RH))/h3 - 
     &  1.0034435154436134*h6*q*T*Log(H2SO4)**2*Log(RH) + 
     &  (1.6021411957310931*h6*q*T*Log(H2SO4)**2*Log(RH))/h3 + 
     &  0.10491613716181539*T**2*Log(H2SO4)**2*Log(RH) - 
     &  (0.1492097929648805*T**2*Log(H2SO4)**2*Log(RH))/h3 - 
     &  4.0934834940398215*h6*T**2*Log(H2SO4)**2*Log(RH) + 
     &  (1.896735315742051*h6*T**2*Log(H2SO4)**2*Log(RH))/h3 - 
     &  0.00022870693934890214*q*T**2*Log(H2SO4)**2*Log(RH) + 
     &  (0.000013175311622036916*q*T**2*Log(H2SO4)**2*Log(RH))/h3- 
     &  0.00231560013471937*h6*q*T**2*Log(H2SO4)**2*Log(RH) - 
     &  0.00015023186988448507*T**3*Log(H2SO4)**2*Log(RH) + 
     &  (0.0003388679191513823*T**3*Log(H2SO4)**2*Log(RH))/h3 + 
     &  0.0015784630518057667*h6*T**3*Log(H2SO4)**2*Log(RH) - 
     &  1.0056424638350849d-7*q*T**3*Log(H2SO4)**2*Log(RH) - 
     &  39.74503165968021*Log(H2SO4)**3*Log(RH) - 
     &  (68.29220045487304*Log(H2SO4)**3*Log(RH))/h3 + 
     &  3236.6794930454566*h6*Log(H2SO4)**3*Log(RH) - 
     &  (2511.4313482742737*h6*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.3044230466506939*q*Log(H2SO4)**3*Log(RH) - 
     &  (0.321694018076018*q*Log(H2SO4)**3*Log(RH))/h3 - 
     &  7.83632554091337*h6*q*Log(H2SO4)**3*Log(RH) + 
     &  (1.2615699211576334*h6*q*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.5071743553774433*T*Log(H2SO4)**3*Log(RH) + 
     &  (0.6993184748948513*T*Log(H2SO4)**3*Log(RH))/h3 - 
     &  28.064969097031586*h6*T*Log(H2SO4)**3*Log(RH) + 
     &  (12.010607450658807*h6*T*Log(H2SO4)**3*Log(RH))/h3 - 
     &  0.0027511166431233884*q*T*Log(H2SO4)**3*Log(RH) + 
     &  (0.0015662272224653852*q*T*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.03815288545975011*h6*q*T*Log(H2SO4)**3*Log(RH) - 
     &  0.0021639947693271114*T**2*Log(H2SO4)**3*Log(RH) - 
     &  (0.0017767382579465334*T**2*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.059762491440998364*h6*T**2*Log(H2SO4)**3*Log(RH) + 
     &  6.130388542051968d-6*q*T**2*Log(H2SO4)**3*Log(RH) + 
     &  3.1017462628108047d-6*T**3*Log(H2SO4)**3*Log(RH) - 
     &  7689.12786121461*Log(RH)**2-(25856.380361656236*Log(RH)**2)/h3- 
     &  1.3915437178042033d6*h6*Log(RH)**2 + 
     &  (1.9928171556602388d6*h6*Log(RH)**2)/h3 + 
     &  53.79216985973448*q*Log(RH)**2 - 
     &  (134.8787674139501*q*Log(RH)**2)/h3 - 
     &  1440.4837667461468*h6*q*Log(RH)**2 - 
     &  (3141.0130344871345*h6*q*Log(RH)**2)/h3 + 
     &  94.7545268444358*T*Log(RH)**2 + 
     &  (342.5857682014369*T*Log(RH)**2)/h3 + 
     &  11341.294070320842*h6*T*Log(RH)**2 - 
     &  (18070.160441620144*h6*T*Log(RH)**2)/h3 - 
     &  0.44696864139255954*q*T*Log(RH)**2 + 
     &  (1.1199451205710647*q*T*Log(RH)**2)/h3 + 
     &  19.01601691757823*h6*q*T*Log(RH)**2 + 
     &  (31.604766584775273*h6*q*T*Log(RH)**2)/h3 - 
     &  0.3470127737474881*T**2*Log(RH)**2 - 
     &  (1.4877879003131458*T**2*Log(RH)**2)/h3 - 
     &  18.727169614346867*h6*T**2*Log(RH)**2 + 
     &  (42.5448759900235*h6*T**2*Log(RH)**2)/h3 + 
     &  0.0007727701515171791*q*T**2*Log(RH)**2 - 
     &  (0.0032041999804810098*q*T**2*Log(RH)**2)/h3 - 
     &  0.055631814243943*h6*q*T**2*Log(RH)**2 - 
     &  (0.07696588905487615*h6*q*T**2*Log(RH)**2)/h3 + 
     &  0.00033600378567117826*T**3*Log(RH)**2 + 
     &  (0.002069337837641118*T**3*Log(RH)**2)/h3 - 
     &  0.015788853514098967*h6*T**3*Log(RH)**2 - 
     &  (0.009301939734036398*h6*T**3*Log(RH)**2)/h3 + 
     &  5.572064421478949d-7*q*T**3*Log(RH)**2 + 
     &  (5.677455521083783d-6*q*T**3*Log(RH)**2)/h3 - 
     &  0.00003507762649727548*h6*q*T**3*Log(RH)**2 + 
     &  1392.9578972108989*Log(H2SO4)*Log(RH)**2 + 
     &  (1911.734041113941*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  179329.36146873524*h6*Log(H2SO4)*Log(RH)**2 - 
     &  (119328.53674704138*h6*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  6.462104618811676*q*Log(H2SO4)*Log(RH)**2 + 
     &  (11.294074621589449*q*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  13.39033080899747*h6*q*Log(H2SO4)*Log(RH)**2 - 
     &  (8.52831107206456*h6*q*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  16.8958062365253*T*Log(H2SO4)*Log(RH)**2 - 
     &  (24.81735534395592*T*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  1464.394638508739*h6*T*Log(H2SO4)*Log(RH)**2 + 
     &  (981.4137808275404*h6*T*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.051952833287610724*q*T*Log(H2SO4)*Log(RH)**2 - 
     &  (0.061982178004231*q*T*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  0.7507183329339466*h6*q*T*Log(H2SO4)*Log(RH)**2 - 
     &  (0.05082083499340442*h6*q*T*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.060317519684020285*T**2*Log(H2SO4)*Log(RH)**2 + 
     &  (0.10238074986938017*T**2*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  2.654746963853928*h6*T**2*Log(H2SO4)*Log(RH)**2 - 
     &  (1.9505367684260284*h6*T**2*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  0.00008828022840243558*q*T**2*Log(H2SO4)*Log(RH)**2 - 
     &  (7.843133775464917d-6*q*T**2*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.00439346543661198*h6*q*T**2*Log(H2SO4)*Log(RH)**2 - 
     &  0.00005488176049657228*T**3*Log(H2SO4)*Log(RH)**2 - 
     &  (0.00012357549613231235*T**3*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.000877040667994904*h6*T**3*Log(H2SO4)*Log(RH)**2 - 
     &  4.394231627146249d-8*q*T**3*Log(H2SO4)*Log(RH)**2 - 
     &  72.11406227014747*Log(H2SO4)**2*Log(RH)**2 - 
     &  (16.86445272085156*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  6535.765336704547*h6*Log(H2SO4)**2*Log(RH)**2 + 
     &  (941.2813201741736*h6*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  0.2428031989628447*q*Log(H2SO4)**2*Log(RH)**2 - 
     &  (0.3050220578412721*q*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  3.3649457326647947*h6*q*Log(H2SO4)**2*Log(RH)**2 + 
     &  (0.9000970248172522*h6*q*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  0.8503624487797692*T*Log(H2SO4)**2*Log(RH)**2 + 
     &  (0.2527834506442179*T*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  50.73348890297342*h6*T*Log(H2SO4)**2*Log(RH)**2 - 
     &  (4.356161470448808*h6*T*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  0.0018342607279861292*q*T*Log(H2SO4)**2*Log(RH)**2 + 
     &  (0.0020365272689109475*q*T*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  0.03276629834613374*h6*q*T*Log(H2SO4)**2*Log(RH)**2 - 
     &  0.0028420227040081517*T**2*Log(H2SO4)**2*Log(RH)**2 - 
     &  (0.0009263932662681723*T**2*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  0.08964774960777555*h6*T**2*Log(H2SO4)**2*Log(RH)**2 + 
     &  3.036637047548649d-6*q*T**2*Log(H2SO4)**2*Log(RH)**2 + 
     &  2.0594310710626316d-6*T**3*Log(H2SO4)**2*Log(RH)**2 + 
     &  0.9950338851031777*Log(H2SO4)**3*Log(RH)**2 - 
     &  (0.622965174549953*Log(H2SO4)**3*Log(RH)**2)/h3 + 
     &  48.41845482837981*h6*Log(H2SO4)**3*Log(RH)**2 - 
     &  (1.7738934087048919*h6*Log(H2SO4)**3*Log(RH)**2)/h3 - 
     &  0.0022234755298572505*q*Log(H2SO4)**3*Log(RH)**2 - 
     &  (0.002437957119890622*q*Log(H2SO4)**3*Log(RH)**2)/h3 + 
     &  0.08069009533680993*h6*q*Log(H2SO4)**3*Log(RH)**2 - 
     &  0.010935011617216193*T*Log(H2SO4)**3*Log(RH)**2 + 
     &  (0.003456846985810427*T*Log(H2SO4)**3*Log(RH)**2)/h3 - 
     &  0.24482068250824904*h6*T*Log(H2SO4)**3*Log(RH)**2 + 
     &  0.000011334503487127534*q*T*Log(H2SO4)**3*Log(RH)**2 + 
     &  0.000029425270779265584*T**2*Log(H2SO4)**3*Log(RH)**2

      h1=EXP(h1)


      h2=-32043.03148295406 + 59725.428570008815/h3 +  
     & 7.128537634261564d6*h6 - (1.3833467233343722d7*h6)/h3+
     & 33.63110252227136*q-(48.61215633992165*q)/h3 - 16602.414377611287
     & *h6*q + (40754.788181739124* h6*q)/h3 -2.3397851800516185*q**2 + 
     & (4.426964073992281*q**2)/h3 + 18.971418767591036*h6*q**2 +
     & (60.446718551038344*h6*q**2)/h3+396.33752593131607*T-
     & (650.8601684277011*T)/h3-83753.25337253512*h6*T+
     & (140932.8771905448*h6*T)/h3 + 1.7779514612590905*q*T -  
     & (4.1201474547289845*q*T)/h3 - 14.23100399324848*h6*q*T - 
     & (199.64146093004214*h6*q*T)/h3 + 0.017586270137944494*q**2*T - 
     & (0.040669111407675304*q**2*T)/h3 - 0.5190592639152767*h6*q**2*T + 
     & (0.4152291779457336*h6*q**2*T)/h3 - 1.6239261240452716*T**2 + 
     & (2.290370025362425*T**2)/h3 + 323.71145981054207*h6*T**2 -  
     & (443.70978129586814*h6*T**2)/h3 - 0.016944490878752067*q*T**2 + 
     & (0.03190748579500371*q*T**2)/h3 + 0.9921294086888162*h6*q*T**2 -  
     & (0.017807643647025154*h6*q*T**2)/h3-0.000026182356192850867*q**2*
     & T**2+(0.00018413948144606674*q**2*T**2)/h3+0.001418259144142961*
     & h6*q**2*T**2 - (0.0005314099952660705*h6*q**2*T**2)/h3 + 
     & 0.002203350630670886*T**3-(0.002628197592847762*T**3)/h3 - 
     & 0.41364399538414304*h6*T**3+(0.4236315079925672*h6*T**3)/h3 + 
     & 0.000036034470518775974*q*T**3-(0.00005072877953687334*q*T**3)/h3
     & -0.0025779163615301994*h6*q*T**3-(0.0001801635475946943*h6*q*
     & T**3)/h3 - 1.0397496813871261d-8*q**2*T**3 -  
     & (4.794510915054545d-7*q**2*T**3)/h3 + 1.0037940606886203d-6*h6*
     & q**2*T**3 +4057.4863055406972*Log(H2SO4) - (4857.4895775507075*
     & Log(H2SO4))/h3 - 873923.0975121224*h6*Log(H2SO4) + 
     & (1.3884698205066123d6*h6*Log(H2SO4))/h3-  
     & 16.567334471680468*q*Log(H2SO4)+(39.06695157867744*q*
     & Log(H2SO4))/h3 + 
     & 4033.3161579197213*h6*q*Log(H2SO4) - (4837.3892389479715*h6*q* 
     & Log(H2SO4))/h3 + 0.32643050327864925*q**2*Log(H2SO4) -  
     & (0.29168894526478634*q**2*Log(H2SO4))/h3 +   3.0197469938934343* 
     & h6*q**2*Log(H2SO4) - (13.620934035683936*h6*q**2*Log(H2SO4))/h3 -  
     & 50.64853826213881*T*Log(H2SO4) + (46.190674056264115*T*
     & Log(H2SO4))/h3 +10309.186386370407*h6*T*Log(H2SO4) - 
     & (13393.112078195805*h6*T* 
     & Log(H2SO4))/h3 - 0.05988575301198189*q*T*Log(H2SO4) + 
     & (0.12991016720908877*q*T*Log(H2SO4))/h3-25.211481437032155*h6*q*
     & T*Log(H2SO4) + (25.76780209283505*h6*q*T*Log(H2SO4))/h3 -  
     & 0.0025433151017881413*q**2*T*Log(H2SO4)-(0.00003896417606619728* 
     & q**2*T*Log(H2SO4))/h3+0.018638367565640454*h6*q**2*T*Log(H2SO4)-  
     & (0.009466280650899394*h6*q**2*T*Log(H2SO4))/h3+
     & 0.20947015213812362*T**2*Log(H2SO4)-(0.1282261301435015*T**2*
     & Log(H2SO4))/h3-40.05118528139879*h6*T**2*Log(H2SO4)+
     & (37.82992227601891* h6*T**2*Log(H2SO4))/h3+0.001426085690256897*
     & q*T**2*Log(H2SO4)-(0.0026110978166799066*q*T**2*Log(H2SO4))/h3 -
     &  0.0017889111275040672*h6*q*T**2*Log(H2SO4) - 
     & (0.0010479732179162718* h6*q*T**2*Log(H2SO4))/h3 + 
     & 4.179579131654065d-6 * q**2*T**2*Log(H2SO4) +
     & (7.970218003902365d-6 * q**2 * 
     & T**2*Log(H2SO4))/h3-0.000116249941623789*h6*q**2*T**2*Log(H2SO4)-  
     & 0.0002870253324673844*T**3*Log(H2SO4)+(0.00009186981667688776*
     & T**3*Log(H2SO4))/h3 + 0.05147918448072947*h6*T**3*Log(H2SO4) -  
     & (0.028081229965757286*h6*T**3*Log(H2SO4))/h3 - 
     & 3.560981861287522d-6* q* 
     & T**3*Log(H2SO4) + (4.961243888005889d-6*q*T**3*Log(H2SO4))/h3 +  
     & 0.0001420993695424454*h6*q*T**3*Log(H2SO4)+1.2719812342415112d-9*
     & q**2*T**3*Log(H2SO4) - 124.19545742393447*Log(H2SO4)**2 +   
     & (41.54331926307027*Log(H2SO4)**2)/h3 + 26490.80092234763*h6*
     & Log(H2SO4)**2- (31462.37707525288*h6*Log(H2SO4)**2)/h3 +  
     & 0.8500524958176608*q* Log(H2SO4)**2 - (3.2674516205233317*q*
     & Log(H2SO4)**2)/h3 - 170.8115272895333*h6*q*Log(H2SO4)**2 + 
     & (120.16514929908757*h6*q* Log(H2SO4)**2)/h3 - 
     & 0.011268864148405246*q**2*Log(H2SO4)**2 +(0.019897271615564007*
     & q**2*Log(H2SO4)**2)/h3 -  0.2574159414899909* 
     & h6*q**2*Log(H2SO4)**2 + (0.5058681955079796*h6*q**2*
     & Log(H2SO4)**2)/h3 +1.5616031047481935*T*Log(H2SO4)**2 +  
     & (0.04762580087464638*T* Log(H2SO4)**2)/h3 - 313.57384684195245*
     & h6*T*Log(H2SO4)**2+(270.50302801880633*h6*T*Log(H2SO4)**2)/h3 -
     &  0.0026826956178692776* q* 
     & T*Log(H2SO4)**2 + (0.02263655181544639*q*T*Log(H2SO4)**2)/h3 +   
     & 1.4230688700541085*h6*q*T*Log(H2SO4)**2 - (0.6649649390384627*
     &  h6* q* T*Log(H2SO4)**2)/h3 + 0.00009215121101630449*q**2*T*
     & Log(H2SO4)**2 -  
     & (0.00010215400234653555*q**2*T*Log(H2SO4)**2)/h3 +  
     & 0.0010136511879910042*h6*q**2*T*Log(H2SO4)**2 -  
     & 0.006504837589925847*T**2*Log(H2SO4)**2 -  (0.002702979699750996* 
     & T**2*Log(H2SO4)**2)/h3 + 1.223462360247441*h6*T**2*Log(H2SO4)**2- 
     & (0.5623561810520002*h6*T**2*Log(H2SO4)**2)/h3 - 
     & 0.000024272415265506407*q*T**2*Log(H2SO4)**2 - 
     & (0.0000317372505770328*q*T**2*Log(H2SO4)**2)/h3 -  
     & 0.0028928533476792*h6*q*T**2*Log(H2SO4)**2-1.7589107450560238d-7* 
     & q**2*T**2*Log(H2SO4)**2 + 8.98009092200891d-6*T**3*Log(H2SO4)**2+   
     & (7.173660190761036d-6*T**3*Log(H2SO4)**2)/h3 -  
     & 0.0015801362599946241*h6*T**3*Log(H2SO4)**2+8.220297383016583d-8*
     & q*T**3*Log(H2SO4)**2 - 12589.220049398413*Log(RH)+
     & (71533.62210328173*Log(RH))/h3 + 1.0678104434003264d6*h6*Log(RH)-
     & (1.6102655953002474d6*h6*Log(RH))/h3+51.82639715490156*q*Log(RH)-
     & (38.71320196661908*q*Log(RH))/h3+3767.8963041701336*h6*q*Log(RH)-
     & (2025.2741146328435*h6*q*Log(RH))/h3+0.43253945689885376*q**2*
     & Log(RH) - (0.5285905724764535*q**2*Log(RH))/h3+2.800088035800146*
     & h6*q**2*Log(RH)+(5.830987444741845*h6*q**2*Log(RH))/h3+
     & 156.76303859222733*T*Log(RH) -  (1050.9788795963977*T* 
     & Log(RH))/h3 + 2968.558570914762*h6*T*Log(RH)-  
     & (13144.694726104603*h6*T*Log(RH))/h3-0.3939836553978611*q*T*
     & Log(RH) - (0.15149542273585964*q*T*Log(RH))/h3-52.00514647221519*
     & h6*q*T*Log(RH)+(33.50499777824741*h6*q*T*Log(RH))/h3-
     & 0.0044939367750510195*q**2*T*Log(RH)+  
     & (0.0037387212699065446*q**2*T*Log(RH))/h3 -  
     & 0.019659867844446746*h6*q**2*T*Log(RH) - (0.013728505197733986*
     & h6*q**2*T*Log(RH))/h3-0.6541703069347068*T**2*Log(RH) + 
     & (5.190500347012039*T**2*Log(RH))/h3-  
     & 87.05017506628032*h6*T**2*Log(RH) +  (179.83498689899525*h6*T**2* 
     & Log(RH))/h3 + 0.00042182566738198824*q*T**2*Log(RH) +  
     & (0.004271845253967774*q*T**2*Log(RH))/h3+0.12589872273911631*h6*
     & q* T**2*Log(RH) - (0.08117178992783239*h6*q*T**2*Log(RH))/h3 +   
     & 0.000010942253625927065*q**2*T**2*Log(RH) - 
     & (3.4136292496834907d-6*q**2* 
     & T**2*Log(RH))/h3 + 3.7003812518449353d-6*h6*q**2*T**2*Log(RH) +  
     & 0.0009227426647911644*T**3*Log(RH) - (0.008600685587907114*T**3* 
     & Log(RH))/h3 + 0.23093344962838908*h6*T**3*Log(RH) -  
     & (0.38790097045592176*h6*T**3*Log(RH))/h3 + 1.2941739001871342d-6*
     & q*T**3*Log(RH)-(0.000013459445557723903*q*T**3*Log(RH))/h3 +  
     & 0.00016881487161048804* h6*q*T**3*Log(RH) + 2.609280671396356d-9*
     & q**2*T**3 *Log(RH)+2130.4592878663534*Log(H2SO4)*Log(RH) -  
     & (4313.3562495734495*Log(H2SO4)*Log(RH))/h3 - 
     & 291848.82747478905*h6*Log(H2SO4)*Log(RH) +  
     & (488820.18312982953*h6*Log(H2SO4)*Log(RH))/h3-9.405300971012823*
     & q*Log(H2SO4)*Log(RH) + (12.85263348658218*q*Log(H2SO4)*
     * Log(RH))/h3+102.49140755096978*h6*q*Log(H2SO4)*Log(RH) - 
     & (289.1347848442807*h6*q*Log(H2SO4)*Log(RH))/h3 - 
     & 0.01500268773941549*q**2*Log(H2SO4)*Log(RH) +  
     & (0.018600134825856034*q**2*Log(H2SO4)*Log(RH))/h3 - 
     & 0.07717355560425733* 
     & h6*q**2*Log(H2SO4)*Log(RH)-(0.161938896012148*h6*q**2*Log(H2SO4)* 
     & Log(RH))/h3 - 27.134659704074295*T*Log(H2SO4)*Log(RH) + 
     & (69.2483692103918*T*Log(H2SO4)*Log(RH))/h3 + 1882.179099085481*
     & h6*T*Log(H2SO4)*Log(RH)-(3332.724995554323*h6*T*Log(H2SO4)*
     & Log(RH))/h3 +   0.0786126279415837*q* 
     & T*Log(H2SO4)*Log(RH) - (0.12353477014379059*q*T*Log(H2SO4)*
     & Log(RH))/h3+2.221848005253306*h6*q*T*Log(H2SO4)*Log(RH) +  
     & (0.532770625490143*h6 *q*T*Log(H2SO4)*Log(RH))/h3 + 
     & 0.00020323051093599427*q**2*T* Log(H2SO4)* 
     & Log(RH) - (0.00014374306638210478*q**2*T*Log(H2SO4)*Log(RH))/h3 +   
     & 0.0006995787862129086*h6*q**2*T*Log(H2SO4)*Log(RH) +   
     & 0.11535250161710112*T**2*Log(H2SO4)*Log(RH)-(0.37254120769567123* 
     & T**2*Log(H2SO4)*Log(RH))/h3 +   0.3210533828417492*h6*T**2* 
     & Log(H2SO4)*Log(RH)+(4.899886058456062*h6*T**2*Log(H2SO4)*
     & Log(RH))/h3-0.00013223802083620908*q*T**2*Log(H2SO4)*Log(RH) +  
     & (0.00031364246871883094*q*T**2*Log(H2SO4)*Log(RH))/h3 - 
     & 0.013017818218354784*h6*q*T**2*Log(H2SO4)*Log(RH) - 
     & 6.285230956812909d-7* 
     & q**2*T**2*Log(H2SO4)*Log(RH) -0.0001644893214543928*T**3* 
     & Log(H2SO4)*Log(RH)+(0.0006683187277334718*T**3*Log(H2SO4)*
     & Log(RH))/h3-0.012564516901538326*h6*T**3*Log(H2SO4)*Log(RH)-
     & 1.2575655284907035d-7*q*T**3*Log(H2SO4)*Log(RH) - 
     & 79.49703918541012*Log(H2SO4)**2*Log(RH) -  
     & (101.22507090731304*Log(H2SO4)**2*Log(RH))/h3 + 
     & 12719.306036991808*h6*Log(H2SO4)**2*Log(RH)-(7977.5339222745315*
     & h6*Log(H2SO4)**2*Log(RH))/h3 
     & +0.3744372178099974*q*Log(H2SO4)**2*Log(RH)-(0.025682479965146*q* 
     & Log(H2SO4)**2*Log(RH))/h3-22.19218619809272*h6*q*Log(H2SO4)**2* 
     & Log(RH) + (5.910834971937343*h6*q*Log(H2SO4)**2*Log(RH))/h3 -   
     & 0.00037981877987097496*q**2*Log(H2SO4)**2*Log(RH) +  
     & (0.00032985668220561373*q**2*Log(H2SO4)**2*Log(RH))/h3 - 
     & 0.0003940881289017262*h6*q**2*Log(H2SO4)**2*Log(RH) + 
     & 1.0255615837685275*T*Log(H2SO4)**2*Log(RH) +  
     & (1.0124457331884233*T* 
     & Log(H2SO4)**2*Log(RH))/h3-111.19643387882614*h6*T*Log(H2SO4)**2* 
     & Log(RH) + (37.058545749337576*h6*T*Log(H2SO4)**2*Log(RH))/h3 - 
     & 0.0034094873321621174*q*T*Log(H2SO4)**2*Log(RH) +  
     & (0.000015597253696421123*q*T*Log(H2SO4)**2*Log(RH))/h3 +   
     & 0.10750304600542902*h6*q*T*Log(H2SO4)**2*Log(RH) + 
     &1.8808923253076581d-6* 
     & q**2*T*Log(H2SO4)**2*Log(RH) - 0.0044066707400722644*T**2*
     & Log(H2SO4)**2*Log(RH) -  (0.0025163539282379434*T**2*
     & Log(H2SO4)**2*Log(RH))/h3 + 0.24020207161395796*h6*T**2*
     & Log(H2SO4)**2*Log(RH)+7.634986243805219d-6*q*T**2*Log(H2SO4)**2*
     & Log(RH) + 6.328401246898811d-6*T**3*Log(H2SO4)**2*Log(RH) + 
     & 3630.6862033225625*Log(RH)**2-(1075.4966438125716*Log(RH)**2)/h3+
     & 54546.69751557024*h6*Log(RH)**2-(49710.530231480734*h6*
     & Log(RH)**2)/h3 - 1.5893360668636096*q*Log(RH)**2 +   
     & (7.970183065343727*q*Log(RH)**2)/h3 + 400.88335315886525*h6*q*
     & Log(RH)**2 -(219.0329650802132*h6*q*Log(RH)**2)/h3 +  
     & 0.015009992625400783*q**2*Log(RH)**2 + (0.0036804085191588995*
     & q**2*Log(RH)**2)/h3+0.4203732847045355*h6*q**2*Log(RH)**2 +  
     & (0.20391771437862574*h6*q**2*Log(RH)**2)/h3 - 44.579109086853684*
     & T*Log(RH)**2 + (8.53054910480424*T* 
     & Log(RH)**2)/h3 + 76.04921125760566*h6*T*Log(RH)**2 -  
     & (329.28669019015933*h6*T*Log(RH)**2)/h3 - 0.019940534315381526*q
     & *T*Log(RH)**2 - (0.050640623913744444*q*T*Log(RH)**2)/h3 - 
     & 1.5765203107981385*h6*q*T*Log(RH)**2 +   (0.3694218342124642*h6*
     & q*T* Log(RH)**2)/h3-0.0001222138921423729*q**2*T*Log(RH)**2 - 
     & (0.000049403831741944986*q**2* 
     & T*Log(RH)**2)/h3 - 0.0025870767829077484*h6*q**2*T*Log(RH)**2 +  
     & 0.17463801499094986*T**2*Log(RH)**2 -  (0.02046987546015939*T**2* 
     & Log(RH)**2)/h3 - 1.3780954635968534*h6*T**2*Log(RH)**2 +  
     & (2.0745682824062985*h6*T**2*Log(RH)**2)/h3 +  
     & 0.00015350888582331981* q*T**2*Log(RH)**2 + 
     & (0.00011285744022151752*q*T**2*Log(RH)**2)/h3 +  
     & 0.00004952110557031714*h6*q*T**2*Log(RH)**2+3.968468979853371d-7*
     & q**2*T**2*Log(RH)**2 - 0.00021374165659919817*T**3*Log(RH)**2 +  
     & (0.000023714077438108623*T**3*Log(RH)**2)/h3 -  
     & 0.0010671420677054659* h6* 
     & T**3*Log(RH)**2 - 1.8562790214558605d-7*q*T**3*Log(RH)**2 - 
     & 279.60556688841575*Log(H2SO4)*Log(RH)**2 +  (44.483779255378025* 
     & Log(H2SO4)*Log(RH)**2)/h3 -11085.947063809159*h6*Log(H2SO4)*
     & Log(RH)**2 +(10910.046235786824*h6*Log(H2SO4)*Log(RH)**2)/h3 + 
     & 0.47924678115681674*q*Log(H2SO4)*Log(RH)**2-(0.3759508982199312*
     & q*Log(H2SO4)*Log(RH)**2)/h3-29.071647637053548*h6*q*Log(H2SO4)*
     & Log(RH)**2 +  (8.887445070958586*h6*q*Log(H2SO4)*Log(RH)**2)/h3-
     & 0.00029636660559253115*q**2 * Log(H2SO4)*Log(RH)**2 + 
     & (0.00046301434118339234*q**2*Log(H2SO4)*Log(RH)**2)/h3 +  
     & 0.005479820108166143*h6*q**2*Log(H2SO4)*Log(RH)**2 +  
     & 3.2937887074139884*T*Log(H2SO4)*Log(RH)**2 -  
     & (0.053358185813312024*T*Log(H2SO4)*Log(RH)**2)/h3 + 
     & 29.117876039650476*h6*T*Log(H2SO4)*Log(RH)**2 -  
     & (34.13451610985067*h6*T*Log(H2SO4)*Log(RH)**2)/h3 -  
     & 0.001150661716592889*q*T*Log(H2SO4)*Log(RH)**2 +  
     & (0.0005497364286653244*q*T* Log(H2SO4)* 
     & Log(RH)**2)/h3+0.09171149584822953*h6*q*T*Log(H2SO4)*Log(RH)**2-  
     & 2.9335804996871176d-6*q**2*T*Log(H2SO4)*Log(RH)**2 -  
     & 0.012024558565450824*T**2*Log(H2SO4)*Log(RH)**2 - 
     & (0.0008282043267863183*T**2*Log(H2SO4)*Log(RH)**2)/h3 + 
     & 0.11350434127907331*h6* T**2* 
     &  Log(H2SO4)*Log(RH)**2 - 3.206307709389464d-6*q* T**2*Log(H2SO4)* 
     & Log(RH)**2 + 0.000012855779740625904*T**3*Log(H2SO4)*Log(RH)**2 +  
     & 3.7577104325072703*Log(H2SO4)**2*Log(RH)**2 -(2.6641547962978636* 
     & Log(H2SO4)**2*Log(RH)**2)/h3 +  476.38347866428205*h6* 
     & Log(H2SO4)**2 *Log(RH)**2 - (127.7453688897455*h6*Log(H2SO4)**2*
     & Log(RH)**2)/h3 -0.01970020856731193*q*Log(H2SO4)**2*Log(RH)**2 + 
     & (0.006983199758249699*q*Log(H2SO4)**2*Log(RH)**2)/h3 +
     &  0.3243000572140283*h6*q*Log(H2SO4)**2* 
     & Log(RH)**2+0.00002903250111944733*q**2*Log(H2SO4)**2*Log(RH)**2-  
     & 0.037500815576593044*T*Log(H2SO4)**2*Log(RH)**2 + 
     & (0.013073560271857765*T*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     & 2.2753830913797373*h6*T* Log(H2SO4)**2* 
     & Log(RH)**2 + 0.00007918460467376976*q*T*Log(H2SO4)**2*Log(RH)**2+  
     & 0.00009291148493939081*T**2*Log(H2SO4)**2*Log(RH)**2

      h2=exp(h2)



      h4=-233.3693139924163 + 3711.127600293859*h6 - 
     &  127375.45943800849*h6**2 - 0.6541599370168311*q - 
     &  8.950348936875036*h6*q + 1420.4060399615116*h6**2*q + 
     &  0.006010885721884837*q**2 - 0.2514391282801529*h6*q**2 + 
     &  11.74107168004114*h6**2*q**2 - 27.242866772851034*RH + 
     &  3230.6550683739456*h6*RH - 7739.349030052802*h6**2*RH + 
     &  0.02657310586465451*q*RH + 0.8072083676135904*h6*q*RH + 
     &  21.19451916114249*h6**2*q*RH - 0.0013789987709190107*q**2*RH - 
     &  0.02985872690339605*h6*q**2*RH + 10.1858919054768*RH**2 + 
     &  1831.5638235525591*h6*RH**2 - 3345.6757256829833*h6**2*RH**2 - 
     &  0.006667965604100408*q*RH**2 - 1.6996949068091805*h6*q*RH**2 - 
     &  0.0002938711589315329*q**2*RH**2 + 20.8476664473087*RH**3 + 
     &  153.46841353011587*h6*RH**3 - 0.023573221099724897*q*RH**3 + 
     &  3.1234153503175706*T - 106.11845318552218*h6*T + 
     &  1616.0925045006472*h6**2*T - 0.000870039433136904*q*T + 
     &  0.07485289069770486*h6*q*T - 2.8492194924162093*h6**2*q*T - 
     &  0.000025041587132070856*q**2*T -
     &  0.00009872033200025474*h6*q**2*T + 0.1139999815738393*RH*T - 
     &  59.3235586131405*h6*RH*T + 188.66428911211182*h6**2*RH*T + 
     &  0.0021232663933975797*q*RH*T + 0.03525906170060478*h6*q*RH*T + 
     &  5.50009952133025d-6*q**2*RH*T - 0.40949236112495135*RH**2*T - 
     &  9.19061496536003*h6*RH**2*T + 
     &  0.00018605332374157369*q*RH**2*T - 
     &  0.15228517482883766*RH**3*T - 0.026779347826859368*T**2 + 
     &  0.7791278560987452*h6*T**2 - 5.635512845808971*h6**2*T**2 - 
     &  0.00004290197829367679*q*T**2 - 
     &  0.00034888935617559023*h6*q*T**2 + 
     &  3.399576950217114d-8*q**2*T**2 + 
     &  0.0047495543899101315*RH*T**2 + 
     &  0.21966980586301704*h6*RH*T**2 - 
     &  5.402900826149558d-6*q*RH*T**2 + 
     &  0.00209481671640675*RH**2*T**2 + 0.00004712913787763615*T**3- 
     &  0.0022395082766246705*h6*T**3 + 7.857921629051029d-8*q*T**3 - 
     &  0.000017637397545586075*RH*T**3 + 
     &  14.487613189612244*Log(H2SO4) + 
     &  834.2954132542228*h6*Log(H2SO4) - 
     &  9818.723117205202*h6**2*Log(H2SO4) + 
     &  0.1293876064755356*q*Log(H2SO4) + 
     &  0.5266286653232314*h6*q*Log(H2SO4) - 
     &  85.4813173795405*h6**2*q*Log(H2SO4) - 
     &  0.0005441269120155648*q**2*Log(H2SO4) + 
     &  0.013203398533271429*h6*q**2*Log(H2SO4) + 
     &  2.9312753776692975*RH*Log(H2SO4) + 
     &  329.1801633051309*h6*RH*Log(H2SO4) - 
     &  2384.7570613779067*h6**2*RH*Log(H2SO4) - 
     &  0.02754902616213419*q*RH*Log(H2SO4) - 
     &  0.39067534860398334*h6*q*RH*Log(H2SO4) + 
     &  0.000059670451934966084*q**2*RH*Log(H2SO4) + 
     &  2.1796062691153963*RH**2*Log(H2SO4) - 
     &  4.33856569910847*h6*RH**2*Log(H2SO4) + 
     &  0.0016947748838741473*q*RH**2*Log(H2SO4) + 
     &  0.7647868504545748*RH**3*Log(H2SO4) + 
     &  0.22920773390043764*T*Log(H2SO4) - 
     &  8.113046295543565*h6*T*Log(H2SO4) + 
     &  68.80414049795816*h6**2*T*Log(H2SO4) + 
     &  0.0014123733231002104*q*T*Log(H2SO4) + 
     &  0.0062552686253774335*h6*q*T*Log(H2SO4) + 
     &  3.138422092020368d-8*q**2*T*Log(H2SO4) - 
     &  0.11854147612732342*RH*T*Log(H2SO4) - 
     &  2.0186043272890317*h6*RH*T*Log(H2SO4) - 
     &  0.000017510448352402354*q*RH*T*Log(H2SO4) - 
     &  0.014662481263341446*RH**2*T*Log(H2SO4) + 
     &  0.001037948515905241*T**2*Log(H2SO4) + 
     &  0.04470222473271311*h6*T**2*Log(H2SO4) - 
     &  7.507637022804968d-7*q*T**2*Log(H2SO4) + 
     &  0.0002763185325508037*RH*T**2*Log(H2SO4) - 
     &  2.046693588795596d-6*T**3*Log(H2SO4) - 
     &  3.3271552306181786*Log(H2SO4)**2 - 
     &  0.6732293679771482*h6*Log(H2SO4)**2 - 
     &  123.7357115914075*h6**2*Log(H2SO4)**2 - 
     &  0.017682151597098173*q*Log(H2SO4)**2 - 
     &  0.06932927383915821*h6*q*Log(H2SO4)**2 + 
     &  0.000028013174553488913*q**2*Log(H2SO4)**2 + 
     &  0.5385357123752617*RH*Log(H2SO4)**2 + 
     &  4.405735442137794*h6*RH*Log(H2SO4)**2 + 
     &  0.0007699683891127761*q*RH*Log(H2SO4)**2 - 
     &  0.01981622019488726*RH**2*Log(H2SO4)**2 - 
     &  0.02863075409562691*T*Log(H2SO4)**2 - 
     &  0.3815368435566718*h6*T*Log(H2SO4)**2 - 
     &  0.000032173458991045656*q*T*Log(H2SO4)**2 + 
     &  0.0006610846557652917*RH*T*Log(H2SO4)**2 + 
     &  5.708155280776268d-6*T**2*Log(H2SO4)**2 + 
     &  0.2841922065650612*Log(H2SO4)**3 + 
     &  1.937916853785423*h6*Log(H2SO4)**3 + 
     &  0.00046645968682534585*q*Log(H2SO4)**3 - 
     &  0.01392240494812053*RH*Log(H2SO4)**3 + 
     &  0.0005604898286672238*T*Log(H2SO4)**3 - 
     &  0.00648946009121241*Log(H2SO4)**4

      H4=EXP(H4)


      h5=68.64045827314231-3277.3575769882523*h6 + 1.0798559249565618*q- 
     &  25.296110707348316*h6*q + 13.398992645698215*RH + 
     &  922.4932305036297*h6*RH - 0.27140107873619296*q*RH + 
     &  20.08312325165439*h6*q*RH + 66.82077511984484*RH**2 + 
     &  1611.1977384351555*h6*RH**2 - 0.02661518788217287*q*RH**2 + 
     &  3.0843227537972138*h6*q*RH**2 - 1.4080258142983926*T - 
     &  1.8568570408634648*h6*T - 0.0037450866352058397*q*T - 
     &  0.2576980690602505*h6*q*T - 1.3810906781490837*RH*T - 
     &  17.730890154257356*h6*RH*T + 0.001076745543266636*q*RH*T - 
     &  0.08909158555166723*h6*q*RH*T - 0.87735596568044*RH**2*T - 
     &  5.603639080904061*h6*RH**2*T + 0.0002224812194282388*q*RH**2*T + 
     &  0.015404311250124585*T**2 + 0.009972776386592903*h6*T**2 - 
     &  7.895537111616416d-7*q*T**2 - 0.0003918985727565079*h6*q*T**2 + 
     &  0.014090613376170114*RH*T**2 + 0.029815313352513736*h6*RH*T**2 - 
     &  3.607688962369381d-6*q*RH*T**2+0.0028248811086430334*RH**2*T**2- 
     &  0.00010105566309790024*T**3 - 0.00008368737713099654*h6*T**3 + 
     &  1.0956121068002908d-8*q*T**3 - 0.000042201374093774165*RH*T**3+ 
     &  2.429718549962505d-7*T**4 - 0.9150072698367441*Log(H2SO4) + 
     &  653.1422380210507*h6*Log(H2SO4) - 0.16169777329907886*q*
     &  Log(H2SO4)+6.860485142132651*h6*q*Log(H2SO4)+13.912279451959678*
     &  RH*Log(H2SO4) - 12.99843208521885*h6*RH*Log(H2SO4) + 
     &  0.021566098556935497*q*RH*Log(H2SO4) - 0.17455209117051437*h6*q*
     &  RH*Log(H2SO4) + 3.4003096116229066*RH**2*Log(H2SO4) - 
     &  36.48723431979357*h6*RH**2*Log(H2SO4) - 0.001787998461490877*q*
     &  RH**2*Log(H2SO4) - 0.11123959974265946*T*Log(H2SO4) +
     &  1.4641341095178422*h6*T*Log(H2SO4) + 0.0005023439404053182*q*T*
     &  Log(H2SO4) + 0.031171816066622334*h6*q*T*Log(H2SO4) - 
     &  0.1666443655044969*RH*T*Log(H2SO4) + 0.9274525071251363*h6*RH*T*
     &  Log(H2SO4) + 0.00003323858217484082*q*RH*T*Log(H2SO4) - 
     &  0.021039758605777333*RH**2*T*Log(H2SO4) + 0.0018827939246916288*
     &  T**2*Log(H2SO4) + 0.0014616406839659635*h6*T**2*Log(H2SO4) - 
     &  2.8887052861955783d-7*q*T**2*Log(H2SO4) + 0.0007036899839828474*
     &  RH*T**2*Log(H2SO4) - 6.077113667352327d-6*T**3*Log(H2SO4) + 
     &  0.8758449734328851*Log(H2SO4)**2 - 61.130559299849466*h6*
     &  Log(H2SO4)**2 + 0.007311434341744119*q*Log(H2SO4)**2 - 
     &  0.4613762652439876*h6*q*Log(H2SO4)**2 + 0.11048996127121676*RH*
     &  Log(H2SO4)**2 - 5.241218416739315*h6*RH*Log(H2SO4)**2 - 
     &  0.0009849834685920686*q*RH*Log(H2SO4)**2 + 0.034417912409269905*
     &  RH**2*Log(H2SO4)**2 - 0.01737690919778605*T*Log(H2SO4)**2 - 
     &  0.10868884345327394*h6*T*Log(H2SO4)**2-0.000013815430650770167*
     &  q*T*Log(H2SO4)**2 - 0.0038903960432723123*RH*T*Log(H2SO4)**2 + 
     &  0.00005879779401379263*T**2*Log(H2SO4)**2 + 0.03393085910585771*
     &  Log(H2SO4)**3 + 2.2221405215918755*h6*Log(H2SO4)**3 - 
     &  0.00009520679794957912*q*Log(H2SO4)**3 + 0.0166181059051036*RH*
     &  Log(H2SO4)**3 - 0.00014027668903237805*T*Log(H2SO4)**3

      if (h1 .gt. 1.d5) then ! take care of weird nuc rate blow up
         h1=0.d0
         h2=0.d0
         h3=0.d0
         h4=0.d0
         h5=0.d0
         h6=0.d0
         return
      endif

      return
      End 

C=======================================================================
C
C *** SUBROUTINE getCCN_kappa
C *** WRITTEN BY Yunha Lee
C *** Compute CCN at 0.1, 0.2, 0.3% 
C
C=======================================================================
C
      SUBROUTINE getCCN_kappa(I,J,L)
C
      USE TOMAS_AEROSOL  
      USE TRACER_COM, only : ntm,trm,tr_mm
     &     ,nbins,xk,trpdens,n_AECIL,
     &       n_AOCIL,n_AOCOB,n_ASO4,n_ANACL,n_ADUST,
     &       n_AECOB
      USE TRDIAG_COM, only: taijls=>taijls_loc,ijlt_ccn_01
     &     ,ijlt_ccn_03,ijlt_ccn_02
      USE CONSTANT, only: pi,gasc
      implicit none 
      REAL*8 SURT,DIAM3(NBINS+1),Tvol,DENS(7),A3
c      REAL*8 Tp,BOXM,BOXV
      integer i, j, l, si, n
      integer k,kk,tracnum

      
!@var constants needed for CCN calculation 
      real*8, parameter :: Mv=18.015d-3
      real*8, parameter :: rhow= 1000.d0
!@var temporal CCN 
      real*8, dimension(nsmax) :: ccn_mod 
!@var temporal Sc 
      real*8 Scnew
!@var Sc at each size boundary
      real*8, dimension (nbins) :: Sc,kappa

!@var CCN at 0.1/0.2/0.3%    (#/m3)
c      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: CCN_TOMAS


C initialize CCN_mod
      CCN_mod(:)=0.
C get density 

      dens(1)=trpdens(n_ASO4(1))
      dens(2)=trpdens(n_ANACL(1))
      dens(3)=trpdens(n_AECOB(1))
      dens(4)=trpdens(n_AECIL(1))
      dens(5)=trpdens(n_AOCOB(1))
      dens(6)=trpdens(n_AOCIL(1))
      dens(7)=trpdens(n_ADUST(1))

C surface tension
      SURT   = 0.0761-1.55E-4*(Temp-273.)

      A3 = (4*Mv*SURT/(gasc*Temp*rhow))**3

      DO N=1,NBINS+1 
C Diameter cubed in each size boundary with assuming density =1800 kg/m3. 
        diam3(n)= xk(n)/1800.d0*6.d0/pi
      ENDDO

      DO N=1,NBINS

        Tvol =Mk(N,1)/dens(1)+Mk(N,2)/dens(2)+Mk(N,3)/dens(3)+
     &       Mk(N,4)/dens(4)+Mk(N,5)/dens(5)+Mk(N,6)/dens(6)
     &       +Mk(N,7)/dens(7)   ! total vol. of species 

        kappa(n)=(0.6*Mk(n,1)/dens(1)+1.28*Mk(n,2)/dens(2)
     &       +0.227*Mk(n,6)/dens(6))/Tvol ! average kappa in a bin 
C note that kappa is hard-coded here. 

        Sc(n) = sqrt(4.d0*A3/27.d0/Diam3(n)/kappa(n))
        Sc(n) = min(Sc(n), 10.d0) ! HACK!!! Yunha must fix this. Chances are that particles of diameter 2.55e-9 in mode 1 are just too small for this calculation?
        Sc(n) = exp(Sc(n))
        Sc(n)=(Sc(n)-1.d0)*100.d0

c        print*,'debug_kappa',n,kappa(n),Sc(n)

        if(Sc(n) .lt. 0.) Sc(n)=1.e-6
      ENDDO

      DO N=1,NBINS

C compute CCN at various Smax
 
        DO SI=1,nsmax 
          IF(SC(N) .LE. SMAX(SI)) CCN_mod(SI)=CCN_mod(SI)
     *         +Nk(n)/boxvol ! unit is cm-3 now 
          
C     INTERPOLATION :
          if(N .LT. NBINS)THEN 
            IF(SC(N+1) .lt. SMAX(si) .and. SC(N) .gt. SMAX(SI) ) THEN 
C     compute new Sc (I+1) using the upper limit Dp to determine the activation fraction  
              Scnew = sqrt(4.d0*A3/27.d0/Diam3(n+1)/kappa(n))
              Scnew = min(Scnew, 10.d0) ! HACK!!! Yunha must fix this. Chances are that particles of diameter 2.55e-9 in mode 1 are just too small for this calculation?
              Scnew = exp(Scnew)
              Scnew=(Scnew-1.d0)*100.d0

              if (Sc(n) .ne. Scnew) ! HACK! This is needed to avoid division by zero. Yunha must verify that this is indeed the correct behavior, as it appears to be the case.
     &        CCN_mod(SI)=CCN_mod(SI)+Nk(n)/boxvol*
     &         (1/(dlog(100.+SMAX(SI)))**(2)-1/(dlog(100.+Scnew))**(2))/
     &         (1/(dlog(100.+Sc(n)))**(2)-1/(dlog(100.+Scnew))**(2))
              
            ENDIF    
          ENDIF
          
        ENDDO ! SMAX

      ENDDO ! size bin
C        PRINT*,'interpolate',i,j,l,si,CCN_mod(si) 
C assing each CCN_mod to 3d array 

      CCN_TOMAS(I,J,L,:)=CCN_mod(:) ! ; unit is not cm-3 yet 

        taijls(i,j,l,ijlt_ccn_01)=taijls(i,j,l,ijlt_ccn_01)+
     &   CCN_mod(1)

       taijls(i,j,l,ijlt_ccn_02)=taijls(i,j,l,ijlt_ccn_02)+
     &   CCN_mod(2)

       taijls(i,j,l,ijlt_ccn_03)=taijls(i,j,l,ijlt_ccn_03)+
     &   CCN_mod(3)



      RETURN
      END
