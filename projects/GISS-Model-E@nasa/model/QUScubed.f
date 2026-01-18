!@sum  Implementation of the GISS Quadratic Upstream Scheme on a
!@+    domain-decomposed cubed sphere grid.
!@+    Main routines (based on modelE conventions):
!@+      QDYNAM: calls AADVQ0 and then calls AADVQ to advect humidity q
!@+      AADVQ0: determines adaptive timestepping parameters
!@+      AADVQ:  driver routine for tracer advection
!@+    For now, AADVQ assumes transported quantities are non-negative.
!@+    Code for the easier case without this constraint will be added.
!@vers 2013/04/02
!@auth M. Kelley

#include "rundeck_opts.h"

      module tracer_adv
      use resolution, only : lm
      use qusdef, only : nmom, mx,my,mz, mxx,myy,mzz, mxy,myz,mzx
      implicit none
      save

      integer :: gid, tile, nxg,nyg,
     &     is,ie,js,je, isd,ied,jsd,jed,
     &     ncyc

      integer, parameter :: ncmax=10
      integer, dimension(lm) :: ncycxy
      real*8, parameter :: mrat_limh=0.25

      real*8 :: xtoy_w,xtoy_e,xtoy_s,xtoy_n

      real*8, dimension(:,:,:), allocatable :: mu,mv,mw

      contains

      subroutine do_edges_and_corners(mma,rm,rmom,mu,mv)
      implicit none
      real*8, dimension(isd:ied,jsd:jed) :: mma,rm
      real*8, dimension(isd:ied+1,jsd:jed) :: mu
      real*8, dimension(isd:ied,jsd:jed+1) :: mv
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j,iup,iam,jup,jbm
      real*8 :: xsigni,xsigno,ysigni,ysigno,am,bm,mnewx,mnewy

      call rotate_edges(rmom)

c
c east-west edges (excluding corners)
c
      do i=1,nxg,nxg-1
        if(i.lt.is .or. i.gt.ie) cycle
        if(i.eq.1) then
          iup=is-1
          iam=is
          xsigno = -1.
        else
          iup=ie+1
          iam=ie+1
          xsigno = +1.
        endif
        xsigni = -xsigno
        do j=max(2,jsd),min(jed,nyg-1)
          am=mu(iam,j)*xsigno
          if(am.ge.0.) then
            call qusout(am,mma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(mxx,i,j),xsigno)
            call lusout(am,mma(i,j),rmom(my,i,j),rmom(mxy,i,j),xsigno)
            call lusout(am,mma(i,j),rmom(mz,i,j),rmom(mzx,i,j),xsigno)
            rmom((/myy,myz,mzz/),i,j) =
     &           rmom((/myy,myz,mzz/),i,j)*(1.-am/mma(i,j))
          else
            am = -am
            call qusin(am,mma(i,j),rm(i,j),rmom(mx,i,j),rmom(mxx,i,j),
     &       mma(iup,j),rm(iup,j),rmom(mx,iup,j),rmom(mxx,iup,j),xsigni)
            call lusin(am,mma(i,j),rmom(my,i,j),rmom(mxy,i,j),
     &           mma(iup,j),rmom(my,iup,j),rmom(mxy,iup,j),xsigni)
            call lusin(am,mma(i,j),rmom(mz,i,j),rmom(mzx,i,j),
     &           mma(iup,j),rmom(mz,iup,j),rmom(mzx,iup,j),xsigni)
            call susin(am,rmom(myy,iup,j),mma(iup,j),rmom(myy,iup,j))
            call susin(am,rmom(myz,iup,j),mma(iup,j),rmom(myz,iup,j))
            call susin(am,rmom(mzz,iup,j),mma(iup,j),rmom(mzz,iup,j))
          endif
          mma(i,j) = mma(i,j) + mu(iam,j)*xsigni
        enddo
      enddo

c
c north-south edges (excluding corners)
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jup=js-1
          jbm=js
          ysigno = -1.
        else
          jup=je+1
          jbm=je+1
          ysigno = +1.
        endif
        ysigni = -ysigno
        do i=max(2,isd),min(ied,nxg-1)
          bm=mv(i,jbm)*ysigno
          if(bm.ge.0.) then
            call qusout(bm,mma(i,j),rm(i,j),
     &           rmom(my,i,j),rmom(myy,i,j),ysigno)
            call lusout(bm,mma(i,j),rmom(mx,i,j),rmom(mxy,i,j),ysigno)
            call lusout(bm,mma(i,j),rmom(mz,i,j),rmom(myz,i,j),ysigno)
            rmom((/mxx,mzx,mzz/),i,j) =
     &           rmom((/mxx,mzx,mzz/),i,j)*(1.-bm/mma(i,j))
          else
            bm = -bm
            call qusin(bm,mma(i,j),rm(i,j),rmom(my,i,j),rmom(myy,i,j),
     &       mma(i,jup),rm(i,jup),rmom(my,i,jup),rmom(myy,i,jup),ysigni)
            call lusin(bm,mma(i,j),rmom(mx,i,j),rmom(mxy,i,j),
     &           mma(i,jup),rmom(mx,i,jup),rmom(mxy,i,jup),ysigni)
            call lusin(bm,mma(i,j),rmom(mz,i,j),rmom(myz,i,j),
     &           mma(i,jup),rmom(mz,i,jup),rmom(myz,i,jup),ysigni)
            call susin(bm,rmom(mxx,i,jup),mma(i,jup),rmom(mxx,i,jup))
            call susin(bm,rmom(mzx,i,jup),mma(i,jup),rmom(mzx,i,jup))
            call susin(bm,rmom(mzz,i,jup),mma(i,jup),rmom(mzz,i,jup))
          endif
          mma(i,j) = mma(i,j) + mv(i,jbm)*ysigni
        enddo
      enddo

c
c corners:
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jup=js-1
          jbm=js
          ysigno = -1.
        else
          jup=je+1
          jbm=je+1
          ysigno = +1.
        endif
        ysigni = -ysigno
        do i=1,nxg,nxg-1
          if(i.lt.is .or. i.gt.ie) cycle
          if(i.eq.1) then
            iup=is-1
            iam=is
            xsigno = -1.
          else
            iup=ie+1
            iam=ie+1
            xsigno = +1.
          endif
          xsigni = -xsigno
          am=mu(iam,j)*xsigno
          bm=mv(i,jbm)*ysigno
          rmom((/mxx,myy,mxy/),i,j) = 0d0
          if(am.ge.0. .and. bm.ge.0.) then ! flow exiting both x and y
            call lusout_2sides(am,bm,mma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(my,i,j),xsigno,ysigno)
            call lusout_2sides(am,bm,mma(i,j),rmom(mz,i,j),
     &           rmom(mzx,i,j),rmom(myz,i,j),xsigno,ysigno)
            rmom(mzz,i,j) = rmom(mzz,i,j)*(1.-(am+bm)/mma(i,j))
          elseif(am.lt.0. .and. bm.lt.0.) then ! flow entering both x and y
            am = -am
            bm = -bm
            call lusin_2sides(am,bm
     &           ,mma(i,j),rm(i,j),rmom(mx,i,j),rmom(my,i,j)
     &           ,mma(iup,j),rm(iup,j),rmom(mx,iup,j),rmom(my,iup,j)
     &           ,mma(i,jup),rm(i,jup),rmom(mx,i,jup),rmom(my,i,jup)
     &           ,xsigni,ysigni)
            call lusin_2sides(am,bm
     &           ,mma(i,j),rmom(mz,i,j),rmom(mzx,i,j),rmom(myz,i,j)
     &        ,mma(iup,j),rmom(mz,iup,j),rmom(mzx,iup,j),rmom(myz,iup,j)
     &        ,mma(i,jup),rmom(mz,i,jup),rmom(mzx,i,jup),rmom(myz,i,jup)
     &           ,xsigni,ysigni)
            call susin_2sides(am,bm
     &           ,mma(i,j),rmom(mzz,i,j)
     &           ,mma(iup,j),rmom(mzz,iup,j)
     &           ,mma(i,jup),rmom(mzz,i,jup)
     &           )
          elseif(am.le.0.) then ! flow exits in y, then enters in x
            call lusout(bm,mma(i,j),rm(i,j),rmom(my,i,j),ysigno)
            call lusout(bm,mma(i,j),rmom(mz,i,j),rmom(myz,i,j),ysigno)
            rmom((/mx,mzx,mzz/),i,j) =
     &           rmom((/mx,mzx,mzz/),i,j)*(1.-bm/mma(i,j))
            mnewy = mma(i,j)-bm
            am = -am
            call lusin(am,mnewy,rm(i,j),rmom(mx,i,j),
     &           mma(iup,j),rm(iup,j),rmom(mx,iup,j),xsigni)
            call lusin(am,mnewy,rmom(mz,i,j),rmom(mzx,i,j),
     &           mma(iup,j),rmom(mz,iup,j),rmom(mzx,iup,j),xsigni)
            call susin(am,rmom(my,i,j),mma(iup,j),rmom(my,iup,j))
            call susin(am,rmom(myz,i,j),mma(iup,j),rmom(myz,iup,j))
            call susin(am,rmom(mzz,i,j),mma(iup,j),rmom(mzz,iup,j))
          elseif(bm.le.0.) then ! flow exits in x, then enters in y
            call lusout(am,mma(i,j),rm(i,j),rmom(mx,i,j),xsigno)
            call lusout(am,mma(i,j),rmom(mz,i,j),rmom(mzx,i,j),xsigno)
            rmom((/my,myz,mzz/),i,j) =
     &           rmom((/my,myz,mzz/),i,j)*(1.-am/mma(i,j))
            mnewx = mma(i,j)-am
            bm = -bm
            call lusin(bm,mnewx,rm(i,j),rmom(my,i,j),
     &           mma(i,jup),rm(i,jup),rmom(my,i,jup),ysigni)
            call lusin(bm,mnewx,rmom(mz,i,j),rmom(myz,i,j),
     &           mma(i,jup),rmom(mz,i,jup),rmom(myz,i,jup),ysigni)
            call susin(bm,rmom(mx,i,j),mma(i,jup),rmom(mx,i,jup))
            call susin(bm,rmom(mzx,i,j),mma(i,jup),rmom(mzx,i,jup))
            call susin(bm,rmom(mzz,i,j),mma(i,jup),rmom(mzz,i,jup))
          endif
          mma(i,j) = mma(i,j) + mu(iam,j)*xsigni + mv(i,jbm)*ysigni
        enddo
      enddo

      return
      end subroutine do_edges_and_corners

      subroutine check_edges_and_corners(mma,rm,rmom,mu,mv)
c HALO ROWS NOT DONE HERE
      implicit none
      real*8, dimension(isd:ied,jsd:jed) :: mma,rm
      real*8, dimension(isd:ied+1,jsd:jed) :: mu ! check these bounds
      real*8, dimension(isd:ied,jsd:jed+1) :: mv ! check these bounds
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j,iam,jbm
      real*8 :: xsigno,ysigno,am,bm
c
c east-west edges (excluding corners)
c
      do i=1,nxg,nxg-1
        if(i.lt.is .or. i.gt.ie) cycle
        if(i.eq.1) then
          iam=is
          xsigno = -1.
        else
          iam=ie+1
          xsigno = +1.
        endif
        do j=max(2,js),min(je,nyg-1)
          am=mu(iam,j)*xsigno
          if(am.gt.0.) then
            call check_qusout(am,mma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(mxx,i,j),xsigno)
          endif
        enddo
      enddo

c
c north-south edges (excluding corners)
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jbm=js
          ysigno = -1.
        else
          jbm=je+1
          ysigno = +1.
        endif
        do i=max(2,is),min(ie,nxg-1)
          bm=mv(i,jbm)*ysigno
          if(bm.gt.0.) then
            call check_qusout(bm,mma(i,j),rm(i,j),
     &           rmom(my,i,j),rmom(myy,i,j),ysigno)
          endif
        enddo
      enddo

c
c corners:
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jbm=js
          ysigno = -1.
        else
          jbm=je+1
          ysigno = +1.
        endif
        do i=1,nxg,nxg-1
          if(i.lt.is .or. i.gt.ie) cycle
          if(i.eq.1) then
            iam=is
            xsigno = -1.
          else
            iam=ie+1
            xsigno = +1.
          endif
          am=mu(iam,j)*xsigno
          bm=mv(i,jbm)*ysigno
          if(am.gt.0. .and. bm.gt.0.) then ! flow exiting both x and y
            call check_lusout_2sides(am,bm,mma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(my,i,j),xsigno,ysigno)
          elseif(am.lt.0. .and. bm.lt.0.) then ! flow entering both x and y
          elseif(am.le.0.) then ! flow exits in y, then enters in x
            call check_lusout(bm,mma(i,j),rm(i,j),rmom(my,i,j),ysigno)
          elseif(bm.le.0.) then ! flow exits in x, then enters in y
            call check_lusout(am,mma(i,j),rm(i,j),rmom(mx,i,j),xsigno)
          endif
        enddo
      enddo
      return
      end subroutine check_edges_and_corners

      subroutine rotate_edges(rmom)
c xtoy =  0 means no change of orientation
c xtoy = +1 or -1 means x' = xtoy*y, y' = -xtoy*x
c xtoy for a direction defines the transformation for data coming from
c       that direction, so this routine is called after a halo_update.
      implicit none
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j
      real*8 :: xtoy,tmp
c
c east-west edges (including corners)
c
      do i=0,nxg+1,nxg+1
        if(i.lt.isd .or. i.gt.ied) cycle
        if(i.eq.0) then
          xtoy = xtoy_w
        else
          xtoy = xtoy_e
        endif
        if(xtoy.eq.0) cycle
        do j=max(1,jsd),min(jed,nyg)
          call rotxy(rmom(mx,i,j),rmom(my,i,j),xtoy)
          call rotxy(rmom(mzx,i,j),rmom(myz,i,j),xtoy)
          rmom(mxy,i,j) = -rmom(mxy,i,j)
          call swapxy(rmom(mxx,i,j),rmom(myy,i,j))
        enddo
      enddo

c
c north-south edges (including corners)
c
      do j=0,nyg+1,nyg+1
        if(j.lt.jsd .or. j.gt.jed) cycle
        if(j.eq.0) then
          xtoy = xtoy_s
        else
          xtoy = xtoy_n
        endif
        if(xtoy.eq.0) cycle
        do i=max(1,isd),min(ied,nxg)
          call rotxy(rmom(mx,i,j),rmom(my,i,j),xtoy)
          call rotxy(rmom(mzx,i,j),rmom(myz,i,j),xtoy)
          rmom(mxy,i,j) = -rmom(mxy,i,j)
          call swapxy(rmom(mxx,i,j),rmom(myy,i,j))
        enddo
      enddo
      return
      end subroutine rotate_edges

      subroutine rotxy(rx,ry,xtoy)
      implicit none
      real*8 :: rx,ry,xtoy
      real*8 :: tmp
      tmp = ry
      ry = -xtoy*rx
      rx = +xtoy*tmp
      return
      end subroutine rotxy

      subroutine swapxy(rx,ry)
      implicit none
      real*8 :: rx,ry
      real*8 :: tmp
      tmp = ry
      ry = rx
      rx = tmp
      return
      end subroutine swapxy

      subroutine checkfobs_x(mma,mu,rm,rmom)
      implicit none
      real*8, dimension(isd:ied,jsd:jed) :: mma,rm
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      real*8, dimension(isd:ied+1,jsd:jed) :: mu
      integer :: i,j
      logical :: changed_mom
      do j=max(1,jsd),min(jed,nyg)
      do i=max(1,isd),min(ied,nxg)
        if(mu(i+1,j).gt.0. .and. mu(i,j).lt.0.) then
          call checkfobs(mu(i,j),mu(i+1,j),mma(i,j),
     &         rm(i,j),rmom(mx,i,j),rmom(mxx,i,j),
     &         changed_mom)
        endif
      enddo
      enddo
      return
      end subroutine checkfobs_x

      subroutine checkfobs_y(mma,mv,rm,rmom)
      implicit none
      real*8, dimension(isd:ied,jsd:jed) :: mma,rm
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      real*8, dimension(isd:ied,jsd:jed+1) :: mv
      integer :: i,j
      logical :: changed_mom
      do j=max(1,jsd),min(jed,nyg)
      do i=is,ie
        if(mv(i,j+1).gt.0. .and. mv(i,j).lt.0.) then
          call checkfobs(mv(i,j),mv(i,j+1),mma(i,j),
     &         rm(i,j),rmom(my,i,j),rmom(myy,i,j),
     &         changed_mom)
        endif
      enddo
      enddo
      return
      end subroutine checkfobs_y

      subroutine checkfobs_z(mma,mw,rm,rmom)
      implicit none
      real*8, dimension(isd:ied,jsd:jed) :: mma,rm
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      real*8, dimension(isd:ied,jsd:jed,2) :: mw
      integer :: i,j
      logical :: changed_mom
      do j=js,je
      do i=is,ie
        if(mw(i,j,2).gt.0. .and. mw(i,j,1).lt.0.) then
          call checkfobs(mw(i,j,1),mw(i,j,2),mma(i,j),
     &         rm(i,j),rmom(mz,i,j),rmom(mzz,i,j),
     &         changed_mom)
        endif
      enddo
      enddo
      return
      end subroutine checkfobs_z

      subroutine checkfobs(aml,amr,m,rm,rxm,rxxm,changed)
      implicit none
      real*8 :: aml,amr,m,rm,rxm,rxxm
      logical :: changed
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
        changed=.true.
      else
        changed=.false.
      endif
      return
      end subroutine checkfobs

      subroutine aadvqx(rm,rmom,mass,mu)
!@sum  AADVQX advection driver for x-direction
!@auth Maxwell Kelley
      implicit none
      REAL*8, dimension(isd:ied,jsd:jed) :: rm,mass
      REAL*8, dimension(isd:ied+1,jsd:jed) :: mu
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,ii,j,ns
      real*8 :: am,frac1,fracm,fw,fe,dm2,mold,mnew,bymnew
      real*8 :: fe0,fe_pass,fw0,rm0,fex_pass,fexx_pass
      real*8, dimension(nmom) :: fmomw,fmome

      call checkfobs_x(mass,mu,rm,rmom)

      do j=max(1,jsd),min(jed,nyg) ! calculate halo rows as well

c-----------------------------------------------------------
      ! calculate tracer mass flux f
c-----------------------------------------------------------
c--------------------------------------------------------------------
      ! calculate tracer fluxes of slopes and curvatures
c--------------------------------------------------------------------
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      i=is-1
      am = mu(i+1,j)
      if(am.le.0.) then      ! air mass flux is negative
        ii=i+1
        frac1=+1.
      else                      ! air mass flux is positive
        ii=i
        frac1=-1.
      endif
      fracm=am/mass(ii,j)
      frac1=fracm+frac1
      fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &     (frac1+fracm)*rmom(mxx,ii,j)))
      fmome(mx)=am*(fracm*fracm*(rmom(mx,ii,j)
     &     -3.*frac1*rmom(mxx,ii,j))-3.*fe)
      fmome(mxx)=am*(am*fracm**3 *rmom(mxx,ii,j)
     &     -5.*(am*fe+fmome(mx)))

      ! cross moments
      fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
      fmome(mxy) = am*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
      fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
      fmome(mzx) = am*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
      fmome(myy) = fracm*rmom(myy,ii,j)
      fmome(mzz) = fracm*rmom(mzz,ii,j)
      fmome(myz) = fracm*rmom(myz,ii,j)

! flux limitations
      fe0 = fe
      fe_pass = fe
      fex_pass = fmome(mx)
      fexx_pass = fmome(mxx)
      if(am.gt.0.) then
        if(fe.lt.0.) then
          fe=0.
          fe_pass=0.
          fex_pass=0
          fexx_pass=0.
        elseif(fe.gt.rm(i,j)) then
          fe=rm(i,j)
          fe_pass = fe
          fex_pass=am*(-3.*fe)
          fexx_pass=am*(-5.*(am*fe+fex_pass))
        endif
      else if(am.lt.0.) then
        if(fe.gt.0.) then
          fe=0.
          fe0=0.
          fmome((/mx,mxx/))=0.
        elseif(fe.lt.-rm(i+1,j)) then
          fe=-rm(i+1,j)
          fe0 = fe
          fmome(mx)=am*(-3.*fe)
          fmome(mxx)=am*(-5.*(am*fe+fmome(mx)))
        endif
      endif
      fw = fe
      fw0 = fe_pass
      fmomw(:) = fmome(:)
      fmomw(mx) = fex_pass
      fmomw(mxx) = fexx_pass

      do i=is,ie
         am = mu(i+1,j)
         if(am.lt.0.) then ! air mass flux is negative
            ii=i+1
            frac1=+1.
         else                 ! air mass flux is positive
            ii=i
            frac1=-1.
         endif
         fracm=am/mass(ii,j)
         frac1=fracm+frac1
         fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &        (frac1+fracm)*rmom(mxx,ii,j)))
         fmome(mx)=am*(fracm*fracm*(rmom(mx,ii,j)
     &        -3.*frac1*rmom(mxx,ii,j))-3.*fe)
         fmome(mxx)=am*(am*fracm**3 *rmom(mxx,ii,j)
     &        -5.*(am*fe+fmome(mx)))
      ! cross moments
         fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
         fmome(mxy) =am*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
         fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
         fmome(mzx) =am*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
         fmome(myy) = fracm*rmom(myy,ii,j)
         fmome(mzz) = fracm*rmom(mzz,ii,j)
         fmome(myz) = fracm*rmom(myz,ii,j)

! flux limitations
         fe0 = fe
         fe_pass = fe
         fex_pass = fmome(mx)
         fexx_pass = fmome(mxx)
         if(am.gt.0.) then
           if(fe.lt.0.) then
             fe=0.
             fe_pass=0.
             fex_pass=0
             fexx_pass=0.
           elseif(fe.gt.rm(i,j)) then
             fe=rm(i,j)
             fe_pass = fe
             fex_pass=am*(-3.*fe)
             fexx_pass=am*(-5.*(am*fe+fex_pass))
           endif
         elseif(am.lt.0.) then
           if(fe.gt.0.) then
             fe=0.
             fe0=0.
             fmome((/mx,mxx/))=0.
           elseif(fe.lt.-rm(i+1,j)) then
             fe=-rm(i+1,j)
             fe0 = fe
             fmome(mx)=am*(-3.*fe)
             fmome(mxx)=am*(-5.*(am*fe+fmome(mx)))
           endif
         endif

         mold=mass(i,j)
         mnew=mold+mu(i,j)-am
         bymnew = 1./mnew
         dm2=mu(i,j)+am
         rm0=rm(i,j)+fw0-fe0
         rm(i,j)=rm(i,j)+fw-fe
      !
         rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm0
     &     +mold*(fw0+fe0))+(fmomw(mx)-fmome(mx)))*bymnew
         rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw0-fe0)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &        mold*(fmomw(my)+fmome(my))) +
     &        (fmomw(mxy)-fmome(mxy)))*bymnew
         rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
         rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmomw(mz)+fmome(mz))) +
     &        (fmomw(mzx)-fmome(mzx)))*bymnew
      !
         rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
         rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

         mass(i,j) = mnew

! clean up roundoff errors
         if(rm(i,j).le.0d0) then
           rm(i,j)=0d0; rmom(:,i,j)=0d0
         endif

         fw = fe
         fmomw(:) = fmome(:)
         fw0 = fe_pass
         fmomw(mx) = fex_pass
         fmomw(mxx) = fexx_pass

      enddo ! i
      enddo ! j

      return
c****
      end subroutine aadvqx

      subroutine aadvqy(rm,rmom,mass,mv)
!@sum  AADVQY advection driver for y-direction
!@auth Maxwell Kelley
      implicit none
      REAL*8, dimension(isd:ied,jsd:jed) :: rm,mass
      REAL*8, dimension(isd:ied,jsd:jed+1) :: mv
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j,jj
      real*8, dimension(isd:ied) :: mvj,fs
      real*8, dimension(nmom,isd:ied) :: fmoms
      real*8, dimension(nmom) :: fmomn
      real*8 :: frac1,fracm,fn,mold,mnew,bymnew,dm2,am
      real*8 :: fn0,fn_pass,rm0,fny_pass,fnyy_pass
      real*8, dimension(isd:ied) :: fs0

      call checkfobs_y(mass,mv,rm,rmom)

c-----------------------------------------------------------
      ! calculate tracer mass flux f
      ! and fluxes of slopes and curvatures fmom
      ! update tracer mass, moments of tracer mass, air mass distribution
c--------------------------------------------------------------------
      j=js-1
      do i=is,ie
        am = mv(i,j+1)
        if(am.le.0.) then ! air mass flux is negative
          jj=j+1
          frac1=+1.
        else                    ! air mass flux is positive
          jj=j
          frac1=-1.
        endif
        fracm=am/mass(i,jj)
        frac1=fracm+frac1
        fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &       (frac1+fracm)*rmom(myy,i,jj)))
        fmomn(my)=am*(fracm*fracm*(rmom(my,i,jj)
     &       -3.*frac1*rmom(myy,i,jj))-3.*fn)
        fmomn(myy)=am*(am*fracm**3 *rmom(myy,i,jj)
     &       -5.*(am*fn+fmomn(my)))
        fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
        fmomn(myz) = am*
     &       (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
        fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
        fmomn(mxy) = am*
     &       (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
        fmomn(mzz) = fracm*rmom(mzz,i,jj)
        fmomn(mxx) = fracm*rmom(mxx,i,jj)
        fmomn(mzx) = fracm*rmom(mzx,i,jj)

! flux limitations
        fn0 = fn
        fn_pass = fn
        fny_pass = fmomn(my)
        fnyy_pass = fmomn(myy)
        if(am.gt.0.) then
          if(fn.lt.0.) then
            fn=0.
            fn_pass=0.
            fny_pass=0
            fnyy_pass=0.
          elseif(fn.gt.rm(i,j)) then
            fn=rm(i,j)
            fn_pass = fn
            fny_pass=am*(-3.*fn)
            fnyy_pass=am*(-5.*(am*fn+fny_pass))
          endif
        elseif(am.lt.0.) then
          if(fn.gt.0.) then
            fn=0.
            fn0=0.
            fmomn((/my,myy/))=0.
          elseif(fn.lt.-rm(i,j+1)) then
            fn=-rm(i,j+1)
            fn0 = fn
            fmomn(my)=am*(-3.*fn)
            fmomn(myy)=am*(-5.*(am*fn+fmomn(my)))
          endif
        endif

        mvj(i) = am
        fs(i) = fn
        fmoms(:,i) = fmomn(:)

        fs0(i) = fn_pass
        fmoms(my,i) = fny_pass
        fmoms(myy,i) = fnyy_pass

      enddo                     ! i

      do j=js,je
      do i=is,ie
         am = mv(i,j+1)
         if(am.lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                   ! air mass flux is positive
            jj=j
            frac1=-1.
         endif
         fracm=am/mass(i,jj)
         frac1=fracm+frac1
         fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &        (frac1+fracm)*rmom(myy,i,jj)))
         fmomn(my)=am*(fracm*fracm*(rmom(my,i,jj)
     &        -3.*frac1*rmom(myy,i,jj))-3.*fn)
         fmomn(myy)=am*(am*fracm**3 *rmom(myy,i,jj)
     &        -5.*(am*fn+fmomn(my)))
         fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
         fmomn(myz) = am*
     &        (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
         fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
         fmomn(mxy) = am*
     &        (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
         fmomn(mzz) = fracm*rmom(mzz,i,jj)
         fmomn(mxx) = fracm*rmom(mxx,i,jj)
         fmomn(mzx) = fracm*rmom(mzx,i,jj)

! flux limitations
         fn0 = fn
         fn_pass = fn
         fny_pass = fmomn(my)
         fnyy_pass = fmomn(myy)
         if(am.gt.0.) then
           if(fn.lt.0.) then
             fn=0.
             fn_pass=0.
             fny_pass=0
             fnyy_pass=0.
           elseif(fn.gt.rm(i,j)) then
             fn=rm(i,j)
             fn_pass = fn
             fny_pass=am*(-3.*fn)
             fnyy_pass=am*(-5.*(am*fn+fny_pass))
           endif
         elseif(am.lt.0.) then
           if(fn.gt.0.) then
             fn=0.
             fn0=0.
             fmomn((/my,myy/))=0.
           elseif(fn.lt.-rm(i,j+1)) then
             fn=-rm(i,j+1)
             fn0 = fn
             fmomn(my)=am*(-3.*fn)
             fmomn(myy)=am*(-5.*(am*fn+fmomn(my)))
           endif
         endif

         mold=mass(i,j)
         mnew=mold+mvj(i)-am
         bymnew = 1./mnew
         dm2=mvj(i)+am
         rm0=rm(i,j)+fs0(i)-fn0
         rm(i,j)=rm(i,j)+fs(i)-fn
      !
         rmom(my,i,j)=(rmom(my,i,j)*mold-3.*(-dm2*rm0
     &     +mold*(fs0(i)+fn0))+(fmoms(my,i)-fmomn(my)))*bymnew
         rmom(myy,i,j) = (rmom(myy,i,j)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fs0(i)-fn0)-fmoms(my,i)
     &     -fmomn(my))+dm2*rmom(my,i,j)*mnew)
     &     +(fmoms(myy,i)-fmomn(myy))) * (bymnew*bymnew)
      ! cross moments
         rmom(mz,i,j)=rmom(mz,i,j)+fmoms(mz,i)-fmomn(mz)
         rmom(myz,i,j)=(rmom(myz,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmoms(mz,i)+fmomn(mz))) +
     &        (fmoms(myz,i)-fmomn(myz)))*bymnew
         rmom(mx,i,j)=rmom(mx,i,j)+fmoms(mx,i)-fmomn(mx)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(mx,i,j) +
     &        mold*(fmoms(mx,i)+fmomn(mx))) +
     &        (fmoms(mxy,i)-fmomn(mxy)))*bymnew
      !
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmoms(mzz,i)-fmomn(mzz)
         rmom(mxx,i,j)=rmom(mxx,i,j)+fmoms(mxx,i)-fmomn(mxx)
         rmom(mzx,i,j)=rmom(mzx,i,j)+fmoms(mzx,i)-fmomn(mzx)

         mass(i,j) = mnew

! clean up roundoff errors
         if(rm(i,j).le.0d0) then
           rm(i,j)=0d0; rmom(:,i,j)=0d0
         endif

         mvj(i) = am
         fs(i) = fn
         fmoms(:,i) = fmomn(:)

         fs0(i) = fn_pass
         fmoms(my,i) = fny_pass
         fmoms(myy,i) = fnyy_pass

      enddo ! i
      enddo ! j

      return
      end subroutine aadvqy

      subroutine aadvqz(rm,rmom,mass,mw,mwdn,fdn,fmomdn,fdn0)
!@sum  AADVQZ advection driver for z-direction
!@auth Maxwell Kelley
      implicit none
      REAL*8, dimension(isd:ied,jsd:jed,2) :: rm,mass
      REAL*8, dimension(nmom,isd:ied,jsd:jed,2) :: rmom
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: fmomdn
      REAL*8, dimension(isd:ied,jsd:jed) :: mw,mwdn,fdn,fdn0
      real*8, dimension(nmom) :: fmomup
      integer :: i,j,l,ll
      real*8 :: frac1,fracm,fup,mold,mnew,bymnew,dm2,am
      real*8 :: fup0,fup_pass,rm0,fupz_pass,fupzz_pass

c-----------------------------------------------------------
      ! calculate tracer mass flux f
      ! and fluxes of slopes and curvatures fmom
      ! update tracer mass, moments of tracer mass, air mass distribution
c--------------------------------------------------------------------
      do j=js,je
      do i=is,ie
         am = mw(i,j)
         if(am.lt.0.) then ! air mass flux is negative
            ll=2
            frac1=+1.
          else                   ! air mass flux is positive
            ll=1
            frac1=-1.
         endif
         fracm=am/mass(i,j,ll)
         frac1=fracm+frac1
         fup=fracm*(rm(i,j,ll)-frac1*(rmom(mz,i,j,ll)-
     &        (frac1+fracm)*rmom(mzz,i,j,ll)))
         fmomup(mz)=am*(fracm*fracm*(rmom(mz,i,j,ll)
     &        -3.*frac1*rmom(mzz,i,j,ll))-3.*fup)
         fmomup(mzz)=am*(am*fracm**3 *rmom(mzz,i,j,ll)
     &        -5.*(am*fup+fmomup(mz)))
         fmomup(my)  = fracm*(rmom(my,i,j,ll)-frac1*rmom(myz,i,j,ll))
         fmomup(myz) = am*
     &        (fracm*fracm*rmom(myz,i,j,ll)-3.*fmomup(my))
         fmomup(mx)  = fracm*(rmom(mx,i,j,ll)-frac1*rmom(mzx,i,j,ll))
         fmomup(mzx) = am*
     &        (fracm*fracm*rmom(mzx,i,j,ll)-3.*fmomup(mx))
         fmomup(myy) = fracm*rmom(myy,i,j,ll)
         fmomup(mxx) = fracm*rmom(mxx,i,j,ll)
         fmomup(mxy) = fracm*rmom(mxy,i,j,ll)

! flux limitations
         fup0 = fup
         fup_pass = fup
         fupz_pass = fmomup(mz)
         fupzz_pass = fmomup(mzz)
         if(am.gt.0.) then
           if(fup.lt.0.) then
             fup=0.
             fup_pass=0.
             fupz_pass=0
             fupzz_pass=0.
           elseif(fup.gt.rm(i,j,1)) then
             fup=rm(i,j,1)
             fup_pass = fup
             fupz_pass=am*(-3.*fup)
             fupzz_pass=am*(-5.*(am*fup+fupz_pass))
           endif
         elseif(am.lt.0.) then
           if(fup.gt.0.) then
             fup=0.
             fup0=0.
             fmomup((/mz,mzz/))=0.
           elseif(fup.lt.-rm(i,j,2)) then
             fup=-rm(i,j,2)
             fup0 = fup
             fmomup(mz)=am*(-3.*fup)
             fmomup(mzz)=am*(-5.*(am*fup+fmomup(mz)))
           endif
         endif

         mold=mass(i,j,1)
         mnew=mold+mwdn(i,j)-am
         bymnew = 1./mnew
         dm2=mwdn(i,j)+am
         rm0=rm(i,j,1)+fdn0(i,j)-fup0
         rm(i,j,1)=rm(i,j,1)+fdn(i,j)-fup
      !
         rmom(mz,i,j,1)=(rmom(mz,i,j,1)*mold-3.*(-dm2*rm0
     &     +mold*(fdn0(i,j)+fup0))+(fmomdn(mz,i,j)-fmomup(mz)))*bymnew
         rmom(mzz,i,j,1) = (rmom(mzz,i,j,1)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fdn0(i,j)-fup0)-fmomdn(mz,i,j)
     &     -fmomup(mz))+dm2*rmom(mz,i,j,1)*mnew)
     &     +(fmomdn(mzz,i,j)-fmomup(mzz))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j,1)=rmom(my,i,j,1)+fmomdn(my,i,j)-fmomup(my)
         rmom(myz,i,j,1)=(rmom(myz,i,j,1)*mold-3.*(-dm2*rmom(my,i,j,1) +
     &        mold*(fmomdn(my,i,j)+fmomup(my))) +
     &        (fmomdn(myz,i,j)-fmomup(myz)))*bymnew
         rmom(mx,i,j,1)=rmom(mx,i,j,1)+fmomdn(mx,i,j)-fmomup(mx)
         rmom(mzx,i,j,1)=(rmom(mzx,i,j,1)*mold-3.*(-dm2*rmom(mx,i,j,1) +
     &        mold*(fmomdn(mx,i,j)+fmomup(mx))) +
     &        (fmomdn(mzx,i,j)-fmomup(mzx)))*bymnew
      !
         rmom(myy,i,j,1)=rmom(myy,i,j,1)+fmomdn(myy,i,j)-fmomup(myy)
         rmom(mxx,i,j,1)=rmom(mxx,i,j,1)+fmomdn(mxx,i,j)-fmomup(mxx)
         rmom(mxy,i,j,1)=rmom(mxy,i,j,1)+fmomdn(mxy,i,j)-fmomup(mxy)

         mass(i,j,1) = mnew

! clean up roundoff errors
         if(rm(i,j,1).le.0d0) then
           rm(i,j,1)=0d0; rmom(:,i,j,1)=0d0
         endif

         mwdn(i,j) = am
         fdn(i,j) = fup
         fmomdn(:,i,j) = fmomup(:)

         fdn0(i,j) = fup_pass
         fmomdn(mz,i,j) = fupz_pass
         fmomdn(mzz,i,j) = fupzz_pass

      enddo ! i
      enddo ! j

      return
c****
      end subroutine aadvqz

      subroutine qusin(am,m,rm,rxm,rxxm,m_up,rm_up,rxm_up,rxxm_up,xsign)
c update rxm,rxxm as if upwind box had zero moments and r = fm/am
      implicit none
      real*8 :: m,m_up,am,rm,rxm,rxxm,rm_up,rxm_up,rxxm_up,xsign
      real*8 :: mnew,a,fm,dr,ddr,mrat
      A    = AM / M_UP
      FM   = A*(RM_UP + (1.-A)*(XSIGN*RXM_UP + (1.-2.*A)*RXXM_UP))
      MNEW   =  M + AM
      rxm = rxm*xsign
      dr  = rm/m - fm/am
      ddr = rxm/m-dr*(1d0-am/m)
      mrat = m/mnew
      RM   =  RM   + FM
      RXM  = (RXM  + 3.*AM*DR )*MRAT
      RXXM = (RXXM + 5.*AM*DDR)*MRAT*MRAT
      rxm = rxm*xsign
      return
      end subroutine qusin

      subroutine lusin(am,m,rm,rxm,m_up,rm_up,rxm_up,xsign)
c update rxm as if upwind box had zero moments and r = fm/am
      implicit none
      real*8 :: m,m_up,am,rm,rxm,rm_up,rxm_up,xsign
      real*8 :: mnew,a,fm,amdr
      A    = AM / M_UP
      FM   = A*(RM_UP + (1.-A)*(XSIGN*RXM_UP))
      MNEW   =  M + AM
      amdr = am*rm/m - fm
      RM  = RM + FM
      RXM = (RXM + 3.*XSIGN*AMDR) * M / MNEW
      return
      end subroutine lusin

      subroutine susin(am,rm,m_up,rm_up)
      implicit none
      real*8 :: m_up,am,rm,rm_up
      RM = RM + AM*RM_UP/M_UP
      return
      end subroutine susin

      subroutine qusout(am,m,rm,rxm,rxxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,rxxm,xsign
      real*8 :: bym,mnew,a,rl,rr,x
      real*8, parameter :: by3=1d0/3d0
      RXM = RXM*XSIGN
      bym = 1./m
      A  = AM*bym
      RL = (RM - RXM + RXXM)*bym
      x = 1.-2.*a
      RR = (RM + RXM*x + RXXM*1.5*(X*X-by3))*bym
      MNEW = M - AM
      RM = RM - A*(RM + (1.-A)*(RXM + (1.-2.*A)*RXXM))
      RXM = MNEW*.5*(RR-RL)
      RXXM = MNEW*.5*(RR+RL)-RM
      RXM = RXM*XSIGN
      return
      end subroutine qusout

      subroutine check_qusout(am,m,rm,rxm,rxxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,rxxm,xsign
      real*8 :: a,fxtra,rm0,arm
      A      = AM / M
      fxtra = A*(1.-A)*(RXM*XSIGN + (1.-2.*A)*RXXM)
      RM0 = RM*(1.-A)
      arm = a*rm
      if(RM0 .lt. fxtra) then
        rxm = rxm*rm0/fxtra
        rxxm = rxxm*rm0/fxtra
      elseif(fxtra.lt.-arm) then
        rxm = -rxm*arm/fxtra
        rxxm = -rxxm*arm/fxtra
      endif
      return
      end subroutine check_qusout

      subroutine lusout(am,m,rm,rxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,xsign
      real*8 :: mnew,a
      A      = AM / M
      MNEW   =  M-AM
      RM  = RM - A*(RM + (1.-A)*RXM*XSIGN)
      rxm = rxm*(mnew/m)**2
      return
      end subroutine lusout

      subroutine check_lusout(am,m,rm,rxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,rxxm,xsign
      real*8 :: a,fxtra,rm0
      A      = AM / M
      fxtra = A*(1.-A)*RXM*XSIGN
      RM0 = RM*(1.-A)
      if(RM0 .lt. fxtra) then
        rxm = rxm*rm0/fxtra
      elseif(fxtra.lt.-a*rm) then
        rxm = -rxm*a*rm/fxtra
      endif
      return
      end subroutine check_lusout

      subroutine lusout_2sides(am,bm,m,rm,rxm,rym,xsign,ysign)
      implicit none
      real*8 :: am,bm,m,rm,rxm,rym,xsign,ysign
      real*8 :: mnew,mnewx,mnewy,a,b
      A      = AM / M
      B      = BM / M
      MNEWX   =  M-AM
      MNEWY   =  M-BM
      MNEW    =  M-AM-BM
      RM  = RM -A*(RM + (1.-A)*RXM*XSIGN) -B*(RM + (1.-B)*RYM*YSIGN)
      rxm = rxm*(mnew*mnewx/(m*m))
      rym = rym*(mnew*mnewy/(m*m))
      return
      end subroutine lusout_2sides

      subroutine check_lusout_2sides(am,bm,m,rm,rxm,rym,xsign,ysign)
      implicit none
      real*8 :: am,bm,m,rm,rxm,rym,xsign,ysign
      real*8 :: a,b,fxtra,rm0
      A      = AM / M
      B      = BM / M
      fxtra = A*(1.-A)*RXM*XSIGN
      if(fxtra.lt.-a*rm) rxm = -rxm*a*rm/fxtra
      fxtra = B*(1.-B)*RYM*YSIGN
      if(fxtra.lt.-b*rm) rym = -rym*b*rm/fxtra
      fxtra = +A*(1.-A)*RXM*XSIGN + B*(1.-B)*RYM*YSIGN
      RM0 = RM*(1.-A-B)
      if(RM0 .lt. fxtra) then
        rxm = rxm*rm0/fxtra
        rym = rym*rm0/fxtra
      endif
      return
      end subroutine check_lusout_2sides

      subroutine lusin_2sides(am,bm
     &     ,m,rm,rxm,rym
     &     ,m_w,rm_w,rxm_w,rym_w
     &     ,m_s,rm_s,rxm_s,rym_s
     &     ,xsign,ysign)
      implicit none
      real*8 :: am,bm
     &     ,m,rm,rxm,rym
     &     ,m_w,rm_w,rxm_w,rym_w
     &     ,m_s,rm_s,rxm_s,rym_s
     &     ,xsign,ysign
      real*8 :: mnew,mnewx,mnewy,a,b,fm_w,fm_s,r_oldx,r_oldy
      RXM = RXM*XSIGN
      RYM = RYM*YSIGN
      A      = AM / M_W
      B      = BM / M_S
      MNEWX   =  M+AM
      MNEWY   =  M+BM
      MNEW    =  M+AM+BM
      FM_W   = A*(RM_W + (1.-A)*RXM_W*XSIGN)
      FM_S   = B*(RM_S + (1.-B)*RYM_S*YSIGN)
      R_OLDX = RM/M
      R_OLDY = RM/M
      R_OLDX = .5*(R_OLDX + (RM+FM_S)/MNEWY)
      R_OLDY = .5*(R_OLDY + (RM+FM_W)/MNEWX)
      RM  = RM + FM_W + FM_S
      RXM = RXM + B*RXM_S*XSIGN
      RYM = RYM + A*RYM_W*YSIGN
      RXM = (RXM + 3.*(AM*R_OLDX - FM_W)) * MNEWY/MNEW
      RYM = (RYM + 3.*(BM*R_OLDY - FM_S)) * MNEWX/MNEW
      RXM = RXM*XSIGN
      RYM = RYM*YSIGN
      return
      end subroutine lusin_2sides

      subroutine susin_2sides(am,bm,m,rm,m_w,rm_w,m_s,rm_s)
      implicit none
      real*8 :: am,bm,m,rm,m_w,rm_w,m_s,rm_s
      real*8 :: a,b,fm_w,fm_s
      A      = AM / M_W
      B      = BM / M_S
      RM  = RM + A*RM_W + B*RM_S
      return
      end subroutine susin_2sides

      end module tracer_adv

      subroutine alloc_tracer_adv(grid)
      use domain_decomp_atm, only : dist_grid
      use tracer_adv
      implicit none
      type(dist_grid), intent(in) :: grid

      gid  = grid%gid
      tile = grid%tile

      nxg = grid%npx
      nyg = grid%npy

      is = grid%is
      ie = grid%ie
      js = grid%js
      je = grid%je

      isd = grid%is-1
      ied = grid%ie+1
      jsd = grid%js-1
      jed = grid%je+1

      allocate(
     &     mu(isd:ied+1,jsd:jed,lm)
     &    ,mv(isd:ied,jsd:jed+1,lm)
     &    ,mw(isd:ied,jsd:jed,lm)
     &     )

      return
      end subroutine alloc_tracer_adv

      subroutine init_qus(grid,im_dum,jm_dum,lm_dum)
      use domain_decomp_atm, only : dist_grid
      use tracer_adv
      implicit none
      type(dist_grid), intent(in) :: grid
      integer, intent(in) :: im_dum,jm_dum,lm_dum

      xtoy_n = 0
      if(je==nyg.and.mod(tile,2).eq.1) xtoy_n = -1
      xtoy_s = 0
      if(js==1  .and.mod(tile,2).eq.0) xtoy_s = -1
      xtoy_w = 0
      if(is==1  .and.mod(tile,2).eq.1) xtoy_w = +1
      xtoy_e = 0
      if(ie==nxg.and.mod(tile,2).eq.0) xtoy_e = +1

      return
      end subroutine init_qus

      subroutine qdynam
      use atm_com, only : q,ps,mb,mma
      use somtq_com, only : qmom
      use tracer_adv
      use domain_decomp_atm, only : grid,halo_update
      implicit none
      integer :: i,j,l
      real*8 :: bymma

      call calc_amp(ps,mb)
      call halo_update(grid, mb) ! for convenience. but calc_amp could fill in the halo.
      call aadvq0

c
c topographic adjustments to moments
c
      call qmom_topo_adjustments


c
c convert from concentration to mass units
c
      do l=1,lm
      do j=js,je
      do i=is,ie
        q(i,j,l)=q(i,j,l)*mb(i,j,l)
        qmom(:,i,j,l)=qmom(:,i,j,l)*mb(i,j,l)
      enddo
      enddo
      enddo

c
c advect
c
      call aadvq (q,qmom, .true. ,'q       ')

c
c convert from mass to concentration units (using updated mma)
c
      do l=1,lm
      do j=js,je
      do i=is,ie
        bymma = 1 / mma(i,j,l)
        q(i,j,l) = q(i,j,l)*bymma
        qmom(:,i,j,l) = qmom(:,i,j,l)*bymma
      enddo
      enddo
      enddo
      return
      end subroutine qdynam

      subroutine aadvq0
      use tracer_adv
      use atm_com, only : MUs,MVs,MWs,mb,mma
      use domain_decomp_atm, only : grid,halo_update,am_i_root
      use domain_decomp_1d, only : globalmax,globalsum
      implicit none
      integer :: i,j,l
      integer :: nc3d,nbad,nc,iam,jbm,ncycxy_loc,nbad_loc
      real*8 :: byn,am,bm,xsigno,ysigno,mubyn,mvbyn,mwbyn
      real*8, dimension(isd:ied,jsd:jed,lm,2) :: uven
      real*8, dimension(isd:ied,jsd:jed) :: ma2d,mb2d

      call halo_update(grid,MUs)
      call halo_update(grid,MVs)
c at the top edge of an odd face, mv is the halo mu
      if(mod(tile,2).eq.1 .and. je.eq.nyg) then
        MVs(isd:ied,nyg+1,:) = MUs(isd:ied,nyg+1,:)
      endif
c at the right edge of an even face, mu is the halo mv
      if(mod(tile,2).eq.0 .and. ie.eq.nxg) then
        MUs(nxg+1,jsd:jed,:) = MVs(nxg+1,jsd:jed,:)
      endif

      do l=1,lm
      do j=jsd,jed
      do i=isd,ied
        mu(i,j,l) = MUs(i,j,l)
        mv(i,j,l) = MVs(i,j,l)
      enddo
      enddo
      enddo
      do l=1,lm-1
        do j=js,je
        do i=is,ie
          mw(i,j,l) = -MWs(i,j,l)
        enddo
        enddo
      enddo
      mw(:,:,lm) = 0d0
      call halo_update(grid,mw)

c      call mpp_update_domains(mu, mv, domain, gridtype=CGRID_SW)

      if(grid%nprocx.gt.1) then
        do l=1,lm
          do j=jsd,jed
            uven(is,j,l,1) = MUs(is+1,j,l)
          enddo
          do i=is,ie
            uven(i,js,l,2) = MVs(i,js+1,l)
          enddo
        enddo
        call halo_update(grid,uven)
        if(ie.lt.nxg) then
          do l=1,lm
            do j=jsd,jed
              mu(ie+2,j,l) = uven(ie+1,j,l,1)
            enddo
          enddo
        endif
        if(je.lt.nyg) then
          do l=1,lm
            do i=is,ie
              mv(i,je+2,l) = uven(i,je+1,l,2)
            enddo
          enddo
        endif
      endif


c
c Determine the horizontal-vertical ncyc
c
      nbad = 1
      ncyc = 0
      do while(nbad.gt.0)
        ncyc = ncyc + 1
        if(ncyc.gt.ncmax) then
          if(am_i_root()) write(6,*) 'stop: ncyc>ncmax in AADVQ0'
          call stop_model('AADVQ0: ncyc>ncmax',255)
        end if
        byn = 1./ncyc
        nbad_loc = 0
        do nc=1,ncyc

C**** check whether horizontal fluxes reduce mass too much
          do l=1,lm
            if(nc.eq.1) then
              mma(is:ie,js:je,l) = mb(is:ie,js:je,l)
            endif
            do j=js,je
            do i=is,ie
              mma(i,j,l) = mma(i,j,l) +
     &             byn*(mu(i,j,l)-mu(i+1,j,l)+mv(i,j,l)-mv(i,j+1,l))
              if (mma(i,j,l).lt.mrat_limh*mb(i,j,l)) then
                nbad_loc = nbad_loc + 1
              endif
            enddo
            enddo
          enddo ! l
          call globalsum(grid, nbad_loc, nbad, all=.true.)

          if(nbad.gt.0) exit    ! nc loop

c check courant numbers in the z direction
          do l=1,lm
            if(l.lt.lm) then
              do j=js,je
              do i=is,ie
                mwbyn = mw(i,j,l)*byn
                if((mma(i,j,l)-mwbyn)*(mma(i,j,l+1)+mwbyn).lt.0.) then
                  nbad_loc = nbad_loc + 1
                endif
              enddo
              enddo
            endif
c update mass from z fluxes
            if(nc.lt.ncyc) then
              if(l.eq.1) then   ! lowest layer
                do j=js,je
                do i=is,ie
                  mma(i,j,l) = mma(i,j,l)-mw(i,j,l)*byn
                enddo
                enddo
              elseif(l.eq.lm) then ! topmost layer
                do j=js,je
                do i=is,ie
                  mma(i,j,l) = mma(i,j,l)+mw(i,j,l-1)*byn
                enddo
                enddo
              else              ! interior layers
                do j=js,je
                do i=is,ie
                  mma(i,j,l) = mma(i,j,l)+(mw(i,j,l-1)-mw(i,j,l))*byn
                enddo
                enddo
              endif
            endif               ! nc.lt.ncyc
          enddo                 ! l
          call globalsum(grid, nbad_loc, nbad, all=.true.)

          if(nbad.gt.0) exit    ! nc loop

        enddo                   ! nc loop
      enddo                     ! nbad .gt. 0

      if(ncyc.gt.1) then
        byn = 1./ncyc
        do l=1,lm
          mu(:,:,l)=mu(:,:,l)*byn
          mv(:,:,l)=mv(:,:,l)*byn
          mw(:,:,l)=mw(:,:,l)*byn
        enddo
      endif

c
c Determine ncycxy for each level
c
      do l=1,lm
        ncycxy(l) = 1
        nc3dloop: do nc3d=1,ncyc
          if(nc3d.eq.1) then
            ma2d(isd:ied,jsd:jed) = mb(isd:ied,jsd:jed,l)
          else
            ma2d(isd:ied,jsd:jed) = mb2d(isd:ied,jsd:jed)
          endif
          nbad = 1
          do while(nbad.gt.0)
            if(ncycxy(l).gt.ncmax) exit nc3dloop
            byn = 1./ncycxy(l)
            nbad = 0
            do nc=1,ncycxy(l)

c
c courant number checks
c

c
c corners
c
              do j=1,nyg,nyg-1
                if(j.lt.js .or. j.gt.je) cycle
                if(j.eq.1) then
                  jbm=js
                  ysigno = -1.
                else
                  jbm=je+1
                  ysigno = +1.
                endif
                do i=1,nxg,nxg-1
                  if(i.lt.is .or. i.gt.ie) cycle
                  if(i.eq.1) then
                    iam=is
                    xsigno = -1.
                  else
                    iam=ie+1
                    xsigno = +1.
                  endif
                  am=mu(iam,j,l)*xsigno*byn
                  bm=mv(i,jbm,l)*ysigno*byn
                  if(am.gt.0. .and. bm.gt.0.) then ! flow exiting both x and y
                    if(am+bm.gt.ma2d(i,j)) then
                      nbad = nbad + 1
                    endif
                  elseif(am.lt.0. .and. bm.lt.0.) then ! flow entering both x and y
                  elseif(am.le.0.) then ! flow exits in y, then enters in x
                    if(bm.gt.ma2d(i,j)) then
                      nbad = nbad + 1
                    endif
                  elseif(bm.le.0.) then ! flow exits in x, then enters in y
                    if(am.gt.ma2d(i,j)) then
                      nbad = nbad + 1
                    endif
                  endif
                  ma2d(i,j) = ma2d(i,j) - am - bm
                enddo
              enddo
              if(nbad.gt.0) then ! exit already.  corners most likely to need ncycxy>1
                ncycxy(l) = ncycxy(l) + 1
                if(nc3d.eq.1) then
                  ma2d(isd:ied,jsd:jed) = mb(isd:ied,jsd:jed,l)
                else
                  ma2d(isd:ied,jsd:jed) = mb2d(isd:ied,jsd:jed)
                endif
                exit            ! nc loop
              endif


c
c east-west edges (excluding corners)
c
              do i=1,nxg,nxg-1
                if(i.lt.is .or. i.gt.ie) cycle
                if(i.eq.1) then
                  iam=is
                  xsigno = -1.
                else
                  iam=ie+1
                  xsigno = +1.
                endif
                do j=max(2,jsd),min(je,nyg-1)
                  am=mu(iam,j,l)*xsigno*byn
                  if(am.gt.ma2d(i,j)) then
                    nbad = nbad + 1
                  endif
                  ma2d(i,j) = ma2d(i,j) - am
                enddo
              enddo

c
c north-south edges (excluding corners)
c
              do j=1,nyg,nyg-1
                if(j.lt.js .or. j.gt.je) cycle
                if(j.eq.1) then
                  jbm=js
                  ysigno = -1.
                else
                  jbm=je+1
                  ysigno = +1.
                endif
                do i=max(2,isd),min(ie,nxg-1)
                  bm=mv(i,jbm,l)*ysigno*byn
                  if(bm.gt.ma2d(i,j)) then
                    nbad = nbad + 1
                  endif
                  ma2d(i,j) = ma2d(i,j) - bm
                enddo
              enddo

              if(nbad.gt.0) then
                ncycxy(l) = ncycxy(l) + 1
                if(nc3d.eq.1) then
                  ma2d(isd:ied,jsd:jed) = mb(isd:ied,jsd:jed,l)
                else
                  ma2d(isd:ied,jsd:jed) = mb2d(isd:ied,jsd:jed)
                endif
                exit            ! nc loop
              endif

c
c interior, x
c
              do j=js,je
                do i=max(2,is),ie
                  mubyn = mu(i,j,l)*byn
                  if((ma2d(i-1,j)-mubyn)*(ma2d(i,j)+mubyn).lt.0.) then
                    nbad = nbad + 1
                  endif
                enddo
                i=is; if(i.eq.1  ) then
                  ma2d(i,j) = ma2d(i,j) -mu(i+1,j,l)*byn
                  if(ma2d(i,j).lt.0.) nbad = nbad + 1
                endif
                i=ie; if(i.eq.nxg) then
                  ma2d(i,j) = ma2d(i,j) +mu(i  ,j,l)*byn
                  if(ma2d(i,j).lt.0.) nbad = nbad + 1
                endif
                do i=max(2,is),min(ie,nxg-1)
                  ma2d(i,j) = ma2d(i,j) + (mu(i,j,l)-mu(i+1,j,l))*byn
                  if(ma2d(i,j).lt.0.) nbad = nbad + 1
                enddo
              enddo
              if(nbad.gt.0) then
                ncycxy(l) = ncycxy(l) + 1
                if(nc3d.eq.1) then
                  ma2d(isd:ied,jsd:jed) = mb(isd:ied,jsd:jed,l)
                else
                  ma2d(isd:ied,jsd:jed) = mb2d(isd:ied,jsd:jed)
                endif
                exit            ! nc loop
              endif

c
c interior, y
c
              j=jsd; if(j.gt.0) then ! update s boundary airmasses
                i=is; if(i.eq.1  ) then
                  ma2d(i,j) = ma2d(i,j) -mu(i+1,j,l)*byn
                endif
                i=ie; if(i.eq.nxg) then
                  ma2d(i,j) = ma2d(i,j) +mu(i  ,j,l)*byn
                endif
                do i=max(2,is),min(ie,nxg-1)
                  ma2d(i,j) = ma2d(i,j) + (mu(i,j,l)-mu(i+1,j,l))*byn
                enddo
              endif
              do j=max(2,js),je
                do i=is,ie
                  mvbyn = mv(i,j,l)*byn
                  if((ma2d(i,j-1)-mvbyn)*(ma2d(i,j)+mvbyn).lt.0.) then
                    nbad = nbad + 1
                  endif
                enddo
              enddo
              if(nbad.gt.0) then
                ncycxy(l) = ncycxy(l) + 1
                if(nc3d.eq.1) then
                  ma2d(isd:ied,jsd:jed) = mb(isd:ied,jsd:jed,l)
                else
                  ma2d(isd:ied,jsd:jed) = mb2d(isd:ied,jsd:jed)
                endif
                exit            ! nc loop
              endif

c add interior mv contribution to airmasses
              if(nc3d.lt.ncyc .or. nc.lt.ncycxy(l)) then
                j=js; if(j.eq.1) then
                  do i=is,ie
                    ma2d(i,j) = ma2d(i,j) -mv(i,j+1,l)*byn
                  enddo
                endif
                j=je; if(j.eq.nyg) then
                  do i=is,ie
                    ma2d(i,j) = ma2d(i,j) +mv(i,j,l)*byn
                  enddo
                endif
                do j=max(2,js),min(je,nyg-1)
                  do i=is,ie
                    ma2d(i,j) = ma2d(i,j) + (mv(i,j,l)-mv(i,j+1,l))*byn
                  enddo
                enddo
              endif

c update w/s interior boundary airmasses to avoid halo updates
              if(nc3d.lt.ncyc .or. nc.lt.ncycxy(l)) then
                i=isd; if(i.gt.0) then
                  j=js; if(j.eq.1) then ! edge mv already done
                    ma2d(i,j) = ma2d(i,j) +
     &                 (mu(i,j,l)-mu(i+1,j,l)-mv(i,j+1,l))*byn
                  endif
                  j=je; if(j.eq.nyg) then ! edge mv already done
                    ma2d(i,j) = ma2d(i,j) +
     &                 (mu(i,j,l)-mu(i+1,j,l)+mv(i,j,l))*byn
                  endif
                  do j=max(2,js),min(je,nyg-1)
                    ma2d(i,j) = ma2d(i,j) +
     &                 (mu(i,j,l)-mu(i+1,j,l)+mv(i,j,l)-mv(i,j+1,l))*byn
                  enddo
                endif
                j=jsd; if(j.gt.0) then ! mu contribution already done
                  do i=is,ie
                    ma2d(i,j) = ma2d(i,j) + (mv(i,j,l)-mv(i,j+1,l))*byn
                  enddo
                endif
              endif

            enddo               ! nc loop
          enddo                 ! nbad.gt.0

c now add the z mass tendency at this level
          if(nc3d.lt.ncyc) then
            if(l.eq.1) then     ! lowest layer
              do j=max(1,jsd),min(nyg,jed)
              do i=max(1,isd),min(nxg,ied)
                mb2d(i,j) = ma2d(i,j)-mw(i,j,l)
              enddo
              enddo
            else if(l.eq.lm) then ! topmost layer
              do j=max(1,jsd),min(nyg,jed)
              do i=max(1,isd),min(nxg,ied)
                mb2d(i,j) = ma2d(i,j)+mw(i,j,l-1)
              enddo
              enddo
            else                ! interior layers
              do j=max(1,jsd),min(nyg,jed)
              do i=max(1,isd),min(nxg,ied)
                mb2d(i,j) = ma2d(i,j)+(mw(i,j,l-1)-mw(i,j,l))
              enddo
              enddo
            endif
          endif

        enddo nc3dloop

c globalmax of ncycxy. DOMAIN_DECOMP needs a vector version of globalmax
        ncycxy_loc = ncycxy(l)
        call globalmax(grid, ncycxy_loc, ncycxy(l))
        if(ncycxy(l).gt.ncmax) then
          if(am_i_root()) write(6,*)
     &         'stop: ncycxy>ncmax in AADVQ0',l,ncycxy(l)
          call stop_model('AADVQ0: ncycxy>ncmax',255)
        endif

      enddo ! l

c If more than 1/2 of the layers require ncycxy>1, increase ncyc
c and make only 1 xy cycle at each level.
      if(count(ncycxy.gt.1).gt.lm/2) then
        ncyc = ncyc*maxval(ncycxy)
        byn = 1./maxval(ncycxy) ! division by first ncyc already done
        do l=1,lm
          mu(:,:,l)=mu(:,:,l)*byn
          mv(:,:,l)=mv(:,:,l)*byn
          mw(:,:,l)=mw(:,:,l)*byn
        enddo
        ncycxy(:) = 1
      endif

c Divide the xy mass fluxes by the number of xy cycles at each level
      do l=1,lm
        if(ncycxy(l).gt.1) then
          byn = 1./ncycxy(l)
          mu(:,:,l)=mu(:,:,l)*byn
          mv(:,:,l)=mv(:,:,l)*byn
c          if(am_i_root()) write(6,*) 'ncycxy ',l,ncycxy(l),ncyc
        endif
      enddo

c      if(am_i_root().and. ncyc.gt.1) write(6,*) 'AADVQ0: ncyc>1',ncyc

      return
      end subroutine aadvq0

      SUBROUTINE AADVQ(RM,RMOM,qlimit,tname)
      use tracer_adv
      use atm_com, only : mma,mb
      use domain_decomp_atm, only : grid,halo_update
      IMPLICIT NONE
      REAL*8, dimension(isd:ied,jsd:jed,lm) :: rm
      REAL*8, dimension(nmom,isd:ied,jsd:jed,lm) :: rmom
      logical, intent(in) :: qlimit
      character*8 tname          !tracer name
c locals
      REAL*8, dimension(isd:ied,jsd:jed) :: mflx,mwdn,fdn,fdn0
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: fmomdn
      real*8, dimension(jsd:jed) :: mu_west,mu_east
      real*8, dimension(isd:ied) :: mv_south,mv_north

      INTEGER :: I,J,L,nc,ncxy

C****
C**** Load mass after advection from mass before advection
C****
      DO L=1,LM
         MMA(:,:,L) = MB(:,:,L) ! fill in halo also
      ENDDO

      do nc=1,ncyc

c horizontal flux limits at edges
      do l=1,lm
        call check_edges_and_corners(mma(isd,jsd,l),rm(isd,jsd,l),
     &       rmom(1,isd,jsd,l),mu(isd,jsd,l),mv(isd,jsd,l))
      enddo

c 3D halo updates
      if(nc.gt.1) call halo_update(grid,mma)
      call halo_update(grid,rm)
      call halo_update(grid,rmom,jdim=3)

      do j=js,je
      do i=is,ie
        mwdn(i,j) = 0.
        fdn(i,j) = 0.
        fdn0(i,j) = 0.
        fmomdn(:,i,j) = 0.
      enddo
      enddo

      DO L=1,LM+1 ! lm+1 is for the uppermost aadvqz call

        if(l.le.lm) then

        do ncxy=1,ncycxy(l) ! loop over horizontal cycles

          if(ncxy.gt.1) then    ! halo update this level if necessary
            call check_edges_and_corners(mma(isd,jsd,l),rm(isd,jsd,l),
     &           rmom(1,isd,jsd,l),mu(isd,jsd,l),mv(isd,jsd,l))
            call halo_update(grid, mma(:,:,l))
            call halo_update(grid, rm(:,:,l))
            call halo_update(grid, rmom(:,:,:,l), jdim=3)
          endif

c horz transport across the edges of faces
          call do_edges_and_corners(mma(isd,jsd,l),rm(isd,jsd,l),
     &         rmom(1,isd,jsd,l),mu(isd,jsd,l),mv(isd,jsd,l))

c horz transport in the interior of faces.
c mu,mv temporarily set to zero at edges.
          if(is.eq.1) then
            do j=jsd,jed
              mu_west(j) = mu(is  ,j,l)
              mu(is  ,j,l) = 0d0
            enddo
          endif
          if(ie.eq.nxg) then
            do j=jsd,jed
              mu_east(j) = mu(ie+1,j,l)
              mu(ie+1,j,l) = 0d0
            enddo
          endif
          CALL AADVQX(RM(isd,jsd,l),RMOM(1,isd,jsd,l),
     &         MMA(isd,jsd,l),mu(isd,jsd,l))
          if(is.eq.1) then
            do j=jsd,jed
              mu(is  ,j,l) = mu_west(j)
            enddo
          endif
          if(ie.eq.nxg) then
            do j=jsd,jed
              mu(ie+1,j,l) = mu_east(j)
            enddo
          endif

          if(js.eq.1) then
            do i=is,ie
              mv_south(i) = mv(i,js  ,l)
              mv(i,js  ,l) = 0d0
            enddo
          endif
          if(je.eq.nyg) then
            do i=is,ie
              mv_north(i) = mv(i,je+1,l)
              mv(i,je+1,l) = 0d0
            enddo
          endif
          CALL AADVQY(RM(isd,jsd,l),RMOM(1,isd,jsd,l),
     &         MMA(isd,jsd,l),mv(isd,jsd,l))
          if(js.eq.1) then
            do i=is,ie
              mv(i,js  ,l) = mv_south(i)
            enddo
          endif
          if(je.eq.nyg) then
            do i=is,ie
              mv(i,je+1,l) = mv_north(i)
            enddo
          endif

        enddo ! end loop over horizontal cycles

        endif

c vertical transport
        if(l.gt.1 .and. l.lt.lm) then
          call checkfobs_z(mma(isd,jsd,l),mw(isd,jsd,l-1),
     &         rm(isd,jsd,l),rmom(1,isd,jsd,l))
        endif
        if(l.gt.1) then
          CALL AADVQZ(RM(isd,jsd,l-1),RMOM(1,isd,jsd,l-1),
     &         MMA(isd,jsd,l-1),mw(isd,jsd,l-1),
     &         mwdn,fdn,fmomdn,fdn0)
        endif

      ENDDO ! end loop over l

      enddo ! nc=1,ncyc

      RETURN
      END SUBROUTINE AADVQ

#ifdef TRACERS_ON
      SUBROUTINE TrDYNAM
      USE MODEL_COM, only: itime
      USE TRACER_COM, only: itime_tr0,trm,trmom,trname,t_qlimit,ntm
      IMPLICIT NONE
      INTEGER N

      DO N=1,NTM
        IF (itime.LT.itime_tr0(N)) cycle
        CALL AADVQ (TRM(:,:,:,n),TrMOM(:,:,:,:,n),t_qlimit(n),trname(n))
      ENDDO

      RETURN
      END SUBROUTINE TrDYNAM
#endif
