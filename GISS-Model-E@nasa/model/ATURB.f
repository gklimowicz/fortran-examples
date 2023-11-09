#include "rundeck_opts.h"

      subroutine atm_diffus(lbase_min,lbase_max,dtime)
!@sum  atm_diffus updates u,v,t,q due to turbulent transport throughout
!@+    all GCM layers using a non-local turbulence model
!@vers 2013/03/27
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@cont atm_diffus,getdz,dout,de_solver_main,de_solver_edge,l_gcm,k_gcm,
!@+    e_gcm,find_pbl_top,zze,apply_fluxes_to_atm
!@var lbase_min/max levels through which to apply turbulence (dummy)
!@var dtime time step

!@var qmin minimum value of specific humidity
!@var itest/jtest longitude/latitude at which dout may be called
!@var call_diag logical variable whether dout is called

      USE CONSTANT, only : grav,deltx,lhe,sha,by3,teeny,mb2kg
      USE RESOLUTION, only : im,jm,lm
      USE RESOLUTION, only : psf,pmtop
      USE MODEL_COM, only : itime
      USE ATM_COM, only : u_3d=>u,v_3d=>v,t_3d=>t,q_3d=>q
cc      USE QUSDEF, only : nmom,zmoms,xymoms
cc      USE SOMTQ_COM, only : tmom,qmom
      USE GEOM, only : imaxj,byaxyp,axyp
      Use ATM_COM,    Only: MA,byMA, PDSIG,PMID,PEDN,PK,PEK
     &     ,u_3d_agrid=>ualij,v_3d_agrid=>valij
      USE DOMAIN_DECOMP_ATM, ONLY : grid, getDomainBounds, halo_update
      USE DIAG_COM, only : jl_trbhr,jl_damdc,jl_trbke,jl_trbdlht
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,itime_tr0,trm,t_qlimit  !,trmom
      USE TRDIAG_COM, only: jlnt_turb
      USE FLUXES, only : trflux1
#endif
      USE SOCPBL, only : b1,b123,prt,kappa,zgs,ustar_min
     &  ,lmonin_min,lmonin_max
      USE PBLCOM, only : dclev,pblht,pblptop
     *     ,e_3d=>egcm,w2_3d=>w2gcm !,t2_3d=>t2gcm
     &     ,t1_after_aturb,u1_after_aturb,v1_after_aturb
      USE FLUXES, only : uflux1,vflux1,tflux1,qflux1,atmsrf


      IMPLICIT NONE

      integer, intent(in) :: lbase_min,lbase_max
      real*8, intent(in) :: dtime

      real*8, parameter :: qmin=1.d-12,pblp_max=400.!mb
      integer, parameter :: itest= 1,jtest=11
      logical, parameter :: call_diag=.false.

      real*8, dimension(lm) :: u,v,t,q,e,u0,v0,t0,q0,e0,p
     &    ,dudz,dvdz,dtdz,dqdz,g_alpha,as2,an2
     &    ,rhoebydz,bydzerho,rho,rhoe,dz,dze
     &    ,km,kh,kq,ke,wt_nl,wq_nl,uw_nl,vw_nl
     &    ,lscale,qturb,p3,p4,rhobydze,bydzrhoe,w2,uw,vw,wt,wq,z
      real*8, dimension(lm+1) :: ze

      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo) ::
     &     rho_3d,rhoe_3d,dz_3d,dze_3d,km_3d,t_3d_virtual
     &     ,uasv,uw_nl_3d,vw_nl_3d ! for wind tendency diagnostic
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     tvsurf,dz0,uflxa,vflxa
      real*8, dimension((1+grid%i_stop_halo-grid%i_strt_halo)*
     &                  (1+grid%j_stop_halo-grid%j_strt_halo)*2) ::
     &     uvflux_vgrid
cc      real*8, dimension(nmom,lm) :: tmomij,qmomij

      real*8 :: uflx,vflx,tvflx,qflx,tvs
     &   ,ustar2,t0ijl,tijl,rak,ustar
     &   ,flux_bot,flux_top,x_surf
     &   ,wstar,dbl,lmonin,tpe0,tpe1,ediff,den,tmp
      integer :: idik,idjk,ldbl,kmax,ldbl_max,
     &    i,j,l,k,n,iter !@i,j,l,k,n,iter loop variable
#ifdef TRACERS_ON
!@var tr0ij initial vertical tracer concentration profile (kg/kg)
!@var trij vertical tracer concentration profile (kg/kg)
!@var trmomij vertical tracer concentration moment profile (kg/kg)
!@var wc_nl non-local fluxes of tracers
      real*8, dimension(lm,ntm) :: tr0ij,trij,wc_nl
      real*8, dimension(lm) :: dtrm,amkg,byMMA
cc      real*8, dimension(nmom,lm,ntm) :: trmomij
!@var trflx surface tracer flux (-w tr) (kg/kg m/s)
      real*8, dimension(ntm) :: trflx
      integer nta,nx,ntix(ntm)
#endif

c vars for velocity diffusion
      real*8, dimension(lm) :: uv,uv0
      integer, dimension(4) :: ilist,jlist
      real*8, dimension(4) :: wts
      integer :: nuv,vpkey,vpkey_last,nnbr

      INTEGER :: I_0, I_1, J_1, J_0, J_0H, J_1H
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

      ! Note that lbase_min/max are here for backwards compatibility
      ! with original drycnv. They are only used to determine where the
      ! routine has been called from.

      if (lbase_min.eq.2) return       ! quit if called from main

      !  convert input T to virtual T

      do j=J_0, J_1
        do i=I_0,imaxj(j)
          !@var tvsurf(i,j) surface virtual temperature
          !@var tsavg(i,j) composite surface air temperature (k)
          tvsurf(i,j)=atmsrf%tsavg(i,j)*(1.d0+deltx*atmsrf%qsavg(i,j))
          do l=1,lm
            ! t_3d_virtual is virtual potential temp. referenced at 1 mb
            t_3d_virtual(l,i,j)=t_3d(i,j,l)*(1.d0+deltx*q_3d(i,j,l))
          end do
        end do
      end do

#ifdef TRACERS_ON
      nx=0
      do n=1,ntm
        if (itime_tr0(n).le.itime) then
          nx=nx+1
          ntix(nx)=n
        end if
      end do
      nta=nx
#endif

      ! integrate equations other than u,v at agrids

      ! get u_3d_agrid and v_3d_agrid
#ifndef SCM
c      call ave_uv_to_agrid(u_3d,v_3d,u_3d_agrid,v_3d_agrid,lm)
#else
      do j=J_0,J_1
         do i=1,imaxj(j)
            do L=1,LM
               u_3d_agrid(L,i,j) = u_3d(i,j,L)
               v_3d_agrid(L,i,j) = v_3d(i,j,L)
            enddo
         enddo
      enddo
#endif

      call getdz(t_3d_virtual,dz_3d,dze_3d,rho_3d,rhoe_3d,tvsurf
     &     ,dz0,im,jm,lm)

      loop_j_tq: do j=J_0, J_1
        loop_i_tq: do i=I_0,imaxj(j)

          do l=1,lm
            u(l)=u_3d_agrid(l,i,j)
            v(l)=v_3d_agrid(l,i,j)
            uasv(l,i,j) = u_3d_agrid(l,i,j)
            ! virtual potential temp. referenced at 1 mb
            t(l)=t_3d_virtual(l,i,j)
            q(l)=q_3d(i,j,l)
cc            qmomij(:,l)=qmom(:,i,j,l)
cc            tmomij(:,l)=tmom(:,i,j,l) ! vert. grad. should virtual ?
c            if(q(l).lt.qmin) q(l)=qmin
            e(l)=e_3d(l,i,j) !e_3d was called egcM
            p(l)=pmid(l,i,j)
            ! t2(l)=t2_3d(l,i,j)  ! not in use
            rho(l)=rho_3d(l,i,j)
            rhoe(l)=rhoe_3d(l,i,j)
            t0(l)=t(l)
            q0(l)=q(l)
            e0(l)=e(l)
            qturb(l)=sqrt(2.d0*e(l))
            dze(l)=dze_3d(l,i,j)
            dz(l)=dz_3d(l,i,j)
            bydzerho(l)=1.d0/(dze(l)*rho(l))
            rhobydze(l)=rho(l)/dze(l)
          end do
          bydzrhoe(1)=0.
          rhoebydz(1)=0.
          do l=1,lm-1
            bydzrhoe(l+1)=1.d0/(dz(l)*rhoe(l+1))
            rhoebydz(l+1)=rhoe(l+1)/dz(l)
          end do

#ifdef TRACERS_ON
          do l=1,lm
            byMMA(l) = byMA(l,i,j)*byaxyp(i,j)
          enddo
          do nx=1,nta
            n=ntix(nx)
            do l=1,lm
              trij(l,nx) = trm(i,j,l,n)*byMMA(l)
cc                trmomij(:,l,nx)=trmom(:,i,j,l,n)
              tr0ij(l,nx)=trij(l,nx)
            enddo
          enddo
#endif

          ! tvs is surface virtual potential temp. referenced at 1 mb
          tvs=tvsurf(i,j)/pek(1,i,j)
          uflx=uflux1(i,j)/rhoe(1)
          vflx=vflux1(i,j)/rhoe(1)
          qflx=qflux1(i,j)/rhoe(1)
          ! tvflx is virtual, potential temp. flux referenced at 1 mb
          tvflx=tflux1(i,j)*(1.d0+deltx*atmsrf%qsavg(i,j))/
     &         (rhoe(1)*pek(1,i,j))
     &         +deltx*atmsrf%tsavg(i,j)/pek(1,i,j)*qflx
          uflxa(i,j)=uflx
          vflxa(i,j)=vflx
#ifdef TRACERS_ON
          do nx=1,nta
            n=ntix(nx)
C**** minus sign needed for ATURB conventions
            trflx(nx)=-trflux1(i,j,n)/rhoe(1)
          end do
#endif

          if(abs(uflx).lt.teeny) uflx=sign(teeny,uflx)
          if(abs(vflx).lt.teeny) vflx=sign(teeny,vflx)
          ustar=(uflx*uflx+vflx*vflx)**(0.25d0)
          ustar=max(ustar,ustar_min)
          ustar2=ustar*ustar

          ! calculate z-derivatives at the surface

          ! @var zgs height of surface layer (m), imported from SOCPBL
          tmp=1d0/(ustar*kappa*zgs)
          dudz(1)=uflx*tmp
          dvdz(1)=vflx*tmp
          dtdz(1)=tvflx*prt*tmp
          dqdz(1)= qflx*prt*tmp

          g_alpha(1)=grav/tvs
          den=kappa*g_alpha(1)*tvflx
          if(den.eq.0.) den=teeny
          lmonin=ustar**3/den
          if(abs(lmonin).lt.lmonin_min) lmonin=sign(lmonin_min,lmonin)
          if(abs(lmonin).gt.lmonin_max) lmonin=sign(lmonin_max,lmonin)

          ! calculate z-derivatives on the edges of the layers

          do l=2,lm
            dudz(l)=(u(l)-u(l-1))/dz(l-1)
            dvdz(l)=(v(l)-v(l-1))/dz(l-1)
            dtdz(l)=(t(l)-t(l-1))/dz(l-1)
            dqdz(l)=(q(l)-q(l-1))/dz(l-1)
            g_alpha(l)=grav*2.d0/(t(l)+t(l-1))
          end do

          !@var an2 brunt-vassala frequency
          !@var as2 shear number squared

          do l=1,lm
            an2(l)=g_alpha(l)*dtdz(l)
            if(abs(an2(l)).lt.teeny) an2(l)=sign(teeny,an2(l))
            as2(l)=max(dudz(l)*dudz(l)+dvdz(l)*dvdz(l),teeny)
          end do

          call zze(dz,dze,dz0(i,j),z,ze,lm)

          !@var ldbl_max the maximum number of layers allowed in the pbl
          ldbl_max=1
          do l=2,lm
            if(p(l).ge.pblp_max) then
              ldbl_max=ldbl_max+1
            else
              exit
            endif
          end do

          ! calculate pbl depth and pbl top level

          call find_pbl_top(z,u,v,t,ustar,ustar2,tvflx,lmonin
     &        ,dbl,ldbl,ldbl_max,lm)


          if(tvflx.lt.0.) then ! convective case
            wstar=(-g_alpha(1)*tvflx*dbl)**by3
          else
            wstar=teeny
          endif

          ! calculate turbulence length scale lscale

          call l_gcm(ze,dbl,lmonin,ustar,qturb,an2,lscale,lm)

          ! calculate turbulent kinetic energy e

          call e_gcm(tvflx,wstar,ustar,dbl,lmonin,ze,g_alpha
     &              ,an2,as2,lscale,e,lm)

          ! calculate turbulent diffusivities and non-local terms

          call k_gcm(tvflx,qflx,uflx,vflx,ustar,wstar,dbl,lmonin
     &        ,ze,lscale,e,qturb,an2,as2,dtdz,dqdz,dudz,dvdz
     &        ,kh,kq,km,ke,wt,wq,w2,uw,vw,wt_nl,wq_nl,uw_nl,vw_nl
#ifdef TRACERS_ON
     &        ,trflx,wc_nl,nta
#endif
     &        ,lm)

c         ! integrate differential eqn for e
c         p3(1)=0. ; p3(lm)=0.
c         p4(1)=0. ; p4(lm)=0.
c         do l=2,lm-1
c             p3(l)=2.d0*(qturb(l)/(b1*lscale(l)))
c             p4(l)=-uw(l)*dudz(l)-vw(l)*dvdz(l)+g_alpha(l)*wt(l)
c             ! p4(l)=km(l)*as2(l)-kh(l)*an2(l)
c         end do
c         x_surf=0.5d0*b123*ustar2
c         call de_solver_edge(e,e0,ke,p3,p4,
c    &        rhobydze,bydzrhoe,x_surf,dtime,lm,.true.)

          do l=1,lm
              qturb(l)=sqrt(2.d0*e(l))
          end do

          ! integrate differential eqn for T
          do l=2,lm-1
              p4(l)=-(rhoe(l+1)*wt_nl(l+1)-rhoe(l)*wt_nl(l))
     &              *bydzerho(l)
          end do
          flux_bot=rhoe(1)*tvflx+rhoe(2)*wt_nl(2)
          flux_top=rhoe(lm)*wt_nl(lm)
          call de_solver_main(t,t0,kh,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)

C**** also diffuse moments
cc        call diff_mom(tmomij)

          ! integrate differential eqn for Q
          flux_bot=rhoe(1)*qflx+rhoe(2)*wq_nl(2)
C**** fix first layer for rare tracer problems
C**** Does this ever happen for q? (put this in just in case)
            if ( q0(1)-dtime*bydzerho(1)*flux_bot.lt.0 ) then
              flux_bot=q0(1)/(dtime*bydzerho(1))
              wq_nl(2)=(flux_bot-rhoe(1)*qflx)/rhoe(2)
            end if
          do l=2,lm-1
              p4(l)=-(rhoe(l+1)*wq_nl(l+1)-rhoe(l)*wq_nl(l))
     &              *bydzerho(l)
C**** check on physicality of non-local fluxes....
              if ( p4(l)*dtime+q0(l).lt.0 ) then
                p4(l)=-q(l)/dtime
                wq_nl(l+1)=(q0(l)/(dtime*bydzerho(l))+rhoe(l)
     *               *wq_nl(l))/rhoe(l+1)
              end if
          end do
          flux_top=rhoe(lm)*wq_nl(lm)
          call de_solver_main(q,q0,kq,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.true.)
c          do l=1,lm
c              if(q(l).lt.qmin) q(l)=qmin
c          end do

C**** also diffuse moments
cc        call diff_mom(qmomij)

#ifdef TRACERS_ON
C**** Use q diffusion coefficient for tracers
C**** Note that non-local effects for tracers can be included
C**** parallel to the case of Q
          do n=1,nta
            flux_bot=rhoe(1)*trflx(n)+rhoe(2)*wc_nl(2,n) !tr0ij(1,n)
C**** fix first layer for rare tracer problems
            if ( t_qlimit(n) .and.
     &           tr0ij(1,n)-dtime*bydzerho(1)*flux_bot.lt.0 ) then
              flux_bot=tr0ij(1,n)/(dtime*bydzerho(1))
              wc_nl(2,n)=(flux_bot-rhoe(1)*trflx(n))/rhoe(2)
            end if
            
            do l=2,lm-1
              p4(l)=-(rhoe(l+1)*wc_nl(l+1,n)-rhoe(l)*wc_nl(l,n))
     &             *bydzerho(l)
C**** check on physicality of non-local fluxes....
              if ( t_qlimit(n) .and. p4(l)*dtime+tr0ij(l,n).lt.0 ) then
                p4(l)=-tr0ij(l,n)/dtime
                wc_nl(l+1,n)=(tr0ij(l,n)/(dtime*bydzerho(l))+rhoe(l)
     *               *wc_nl(l,n))/rhoe(l+1)
              end if
            end do
            flux_top=rhoe(lm)*wc_nl(lm,n)
            call de_solver_main(trij(1,n),tr0ij(1,n),kq,p4,
     &        rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,t_qlimit(n))
cc          call diff_mom(trmomij)

         end do
#endif
          dclev(i,j)=real(ldbl)
          pblht(i,j)=dbl

          if(ldbl.le.1) then
            pblptop(i,j)=p(1)
          else
            pblptop(i,j) = p(ldbl-1)+
     &           (p(ldbl)-p(ldbl-1))*(dbl-z(ldbl-1))/(z(ldbl)-z(ldbl-1))
          endif
C**** calculate possible energy loss
          tpe0=-tflux1(i,j)*dtime*sha
          tpe1=0.
          do l=1,lm
            tpe0=tpe0+t_3d(i,j,l)*pk(l,i,j)*pdsig(l,i,j)*sha*mb2kg
            tpe1=tpe1+t(l)*pk(l,i,j)*pdsig(l,i,j)*sha*mb2kg/(1.d0+deltx
     *           *q(l))
          end do
          ediff=(tpe1-tpe0)/((psf-pmtop)*sha*mb2kg)        ! C

          do l=1,lm
C**** correct virt pot t for energy error
            t(l) = t(l) - ediff*(1.d0+deltx*q(l))/pk(l,i,j)
            ! update 3-d t,q,e and km
            t0ijl=t_3d(i,j,l)
            tijl=t(l)/(1.d0+deltx*q(l))
            t_3d(i,j,l)=tijl
C**** moment variation to be added
cc            qmom(:,i,j,l)=qmomij(:,l)
cc            tmom(:,i,j,l)=tmomij(:,l)

            q_3d(i,j,l)=q(l)
            e_3d(l,i,j)=e(l)
            w2_3d(l,i,j)=w2(l)
            ! t2_3d(l,i,j)=t2(l)  ! not in use
            km_3d(l,i,j)=km(l)
            uw_nl_3d(l,i,j)=uw_nl(l)
            vw_nl_3d(l,i,j)=vw_nl(l)
            ! ACCUMULATE DIAGNOSTICS for t and q
            call inc_ajl(i,j,l,JL_TRBHR,
     &           (tijl-t0ijl)*PK(L,I,J)*PDSIG(L,I,J))
            call inc_ajl(i,j,l,JL_TRBDLHT,
     &           (q(l)-q0(l))*PDSIG(L,I,J)*LHE/SHA)
            call inc_ajl(i,j,l,JL_TRBKE,e(l))
          end do

          t1_after_aturb(i,j) = t_3d(i,j,1)

#ifdef TRACERS_ON
          do l=1,lm
            amkg(l) = MA(l,i,j)*axyp(i,j)
          enddo
          do nx=1,nta
            n=ntix(nx)
            do l=1,lm
              trm(i,j,l,n)=trij(l,nx)*amkg(l)
cc            trmom(:,i,j,l,n)=trmomij(:,l,nx)
#ifndef SKIP_TRACER_DIAGS
              dtrm(l) = (trij(l,nx)-tr0ij(l,nx))*amkg(l)
#endif
            enddo
#ifndef SKIP_TRACER_DIAGS
            call inc_tajln_column(i,j,1,lm,lm,jlnt_turb,n,dtrm)
#endif
          enddo
#endif

          ! Write out diagnostics if at selected grid point:

          if (call_diag.and.(i.eq.itest).and.(j.eq.jtest)) then
            call dout(z,ze,dz,u,v,t,q,ke,dtdz,dqdz,an2,as2
     &         ,wt_nl,wq_nl,kh,km,e,lscale
     &         ,uflx,vflx,tvflx,qflx,dbl,ldbl,i,j,lm)
          endif

#ifdef SCM
c diffuse velocities on the primary grid in single-column model
          u0(:) = u(:)
          v0(:) = v(:)
          flux_top=0.
          do l=2,lm-1
            p4(l)=-(rhoe(l+1)*uw_nl(l+1)-rhoe(l)*uw_nl(l))
     &            *bydzerho(l)
          end do
          flux_bot=uflux1(i,j)+rhoe(2)*uw_nl(2)
          call de_solver_main(u,u0,km,p4,
     &       rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)
          do l=2,lm-1
            p4(l)=-(rhoe(l+1)*vw_nl(l+1)-rhoe(l)*vw_nl(l))
     &            *bydzerho(l)
          end do
          flux_bot=vflux1(i,j)+rhoe(2)*vw_nl(2)
          call de_solver_main(v,v0,km,p4,
     &       rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)
          u_3d(i,j,:) = u(:)
          v_3d(i,j,:) = v(:)
#endif

        end do loop_i_tq
      end do loop_j_tq

#ifndef SCM
c
c integrate differential eqns for U and V on velocity grid
c

c
c first, fill poles for km,rho,dz to make life easier.
c
      if(have_south_pole) then
        j=1
        do i=2,im
        do l=1,lm
          km_3d(l,i,j)   = km_3d(l,1,j)
          uw_nl_3d(l,i,j)= uw_nl_3d(l,1,j)
          vw_nl_3d(l,i,j)= vw_nl_3d(l,1,j)
          dz_3d(l,i,j)   = dz_3d(l,1,j)
          dze_3d(l,i,j)  = dze_3d(l,1,j)
          rho_3d(l,i,j)  = rho_3d(l,1,j)
          rhoe_3d(l,i,j) = rhoe_3d(l,1,j)
        enddo
        enddo
      endif
      if(have_north_pole) then
        j=jm
        do i=2,im
        do l=1,lm
          km_3d(l,i,j)   = km_3d(l,1,j)
          uw_nl_3d(l,i,j)= uw_nl_3d(l,1,j)
          vw_nl_3d(l,i,j)= vw_nl_3d(l,1,j)
          dz_3d(l,i,j)   = dz_3d(l,1,j)
          dze_3d(l,i,j)  = dze_3d(l,1,j)
          rho_3d(l,i,j)  = rho_3d(l,1,j)
          rhoe_3d(l,i,j) = rhoe_3d(l,1,j)
        enddo
        enddo
      endif


      CALL HALO_UPDATE(grid, km_3d, jdim=3)
      CALL HALO_UPDATE(grid, uw_nl_3d, jdim=3)
      CALL HALO_UPDATE(grid, vw_nl_3d, jdim=3)
      CALL HALO_UPDATE(grid, dz_3d, jdim=3)
      CALL HALO_UPDATE(grid, dze_3d, jdim=3)
      CALL HALO_UPDATE(grid, rho_3d, jdim=3)
      CALL HALO_UPDATE(grid, rhoe_3d, jdim=3)
c
c put A-grid surface stresses on the velocity grid
c
      call regrid_atov_1d(uflxa,vflxa,uvflux_vgrid)
c
c loop over U,V data
c
      flux_top=0.d0
      vpkey_last = -1000
      call get_nuv(nuv)
      do n=1,nuv
c copy 3D velocity to 1D array
        call get_uv_of_n(n,uv)
        uv0(:) = uv(:)
        call get_vpkey_of_n(n,vpkey)
        if(vpkey.ne.vpkey_last) then
c we have moved to a new velocity point.
c interpolate a-grid rho,km,dz to the new point.
          call get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
          do l=1,lm
            rho(l)=0d0
            rhoe(l)=0d0
            km(l)=0d0
            uw_nl(l)=0d0
            vw_nl(l)=0d0
            dze(l)=0d0
            dz(l)=0d0
          enddo
          do k=1,nnbr
            i = ilist(k)
            j = jlist(k)
            do l=1,lm
              rho(l) =rho(l) +wts(k)*rho_3d(l,i,j)
              rhoe(l)=rhoe(l)+wts(k)*rhoe_3d(l,i,j)
              km(l)  =km(l)  +wts(k)*km_3d(l,i,j)
              uw_nl(l)=uw_nl(l)+wts(k)*uw_nl_3d(l,i,j)
              vw_nl(l)=vw_nl(l)+wts(k)*vw_nl_3d(l,i,j)
              dze(l) =dze(l) +wts(k)*dze_3d(l,i,j)
              dz(l)  =dz(l)  +wts(k)*dz_3d(l,i,j)
            enddo
          enddo
          do l=1,lm
            bydzerho(l)=1.d0/(dze(l)*rho(l))
          enddo
          rhoebydz(1)=0.
          do l=1,lm-1
            rhoebydz(l+1)=rhoe(l+1)/dz(l)
          enddo
        endif
c perform the diffusion
        flux_bot=rhoe(1)*uvflux_vgrid(n)
        if(mod(n,2).eq.1) then
          flux_bot=flux_bot+rhoe(2)*uw_nl(2)
          do l=2,lm-1
            p4(l)=-(rhoe(l+1)*uw_nl(l+1)-rhoe(l)*uw_nl(l))
     &            *bydzerho(l)
          end do
        else
          flux_bot=flux_bot+rhoe(2)*vw_nl(2)
          do l=2,lm-1
            p4(l)=-(rhoe(l+1)*vw_nl(l+1)-rhoe(l)*vw_nl(l))
     &            *bydzerho(l)
          end do
        endif
        call de_solver_main(uv,uv0,km,p4,
     &       rhoebydz,bydzerho,flux_bot,flux_top,dtime,lm,.false.)
c store the updated velocity
        call store_uv_of_n(n,uv)
        vpkey_last = vpkey
      enddo
#endif

c
c wind tendency diagnostic on the A grid
c
      call recalc_agrid_uv ! add option for tendency computation?
      DO J=J_0S,J_1S
      DO I=I_0,I_1
      DO L=1,LM
        call inc_ajl(i,j,l,JL_DAMDC,
     &       (u_3d_agrid(L,I,J)-uasv(L,I,J))*PDSIG(L,I,J))
      END DO
      END DO
      END DO

      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          u1_after_aturb(i,j) = u_3d_agrid(1,i,j)
          v1_after_aturb(i,j) = v_3d_agrid(1,i,j)
        END DO
      END DO

      return
      end subroutine atm_diffus

      subroutine getdz(tv,dz,dze,rho,rhoe,tvsurf,dz0,im,jm,lm)
!@sum  getdz computes the 3-d finite difference dz and dze
!@+    as well as the 3-d density rho and rhoe
!@+    called at the primary grid (A-grid)
!@auth Ye Cheng/G. Hartke
!@var  tv virtual potential temp. referenced at 1 mb
!@var  dz0 z(1)-ze(1)
!@var  dz main grid spacing
!@var  dze edge grid spacing
!@var  rho,rhoe air density at the main/edge grids
!@var  tvsurf surface virtual temperature
!@var  im,jm,lm 3-d grids

      !
      !     Grids:
      !
      !                   -------------------------  lm+1
      !                lm  - - - - - - - - - - - - -
      !                   -------------------------  lm
      !                l+1 - - - - - - - - - - - - -
      !                    -------------------------  l+1
      !     (main)     l   - - - - - - - - - - - - -
      !                    -------------------------  l     (edge)
      !                l-1 - - - - - - - - - - - - -
      !                    -------------------------  l-1
      !                2   - - - - - - - - - - - - -
      !                    -------------------------    2
      !                1   - - - - - - - - - - - - -
      !                    -------------------------    1
      !
      !           dz0(i,j) z(1,i,j) - ze(1,i,j)
      !           dz(l,i,j) z(l+1,i,j) - z(l,i,j)
      !           dze(l,i,j) ze(l+1,i,j) - ze(l,i,j)
      !           rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
      !           rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dze(l,i,j))
      !
      !     at main: u,v,tv,q,ke
      !     at edge: e,lscale,km,kh,gm,gh
      !
      USE CONSTANT, only : grav,rgas
      USE GEOM, only : imaxj
      USE ATM_COM, only : pmid,pk,pedn
      USE DOMAIN_DECOMP_ATM, ONLY : grid

      implicit none

      integer, intent(in) :: im,jm,lm
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo),
     &        intent(in) :: tvsurf
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo),
     &        intent(out) :: dz0
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo),
     &        intent(in) :: tv
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo),
     &        intent(out) :: rho,rhoe,dz,dze

      real*8 :: temp0,temp1,temp1e,pl1,pl,pl1e,ple,plm1e
      integer :: i,j,l  !@var i,j,l loop variable

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S

C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP

      !@ temp0 virtual temperature (K) at (i,j) and SIG(l)
      !@ temp1 virtual temperature (K) at (i,j) and SIG(l+1)
      !@ temp1e average of temp0 and temp1
      do j=J_0, J_1
        do i=I_0,imaxj(j)
          do l=1,lm-1

            pl1 =pmid(l+1,i,j)
            pl  =pmid(l,i,j)
            pl1e=pedn(l+1,i,j)
            ple =pedn(l,i,j)
            temp0 =tv(l,i,j)*pk(l,i,j)
            temp1 =tv(l+1,i,j)*pk(l+1,i,j)
            temp1e=0.5d0*(temp0+temp1)
            dz(l,i,j)    =-(rgas/grav)*temp1e*log(pl1/pl)
            dze(l,i,j)=-(rgas/grav)*temp0 *log(pl1e/ple)
            rhoe(l+1,i,j)=100.d0*(pl-pl1)/(grav*dz(l,i,j))
            rho(l,i,j)=100.d0*(ple-pl1e)/(grav*dze(l,i,j))
            if(l.eq.1) then
              dz0(i,j)=-(rgas/grav)*.5d0*(temp0+tvsurf(i,j))
     &                 *log(pl/ple)
              rhoe(1,i,j)=100.d0*ple/(tvsurf(i,j)*rgas)
            endif
            if(l.eq.lm-1) then
              plm1e=pedn(lm+1,i,j)
              dze(lm,i,j)=-(rgas/grav)*temp1 *log(plm1e/pl1e)
              dz(lm,i,j)=0.
              rho(lm,i,j)=100.d0*(pl1e-plm1e)/(grav*dze(lm,i,j))
            endif

          end do
        end do
      end do

      return
      end subroutine getdz

      subroutine dout(z,ze,dz,u,v,t,q,ke,dtdz,dqdz,an2,as2
     &      ,wt_nl,wq_nl,kh,km,e,lscale
     &      ,uflx,vflx,tvflx,qflx,dbl,ldbl,i,j,n)
!@sum dout writes out diagnostics at (i,j)
!@auth  Ye Cheng/G. Hartke
!@var p  pressure at main grid z
!@var pe  pressure at secondary grid ze
!@var u  west-east   velocity component
!@var v  south-north velocity component
!@var t  virt. pot. temperature (referenced to 1mb)
!@var q  specific humidity
!@var e  turbulent kinetic energy
!@var ze hight at edge lyer
!@var dz(j) z(j+1)-z(j)
!@var dze(j) ze(j+1)-ze(j)
!@var as2 dudz^2+dvdz^2
!@var dtdz z-derivative of t at edge grid ze
!@var dqdz z-derivative of q at edge grid ze
!@var an2 g_alpha*dtdz, brunt-vasala frequency
!@var km turbulent viscosity for u and v equations
!@var kh turbulent diffusivity for t,q
!@var ke turbulent diffusivity for e
!@var lscale  turbulent length scale
!@var uflx momentun flux -uw at surface, ze(1)
!@var vflx momentun flux -vw at surface, ze(1)
!@var tvflx heat flux -wt at surface, ze(1)
!@var qflx moisture flux -wq at surface, ze(1)
!@var i/j horizontal location at which the output is written
!@var n number of vertical main layers

      USE ATM_COM, only : pmid,pedn,pk,pek
      USE SOCPBL, only : rimax

      implicit none

      integer, intent(in) :: ldbl,i,j,n
      real*8, dimension(n), intent(in) ::
     &   dz,u,v,t,q,ke,dtdz,dqdz,an2,as2
     &   ,wt_nl,wq_nl,kh,km,e,lscale,z
      real*8, dimension(n+1), intent(in) :: ze
      real*8, intent(in) :: uflx,vflx,tvflx,qflx,dbl

      real*8, dimension(n) :: p,pe
      real*8 :: ri,wt_lcl,wq_lcl
      integer :: l  !@var l loop variable

      Write (67,1100) "i=",i,"j=",j
      Write (67,1200) "uflx=",uflx,"vflx=",vflx
      Write (67,1200) "tvflx=",tvflx*pek(1,i,j),"qflx=",qflx
      Write (67,1250) "ldbl=",ldbl,"dbl=",dbl,"rimax=",rimax

      ! pressure at main and edge layers

      do l=1,n
          p(l)=100.*pmid(l,i,j)
          pe(l)=100.*pedn(l,i,j)
      end do

      ! fields at layer center

      write (67,1300)
      do l=1,n
        write (67,2000) l,p(l),z(l),dz(l),u(l),v(l),t(l)*pk(l,i,j),q(l)
     &                  ,ke(l)
      end do
      write (67,*)

      ! fields at layer edge

      write (67,1400)
      do l=1,n
        wt_lcl=-kh(l)*dtdz(l)*pek(l,i,j)
        wq_lcl=-kh(l)*dqdz(l)
        ri=an2(l)/as2(l)
        write (67,2000) l,pe(l),ze(l),wt_lcl,wt_nl(l)*pek(l,i,j)
     &    ,wq_lcl,wq_nl(l),kh(l),km(l),lscale(l),ri,e(l)
      end do
      write (67,1500)

      return
1100  format(3(2x,a,i6))
1200  format(2(2x,a,1pe14.4))
1250  format(2x,a,i6,2x,a,1pe14.4,2x,a,1pe14.4)
1300  format (' ',' l',1x,
     1            '     p     ',1x,
     2            '     z     ',1x,'     dz    ',1x,'     u     ',1x,
     3            '     v     ',1x,'     t     ',1x,'     q     ',1x,
     4            '     ke    ',1x,'           ',1x,'           ',1x,
     5            '           ',/)
        write (67,2000) l,pe(l),ze(l),wt_lcl,wt_nl(l),wq_lcl,wq_nl(l)
     2    ,kh(l),km(l),lscale(l),ri,e(l)
1400  format (' ',' l',1x,
     2            '  p (edge) ',1x,
     2            '  z (edge) ',1x,'  wt_lcl   ',1x,'   wt_nl   ',1x,
     3            '   wq_lcl  ',1x,'  wq_nl    ',1x,'     kh    ',1x,
     4            '     km    ',1x,'  lscale   ',1x,'     ri    ',1x,
     5            '     e     ',/)
1500  format(132('-'))
2000  format (1x,i2,1x,1pe11.4,9(1x,1pe11.4),1x,1pe10.3)
      end subroutine dout

      subroutine de_solver_main(x,x0,p1,p4,
     &    rhoebydz,bydzerho,flux_bot,flux_top,dtime,n,qlimit)
!@sum differential equation solver using tridiagonal method.
!@+  The differential equation is expressed as
!@+  d/dt x = d/dz (P1 d/dz x) + P4
!@+  where x is the unknown to be solved, 
!@+  x and P4 are at the layer middle z, while
!@+  P1 is at the layer edge ze.
!@auth  Ye Cheng/G. Hartke
!@var x the unknown to be solved (at main drid)
!@var x0 x at previous time step
!@var p1,p4 coeff. of the d.e.
!@var rhoebydz(j) rhoe(j+1)/dz(j)
!@var bydzerho(j) 1/(dze(j)*rho(j))
!@var flux_bot flux from the bottom
!@var flux_top flux at the top
!@var dtime time step
!@var n number of main grid
!@var qlimit true if tracer must be positive definite
      use TRIDIAG_MOD, only : TRIDIAG
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) ::
     &        x0,p1,p4,rhoebydz,bydzerho
      real*8, intent(in) :: flux_bot,flux_top,dtime
      real*8, dimension(n), intent(out) :: x

      real*8, dimension(n) :: sub,dia,sup,rhs
      real*8 :: alpha
      integer :: j  !@var j loop variable
      logical, intent(in) :: qlimit

      !sub(j)*x_jm1_kp1+dia(j)*x_j_kp1+sup(j)*x_jp1_kp1 = rhs(j)
      !k refers to time step, j refers to main grid
      !except p1(j) which is defined on the edge grid

      do j=2,n-1
          sub(j)=-dtime*p1(j)*rhoebydz(j)*bydzerho(j)
          sup(j)=-dtime*p1(j+1)*rhoebydz(j+1)*bydzerho(j)
          dia(j)=1.d0-(sub(j)+sup(j))
          rhs(j)=x0(j)+dtime*p4(j)
          if (qlimit.and. rhs(j).lt.0) rhs(j)=0.  ! prevent roundoff error
      end do

      ! Lower boundary conditions(x=T) :
      ! at main grid 1
      ! (T(1)-T0(1))/dtime = -(1/rho)*d/dz(rho*wt)
      ! -d/dz(rho*wt)=-(rhoe(2)*wt(2)-rhoe(1)*wt(1))/dze(1)
      ! wt(2)=-p1(2)*(T(2)-T(1))/dz(1)+wt_nl(2)
      ! the above together yield
      ! (1+alpha)*T(1)-alpha*T(2)=T0(1)-dtime*bydzerho(1)*flux_bot
      ! where flux_bot=-rhoe(1)*wt(1)+rhoe(2)*wt_nl(2)
      ! and specify wt(1)=-tvflx

      alpha=dtime*p1(2)*rhoebydz(2)*bydzerho(1)
      dia(1)=1.d0+alpha
      sup(1)=-alpha
      rhs(1)=x0(1)-dtime*bydzerho(1)*flux_bot

      ! Upper boundary conditions(x=T) :
      ! at main grid n(=lm)
      ! (T(n)-T0(n))/dtime = -(1/rho)*d/dz(rho*wt)
      ! -d/dz(rho*wt)=-(rhoe(n+1)*wt(n+1)-rhoe(n)*wt(n))/dze(n)
      ! wt(n)=-p1(n)*(T(n)-T(n-1))/dz(n-1)+wt_nl(n)
      ! the above together yield
      ! -alpha*T(n-1)+(1+alpha)*T(n)=T0(n)+dtime*bydzerho(n)*flux_top
      ! where flux_top=-rhoe(n+1)*wt(n+1)+rhoe(n)*wt_nl(n)
      ! and specify wt(n+1)=0

      alpha=dtime*p1(n)*rhoebydz(n)*bydzerho(n)
      sub(n)=-alpha
      dia(n)=1.d0+alpha
      rhs(n)=x0(n)+dtime*bydzerho(n)*flux_top

      call tridiag(sub,dia,sup,rhs,x,n)

      return
      end subroutine de_solver_main

      subroutine de_solver_edge(x,x0,p1,p3,p4,
     &    rhobydze,bydzrhoe,x_surf,dtime,n)
!@sum differential equation solver using tridiagonal method.
!@+  The differential equation is expressed as
!@+  d/dt x = d/dz (P1 d/dz x) - P3 x + P4
!@+  where x is the unknown to be solved,
!@+  x, P3 and P4 are at the layer edge ze, while
!@+  P1 is at the layer middle z.
!@auth  Ye Cheng/G. Hartke
!@var x the unknown to be solved (at edge drid)
!@var x0 x at previous time step
!@var p1,p3,p4 coeff. of the d.e.
!@var rhobydze(j) rho(j)/dze(j)
!@var bydzrhoe(j) 1/(dz(j-1)*rhoe(j))
!@var x_surf surface value of x
!@var dtime time step
!@var n number of vertical edge grid

      use TRIDIAG_MOD, only : TRIDIAG
      use CONSTANT, only : teeny
      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) ::
     &        x0,p1,p3,p4,rhobydze,bydzrhoe
      real*8, intent(in) :: x_surf,dtime
      real*8, dimension(n), intent(out) :: x

      real*8, dimension(n) :: sub,dia,sup,rhs
      integer :: j  !@var j loop variable

      ! sub(j)*x_jm1_kp1+dia(j)*x_j_kp1+sup(j)*x_jp1_kp1 = rhs(j)
      ! k refers to time step, j refers to edge grid
      ! except p1(j) which is defined on the main grid

      do j=2,n-1
          sub(j)=-dtime*p1(j-1)*rhobydze(j-1)*bydzrhoe(j)
          sup(j)=-dtime*p1(j)*rhobydze(j)*bydzrhoe(j)
          dia(j)=1.d0-(sub(j)+sup(j))+dtime*p3(j)
          rhs(j)=x0(j)+dtime*p4(j)
      end do

      ! Boundary conditions:

      dia(1)=1.d0
      sup(1)=0.d0
      rhs(1)=x_surf

      sub(n)=-1.d0
      dia(n)=1.d0
      rhs(n)=0.d0

      call tridiag(sub,dia,sup,rhs,x,n)
      ! assuming x is a non-negative scalar:
      do j=1,n
         if(x(j).lt.teeny) x(j)=teeny
      end do

      return
      end subroutine de_solver_edge

      subroutine apply_fluxes_to_atm
!@sum a dummy subroutine that replaces the real one needed by DRYCNV.
!@auth I. Aleinov
      return
      end subroutine apply_fluxes_to_atm

      subroutine zze(dz,dze,dz0,z,ze,n)
!@sum finds the layer middle and edge heights, z and ze
!@+  Note that z(L) is between ze(L) and ze(L+1).
!@auth  Ye Cheng/G. Hartke
!@var  z vertical coordinate of mid points
!@var  ze vertical coordinate at edges
!@var  dz(l) z(l+1) - z(l)
!@var  dze(l) ze(l+1) - ze(l)
!@var  dz0 z(1) - ze(1)
!@var  n number of layers

      USE SOCPBL, only : zgs

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n), intent(in) :: dz,dze
      real*8, intent(in) :: dz0
      real*8, dimension(n), intent(out) :: z
      real*8, dimension(n+1), intent(out) :: ze
      integer :: j  !@var j loop variable

      ze(1)=zgs
      z(1)=ze(1)+dz0
      do j=2,n
          ze(j)=ze(j-1)+dze(j-1)
          z(j)=z(j-1)+dz(j-1)
      end do
      ze(n+1)=ze(n)+dze(n)

      return
      end subroutine zze

      subroutine l_gcm(ze,dbl,lmonin,ustar,qturb,an2,lscale,n)
!@sum calculates the turbulent length scale
!@+  (lscale, in meters). Within the PBL, it is according
!@+  to Nakanishi(2001); above the PBL, we generalized and 
!@+  employed a formula by Holtslag and Boville (1993).
!@auth Ye Cheng/G. Hartke
!@var ze height (meters) of layer edge
!@var dbl pbl depth (meters)
!@var lscale turbulent length scale
!@var n number of layers

      USE CONSTANT, only : teeny,by3
      USE SOCPBL, only : kappa

      implicit none

      integer, intent(in) :: n
      real*8, dimension(n+1), intent(in) :: ze
      real*8, dimension(n), intent(in) :: qturb,an2
      real*8, intent(in) :: dbl,lmonin,ustar
      real*8, dimension(n), intent(out) :: lscale

      real*8, parameter :: fac=0.3d0,l0min=30.d0
      real*8 :: z,zeta,kz,l0,ls,lb,l1,an,qty
      integer :: j  !@var j loop variable

      do j=1,n
         z=ze(j)
         kz=kappa*z
         if(z.lt.dbl) then         ! within pbl
            zeta=z/lmonin
            l0=.3d0*dbl
            if(zeta.ge.1.) then
               ls=kz/3.7d0
            elseif(zeta.ge.0.) then
               ls=kz/(1.+2.7d0*zeta)
            else
c              ls=kz*(1.-100.*zeta)**0.2d0
               ls=kz
            endif
            if (an2(j).gt.0.) then
               an=sqrt(an2(j))
               if(zeta.ge.0.) then
                  lb=qturb(j)/an
               else
                  qty=(ustar/((-kappa*lmonin*l0*l0)**by3*an))**0.5d0
                  lb=qturb(j)*(1.+5.*qty)/an
               endif
            else
               lb=1.d30
            endif
            lscale(j)=l0*ls*lb/(l0*ls+l0*lb+ls*lb)
        else                     ! above pbl
            l1=l0min+max(fac*dbl-l0min,0.d0)*exp(1.-z/dbl)
            lscale(j)=l1*kz/(l1+kz)
        endif
      end do

      return
      end subroutine l_gcm

      subroutine k_gcm(tvflx,qflx,uflx,vflx,ustar,wstar,dbl,lmonin
     &  ,ze,lscale,e,qturb,an2,as2,dtdz,dqdz,dudz,dvdz
     &  ,kh,kq,km,ke,wt,wq,w2,uw,vw,wt_nl,wq_nl,uw_nl,vw_nl
#ifdef TRACERS_ON
     &  ,trflx,wc_nl,nta
#endif
     &  ,n)
!@sum computes the turbulent stability functions Km (for momentum)
!@+  and Kh (for heat and moisture), as well as the fluxes
!@+  (local and non-local).
!@+  Within the convective PBL, it is according to
!@+  Holtslag and Boville (1993); within the stable PBL or above
!@+  the PBL, it is according to Cheng et al. (2002).
!@auth  Ye Cheng/G. Hartke
!@var tvflx virtual potential temperature flux at surface
!@var qflx moisture flux at surface
!@var ustar friction velocity
!@var wstar Deardorff convective velocity scale
!@var dbl the height of the PBL (real*8, in meters)
!@var kh turbulent diffusivity for scalars (heat, moisture,...)
!@var ze height (in meters) of layer edges
!@var lscale turbulent length scale
!@var e turbulent kinetic energy
!@var qturb sqrt(2*e)
!@var an2 g_alpha*dTdz
!@var as2 (dUdz)**2+(dVdz)**2
!@var dtdz dt/dz
!@var dqdz dq/dz
!@var dudz du/dz
!@var dvdz dv/dz
!@var km turbulent diffusivity for momentun
!@var ke turbulent diffusivity for e
!@var w2 vertical component of 2*e
!@var wt,wq,uw,vw turbulent fluxes
!@var wt_nl non-local part of heat flux wt
!@var wq_nl non-local part of moisture flux wq
!@var wc_nl non-local part of tracer flux
!@var n number of layers

      USE CONSTANT, only : teeny,by3,sha
      USE SOCPBL, only : kappa,prt,ghmin,ghmax,d1,d2,d3,d4,d5
     &   ,s0,s1,s2,s4,s5,s6,s7,s8,b1,g5
     &   ,k_max,kmmin,khmin,find_phim0,find_phih

      implicit none

      integer, intent(in) :: n
#ifdef TRACERS_ON
      integer, intent(in) :: nta
      real*8, dimension(n,nta),intent(out) :: wc_nl
      real*8, dimension(nta),intent(in) :: trflx
      real*8, dimension(nta) :: cgtr
      integer nt
#endif

      real*8, intent(in) :: tvflx,qflx,uflx,vflx,ustar,wstar,dbl,lmonin
      real*8, dimension(n), intent(in) :: lscale,e,qturb,an2,as2
     &        ,dtdz,dqdz,dudz,dvdz
      real*8, dimension(n+1), intent(in) :: ze
      real*8, dimension(n), intent(out) ::
     &  kh,kq,km,ke,wt,wq,w2,uw,vw,wt_nl,wq_nl,uw_nl,vw_nl

      real*8 :: tmp,tau,gm,gh,gmmax,byden,sm,sh
     &    ,ustar2,wstar3,zzi,tau_pt,w2j,phih1,by_phim1,wm1,pr1
     &    ,cgh1,km_n,kh_n,pr,cgh,cgq,zet,phih,phim,by_phim,wm,kz
     &    ,cgu1,cgv1
      integer :: j  !@var j loop variable

      !@ Non-local model: Holtslag and Boville, 1993.
      !@ Local model: Cheng et al. 2002.
      ustar2=ustar*ustar
      wstar3=wstar*wstar*wstar

      !@ some quantities independent of z:
      zet=.1d0*dbl/lmonin
      if(zet.lt.0.d0) then
        call find_phim0(zet,phim)
        call find_phih(zet,phih1)
        by_phim1=1./phim
        wm1=ustar*by_phim1
        pr1=phih1*by_phim1+.72d0*kappa*wstar/wm1
        cgh1=7.2d0*wstar*(-tvflx)/(wm1**2*dbl)
c       cgu1=7.2d0*wstar*(-uflx)/(wm1**2*dbl)
c       cgv1=7.2d0*wstar*(-vflx)/(wm1**2*dbl)
        ! counter-gradient uw,vw turned off
        cgu1=0.;cgv1=0.
#ifdef TRACERS_ON
        do nt=1,nta
          cgtr(nt)=7.2*wstar*max(-trflx(nt),0.d0)/(wm1**2*dbl)
        end do
#endif
      else
        phih1=0.;by_phim1=0.;wm1=0.;pr1=0.;cgh1=0.
        cgu1=0.;cgv1=0.
#ifdef TRACERS_ON
        do nt=1,nta
          cgtr(nt)=0.
        end do
#endif
      endif
      do j=1,n
          kz=kappa*ze(j)
          tau=b1*lscale(j)/(qturb(j)+teeny)
          gh=tau*tau*an2(j)
          gm=tau*tau*as2(j)
          if(gh.lt.ghmin) gh=ghmin
          if(gh.gt.ghmax) gh=ghmax
          gmmax=(1.+d1*gh+d3*gh*gh)/(d2+d4*gh)
          if(gm.gt.gmmax) gm=gmmax
          byden=1./(1.+d1*gh+d2*gm+d3*gh*gh+d4*gh*gm+d5*gm*gm)
          sm=(s0+s1*gh+s2*gm)*byden
          sh=(s4+s5*gh+s6*gm)*byden
          km(j)=min(max(tau*e(j)*sm,kmmin),k_max)
          kh(j)=min(max(tau*e(j)*sh,khmin),k_max)
          kq(j)=kh(j)
          wt_nl(j)=0.
          wq_nl(j)=0.
          uw_nl(j)=0.
          vw_nl(j)=0.
#ifdef TRACERS_ON
          do nt=1,nta
            wc_nl(j,nt)=0.
          end do
#endif
          zzi=ze(j)/dbl
          zet=ze(j)/lmonin
          if((zzi.le.1.d0).and.(zet.lt.0.d0)) then ! within unstable pbl
             if(zzi.lt.(.1d0)) then     !!! in surface layer
               call find_phim0(zet,phim)
               call find_phih(zet,phih)
               by_phim=1./phim
               wm=ustar*by_phim
               km_n=kz*wm*(1.-zzi)**2
               pr=phih*by_phim
               kh_n=km_n/pr
             else                     !!! in outer layer of pbl
               km_n=kz*wm1*(1.-zzi)**2
               kh_n=km_n/pr1
               wt_nl(j)=kh_n*cgh1
               uw_nl(j)=km_n*cgu1
               vw_nl(j)=km_n*cgv1
!#ifdef TRACERS_ON
!               do nt=1,nta
!                 wc_nl(j,nt)=kh_n*cgtr(nt)
!               end do
!#endif
             endif
             tmp=(1.6d0*ustar2*(1.-zzi)+teeny)**1.5d0
     &           +1.2d0*wstar3*zzi*(1.-.9d0*zzi)**1.5d0
             w2j=tmp**(2.*by3)
             km(j)=min(max(km_n,kmmin),k_max)
             kh(j)=min(max(kh_n,khmin),k_max)
             kq(j)=kh(j)
          else                 ! above the pbl
             w2j=by3*(2.*e(j)-tau*(s7*km(j)*as2(j)+s8*kh(j)*an2(j)))
          endif

          ke(j)=5.*km(j)
          wt(j) = -kh(j)*dtdz(j)+wt_nl(j)
          wq(j) = -kq(j)*dqdz(j)
          w2(j) = min(max(0.24d0*e(j),w2j),2.*e(j))
          uw(j) = -km(j)*dudz(j)+uw_nl(j)
          vw(j) = -km(j)*dvdz(j)+vw_nl(j)
      end do

      return
      end subroutine k_gcm

      subroutine e_gcm(tvflx,wstar,ustar,dbl,lmonin,ze,g_alpha
     &    ,an2,as2,lscale,e,n)
!@sum finds the turbulent kinetic energy (e, in m^2/s^2). 
!@+  Within the PBL, e is determined according to the
!@+  parameterization of the Large Eddy Simulation (LES) data
!@+  (Moeng and Sullivan, 1994), above the PBL, e is calculated
!@+  by the second order closure model of Cheng et al. (2002).
!@auth  Ye Cheng
!@var (see subroutine k_gcm)
!@var lmonin Monin-Obukov length
!@var g_alpha grav*alpha
      USE CONSTANT, only : teeny,by3
      USE SOCPBL, only : kappa,emax,rimax,b1,c1,c2,c3,c4,c5,gm_at_rimax
     &                  ,emin,find_phim0

      implicit none

      integer, intent(in) :: n   !@var n  array dimension
      real*8, intent(in) :: tvflx,wstar,ustar,dbl,lmonin
      real*8, dimension(n), intent(in)   :: g_alpha,an2,as2,lscale
      real*8, dimension(n+1), intent(in) :: ze
      real*8, dimension(n), intent(out) :: e
      integer :: j !@var j loop variable
      real*8 :: ri,gm,aa,bb,cc,phim,tmp,wstar3,ustar3,zj,zet,eps,ej
     &         ,kz

      ustar3=ustar*ustar*ustar
      wstar3=wstar*wstar*wstar
      do j=1,n
        zj=ze(j)
        kz=kappa*zj
        if(zj.le.dbl) then !Hogstrom 1988,1996
          zet=zj/lmonin
          call find_phim0(zet,phim)
          eps=.4d0*wstar3/dbl+ustar3*(1.-zj/dbl)*phim/kz
          ej=.5d0*(19.3d0*lscale(j)*eps)**(2.*by3)
          e(j)=min(max(ej,emin),emax)
        else
          ri=an2(j)/max(as2(j),teeny)
          if(ri.lt.rimax) then
            aa=c1*ri*ri-c2*ri+c3
            bb=c4*ri+c5
            cc=2.d0
            if(abs(aa).lt.1d-8) then
              gm= -cc/bb
            else
              tmp=bb*bb-4.*aa*cc
              gm=(-bb-sqrt(tmp))/(2.*aa)
            endif
          else
            ri=rimax
            gm=gm_at_rimax
          endif
          tmp=0.5d0*(b1*lscale(j))**2*as2(j)/max(gm,teeny)
          e(j)=min(max(tmp,emin),emax)
        endif
      end do
      return
      end subroutine e_gcm

      subroutine find_pbl_top(z,u,v,t,ustar,ustar2,tvflx,lmonin
     &   ,dbl,ldbl,ldbl_max,n)
!@sum  finds the PBL height (dbl, in meters)
!@+  and the main level index immediately above (ldbl),
!@+  using the bulk Richardson number criterion
!@+  (Holtslag and Boville, 1993).
!@auth Ye Cheng
!@var z main layer height (in meters)
!@var dbl the pbl height (in meters)
!@var ldbl the main layer immediately above the pbl height dbl
!@var tvflx minus virtual heat flux at the surface
!@var ldbl_max the maximum allowable number of layers in the pbl
!@var n total number of layers
c***  this dbl is different from the dbl in module socpbl,
c***  the latter is itype dependent

      USE SOCPBL, only : find_phim0
      USE CONSTANT, only : grav,teeny,by3

      implicit none

      integer, intent(in) :: ldbl_max,n
      real*8, dimension(n), intent(in) :: z,u,v,t
      real*8, intent(in) :: ustar,ustar2,tvflx,lmonin

      real*8, intent(out) :: dbl
      integer, intent(out) :: ldbl

      real*8, parameter :: fac=100.,ri_cr=0.50d0,b=8.5d0
      REAL*8, parameter :: dbl_max=4000.d0 ! meters

      real*8, dimension(n) :: ri
      real*8 :: v2l,wtvs,wm,t1_w_excess,den,dbls,zet,phim
      integer :: l

      dbl=z(1)
      ri(1)=0.
      do l=2,ldbl_max
        v2l=max((u(l)-u(1))**2+(v(l)-v(1))**2+fac*ustar2,teeny)
        ri(l)=(z(l)-z(1))*grav*(t(l)-t(1))/(t(1)*v2l)
        if(ri(l).ge.ri_cr) then
          den=ri(l)-ri(l-1)
          if(den.eq.0.) den=teeny
          dbl=z(l-1)+(z(l)-z(l-1))*(ri_cr-ri(l-1))/den
          exit
        endif
        if(l.eq.ldbl_max) dbl=z(l)
      end do
      ! tvflx = - <w*tv> at surface
      wtvs=-tvflx
      if(wtvs.gt.0.) then
        zet=.1d0*dbl/lmonin
        call find_phim0(zet,phim)
        wm=ustar/phim
        t1_w_excess=t(1)+b*wtvs/(wm+teeny)
        do l=2,ldbl_max
          v2l=max((u(l)-u(1))**2+(v(l)-v(1))**2+fac*ustar2,teeny)
          ri(l)=(z(l)-z(1))*grav*(t(l)-t1_w_excess)/(t(1)*v2l)
          if(ri(l).ge.ri_cr) then
            den=ri(l)-ri(l-1)+teeny
            dbl=z(l-1)+(z(l)-z(l-1))*(ri_cr-ri(l-1))/den
            exit
          endif
          if(l.eq.ldbl_max) dbl=z(l)
        end do
      endif
      dbl=min(dbl,dbl_max)

      !@var ldbl the level immediately above the pbl height 
      if (dbl.le.z(1)) then
        ldbl=1
      else
        do l=2,ldbl_max
          if(dbl.gt.z(l-1).and.dbl.le.z(l)) then
            ldbl=l
            exit
          endif
          if(l.eq.ldbl_max) ldbl=ldbl_max
        end do
      endif
      
      return
      end subroutine find_pbl_top

