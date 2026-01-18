#include "rundeck_opts.h"

#ifdef TRACERS_ON

      subroutine diag_trac_prep
      implicit none
      call gather_zonal_trdiag
      call diagjlt_prep
      call diagijt_prep
      call diagijlt_prep
      return
      end subroutine diag_trac_prep

      SUBROUTINE DIAGJLT_prep
! comments to be added
      use OldTracer_mod, only:  ntm_power
      use OldTracer_mod, only:  trname, src_dist_index
      use constant, only: teeny, grav
      use resolution, only: lm
      USE MODEL_COM, only: idacc
      USE TRACER_COM, only: ntm, n_Water, n_CH4, n_O3
#ifdef TRACERS_WATER
      use OldTracer_mod, only:trw0, dowetdep
#endif
#ifdef TRACERS_SPECIAL_O18
      USE TRACER_COM, only: n_H2O18,n_HDO,n_H2O17
#endif
      USE GEOM, only : areag
      USE DIAG_COM, only: jm=>jm_budg,dxyp=>dxyp_budg,
     &     ia_dga,ia_src,ajl,cdl_jl_template
     &     ,jl_dpa,jl_dpasrc,jl_dwasrc
#ifndef CUBED_SPHERE
     &     ,fim
#endif
      USE MDIAG_COM, only : make_timeaxis
      USE TRDIAG_COM, only : tajln, tajls, lname_jln, sname_jln,
     *     units_jln,  scale_jln, lname_jls, sname_jls, units_jls,
     *     scale_jls, jls_power, jls_ltop, ia_jls, jwt_jls, jgrid_jls,
     *     jls_3Dsource, jlnt_conc, jlnt_mass, jlnt_nt_tot, jlnt_nt_mm,
     *     jlnt_lscond,  jlnt_turb,  jlnt_vt_tot, jlnt_vt_mm, jlnt_mc,
     *     jgrid_jlq, ia_jlq, scale_jlq, jlq_power, ktajls, jls_source
#ifdef TRACERS_WATER
     *     ,jlnt_cldh2o
#endif
     &     ,tajl=>tajl_out, ktajl_, ktajl_out, cdl_tajl, hemis_tajl
     &     ,vmean_tajl, denom_tajl,ia_tajl,sname_tajl,lname_tajl
     &     ,ltop_tajl,units_tajl,scale_tajl,pow_tajl
     &     ,jgrid_tajl,lgrid_tajl
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      USE TRDIAG_COM, only : to_per_mil
#endif
#ifndef CUBED_SPHERE
      USE BDJLT
#endif
      use domain_decomp_atm, only : am_i_root
      use cdl_mod
      IMPLICIT NONE
      real*8 :: bydxyp(jm),byapo(jm),onespo(jm),fj(jm)
      INTEGER :: J,L,N,K,KK,KKK,jtpow,n1,n2,k_dpa,k_dwa,k_vap,k_cnd,
     &     j1,j2
      REAL*8 :: dD, d18O, d17O, byiacc, hemfac
      real*8, dimension(:,:,:), allocatable :: tajl_tmp
      character(len=10) :: zstr
      character(len=3) :: ltopstr,powstr
      logical, dimension(ktajl_) :: per_area,output_vsum
      logical :: set_miss

      if(.not. am_i_root()) return

      onespo = 1.
      bydxyp = 1d0/dxyp
      byapo = bydxyp

#ifndef CUBED_SPHERE
      call JLt_TITLEX ! needed for some extra titles
      onespo(1)  = fim
      onespo(jm) = fim
      byapo = bydxyp*onespo/fim
#endif


      do k=1,ktajl_
        denom_tajl(k) = 0
        ia_tajl(k) = ia_src
        sname_tajl(k) = 'unused'
        lname_tajl(k) = 'unused'
        units_tajl(k) = 'unused'
        scale_tajl(k) = 1.
        pow_tajl(k) = 0
        jgrid_tajl(k) = 1
        lgrid_tajl(k) = 1
        ltop_tajl(k) = lm
        per_area(k) = .true.
        output_vsum(k) = .false. ! whether to output vertical sums
      enddo

      k = 0

c
      k = k + 1
      k_dpa = k
      do l=1,lm
        tajl(:,l,k) = ajl(:,l,jl_dpasrc)*bydxyp(:)
      enddo
c
      k = k + 1
      k_dwa = k
      do l=1,lm
        tajl(:,l,k) = ajl(:,l,jl_dwasrc)*bydxyp(:)
      enddo

#ifndef TRACERS_AMP
#ifdef TRACERS_WATER
c
      k = k + 1
      k_vap = k
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,jlnt_mass,n_water)*bydxyp(:)
      enddo
c
      k = k + 1
      k_cnd = k
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,jlnt_cldh2o,n_water)*bydxyp(:)
      enddo

#endif
#endif

C****
C**** LOOP OVER TRACERS
C****

      DO N=1,NTM
      if (src_dist_index(n)>1) cycle
C****
C**** TRACER CONCENTRATION
C****
      k = k + 1
      kk = jlnt_conc
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,jlnt_mass,n)*bydxyp(:)
      enddo

#ifdef TRACERS_WATER
      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
        denom_tajl(k) = k_vap
        tajl(:,:,k) = 1d3*(tajl(:,:,k)/trw0(n)-tajl(:,:,k_vap))
      else
#endif

      denom_tajl(k) = k_dpa
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jln(n)*scale_jlq(kk)*10.**(-jtpow)

#ifdef TRACERS_WATER
      end if
#endif

#ifdef TRACERS_COSMO
      if(n.eq.n_Be7) k_Be7 = k
      if(n.eq.n_Pb210) k_Pb210 = k
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_SPECIAL_Shindell) ||\
    (defined TRACERS_AEROSOLS_SEASALT)
C****
C**** Mass diagnostic (this is saved for everyone, but only output
C**** for Dorothy and Drew for the time being)
C****
      k = k + 1
      kk=jlnt_mass
      per_area(k) = .false.
      output_vsum(k) = .true.
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      tajl(:,:,k) = tajln(:,:,kk,n)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AEROSOLS_SEASALT)
      jtpow = jtpow+13
#else
      denom_tajl(k) = k_dpa
      byiacc = 1d0/(idacc(ia_tajl(k))+teeny)
      do l=1,lm
        tajl(:,l,k) = tajl(:,l,k)*byapo(:)*byiacc*tajl(:,l,k_dpa)
      enddo
#endif
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
#endif

#ifdef TRACERS_WATER
C****
C**** TRACER CLOUD WATER CONCENTRATION
C****
      if (dowetdep(n)) then
      k = k + 1
      kk = jlnt_cldh2o

      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      do l=1,lm
        tajl(:,l,k) = tajln(:,l,kk,n)*bydxyp(:)
      enddo

      if (to_per_mil(n).gt.0) then
C**** Note permil concentrations REQUIRE trw0 and n_water to be defined!
        denom_tajl(k) = k_cnd
        tajl(:,:,k) = 1d3*(tajl(:,:,k)/trw0(n)-tajl(:,:,k_cnd))
      else
        denom_tajl(k) = k_dwa
        jtpow = ntm_power(n)+jlq_power(kk)
        scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      end if
      end if
#endif
C****
C**** NORTHWARD TRANSPORTS: Total and eddies
C****
#ifndef CUBED_SPHERE
      k = k + 1
      kk = jlnt_nt_tot
      per_area(k) = .false.
      output_vsum(k) = .true.
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      jgrid_tajl(k) = 2
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(2:jm,:,k) = tajln(1:jm-1,:,kk,n)
c
      k = k + 1
      kk = jlnt_nt_eddy
      per_area(k) = .false.
      output_vsum(k) = .true.
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      jgrid_tajl(k) = 2
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(2:jm,:,k) = tajln(1:jm-1,:,jlnt_nt_tot,n)
     &                -tajln(1:jm-1,:,jlnt_nt_mm ,n)
#endif
C****
C**** VERTICAL TRANSPORTS: Total and eddies
C****
#ifndef CUBED_SPHERE
      k = k + 1
      kk = jlnt_vt_tot
      per_area(k) = .false.
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      lgrid_tajl(k) = 2
      ltop_tajl(k) = lm-1
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(:,:,k) = tajln(:,:,kk,n)
c
      k = k + 1
      kk = jlnt_vt_eddy
      per_area(k) = .false.
      sname_tajl(k) = sname_jln(kk,n)
      lname_tajl(k) = lname_jln(kk,n)
      units_tajl(k) = units_jln(kk,n)
      lgrid_tajl(k) = 2
      ltop_tajl(k) = lm-1
      ia_tajl(k) = ia_jlq(kk)
      jtpow = ntm_power(n)+jlq_power(kk)
      scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
      tajl(:,:,k) = tajln(:,:,jlnt_vt_tot,n)-tajln(:,:,jlnt_vt_mm,n)
#endif
c
c tendencies from various processes
c
      do kkk=1,3
        k = k + 1
        if(kkk.eq.1) kk = jlnt_mc
        if(kkk.eq.2) kk = jlnt_lscond
        if(kkk.eq.3) kk = jlnt_turb
        per_area(k) = .false.
        output_vsum(k) = .true.
        sname_tajl(k) = sname_jln(kk,n)
        lname_tajl(k) = lname_jln(kk,n)
        units_tajl(k) = units_jln(kk,n)
        ia_tajl(k) = ia_jlq(kk)
        jtpow = ntm_power(n)+jlq_power(kk)
        scale_tajl(k) = scale_jlq(kk)*10.**(-jtpow)
        do l=1,lm
          tajl(:,l,k) = tajln(:,l,kk,n)*onespo(:)
        enddo
      enddo

      enddo ! end loop over tracers

C****
C**** JL Specials (incl. Sources and sinks)
C**** Partial move towards correct units (kg/(mb m^2 s)).
C**** Plot depends on jwt_jls.
C**** Note that only jwt_jls=3 is resolution independent.
C****
      do kk=1,ktajls
        if (sname_jls(kk).eq."daylight" .or. sname_jls(kk).eq."H2O_mr"
     *       .or. lname_jls(kk).eq."unused") cycle

        k = k+1
        sname_tajl(k) = sname_jls(kk)
        lname_tajl(k) = lname_jls(kk)
        units_tajl(k) = units_jls(kk)
        ia_tajl(k) = ia_jls(kk)
        ltop_tajl(k) = jls_ltop(kk)
        jtpow = jls_power(kk)
        scale_tajl(k) = scale_jls(kk)*10.**(-jtpow)
        tajl(:,:,k) = tajls(:,:,kk)
        byiacc = 1d0/(idacc(ia_tajl(k))+teeny)
        if(jwt_jls(kk).eq.1) then
          per_area(k) = .false.
          output_vsum(k) = .true.
          do l=1,lm
            tajl(:,l,k) = tajl(:,l,k)*onespo(:)
          enddo
        elseif(jwt_jls(kk).eq.2) then
          denom_tajl(k) = k_dpa
          do l=1,lm
            tajl(:,l,k) = tajl(:,l,k)*bydxyp(:)*byiacc*tajl(:,l,k_dpa)
          enddo
        elseif(jwt_jls(kk).eq.3) then
          denom_tajl(k) = k_dpa
          do l=1,lm
            tajl(:,l,k) = tajl(:,l,k)*bydxyp(:)*100./grav
          enddo
        endif

c        select case (jwt_jls(kk))
c        case (1)   !  simple sum (like kg/s),
c          CALL JLMAP_t (lname_jls(kk),sname_jls(kk),units_jls(kk),plm,
c     *         tajls(1,1,kk),scalet,onespo,ones,jls_ltop(kk),jwt_jls(kk)
c     *         ,jgrid_jls(kk))
c        case (2)   !  area weighting (like kg/m^2 s)
c          CALL JLMAP_t (lname_jls(kk),sname_jls(kk),units_jls(kk),plm
c     *         ,tajls(1,1,kk),scalet,bydxyp,ones,jls_ltop(kk),jwt_jls(kk)
c     *         ,jgrid_jls(kk))
c        case (3)   !  area + pressure weighting (like kg/mb m^2 s)
c          CALL JLMAP_t (lname_jls(kk),sname_jls(kk),units_jls(kk),plm
c     *         ,tajls(1,1,kk),scalet,byapo,ones,jls_ltop(kk),jwt_jls(kk)
c     *         ,jgrid_jls(kk))
c        end select

        end do

#ifdef TRACERS_SPECIAL_Lerner
C**** some special combination diagnostics

C**** total chemical change for CH4
      if (n_CH4.gt.0) then
        k = k + 1
        kk=jls_3Dsource(1,n_CH4)
        sname_tajl(k) = 'Total_Chem_change_'//trim(trname(n_CH4))
        lname_tajl(k) = 'TOTAL CHANGE OF '//trim(trname(n_CH4))//
     &                  ' BY CHEMISTRY'
        units_tajl(k) = units_jls(kk)
        ia_tajl(k) = ia_jls(kk)
        ltop_tajl(k) = jls_ltop(kk)
        scale_tajl(k) = scale_jls(kk)*10.**(-jls_power(kk))
        tajl(:,:,k) = tajls(:,:,jls_3Dsource(1,n_CH4))
     &              + tajls(:,:,jls_3Dsource(2,n_CH4))
        do l=1,lm
          tajl(:,l,k) = tajl(:,l,k)*onespo(:)
        enddo
        per_area(k) = .false.
        output_vsum(k) = .true.
      end if
C**** total chemical change for O3
      if (n_O3.gt.0) then
        k = k + 1
        kk=jls_3Dsource(1,n_O3)
        sname_tajl(k) = 'Total_change_chem+depo_'//trim(trname(n_O3))
        lname_tajl(k) = 'Total Change of '//trim(trname(n_O3))//
     &                  ' by Chemistry and deposition'
        units_tajl(k) = units_jls(kk)
        ia_tajl(k) = ia_jls(kk)
        ltop_tajl(k) = jls_ltop(kk)
        scale_tajl(k) = scale_jls(kk)*10.**(-jls_power(kk))
        tajl(:,:,k) = tajls(:,:,jls_3Dsource(1,n_O3))
     &              + tajls(:,:,jls_3Dsource(2,n_O3))
     &              + tajls(:,:,jls_3Dsource(3,n_O3))
     &              + tajls(:,:,jls_source(1,n_O3))
        do l=1,lm
          tajl(:,l,k) = tajl(:,l,k)*onespo(:)
        enddo
        per_area(k) = .false.
        output_vsum(k) = .true.
      end if
#endif

#ifdef TRACERS_COSMO
C**** ratios : Be7/Pb210 and Be10/Be7
      if (n_Be7.gt.0 .and. n_Be10.gt.0 .and. n_Pb210.gt.0) then
C*** ratio Be10/Be7
        k = k + 1
        lname_tajl(k)="Be10 to Be7 ratio"
        sname_tajl(k)="be10be7"
        units_tajl(k)=" "
        denom_tajl(k)=k_Be7
        tajl(:,:,k) = tajln(:,:,jlnt_mass,n_Be10)

C*** ratio Be7/Pb210
        k = k + 1

        lname_tajl(k)="Be7 to Pb210 ratio"
        sname_tajl(k)="be7pb210"
        units_tajl(k)="mBq/mBq"        !be sure this is 1/scalet
        denom_tajl(k)=k_Pb210
C*** scale by (Be7decay/mm_Be7)/(Pb210decay/mm_Pb210) to convert to mBq
        scale_tajl(k)=
     &       trdecay(n_Be7)*tr_mm(n_Pb210)/trdecay(n_Pb210)/tr_mm(n_Be7)
        tajl(:,:,k) = tajln(:,:,jlnt_mass,n_Be7)

      end if

#endif

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-8*d18O)
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** Concentration in water vapour
        k = k + 1
        kk=jlnt_mass
        lname_tajl(k)="Deuterium excess"
        sname_tajl(k)="dexcess"
        units_tajl(k)="per mil"
        denom_tajl(k)=k_vap
        do l=1,lm
          do j=1,jm
            d18O=tajln(j,l,kk,n1)/trw0(n1)-tajln(j,l,kk,n_water)
            dD=tajln(j,l,kk,n2)/trw0(n2)-tajln(j,l,kk,n_water)
            tajl(j,l,k)=1d3*(dD-8.*d18O)*bydxyp(j)
          enddo
        enddo

C**** Concentration in cloud water
        k = k + 1
        kk=jlnt_cldh2o
        lname_tajl(k)="Deuterium excess in cloud water"
        sname_tajl(k)="dexcess_cldh2o"
        units_tajl(k)="per mil"
        denom_tajl(k)=k_cnd
        do l=1,lm
          do j=1,jm
            d18O=tajln(j,l,kk,n1)/trw0(n1)-tajln(j,l,kk,n_water)
            dD=tajln(j,l,kk,n2)/trw0(n2)-tajln(j,l,kk,n_water)
            tajl(j,l,k)=1d3*(dD-8.*d18O)*bydxyp(j)
          enddo
        enddo

      end if

C****
C**** Calculations of 17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C**** Note: some of these definitions probably belong in JLt_TITLEX
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17

C**** Concentration in water vapour
        k = k + 1
        kk=jlnt_mass

        lname_tajl(k)="D17O excess"
        sname_tajl(k)="D17O_excess"
        units_tajl(k)="per meg"
        denom_tajl(k)=k_vap

        do l=1,lm
          do j=1,jm
            if (tajln(j,l,kk,n_water).gt.0) then
              d18O=tajln(j,l,kk,n1)/(trw0(n1)*tajln(j,l,kk,n_water))
              d17O=tajln(j,l,kk,n2)/(trw0(n2)*tajln(j,l,kk,n_water))
              tajl(j,l,k)=1d6*(log(d17O)-0.529d0*log(d18O))*
     &             tajl(j,l,k_vap)
            else
              tajl(j,l,k)=0.
            end if
          end do
        end do

C**** Concentration in cloud water
        k = k + 1
        kk=jlnt_cldh2o
        lname_tajl(k)="D17O excess in cloud water"
        sname_tajl(k)="D17O_excess_cldh2o"
        units_tajl(k)="per meg"
        denom_tajl(k)=k_cnd

        do l=1,lm
          do j=1,jm
            if (tajln(j,l,kk,n_water).gt.0) then
              d18O=tajln(j,l,kk,n1)/(trw0(n1)*tajln(j,l,kk,n_water))
              d17O=tajln(j,l,kk,n2)/(trw0(n2)*tajln(j,l,kk,n_water))
              tajl(j,l,k)=1d6*(log(d17O)-0.529d0*log(d18O))*
     &             tajl(j,l,k_cnd)
            else
              tajl(j,l,k)=0.
            end if
          end do
        end do

      end if
#endif

      ktajl_out = k

c
c if necessary, reallocate tajl to be the right size
c
      if(ktajl_out.lt.size(tajl,3)) then
        allocate(tajl_tmp(jm,lm,ktajl_out))
        tajl_tmp(:,:,1:ktajl_out) = tajl(:,:,1:ktajl_out)
        deallocate(tajl)
        allocate(tajl(jm,lm,ktajl_out))
        tajl = tajl_tmp
        deallocate(tajl_tmp)
      elseif(ktajl_out.gt.size(tajl,3)) then
        call stop_model('error: ktajl_out > size(tajl,3)',255)
      endif

#ifdef TRACERS_SPECIAL_Lerner
c adjust sname elements beginning with 14C
      do k=1,ktajl_out
        if(sname_tajl(k)(1:3).eq.'14C') then
          sname_tajl(k) = 'a'//trim(sname_tajl(k))
        endif
      enddo
#endif

c
c compute hemispheric, global, vertical means/sums
c
      if(.not.allocated(hemis_tajl)) then
        allocate(hemis_tajl(3,lm,ktajl_out))
        allocate(vmean_tajl(jm+3,1,ktajl_out))
      endif
      hemis_tajl = 0.
      vmean_tajl = 0.
      do k=1,ktajl_out
        if(per_area(k)) then
          hemfac = 2./areag
        else
          hemfac = 1.
        endif
        do l=1,lm
          fj(:) = tajl(:,l,k)
          if(.not. per_area(k)) fj(:) = fj(:)/dxyp(:)
          j1 = 1; j2 = jm/2
          hemis_tajl(1,l,k) = hemfac*sum(fj(j1:j2)*dxyp(j1:j2))
          j1 = jm/2+1; j2 = jm
          hemis_tajl(2,l,k) = hemfac*sum(fj(j1:j2)*dxyp(j1:j2))
          hemis_tajl(3,l,k) = hemis_tajl(1,l,k)+hemis_tajl(2,l,k)
        enddo
        if(per_area(k)) hemis_tajl(3,:,k) = .5*hemis_tajl(3,:,k)
        vmean_tajl(1:jm,1,k) = sum(tajl(:,:,k),dim=2)
        vmean_tajl(jm+1:jm+3,1,k) = sum(hemis_tajl(:,:,k),dim=2)
      enddo

c
c Declare the dimensions and metadata of TAJL output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      cdl_tajl = cdl_jl_template ! invoke a copy method later
      do k=1,ktajl_out
        if(trim(units_tajl(k)).eq.'unused') cycle
        if(lgrid_tajl(k).eq.1) then
          zstr='plm'
        else
          zstr='ple'
        endif
        set_miss = denom_tajl(k).ne.0
        call add_var(cdl_tajl,'float '//trim(sname_tajl(k))//'('//
     &       trim(zstr)//',lat_budg) ;',
     &       long_name=trim(lname_tajl(k)),
     &       units=trim(units_tajl(k)),
     &       auxvar_string=
     &          'float '//trim(sname_tajl(k))//'_hemis('//
     &          trim(zstr)//',shnhgm) ;',
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
        if(pow_tajl(k).ne.0) then
          write(powstr,'(i2)') pow_tajl(k)
          call add_varline(cdl_tajl,
     &         trim(sname_tajl(k))//':prtpow = '//trim(powstr)//' ;')
        endif
        if(ltop_tajl(k).ne.lm) then
          write(ltopstr,'(i3)') ltop_tajl(k)
          call add_varline(cdl_tajl,
     &         trim(sname_tajl(k))//':ltop = '//trim(ltopstr)//' ;')
        endif
        if(denom_tajl(k).gt.0 .or. output_vsum(k)) then
          if(make_timeaxis) then ! hacky logic
            call add_var(cdl_tajl, 'float '//trim(sname_tajl(k))//
     &           '_vmean(time,lat_budg_plus3) ;')
          else
            call add_var(cdl_tajl, 'float '//trim(sname_tajl(k))//
     &           '_vmean(lat_budg_plus3) ;')
          endif
        endif
      enddo

      RETURN
      END SUBROUTINE DIAGJLT_prep

      subroutine diagijt_prep
! comments to be added
      use resolution, only : im,jm
      use model_com, only: idacc
      use OldTracer_mod, only: dodrydep, dowetdep
      use OldTracer_mod, only: trname, trw0, src_dist_index
      use tracer_com, only: ntm, n_water
#ifdef TRACERS_SPECIAL_O18
      use tracer_com, only: n_h2o18, n_hdo, n_h2o17
#endif
      use diag_com
      use diag_com_rad, only : ij_cldcv
      use mdiag_com, only : sname_strlen,make_timeaxis
      use trdiag_com, only : taijn=>taijn_loc, taijs=>taijs_loc,
     &     ktaij_,ktaij_out,taij=>taij_out,
     &     scale_taij,cdl_taij,cdl_taij_latlon,hemis_taij,
     &     ir_taij,ia_taij,denom_taij,lname_taij,sname_taij,units_taij,
     &     sname_tij, lname_tij, ijts_HasArea, ijts_clrsky, denom_ijts,
     &     units_tij, scale_tij, lname_ijts,  sname_ijts,
     &     units_ijts,  scale_ijts,  ia_ijts, ktaij, ktaijs, dname_ijts,
     &     tij_drydep, tij_gsdep, tij_surf, tij_grnd, tij_prec,
     &     dname_tij, ijts_pocean
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
     &     ,to_per_mil
#endif
      use constant, only : teeny
      use domain_decomp_atm, only : grid
      use domain_decomp_atm, only : getDomainBounds,am_i_root,sumxpe
      use geom, only : byaxyp,axyp,lat2d,areag
      use cdl_mod
      implicit none
      integer ::  i,j,k,kx,k1,n,n1,n2,khem
      integer :: k_Be7,k_Pb210
      integer :: i_0,i_1,j_0,j_1, i_0h,i_1h,j_0h,j_1h
      real*8, dimension(:,:,:), allocatable :: taij_tmp
      character(len=16) :: ijstr
      real*8, dimension(:,:), allocatable :: shnh_loc,shnh
      character(len=sname_strlen), dimension(ktaij_) :: dname_taij
      logical :: set_miss

      logical :: have_south_pole, have_north_pole
      call getDomainBounds(grid, have_south_pole = have_south_pole,
     &               have_north_pole = have_north_pole)

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

C**** Fill in the undefined pole box duplicates
      if(have_south_pole) then
        do i=2,im
          taijn(i,1,:,:) = taijn(1,1,:,:)
          taijs(i,1,:) = taijs(1,1,:)
        enddo
      endif
      if(have_north_pole) then
        do i=2,im
          taijn(i,jm,:,:) = taijn(1,jm,:,:)
          taijs(i,jm,:) = taijs(1,jm,:)
        enddo
      endif

      do k=1,ktaij_
        denom_taij(k) = 0
        ia_taij(k) = ia_src
        sname_taij(k) = 'unused'
        dname_taij(k) = ''
        lname_taij(k) = 'unused'
        units_taij(k) = 'unused'
        scale_taij(k) = 1.
      enddo

      if(ijts_clrsky.gt.0) then
        taijs(:,:,ijts_clrsky) =
     &       real(idacc(ia_rad))-aij_loc(:,:,ij_cldcv)
      endif
      if(ijts_pocean.gt.0) then
        taijs(:,:,ijts_pocean) = aij_loc(:,:,ij_pocean)
      endif

      k = 0

c
c Tracer sums/means and ground conc
c
      do n=1,NTM
      if (src_dist_index(n)>1) cycle
      do kx=1,ktaij
        if (index(lname_tij(kx,n),'unused').gt.0) cycle
        k = k+1

        taij(i_0:i_1,j_0:j_1,k) = taijn(i_0:i_1,j_0:j_1,kx,n)

        sname_taij(k) = sname_tij(kx,n)
        dname_taij(k) = dname_tij(kx,n)
        lname_taij(k) = lname_tij(kx,n)
        units_taij(k) = units_tij(kx,n)
        scale_taij(k) = scale_tij(kx,n)
        ir_taij(k) = ir_log2

#ifdef TRACERS_COSMO
        if(kx.eq.tij_surf .and. n.eq.n_Be7) k_Be7 = k
        if(kx.eq.tij_surf .and. n.eq.n_Pb210) k_Pb210 = k
#endif

#ifdef TRACERS_WATER
        if(index(units_taij(k),'er mil').gt.0 .and. n.ne.n_water) then
          do j=j_0,j_1; do i=i_0,i_1
            taij(i,j,k) =
     &           1d3*(taij(i,j,k)/trw0(n)-taijn(i,j,kx,n_water))
          enddo       ; enddo
        endif
#endif

      enddo ! end loop over kx

#if (defined TRACERS_WATER) && (defined TRACERS_DRYDEP)
c
c dry/(dry+wet) deposition fraction
c
      if (dodrydep(n).and.dowetdep(n)) then
        k=k+1
        sname_taij(k) = trim(trname(n))//"_pc_dry_dep"
        lname_taij(k) = trim(trname(n))//" Percent Dry Deposition"
        units_taij(k) = "%"
        ir_taij(k) = ir_pct
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = taijn(i,j,tij_drydep,n)+taijn(i,j,tij_gsdep,n)
          taij(i,j,k+1) = taij(i,j,k) + taijn(i,j,tij_prec,n)
          taij(i,j,k) = taij(i,j,k)*100.
        enddo
        enddo
        denom_taij(k) = k+1
        k=k+1
      endif
#endif

      enddo ! end loop over tracers

#ifdef TRACERS_COSMO
      if (n_Be7.gt.0 .and. n_Be10.gt.0 .and. n_Pb210.gt.0) then
C**** Be10/Be7
        k=k+1
        sname_taij(k) = "be10be7_ij"
        lname_taij(k) = "surface ratio Be10 to Be7"
        units_taij(k) = " "
        ir_taij(k) = ir_0_18
        ia_taij(k) = ia_srf
        denom_taij(k) = k_Be7
        taij(i_0:i_1,j_0:j_1,k) = taijn(i_0:i_1,j_0:j_1,tij_surf,n_Be10)
C**** Be7/Pb210
        k=k+1
        sname_taij(k) = "be7pb210_ij"
        lname_taij(k) = "surface ratio Be7 to Pb210"
        units_taij(k) = "mBq/mBq "
        ir_taij(k) = ir_0_180
        ia_taij(k) = ia_srf
        denom_taij(k) = k_Pb210
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = taij(i,j,tij_surf,n_Be7)
C*** scale by (Be7decay/mm_Be7)/(Pb210decay/mm_Pb210) to convert to mBq
          taij(i,j,k)=taij(i,j,k)*trdecay(n_Be7)*tr_mm(n_Pb210)
     &             /trdecay(n_Pb210)/tr_mm(n_Be7)
        enddo
        enddo
      endif
#endif

#ifdef TRACERS_SPECIAL_O18
C****
C**** Calculations of deuterium excess (d=dD-d18O)
C****
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** precipitation
        k=k+1
        sname_taij(k) = "prec_ij_dex"
        dname_taij(k) = sname_tij(tij_prec,n_Water)
        lname_taij(k) = "Deuterium excess in precip"
        units_taij(k) = "per mil"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = 1d3*(taijn(i,j,tij_prec,n2)/trw0(n2)-
     &                  8.*taijn(i,j,tij_prec,n1)/trw0(n1)+
     &                  7.*taijn(i,j,tij_prec,n_water))
        enddo
        enddo
C**** ground concentration
        k=k+1
        sname_taij(k) = "grnd_ij_dex"
        dname_taij(k) = sname_tij(tij_grnd,n_Water)
        lname_taij(k) = "Deuterium excess at Ground"
        units_taij(k) = "per mil"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
        do i=i_0,i_1
          taij(i,j,k) = 1d3*(taijn(i,j,tij_grnd,n2)/trw0(n2)-
     &                  8.*taijn(i,j,tij_grnd,n1)/trw0(n1)+
     &                  7.*taijn(i,j,tij_grnd,n_water))
        enddo
        enddo
      end if

C****
C**** Calculations of D17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C****
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
C**** precipitation
        k=k+1
        sname_taij(k) = "prec_ij_D17O"
        dname_taij(k) = sname_tij(tij_prec,n_Water)
        lname_taij(k) = "D17O excess in precip"
        units_taij(k) = "per meg"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
          do i=i_0,i_1
            if (taijn(i,j,tij_prec,n_water).gt.0) then
              taij(i,j,k) = 1d6*taijn(i,j,tij_prec,n_water)*
     &             (log(taijn(i,j,tij_prec,n2)/trw0(n2))-
     &             0.529d0*log(taijn(i,j,tij_prec,n1)/trw0(n1))-
     &             0.471d0*log(taijn(i,j,tij_prec,n_water)))
            else
              taij(i,j,k)=0.
            end if
          end do
        end do
C**** ground concentration
        k=k+1
        sname_taij(k) = "grnd_ij_D17O"
        dname_taij(k) = sname_tij(tij_grnd,n_Water)
        lname_taij(k) = "D17O excess at Ground"
        units_taij(k) = "per meg"
        ir_taij(k) = ir_m45_130
        ia_taij(k) = ia_src
        do j=j_0,j_1
          do i=i_0,i_1
            if (taijn(i,j,tij_grnd,n_water).gt.0) then
              taij(i,j,k) = 1d6*taijn(i,j,tij_grnd,n_water)*
     &             (log(taijn(i,j,tij_grnd,n2)/trw0(n2))-
     &             0.529d0*log(taijn(i,j,tij_grnd,n1)/trw0(n1))-
     &             0.471d0*log(taijn(i,j,tij_grnd,n_water)))
            else
              taij(i,j,k)=0.
            end if
          end do
        end do
      end if
#endif

c
c tracer sources and sinks
c
      do kx=1,ktaijs
        if (index(lname_ijts(kx),'unused').gt.0) cycle

        k = k+1

        sname_taij(k) = sname_ijts(kx)
        dname_taij(k) = dname_ijts(kx)
        lname_taij(k) = lname_ijts(kx)
        units_taij(k) = units_ijts(kx)
        ir_taij(k) = ir_log2    ! should be the correct default
        ia_taij(k) = ia_ijts(kx)
        scale_taij(k) = scale_ijts(kx)

        taij(i_0:i_1,j_0:j_1,k) = taijs(i_0:i_1,j_0:j_1,kx)

        if(ijts_HasArea(kx)) then
          do j=j_0,j_1
          do i=i_0,i_1
            taij(i,j,k) = taij(i,j,k)*byaxyp(i,j)
          enddo
          enddo
        endif

      enddo

      if(any(dname_taij(1:k).eq.'oicefr')) then      
      k=k+1
      sname_taij(k) = 'oicefr'
      lname_taij(k) = lname_ij(ij_rsoi)
      units_taij(k) = units_ij(ij_rsoi)
      ia_taij(k) = ia_ij(ij_rsoi)
      scale_taij(k) = scale_ij(ij_rsoi)
      taij(:,:,k) = aij_loc(:,:,ij_rsoi)
      endif

#ifdef TRACERS_SPECIAL_O18
      !Add ocean fraction to taij array, in
      !order to use it as a denominator:
      if(any(dname_taij(1:k).eq.'ocnfrac')) then
        k=k+1
        sname_taij(k) = 'ocnfrac'
        lname_taij(k) = lname_ij(ij_pocean)
        units_taij(k) = units_ij(ij_pocean)
        ia_taij(k) = ia_ij(ij_pocean)
        scale_taij(k) = scale_ij(ij_pocean)
        taij(:,:,k) = aij_loc(:,:,ij_pocean)
      endif
#endif

      ktaij_out = k

c
c find the indices of string-specified denominators 
c
      call FindStrings(dname_taij,sname_taij,denom_taij,ktaij_out)

c      do k=1,ktaij_out
c        if(len_trim(dname_taij(k)).gt.0) then
c          do kk=1,ktaij_out
c            if(trim(sname_taij(kk)).eq.trim(dname_taij(k))) then
c              denom_taij(k) = kk
c              exit
c            endif
c          enddo
c        endif
c      enddo

c
c if necessary, reallocate taij to be the right size
c
      if(ktaij_out.lt.size(taij,3)) then
        allocate(taij_tmp(i_0:i_1,j_0:j_1,ktaij_out))
        taij_tmp(i_0:i_1,j_0:j_1,1:ktaij_out) =
     &      taij(i_0:i_1,j_0:j_1,1:ktaij_out)
        deallocate(taij)
        allocate(taij(i_0h:i_1h,j_0h:j_1h,ktaij_out))
        taij(i_0:i_1,j_0:j_1,1:ktaij_out) =
     &       taij_tmp(i_0:i_1,j_0:j_1,1:ktaij_out)
        deallocate(taij_tmp)
      elseif(ktaij_out.gt.size(taij,3)) then
        call stop_model('error: ktaij_out > size(taij,3)',255)
      endif

#ifdef TRACERS_SPECIAL_Lerner
c adjust sname elements beginning with 14C
      do k=1,ktaij_out
        if(sname_taij(k)(1:3).eq.'14C') then
          sname_taij(k) = 'a'//trim(sname_taij(k))
        endif
      enddo
#endif

c
c hemispheric and global averages
c
      allocate(shnh_loc(2,ktaij_out),shnh(2,ktaij_out))
      shnh_loc = 0.
      do k=1,ktaij_out
        do j=j_0,j_1
        do i=i_0,i_1
          khem = 1
          if(lat2d(i,j).ge.0.) khem = 2
          shnh_loc(khem,k) = shnh_loc(khem,k) + axyp(i,j)*taij(i,j,k)
        enddo
        enddo
      enddo
      call sumxpe(shnh_loc,shnh)
      if(am_i_root()) then
        if(.not. allocated(hemis_taij)) then
          allocate(hemis_taij(1,3,ktaij_out))
        endif
        do k=1,ktaij_out
          shnh(:,k) = 2.*shnh(:,k)/areag
          hemis_taij(1,1:2,k) = shnh(1:2,k)
          hemis_taij(1,3,k) = .5*(shnh(1,k)+shnh(2,k))
        enddo
      endif

c
c Declare the dimensions and metadata of TAIJ output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      cdl_taij = cdl_ij_template ! invoke a copy method later

#ifdef CUBED_SPHERE
      cdl_taij_latlon = cdl_ij_latlon_template ! invoke a copy method later
      ijstr='(tile,y,x) ;'
#else
      ijstr='(lat,lon) ;'
#endif
      do k=1,ktaij_out
        if(trim(sname_taij(k)).eq.'unused') cycle
        set_miss = denom_taij(k).ne.0
        call add_var(cdl_taij,
     &       'float '//trim(sname_taij(k))//trim(ijstr),
     &       long_name=trim(lname_taij(k)),
     &       auxvar_string=
     &           'float '//trim(sname_taij(k))//'_hemis(shnhgm);',
     &       units=trim(units_taij(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#ifdef CUBED_SPHERE
        call add_var(cdl_taij_latlon,
     &       'float '//trim(sname_taij(k))//'(lat,lon);',
     &       long_name=trim(lname_taij(k)),
     &       auxvar_string=
     &           'float '//trim(sname_taij(k))//'_hemis(shnhgm);',
     &       units=trim(units_taij(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#endif
      enddo

      return
      end subroutine diagijt_prep

      SUBROUTINE DIAGIJLt_prep
!     comments to be added
      use resolution, only : im,jm,lm
      use model_com, only: idacc
      use OldTracer_mod, only: trw0, src_dist_index
      use tracer_com, only: ntm, n_water
#ifdef TRACERS_SPECIAL_O18
      use tracer_com, only: n_h2o18, n_hdo, n_h2o17
#endif
      use diag_com
      use diag_com_rad, only : ij_cldcv,ijl_cf
      use mdiag_com, only : sname_strlen,make_timeaxis 
      use trdiag_com, only : taijln=>taijln_loc, taijls=>taijls_loc,
     &     ktaijl_,ktaijl_out,taijl=>taijl_out,scale_taijl,ir_taijl,
     &     ia_taijl,denom_taijl,lname_taijl,sname_taijl,units_taijl,
     &     cdl_taijl, cdl_taijl_latlon,sname_ijlt,lname_ijlt,dname_ijlt,
     &     units_ijlt, sname_ijt, lname_ijt, units_ijt, scale_ijt,
     &     ir_ijlt, ia_ijlt, scale_ijlt, ktaijl,
     &     ijlt_clrsky2d
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
     &     ,to_per_mil
#endif
      use constant, only : teeny
      use domain_decomp_atm, only : grid,getDomainBounds,am_i_root
      use geom, only : byaxyp
      use cdl_mod
      implicit none
      integer i,j,l,k,kx,kk,n,n1,n2
      real*8 :: r1,r2
      integer :: k_water
      integer :: i_0,i_1,j_0,j_1, i_0h,i_1h,j_0h,j_1h
      real*8, dimension(:,:,:,:), allocatable :: taijl_tmp
      character(len=16) :: zstr,hstr,tstr
      character(len=sname_strlen), dimension(ktaijl_) :: dname_taijl
      logical :: set_miss
      logical :: have_south_pole, have_north_pole
      call getDomainBounds(grid, have_south_pole = have_south_pole,
     &     have_north_pole = have_north_pole)

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

C**** Fill in the undefined pole box duplicates
      if(have_south_pole) then
        do i=2,im
          taijln(i,1,:,:) = taijln(1,1,:,:)
          taijls(i,1,:,:) = taijls(1,1,:,:)
        enddo
      endif
      if(have_north_pole) then
        do i=2,im
          taijln(i,jm,:,:) = taijln(1,jm,:,:)
          taijls(i,jm,:,:) = taijls(1,jm,:,:)
        enddo
      endif

      do k=1,ktaijl_
        denom_taijl(k) = 0
        ia_taijl(k) = ia_src
        sname_taijl(k) = 'unused'
        dname_taijl(k) = ''
        lname_taijl(k) = 'unused'
        units_taijl(k) = 'unused'
        scale_taijl(k) = 1.
      enddo

      if(ijlt_clrsky2d.gt.0) then
        ! following two lines correspond to 3D weighting
!        taijls(:,:,:,ijlt_clrsky) =
!     &       real(idacc(ia_rad))-aijl_loc(:,:,:,ijl_cf)
  ! but RADIA 2D variable OPNSKY is current weight for per-layer opt depths
        do l=1,lm; do j=j_0,j_1; do i=i_0,i_1
          taijls(i,j,l,ijlt_clrsky2d) =
     &         real(idacc(ia_rad))-aij_loc(i,j,ij_cldcv)
        enddo    ; enddo       ; enddo
      endif

      k = 0

C**** Tracer concentrations
      do n=1,NTM
        if (src_dist_index(n)>1) cycle
        k = k+1
        sname_taijl(k) = sname_ijt(n)
        lname_taijl(k) = lname_ijt(n)
        units_taijl(k) = units_ijt(n)
        ir_taijl(k) = ir_log2
        ia_taijl(k) = ia_src
        scale_taijl(k) = scale_ijt(n)
        do l=1,lm
          do j=j_0,j_1; do i=i_0,i_1
            taijl(i,j,l,k) = taijln(i,j,l,n)*byaxyp(i,j)
          enddo       ; enddo
        enddo
#ifdef TRACERS_WATER
        if(n.eq.n_water) k_water = k
        if(to_per_mil(n).gt.0) then
          do l=1,lm
            do j=j_0,j_1; do i=i_0,i_1
              taijl(i,j,l,k) = 1d3*(taijl(i,j,l,k)/trw0(n)
     &             -byaxyp(i,j)*taijln(i,j,l,n_water))
            enddo       ; enddo
          enddo
          denom_taijl(k) = k_water
          if(n.eq.n_water) then ! save a non-per-mil denom
            k = k+1
            k_water = k
            denom_taijl(k-1) = k_water
            do l=1,lm
              do j=j_0,j_1; do i=i_0,i_1
                taijl(i,j,l,k) =  taijln(i,j,l,n)*byaxyp(i,j)
              enddo       ; enddo
            enddo
            ia_taijl(k) = ia_taijl(k-1)
          endif
        endif
#endif
      enddo

C**** Tracer specials 
      do kx=1,ktaijl
        if (index(lname_ijlt(kx),'unused').gt.0) cycle
        k = k+1
        sname_taijl(k) = sname_ijlt(kx)
        dname_taijl(k) = dname_ijlt(kx)
        lname_taijl(k) = lname_ijlt(kx)
        units_taijl(k) = units_ijlt(kx)
        ir_taijl(k) = ir_ijlt(kx)
        ia_taijl(k) = ia_ijlt(kx)
        scale_taijl(k) = scale_ijlt(kx)
        do l=1,lm
          do j=j_0,j_1; do i=i_0,i_1
            taijl(i,j,l,k) = taijls(i,j,l,kx)
          enddo       ; enddo
        enddo
      enddo

#ifdef TRACERS_SPECIAL_O18
C**** 
C**** Calculations of deuterium excess (d=dD-d18O)
C**** 
      if (n_H2O18.gt.0 .and. n_HDO.gt.0) then
        n1=n_H2O18
        n2=n_HDO
C**** water vapour
        k=k+1
        denom_taijl(k) = k_water
        lname_taijl(k) = "Deuterium excess water vapour"
        sname_taijl(k) = "wvap_ij_dex"
        units_taijl(k) = "per mil"
        ir_taijl(k) = ir_m45_130
        ia_taijl(k) = ia_src
        do l=1,lm
          do j=j_0,j_1
            do i=i_0,i_1
              taijl(i,j,l,k) = 1d3*(taijln(i,j,l,n2)/trw0(n2)-
     &             8.*taijln(i,j,l,n1)/trw0(n1)+
     &             7.*taijln(i,j,l,n_water))*byaxyp(i,j)
            enddo
          enddo
        enddo
      end if

C**** 
C**** Calculations of D17O excess (D17O=ln(d17O+1)-0.529*ln(d18O+1))
C**** 
      if (n_H2O18.gt.0 .and. n_H2O17.gt.0) then
        n1=n_H2O18
        n2=n_H2O17
C**** water vapour
        k=k+1
        denom_taijl(k) = k_water
        lname_taijl(k) = "D17O excess water vapour"
        sname_taijl(k) = "wvap_ij_D17O"
        units_taijl(k) = "per meg"
        ir_taijl(k) = ir_m45_130
        ia_taijl(k) = ia_src
        do l=1,lm
          do j=j_0,j_1
            do i=i_0,i_1
              if (taijln(i,j,l,n_water).gt.0) then
                r1 = taijln(i,j,l,n1)/(trw0(n1)*taijln(i,j,l,n_water))
                r2 = taijln(i,j,l,n2)/(trw0(n2)*taijln(i,j,l,n_water))
                taijl(i,j,l,k) = 1d6*taijl(i,j,l,k_water)*
     &               (log(r2)-0.529d0*log(r1))
              else
                taijl(i,j,l,k)=0.
              end if
            end do
          end do             
        end do
      end if
#endif

      ktaijl_out = k

c
c find the indices of string-specified denominators 
c
      call FindStrings(dname_taijl,sname_taijl,denom_taijl,ktaijl_out)

c     
c     if necessary, reallocate taijl to be the right size
c     
      if(ktaijl_out.lt.size(taijl,4)) then
        allocate(taijl_tmp(i_0:i_1,j_0:j_1,lm,ktaijl_out))
        taijl_tmp(i_0:i_1,j_0:j_1,:,1:ktaijl_out) =
     &       taijl(i_0:i_1,j_0:j_1,:,1:ktaijl_out)
        deallocate(taijl)
        allocate(taijl(i_0h:i_1h,j_0h:j_1h,lm,ktaijl_out))
        taijl(i_0:i_1,j_0:j_1,:,1:ktaijl_out) =
     &       taijl_tmp(i_0:i_1,j_0:j_1,:,1:ktaijl_out)
        deallocate(taijl_tmp)
      elseif(ktaijl_out.gt.size(taijl,4)) then
        call stop_model('error: ktaijl_out > size(taijl,4)',255)
      endif

#ifdef TRACERS_SPECIAL_Lerner
c     adjust sname elements beginning with 14C
      do k=1,ktaijl_out
        if(sname_taijl(k)(1:3).eq.'14C') then
          sname_taijl(k) = 'a'//trim(sname_taijl(k))
        endif
      enddo
#endif

c     
c     Declare the dimensions and metadata of TAIJL output fields using
c     netcdf CDL notation.  The C convention for dimension ordering
c     must be used (reversed wrt Fortran).
c     
      cdl_taijl = cdl_ijl_template ! invoke a copy method later

#ifdef CUBED_SPHERE
      cdl_taijl_latlon = cdl_ijl_latlon_template ! invoke a copy method later
      tstr='(tile,'
      hstr=',y,x) ;'
#else
      tstr='('
      hstr=',lat,lon) ;'
#endif
      do k=1,ktaijl_out
        if(trim(sname_taijl(k)).eq.'unused') cycle
        zstr='level'
        set_miss = denom_taijl(k).ne.0
        call add_var(cdl_taijl,
     &       'float '//trim(sname_taijl(k))//
     &       trim(tstr)//trim(zstr)//trim(hstr),
     &       long_name=trim(lname_taijl(k)),
     &       units=trim(units_taijl(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#ifdef CUBED_SPHERE
        call add_var(cdl_taijl_latlon,
     &       'float '//trim(sname_taijl(k))//
     &       '('//trim(zstr)//',lat,lon);',
     &       long_name=trim(lname_taijl(k)),
     &       units=trim(units_taijl(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#endif
      enddo

      return
      end subroutine diagijlt_prep

#endif /* TRACERS_ON */

      SUBROUTINE DIAGTCP_prep
! comments to be added
      USE CONSTANT, only: teeny
      USE MODEL_COM, only : idacc
      USE GEOM, only: areag
      USE DIAG_COM, only: jm=>jm_budg,dxyp=>dxyp_budg,cdl_latbudg
      USE TRDIAG_COM, only : natmtrcons,nocntrcons,tconsrv,
     &     ktcon,ktcon_out,nsum_tcon,scale_tcon,
     &     ia_tcon,title_tcon,hemis_tconsrv,name_tconsrv,tconsrv_out,
     &     ia_tcon_out,scale_tcon_out,sname_tconsrv_out,cdl_tconsrv
      use domain_decomp_atm, only : am_i_root
      use cdl_mod
      IMPLICIT NONE
      INTEGER :: j,n,k,kk,KTCON_max,j1,j2,ntmoa
      real*8 :: hemfac
      character*80 :: sname
      character(len=80), dimension(:), allocatable :: titles

      ntmoa = natmtrcons+nocntrcons

      call gather_zonal_tcons

      if(.not.am_i_root()) return

      if(.not.allocated(tconsrv_out)) then

      allocate(titles(ktcon*ntmoa))

c
c determine how many actual output fields there are
c
        ktcon_out = 0
        do n=1,ntmoa
          do k=ktcon,1,-1
            if(nsum_tcon(k,n).gt.0) ktcon_max = nsum_tcon(k,n)
          enddo
          do k=1,ktcon_max
            ktcon_out = ktcon_out + 1
          enddo
        enddo
c
c allocate space for actual number of outputs
c
        allocate(tconsrv_out(jm,ktcon_out))
        allocate(hemis_tconsrv(3,ktcon_out))
        allocate(ia_tcon_out(ktcon_out))
        allocate(scale_tcon_out(ktcon_out))
        allocate(sname_tconsrv_out(ktcon_out))
c
c copy metadata
c
        kk = 0
        do n=1,ntmoa
          do k=ktcon,1,-1
            if(nsum_tcon(k,n).gt.0) ktcon_max = nsum_tcon(k,n)
          enddo
          do k=1,ktcon_max
            kk = kk + 1
            ia_tcon_out(kk) = ia_tcon(k,n)
            scale_tcon_out(kk) = scale_tcon(k,n)
            sname_tconsrv_out(kk) = name_tconsrv(k,n)
            titles(kk) = title_tcon(k,n)
          enddo
        enddo

c
c Declare the dimensions and metadata of TCONSRV output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
        cdl_tconsrv = cdl_latbudg ! invoke a copy method later
        do k=1,ktcon_out
          sname = trim(sname_tconsrv_out(k))
          call add_var(cdl_tconsrv,
     &         'float '//trim(sname)//'(lat_budg) ;',
     &         long_name=trim(titles(k)))
          call add_var(cdl_tconsrv,
     &         'float '//trim(sname)//'_hemis(shnhgm) ;')
        enddo

        deallocate(titles)

      endif ! memory allocation and setup

c
c copy the nonzero contents of tconsrv into tconsrv_out
c also calculate sums of changes, hemispheric/global avgs
c
      hemfac = 2./sum(dxyp)
      kk = 0
      do n=1,ntmoa
        do k=ktcon,1,-1
C**** LOOP BACKWARDS SO THAT INITIALIZATION IS DONE BEFORE SUMMATION!
          if(nsum_tcon(k,n).eq.0) then
            tconsrv(:,k,n)=0.
          elseif(nsum_tcon(k,n).gt.0) then
            tconsrv(:,nsum_tcon(k,n),n)=
     &      tconsrv(:,nsum_tcon(k,n),n)+tconsrv(:,k,n)
     &           *scale_tcon(k,n)*idacc(12)/(idacc(ia_tcon(k,n))+teeny)
            ktcon_max = nsum_tcon(k,n)
          endif
        enddo
        do k=1,ktcon_max
          kk = kk + 1
          tconsrv_out(:,kk) = tconsrv(:,k,n)
          j1 = 1; j2 = jm/2
          hemis_tconsrv(1,kk) = hemfac*sum(tconsrv_out(j1:j2,kk))
          j1 = jm/2+1; j2 = jm
          hemis_tconsrv(2,kk) = hemfac*sum(tconsrv_out(j1:j2,kk))
          hemis_tconsrv(3,kk) = .5*sum(hemis_tconsrv(1:2,kk))
          tconsrv_out(:,kk) = tconsrv_out(:,kk)/dxyp
        enddo
      enddo

      RETURN
      END SUBROUTINE DIAGTCP_prep
