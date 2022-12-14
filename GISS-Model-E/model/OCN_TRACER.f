#include "rundeck_opts.h"

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@sum  OCN_TRACER: tracer-dependent routines for GISS Ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Tracer initialisation + sources: tracer_ic_ocean
!@+
!@auth Jean Lerner/Gavin Schmidt

      subroutine tracer_ic_ocean(atmocn)
!@sum tracer_ic_ocean initialise ocean tracers
!@auth Gavin Schmidt
!@ver 1.0
      USE MODEL_COM, only: itime,itimei
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#ifdef TRACERS_SPECIAL_O18
      USE OCN_TRACER_COM, only : water_tracer_ic
#endif

      USE OCEAN, only : im,jm,lmo,dxypo,mo,lmm,imaxj,oXYP,use_qus
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo,tzmo,mo,s0m,sxmo,symo,szmo,oc_tracer_mean
#endif
      use ocean, only : nbyzm,i1yzm,i2yzm
      USE SEAICE, only : xsi,lmi
      USE STRAITS, only : nmst!,msist,ssist
#ifdef TRACERS_OCEAN
     *     ,lmst,ist,jst,xst,yst,mmst,s0mst,sxmst,szmst,trmst,txmst
     *     ,tzmst
#endif
!#ifdef TRACERS_WATER
!     *     ,trsist
!#endif
      USE FILEMANAGER, only : openunit,closeunit
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, haveLatitude,
     *     broadcast, GLOBALSUM
      USE OCEANR_DIM, only : grid=>ogrid
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, pack_data, unpack_data
      !USE OCEAN, only : gather_ocean
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      USE Dictionary_mod
      use pario, only : read_dist_data,par_open,par_close
      IMPLICIT NONE
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     tr_ic
      type(atmocn_xchng_vars) :: atmocn
c
      integer n,i,j,l,nst,i1,j1,i2,j2,ll
#ifdef TRACERS_SPECIAL_O18
      integer iu_O18ic,ip1,im1
      character*80 title
      real*4, dimension(im,jm,lmo) :: t0m4,tzm4
      real*8 fwm,afac
#endif
      real*8 t01,t02,trsum,wsum,tratio,frac_tr
      real*8, dimension(im,jm,lmo) :: mo_glob,s0m_glob,trmo_glob
      real*8 :: OTRACJ(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) 
      INTEGER :: J_0S, J_1S, J_0, J_1, J_0H, J_1H
      type(ocn_tracer_entry), pointer :: entry
      logical :: glob_used=.false.
      integer :: fid=-1, nn

      call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     *     J_STRT = J_0, J_STOP = J_1, 
     *     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

#ifdef TRACERS_SPECIAL_O18
      if (is_set_param('water_tracer_ic')) then
        call get_param('water_tracer_ic', water_tracer_ic)
      endif
      !call sync_param ( "water_tracer_ic",water_tracer_ic  )
#endif

C**** Note that only sea ice related arrays are initialised if
C**** only TRACERS_WATER is true.

      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        if (.not.entry%need_ic) cycle
        if (entry%from_file) then
          if (itime.ne.itimei) cycle
          if (fid<0) fid=par_open(grid,'OCN_TRACER_IC','read')
          call read_dist_data(grid,fid,trim(entry%trname),tr_ic)
          do l=1,lmo
            do j=j_0,j_1
              trmo(:,j,l,n) = 0.
              do nn=1,nbyzm(j,l)
                do i=i1yzm(nn,j,l),i2yzm(nn,j,l)
                  trmo(i,j,l,n) = tr_ic(i,j,l)*mo(i,j,l)*dxypo(j)
                enddo
              enddo
            enddo
          enddo
          cycle
        endif
        if (itime.eq.entry%itime_tr0) then
        select case (entry%trname(1:6))

        case default
#ifdef TRACERS_OCEAN
          trmo(:,:,:,n)=0.
          txmo(:,:,:,n)=0.
          tymo(:,:,:,n)=0.
          tzmo(:,:,:,n)=0.
C**** straits
          if (am_i_root()) then
            trmst(:,:,n)=0.
            txmst(:,:,n)=0.
            tzmst(:,:,n)=0.
          end if
          CALL broadcast(grid, trmst)
          CALL broadcast(grid, txmst)
          CALL broadcast(grid, tzmst)

#endif

!#ifdef TRACERS_WATER
!          if (am_i_root()) trsist(:,:,n)=0.
!          CALL broadcast(grid, trsist)
!#endif

#if (defined TRACERS_OCEAN) && (defined TRACERS_ZEBRA)
        case ('zebraL')

           read(entry%trname(7:8),'(I2)') ll ! gets level from name
           do l=1,lmo
             do j=J_0,J_1
               do i=1,im
                 if (l.eq.ll .and. l.le.lmm(i,j)) then
                   trmo(i,j,l,n) = mo(i,j,l)*dxypo(j) ! set conc=1 for l=ll
                 else
                   trmo(i,j,l,n) = 0.   ! zero otherwise
                 end if
               end do
             end do
           end do

           if (am_i_root()) then
             do nst=1,nmst
               do l=1,lmst(nst)
                 if (l.eq.ll) then
                   trmst(nst,l,n)= mmst(nst,l)
                 else
                   trmst(nst,l,n)= 0.
                 endif
               end do
             end do
           end if

! set all gradients to zero initially
           txmo(:,:,:,n) = 0; tymo(:,:,:,n)= 0. ; tzmo(:,:,:,n)=0.
           if(use_qus==1) then
             txxmo(:,:,:,n) = 0; txymo(:,:,:,n)= 0. ; tzxmo(:,:,:,n)=0.
             tyymo(:,:,:,n) = 0; tyzmo(:,:,:,n)= 0. ; tzzmo(:,:,:,n)=0.
           endif
           if (am_i_root()) then
             txmst(:,:,n)=0. ; tzmst(:,:,n)=0.
          endif
          CALL broadcast(grid, trmst)
          CALL broadcast(grid, txmst)
          CALL broadcast(grid, tzmst)
#endif

        case ('Water', 'H2O18', 'HDO', 'H2O17' )
#if (defined TRACERS_OCEAN) && (defined TRACERS_SPECIAL_O18)
C**** Open ic file for isotope tracers
          if(water_tracer_ic.eq.1)
     *         call openunit("H2O18ic",iu_O18ic,.true.,.true.)
#endif

#ifdef TRACERS_OCEAN
#ifndef TRACERS_SPECIAL_O18
C**** main ocean variabiles and gradients
          do j=J_0,J_1
          do i=1,im
            do l=1,lmm(i,j)
              trmo(i,j,l,n)=entry%trw0*(mo(i,j,l)*dxypo(j)-s0m(i,j,l))
              txmo(i,j,l,n)=-entry%trw0*sxmo(i,j,l)
              tymo(i,j,l,n)=-entry%trw0*symo(i,j,l)
              tzmo(i,j,l,n)=-entry%trw0*szmo(i,j,l)
            end do
          end do
          end do
#else
C**** read in initial conditions for isotopes
C**** search through for correct tracer (since there is no guarantee
C**** that they will be in same order as tracer indices).
C**** data are now in 'per mil'
          if(water_tracer_ic.eq.1) then
            rewind (iu_O18ic)
 10         read  (iu_O18ic,err=800,end=810) title,t0m4,tzm4
            if (index(title,trim(entry%trname)).eq.0) goto 10
            write (6,*) 'Read from H2O18ic: ',title
            call closeunit(iu_O18ic)
          else
            t0m4(:,:,:)=0.
            tzm4(:,:,:)=0.
c            if(n.eq.n_Water) t0m4(:,:,:)=1.
          endif

C**** Turn per mil data into mass ratios (using current standard)
          t0m4(:,:,:)=(t0m4(:,:,:)*1d-3+1.)*entry%trw0
          tzm4(:,:,:)=0.    ! tzm4(:,:,:)*1d-3*entry%trw0 corrupted?
C****
          do l=1,lmo
            txmo(:,J_0:J_1,l,n) = 0.
            tymo(:,J_0:J_1,l,n) = 0.
C**** Define East-West horizontal gradients
            im1=im-1
            i=im
            do j=J_0S,J_1S
              do ip1=1,im
                if (lmm(i,j).ge.l) then
                  if (lmm(im1,j).ge.l) then
                    if (lmm(ip1,j).ge.l) then
                      txmo(i,j,l,n)=2.5d-1*(t0m4(ip1,j,l)-t0m4(im1,j,l))
                    else
                      txmo(i,j,l,n)=  5d-1*(t0m4(  i,j,l)-t0m4(im1,j,l))
                    end if
                  else
                    if (lmm(ip1,j).ge.l)
     *                   txmo(i,j,l,n)=5d-1*(t0m4(ip1,j,l)-t0m4(i,j,l))
                  end if
                end if
                im1=i
                i=ip1
              end do
            end do
C**** Define North-South horizontal gradients
            do j=J_0S,J_1S
              do i=1,im
                if (lmm(i,j).ge.l)  then
                  if (lmm(i,j-1).ge.l)  then
                    if (lmm(i,j+1).ge.l)  then
                      tymo(i,j,l,n)=2.5d-1*(t0m4(i,j+1,l)-t0m4(i,j-1,l))
                    else
                      tymo(i,j,l,n)=  5d-1*(t0m4(i,  j,l)-t0m4(i,j-1,l))
                    end if
                  else
                    if (lmm(i,j+1).ge.l)
     *                   tymo(i,j,l,n)=5d-1*(t0m4(i,j+1,l)-t0m4(i,j,l))
                  end if
                end if
              end do
            end do
C**** Multiply ratios by freshwater mass
            do j=J_0,J_1
            do i=1,im
              if (l.le.lmm(i,j)) then
                fwm = mo(i,j,l)*dxypo(j)-s0m(i,j,l)
                trmo(i,j,l,n)=t0m4(i,j,l)*fwm
                txmo(i,j,l,n)=txmo(i,j,l,n)*fwm-sxmo(i,j,l)*t0m4(i,j,l)
                tymo(i,j,l,n)=tymo(i,j,l,n)*fwm-symo(i,j,l)*t0m4(i,j,l)
                tzmo(i,j,l,n)=tzm4(i,j,l)  *fwm-szmo(i,j,l)*t0m4(i,j,l)
              else
                trmo(i,j,l,n)=0.
                txmo(i,j,l,n)=0.
                tymo(i,j,l,n)=0.
                tzmo(i,j,l,n)=0.
              end if
            end do
            end do
          end do
#endif

C**** Initiallise strait values based on adjacent ocean boxes
          !call gather_ocean(1)  ! mo,g0m,gx-zmo,s0m,sx-zmo,trmo,tx-zmo

          call pack_data(grid,trmo(:,:,:,n),trmo_glob)
          if (.not.glob_used) then
            call pack_data(grid,mo,mo_glob)
            call pack_data(grid,s0m,s0m_glob)
            glob_used=.true.
          endif
          if(am_I_root()) then 
          do nst=1,nmst
            i1=ist(nst,1)
            j1=jst(nst,1)
            i2=ist(nst,2)
            j2=jst(nst,2)
            do l=1,lmst(nst)
              t01=trmo_glob(i1,j1,l)/(mo_glob(i1,j1,l)*dxypo(j1)
     *             -s0m_glob(i1,j1,l))
              t02=trmo_glob(i2,j2,l)/(mo_glob(i2,j2,l)*dxypo(j2)
     *             -s0m_glob(i2,j2,l))
              trmst(l,nst,n) = 5d-1*(mmst(l,nst)-s0mst(l,nst))*(t01+t02)
              txmst(l,nst,n) = -5d-1*sxmst(l,nst)*(t01+t02)
              tzmst(l,nst,n) = -5d-1*szmst(l,nst)*(t01+t02)
            end do
            do l=lmst(nst)+1,lmo
              trmst(l,nst,n) = 0.
              txmst(l,nst,n) = 0.
              tzmst(l,nst,n) = 0.
            end do
!#ifdef TRACERS_WATER
!            trsist(n,1:2,nst) = entry%trw0*(msist(1,nst)*xsi(1:2)
!     *           -ssist(1:2,nst))
!            trsist(n,3:lmi,nst)=entry%trw0*(msist(2,nst)*xsi(3:lmi)
!     *           -ssist(3:lmi,nst))
!#endif
          end do
          end if

          call bcast_straits(.false.) ! bcst tracers

C**** Balance tracers so that average concentration is TRW0 
C**** or oc_tracer_mean
          
          CALL CONSERV_OTR(OTRACJ,N)
          CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

          if (AM_I_ROOT()) then
            if (oc_tracer_mean(n).ne.-999.) then
              tratio=entry%trw0*(oc_tracer_mean(n)*1d-3+1.)
            else
              tratio=entry%trw0
            end if
            
            if (entry%trname.eq.'Water') then
              wsum = trsum
              frac_tr = -999.
            else
              frac_tr = tratio/(trsum/wsum)
              write(6,*) "Average oceanic tracer concentration ",
     *           entry%trname,(trsum/(wsum*entry%trw0)-1d0)*1d3,frac_tr
            end if
          end if
          call broadcast(grid,  frac_tr)
         
          if (frac_tr.ne.-999.) then    
            do l=1,lmo
              do j=J_0,J_1
                do i=1,imaxj(j)
                  trmo(i,j,l,n) = trmo(i,j,l,n) * frac_tr
                  if (frac_tr.lt.1.) then
                    TXMO(I,J,L,n)=frac_tr*TXMO(I,J,L,n)
                    TYMO(I,J,L,n)=frac_tr*TYMO(I,J,L,n)
                    TZMO(I,J,L,n)=frac_tr*TZMO(I,J,L,n)
                  end if
                end do
              end do
            end do
C**** straits
            if (am_i_root()) then
              DO NST=1,NMST
                DO L=1,LMST(NST)
                  TRMST(L,NST,N)=frac_tr*TRMST(L,NST,N)
                  if (frac_tr.lt.1) then
                    TXMST(L,NST,N)=frac_tr*TXMST(L,NST,N)
                    TZMST(L,NST,N)=frac_tr*TZMST(L,NST,N)
                  end if
                END DO
              END DO
            end if
            CALL broadcast(grid, TRMST)
            CALL broadcast(grid, TXMST)
            CALL broadcast(grid, TZMST)

C**** Check
            CALL CONSERV_OTR(OTRACJ,N)
            CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

            if (AM_I_ROOT()) then
              tratio=(trsum/(wsum*entry%trw0)-1.)*1000.
              write(6,*) "New ocean tracer mean: ",tratio
     *             ,oc_tracer_mean(n)
            end if
          end if
#endif
C****
        end select
        write(6,*) entry%trname," tracer initialised in ocean"
        end if
      end do

C**** ensure that atmospheric arrays are properly updated (i.e. gtracer)
      ! call now made at higher level
      !CALL TOC2SST(atmocn)

      if (fid>=0) call par_close(grid,fid)
      return
 800  write(6,*) "Error reading input file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
 810  write(6,*) "Tracer ",entry%trname," not found in file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
      end subroutine tracer_ic_ocean
#endif

#ifdef TRACERS_OCEAN
      SUBROUTINE OC_TDECAY(DTS)
!@sum OC_TDECAY decays radioactive tracers in ocean
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only:tracerlist, ocn_tracer_entry,expDecayRate
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      IMPLICIT NONE
      real*8, intent(in) :: dts
      integer n
      type(ocn_tracer_entry), pointer :: entry

      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        if (entry%from_file) cycle
        if (entry%trdecay.gt.0. .and. itime.ge.entry%itime_tr0) then
C**** Oceanic decay
          trmo(:,:,:,n)   = expDecayRate(n)*trmo(:,:,:,n)
          txmo(:,:,:,n)   = expDecayRate(n)*txmo(:,:,:,n)
          tymo(:,:,:,n)   = expDecayRate(n)*tymo(:,:,:,n)
          tzmo(:,:,:,n)   = expDecayRate(n)*tzmo(:,:,:,n)
        end if
      end do
C****
      return
      end subroutine oc_tdecay

      SUBROUTINE OCN_TR_AGE(DTS)
!@sum OCN_TR_AGE age tracers in ocean
!@auth Gavin Schmidt/Natassa Romanou
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_age
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo=>motr, imaxj,
     *     lmm, lmo,use_qus,txxmo,txymo,tzxmo,tyymo,tyzmo,tzzmo

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 age_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=0 and add 1 below
C**** this is mass*age (kg*year)
C**** age=1/(INT_DAYS_PER_YEAR*24*3600) in years
      age_inc=dts/(INT_DAYS_PER_YEAR*SECONDS_PER_DAY)
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,1,n_age)=0 ; TXMO(I,J,1,n_age)=0 
                TYMO(I,J,1,n_age)=0 ; TZMO(I,J,1,n_age)=0
                if(use_qus==1) then
                   TXXMO(I,J,1,n_age)=0 ; TXYMO(I,J,1,n_age)=0
                   TZXMO(I,J,1,n_age)=0 ; TYYMO(I,J,1,n_age)=0
                   TYZMO(I,J,1,n_age)=0 ; TZZMO(I,J,1,n_age)=0
                endif
!all nine set to zero
              else
                TRMO(I,J,L,n_age)= TRMO(I,J,L,n_age) +
     +                age_inc * MO(I,J,L) * oXYP(I,J)
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_age
#endif

      SUBROUTINE OCN_TR_CFC(DTS,icfc)
!@sum OCN_TR_CFC tracer in ocean
!@auth Natassa Romanou
      USE Dictionary_mod, only : get_param
      USE MODEL_COM, only : itime,modelEclock
      USE CONSTANT,   only : grav
      USE OCN_TRACER_COM, only : n_ocfc11,n_ocfc12,n_sf6
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo=>motr,
     &     t3d,s3d,r3d, imaxj, focean,
     *     lmm, lmo,dxypo,olat=>olat2d_dg ! 2D array containing lat at each i,j
     *    ,use_qus,txxmo,txymo,tzxmo,tyymo,tyzmo,tzzmo
      USE OFLUXES,    only : oRSI,oAPRESS,ocnatm
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid
      USE ODIAG, only : ij_cfcair,ij_kw,ij_csat,ij_cfcflux,ij_cfcsolub
     .        ,oij=>oij_loc
     .        ,ij_cfc12air,ij_kw12,ij_csat12,ij_cfc12flux,ij_cfc12solub
     .        ,ij_sf6air,ij_kw_sf6,ij_csat_sf6,ij_sf6flux,ij_sf6solub
      use model_com, only: modeleclock
      use runtimecontrols_mod, only: ocn_cfc

      IMPLICIT NONE
      interface
        subroutine read_atmcfc(iyear,nh1,sh1,nh2,sh2,nh3,sh3)
          real*8, allocatable, dimension(:), intent(out) :: 
     .            nh1,sh1,nh2,sh2,nh3,sh3
          integer, allocatable, dimension(:), intent(out) :: iyear
        end subroutine read_atmcfc
      end interface
      real*8, allocatable, dimension(:), save :: cfc11nh,cfc11sh
      real*8, allocatable, dimension(:), save :: cfc12nh,cfc12sh
      real*8, allocatable, dimension(:), save :: sf6nh,sf6sh
      integer, allocatable, dimension(:), save :: icfcyear
      real*8, intent(in) :: dts
      real*8 :: cfc_inc, Xconv,a,pres,sst,sss,temgs,wind,pnoice,Xkw
     .              ,solub,solub_cfc_sf6,schmidtno_cfc,Sc,kw,cfcair,csat
     .              ,fluxa,flux,flux_tendency,rho_water,dp1d
     .              ,Pnorth,Psouth,trmopro,fluxb,scsf6,dtr,ftr
      real*8 :: ys ! northern boundary of SH constant-value domain (deg N)
      real*8 :: yn ! southern boundary of NH constant-value domain (deg N)
      real*8 :: wt_sh ! weight for SH constant-value domain
      integer i,j,l,k,icfc,trac_ind
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1,year, month, dayOfYear, date
      real*8 :: cfc_conc_const,MW_gas

      if (.not.allocated(icfcyear)) then
        call read_atmcfc(icfcyear,cfc11nh,cfc11sh
     .                           ,cfc12nh,cfc12sh
     .                           ,sf6nh,sf6sh)
      end if

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** interpolate atmospheric values
      call modelEclock%get(year=year, month=month, date=date,
     &     dayOfYear=dayOfYear)

      if (year<icfcyear(1).or.year>icfcyear(105)) return

! set the right tracer index
      if (icfc==11) trac_ind=n_ocfc11
      if (icfc==12) trac_ind=n_ocfc12
      if (icfc==6)  trac_ind=n_sf6
      
      !pick appropriate values for each year
      do i=1,105
        if (year == icfcyear(i)) then
         yn=10.d0
         ys=-10.d0
         if (icfc==11) then
             Pnorth=cfc11nh(i)
             Psouth=cfc11sh(i)
         endif
         if (icfc==12) then
             Pnorth=cfc12nh(i)
             Psouth=cfc12sh(i)
         endif
         if (icfc==6) then
             Pnorth=sf6nh(i)
             Psouth=sf6sh(i)
         endif
        endif
      enddo


C**** at each time step set surface tracer conc=1+flux from atmos
      Xconv = 1/3.6d5
      a = 0.337
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
      IF (focean(i,j)>0) then
      wind=ocnatm%wsavg(i,j)        !owind(i,j)
      pnoice = 1.d0 - oRSI(i,j)     !1-fice
      k = 1    !surface only
        sst = t3d(1,i,j)     !in situ   temperature
        sss = s3d(1,i,j)*1000.d0      !convert to psu (eg. ocean mean salinity=35psu)
        rho_water = r3d(1,i,j)
        dp1d = MO(I,J,K)/rho_water   !local thickenss of each layer in meters

      Xkw = Xconv * a * wind**2       ! units in m/s
      if (icfc==11) solub = solub_cfc_sf6(sst,sss,11)   !mol/m3/pptv
      if (icfc==12) solub = solub_cfc_sf6(sst,sss,12)   !mol/m3/pptv
      if (icfc==6)  solub = solub_cfc_sf6(sst,sss,6)    !mol/m3/pptv

      if (icfc==11) Sc = schmidtno_cfc(sst,11)
      if (icfc==12) Sc = schmidtno_cfc(sst,12)
      if (icfc==6) Sc = scsf6(sst)
      kw = Xkw/sqrt(Sc/660)

#ifdef OCN_CFCconst
!     cfcair = 1.          !pptv to derive greens functions -- corresponds to year=1951
      call get_param('cfc_conc_const',cfc_conc_const)
      cfcair=cfc_conc_const
#else
      if(olat(i,j) <= ys) then
          wt_sh = 1d0
      elseif(olat(i,j) >= yn) then
          wt_sh = 0d0
      else
          wt_sh = (yn-olat(i,j))/(yn-ys)
      endif
      cfcair = (1.-wt_sh) * Pnorth + wt_sh * Psouth
#endif
      !molecular weights of gases
      if (icfc==11) MW_gas = 137.37d0    !Use values from Sarmiento book
      if (icfc==12) MW_gas = 120.90d0
      if (icfc==6)  MW_gas = 146.10d0

!mo units: kg/m2
!     pres = 1. !atm
      pres = oAPRESS(i,j)    !surface atm. pressure ANOMALY in Pa
      pres = (pres + 1013.25d0*100.d0)  !in Pa
     .     * 9.86923266716e-6             !atm
      csat = solub * cfcair * pres/1.d0       ! mol/m3
      fluxa = kw * csat                     ! mol/m2/s
      fluxb = kw * trmo(i,j,1,trac_ind) *1000.d0/MW_gas
     .           * rho_water/mo(i,j,1)/dxypo(j)            !mol/m2/s
      flux = fluxa - fluxb     !mol/m2/s

      flux_tendency= flux /dp1d * 1000.d0*pnoice  !mili-mol/m3/s

      trmopro=trmo(i,j,1,trac_ind)   !save for printout

      dtr = (flux_tendency * DTS)*1e-6*MW_gas
     .                            *mo(i,j,1)*dxypo(j)/rho_water

! make sure dtr+trmo should never be negative
      dtr = max(dtr,-1.d0*trmo(i,j,k,trac_ind))

        if (dtr.lt.0) then
          ftr = -dtr/trmo(i,j,k,trac_ind)
          txmo(i,j,k,trac_ind)=txmo(i,j,k,trac_ind)*(1.-ftr)
          tymo(i,j,k,trac_ind)=tymo(i,j,k,trac_ind)*(1.-ftr)
          tzmo(i,j,k,trac_ind)=tzmo(i,j,k,trac_ind)*(1.-ftr)
          if(use_qus==1) then
            txxmo(i,j,k,trac_ind)=txxmo(i,j,k,trac_ind)*(1.-ftr)
            txymo(i,j,k,trac_ind)=txymo(i,j,k,trac_ind)*(1.-ftr)
            tzxmo(i,j,k,trac_ind)=tzxmo(i,j,k,trac_ind)*(1.-ftr)
            tyymo(i,j,k,trac_ind)=tyymo(i,j,k,trac_ind)*(1.-ftr)
            tyzmo(i,j,k,trac_ind)=tyzmo(i,j,k,trac_ind)*(1.-ftr)
            tzzmo(i,j,k,trac_ind)=tzzmo(i,j,k,trac_ind)*(1.-ftr)
          endif
        endif

      trmo(i,j,1,trac_ind) =  trmo(i,j,1,trac_ind) +dtr   !kg


      if (i.eq.50.and.j.eq.90) then
         write(*,'(a,6i5,12e12.4)')'CFC OUTPUT 1:',
     . date,month,year,i,j,icfc,pnoice,cfcair,pres,solub,csat,kw,fluxa,
     . fluxb,flux,dp1d,trmopro,trmo(i,j,1,trac_ind)
      endif


       if (ocn_cfc) then
        if (icfc==11) then
         OIJ(I,J,IJ_cfcair) = OIJ(I,J,IJ_cfcair) + cfcair
         OIJ(I,J,IJ_csat) = OIJ(I,J,IJ_csat) +  csat
         OIJ(I,J,IJ_kw) = OIJ(I,J,IJ_kw) +  kw
         OIJ(I,J,IJ_cfcsolub) = OIJ(I,J,IJ_cfcsolub) +  solub
         OIJ(I,J,IJ_cfcflux) = OIJ(I,J,IJ_cfcflux) +  flux
        endif
        if (icfc==12) then
         OIJ(I,J,IJ_cfc12air) = OIJ(I,J,IJ_cfc12air) + cfcair
         OIJ(I,J,IJ_csat12) = OIJ(I,J,IJ_csat) +  csat
         OIJ(I,J,IJ_kw12) = OIJ(I,J,IJ_kw) +  kw
         OIJ(I,J,IJ_cfc12solub) = OIJ(I,J,IJ_cfc12solub) +  solub
         OIJ(I,J,IJ_cfc12flux) = OIJ(I,J,IJ_cfc12flux) +  flux
        endif
        if (icfc==6) then
         OIJ(I,J,IJ_sf6air) = OIJ(I,J,IJ_sf6air) + cfcair
         OIJ(I,J,IJ_csat_sf6) = OIJ(I,J,IJ_csat_sf6) +  csat
         OIJ(I,J,IJ_kw_sf6) = OIJ(I,J,IJ_kw_sf6) +  kw
         OIJ(I,J,IJ_sf6solub) = OIJ(I,J,IJ_sf6solub) +  solub
         OIJ(I,J,IJ_sf6flux) = OIJ(I,J,IJ_sf6flux) +  flux
        endif
      endif
      ENDIF
      ENDDO
      ENDDO

      end SUBROUTINE OCN_TR_CFC


      SUBROUTINE read_atmcfc(icfcyear,cfc11nh,cfc11sh
     .                               ,cfc12nh,cfc12sh
     .                               ,sf6nh,sf6sh)
  
      USE FILEMANAGER, only: openunit,closeunit

      implicit none
      real*8, allocatable, dimension(:), intent(out) :: cfc11nh, cfc11sh
      real*8, allocatable, dimension(:), intent(out) :: cfc12nh, cfc12sh
      real*8, allocatable, dimension(:), intent(out) :: sf6nh,  sf6sh
      real*8 :: dummy
      integer, allocatable, dimension(:), intent(out) :: icfcyear
      integer iu_file,i
      character(len=80) :: first_line_dummy

      allocate(icfcyear(105),cfc11nh(105), cfc11sh(105))
      allocate(cfc12nh(105), cfc12sh(105))
      allocate(sf6nh(105), sf6sh(105))

      !read in the values from cfc1112.atm
      call openunit('cfcatm_data',iu_file,.false.,.true.)
      read(iu_file,*)
      do i=1,105
      read(iu_file,*)icfcyear(i), cfc11nh(i), cfc11sh(i)
     .                          , cfc12nh(i), cfc12sh(i)
     .                          , dummy, dummy
     .                          , dummy, dummy
     .                          , sf6nh(i), sf6sh(i)
      enddo
      call closeunit(iu_file)

      end SUBROUTINE read_atmcfc

      real*8 FUNCTION scsf6(temp)
!>    Compute Schmidt number for SF6 in seawater from temperature
!     (from Orr et al, 2016 mocsy codes)

!  Compute Schmidt number of SF6 in seawater w/ formulation from Wanninkhof (Limnol. Oceanogr.: Methods 12, 2014, 351–362)
!  Input is temperature in deg C.

      IMPLICIT NONE

      !  Input & output variables:
      real*8, INTENT(in) :: temp
!     real*8 :: scsf6

      scsf6 = 3177.5 - 200.57*temp + 6.8865*temp**2 - 
     .                 0.13335*temp**3 + 0.0010877*temp**4

      RETURN
      END FUNCTION scsf6


      real*8 function solub_cfc_sf6(pt,ps,kn)

!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /www/html/ipsl/OCMIP/phase2/simulations/CFC/boundcond/RCS/sol_cfc.f$
!_ $Revision: 1.2 $
!_ $Date: 1998/07/17 07:37:02 $   ;  $State: Exp $
!_ $Author: jomce $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: sol_cfc.f,v $
!_ Revision 1.2  1998/07/17 07:37:02  jomce
!_ Fixed slight bug in units conversion: converted 1.0*e-12 to 1.0e-12
!_ following warning from Matthew Hecht at NCAR.
!_
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!     function sol_cfc_sf6=solub_cfc_sf6(pt,ps,kn)
!-------------------------------------------------------------------
!
!     CFC 11 and 12 Solubilities in seawater
!     SF6           Solubility in seawater
!     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!
!     pt:       temperature (degree Celcius)
!     ps:       salinity    (o/oo)
!     kn:       11 = CFC-11, 12 = CFC-12, 6 = SF6
!     sol_cfc:  in mol/m3/pptv
!               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm

!
!     J-C Dutay - LSCE
!-------------------------------------------------------------------
      use constant, only : tf
      implicit none
      real*8 :: a1,a2,a3,a4,b1,b2,b3,ta,d,pt,ps
      integer kn

!c coefficient for solubility in  mol/l/atm
!c ----------------------------------------
!
!     for CFC 11
!     ----------
      if (kn.eq.11) then
            a1 = -229.9261d0
            a2 =  319.6552d0
            a3 =  119.4471d0
            a4 =   -1.39165d0
            b1 =   -0.142382d0
            b2 =    0.091459d0
            b3 =   -0.0157274d0
      endif
!    
!     for CFC/12
!     ----------
      if (kn.eq.12) then
             a1 = -218.0971d0
             a2 =  298.9702d0
             a3 =  113.8049d0
             a4 =   -1.39165d0
             b1 =   -0.143566d0
             b2 =    0.091015d0
             b3 =   -0.0153924d0
       endif
 
!    
!     for SF6   OCMIP2016 protocol
!     ----------
      if (kn.eq.6) then
             a1 = -80.0343d0
             a2 =  117.232d0
             a3 =  29.5817d0
             a4 = -0.0d0
             b1 =  0.0335183d0
             b2 = -0.0373942d0
             b3 =  0.0048472d0
       endif


      ta       = ( pt + tf)* 0.01d0
      d    = (b3*ta + b2)*ta + b1
 
 
      solub_cfc_sf6 = exp(a1 + a2/ta + a3*log(ta) + a4*ta*ta + ps*d)

!     conversion from mol/(l * atm) to mol/(m^3 * atm) 
!     ------------------------------------------------
      solub_cfc_sf6 = 1000.d0 * solub_cfc_sf6
 
!     conversion from mol/(m^3 * atm) to mol/(m3 * pptv) 
!     --------------------------------------------------
      solub_cfc_sf6 = 1.0d-12 * solub_cfc_sf6

      end function solub_cfc_sf6


      real*8 function schmidtno_cfc(t,kn)

!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /home/orr/ocmip/web/OCMIP/phase2/simulations/CFC/boundcond/RCS/sc_c$
!_ $Revision: 1.1 $
!_ $Date: 1998/07/07 15:22:00 $   ;  $State: Exp $
!_ $Author: orr $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: sc_cfc.f,v $
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!---------------------------------------------------
!     CFC 11 and 12 schmidt number 
!     as a function of temperature. 
!
!     ref: Zheng et al (1998), JGR, vol 103,No C1 
!
!     t: temperature (degree Celcius)
!     kn: = 11 for CFC-11,  12 for CFC-12
!
!     J-C Dutay - LSCE
!---------------------------------------------------
!
!   coefficients with t in degree Celcius
!   ------------------------------------
      implicit none
      real*8 :: a1,a2,a3,a4,a5,t
      integer kn

      if (kn==11) then
      !a1 = 3501.8d0
      !a2 = -210.31d0
      !a3 =    6.1851d0
      !a4 =   -0.07513d0
      !new parameterizations based on Wanninkhof (2014)
      a1 = 3579.2d0
      a2 = -222.63d0
      a3 = 7.5749d0
      a4 = -0.14595d0
      a5 = 0.0011874d0
      endif
 
      if (kn==12) then
      !a1 = 3845.4d0
      !a2 = -228.95d0
      !a3 =    6.1908d0
      !a4 =   -0.067430d0
      a1 = 3828.1
      a2 = -249.86
      a3 = 8.7603
      a4 = -0.1716
      a5 = 0.001408
      endif

      schmidtno_cfc = a1 + a2 * t + a3 *t*t + a4 *t*t*t + a5 *t*t*t*t;

      end function schmidtno_cfc
!---------------------------------------------------------------------------------------

      SUBROUTINE DIAGTCO (M,NT0,atmocn)
!@sum  DIAGTCO Keeps track of the conservation properties of ocean tracers
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
      USE OCEAN, only : IMO=>IM, oJ_BUDG, oJ_0B, oJ_1B
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds
      USE OCEANR_DIM, only : oGRID
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT0 index denoting tracer number, NT is index in tconsrv 
      INTEGER, INTENT(IN) :: nt0
!@var atmocn
      type(atmocn_xchng_vars) :: atmocn
c
      INTEGER :: nt
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(IMO,
     &                  oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: TOTAL
      REAL*8, DIMENSION(oJ_0B:oJ_1B) :: TOTALJ
      INTEGER :: nm,ni
      INTEGER :: I, J, J_0, J_1, I_0, I_1

      call getDomainBounds(ogrid, J_STRT=J_0, J_STOP=J_1)
      I_0 = oGRID%I_STRT
      I_1 = oGRID%I_STOP
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.
C**** NOFMT(1,NT) is the index for the instantaneous value.
      nt=nt0+atmocn%natmtrcons
      if (atmocn%nofmt(m,nt).gt.0) then
C**** Calculate current value TOTAL (kg)
        call conserv_otr(total,nt0)
        nm=atmocn%nofmt(m,nt)
        ni=atmocn%nofmt(1,nt)

C**** Calculate zonal sums
        totalj(oj_0b:oj_1b)=0.
        do j=j_0,j_1
        do i=i_0,i_1
          totalj(oj_budg(i,j)) = totalj(oj_budg(i,j)) + total(i,j)
        end do
        end do

c**** Accumulate difference from last time in TCONSRV(NM)
        if (m.gt.1) then
          atmocn%tconsrv(oJ_0b:oJ_1b,nm,nt) =
     &         atmocn%tconsrv(oJ_0b:oJ_1b,nm,nt)
     &         +(totalj(oJ_0b:oJ_1b)-atmocn%tconsrv(oJ_0b:oJ_1b,ni,nt)) 
        end if
C**** Save current value in TCONSRV(NI)
        atmocn%tconsrv(oJ_0b:oJ_1b,ni,nt)=totalj(oJ_0b:oJ_1b)
      end if
      return
      end subroutine diagtco

      SUBROUTINE conserv_OTR(OTR,ITR)
!@sum  conserv_OTR calculates zonal ocean tracer (kg) on ocean grid
!@auth Gavin Schmidt
      USE OCEAN, only : IMO=>IM,JMO=>JM,oXYP,LMM, imaxj,trmo
      USE STRAITS, only : nmst,jst,ist,lmst,trmst
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds
      USE OCEANR_DIM, only : ogrid

      IMPLICIT NONE
!@var OTR zonal ocean tracer (kg)
      REAL*8, DIMENSION(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) :: OTR 
      INTEGER, INTENT(IN) :: ITR
      INTEGER I,J,L,N

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds(ogrid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          OTR(I,J) = SUM(TRMO(I,J,:LMM(I,J),ITR))
        END DO
      END DO
      if (HAVE_SOUTH_POLE) OTR(2:IMO,1)   = OTR(1,1)
      if (HAVE_NORTH_POLE) OTR(2:IMO,JMO) = OTR(1,JMO)
C**** include straits variables
      DO N=1,NMST
        I=IST(N,1)
        J=JST(N,1)
        If (J >= J_0 .and. J <= J_1)
     &       OTR(I,J)=OTR(I,J)+SUM(TRMST(:LMST(N),N,ITR))
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OTR

      SUBROUTINE OCN_TR_VENT(DTS)
!@sum OCN_VENT tracer in ocean
!@auth Natassa Romanou
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_vent
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo=>motr,
     &     imaxj, focean,
     *     lmm, lmo,use_qus,txxmo,txymo,tzxmo,tyymo,tyzmo,tzzmo

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 vent_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=1
      vent_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,L,n_vent)= vent_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_vent)= 0.
                TYMO(I,J,L,n_vent)= 0.
                TZMO(I,J,L,n_vent)= 0.
                if(use_qus==1) then
                   TXXMO(I,J,L,n_vent)= 0.
                   TXYMO(I,J,L,n_vent)= 0.
                   TZXMO(I,J,L,n_vent)= 0.
                   TYYMO(I,J,L,n_vent)= 0.
                   TYZMO(I,J,L,n_vent)= 0.
                   TZZMO(I,J,L,n_vent)= 0.
                endif
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_vent

      SUBROUTINE OCN_TR_GASX(DTS)
!@sum OCN_GASX tracer in ocean
!@auth Natassa Romanou
! same as ventilation tracer but take ice into account
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_gasx
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo=>motr,
     &     imaxj, focean,
     *     lmm, lmo,use_qus,txxmo,tyymo,tzzmo,txymo,tyzmo,tzxmo
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid
      USE OFLUXES, only : oRSI

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 gasx_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

      gasx_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
#ifdef OTRAC_unspecified
                !under ice do not specify anything
                if (oRSI(i,j).ne.1.) then
                TRMO(I,J,L,n_gasx)= gasx_inc * MO(I,J,L) * oXYP(I,J)
     .                            * (1-oRSI(i,j))
                else
                write(*,'(a,3i5,2e12.4)')'ice frac=1',
     .          itime,i,j,orsi(i,j),trmo(i,j,l,n_gasx)
                endif
#else
        !in this formulation we specify 0 under the ice
                TRMO(I,J,L,n_gasx)= gasx_inc * MO(I,J,L) * oXYP(I,J)
     .                            * (1-oRSI(i,j))
                TXMO(I,J,L,n_gasx)= 0.
                TYMO(I,J,L,n_gasx)= 0.
                TZMO(I,J,L,n_gasx)= 0.
                if(use_qus==1) then
                   TXXMO(I,J,L,n_gasx)= 0.
                   TXYMO(I,J,L,n_gasx)= 0.
                   TZXMO(I,J,L,n_gasx)= 0.
                   TYYMO(I,J,L,n_gasx)= 0.
                   TYZMO(I,J,L,n_gasx)= 0.
                   TZZMO(I,J,L,n_gasx)= 0.
                endif
#endif
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
C****
      return
      end subroutine ocn_tr_gasx

      SUBROUTINE OCN_TR_WaterMass(DTS)
!@sum OCN_WaterMass tracer in ocean
!@auth Natassa Romanou
      USE MODEL_COM, only : itime
      USE OCN_TRACER_COM, only : n_wms1,n_wms2,n_wms3
      USE OCEAN, only : trmo,txmo,tymo,tzmo, oxyp, mo=>motr,
     &     imaxj, focean,
     *     lmm, lmo, oLON_DG,oLAT_DG,ZOE=>ZE,use_qus
     *     ,txxmo,txymo,tzxmo,tyymo,tyzmo,tzzmo


      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      real*8, intent(in) :: dts
      real*8 wms1_inc,wms2_inc,wms3_inc
      integer i,j,l
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)

C**** at each time step set surface tracer conc=1
C**** this is mass*age (kg*year)
C**** age=1/(JDperY*24*3600) in years
      if (n_wms1 /= 0) then
      wms1_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          if (oLAT_DG(j,2) .le. -55.d0) then    !Antarctic
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,L,n_wms1)= wms1_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_wms1)= 0.
                TYMO(I,J,L,n_wms1)= 0.
                TZMO(I,J,L,n_wms1)= 0.
                if(use_qus==1) then
                   TXXMO(I,J,L,n_wms1)= 0.
                   TXYMO(I,J,L,n_wms1)= 0.
                   TZXMO(I,J,L,n_wms1)= 0.
                   TYYMO(I,J,L,n_wms1)= 0.
                   TYZMO(I,J,L,n_wms1)= 0.
                   TZZMO(I,J,L,n_wms1)= 0.
                endif
              end if
            end if
          ENDDO
          endif
        ENDDO
      ENDDO
      endif

      if (n_wms2 /= 0) then
      wms2_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          if (oLAT_DG(j,2) .ge. 40.d0) then    !North Atlantic
          DO I=1,IMAXJ(J)
          if (oLON_DG(i,2) .ge. -90.d0 .and. oLON_DG(i,2) .le. 60.d0)
     .      then    !North Atlantic
            if (l.le.lmm(i,j)) then
              if (L.eq.1) then
                TRMO(I,J,L,n_wms2)= wms2_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_wms2)= 0.
                TYMO(I,J,L,n_wms2)= 0.
                TZMO(I,J,L,n_wms2)= 0.
                if(use_qus==1) then
                   TXXMO(I,J,L,n_wms2)= 0.
                   TXYMO(I,J,L,n_wms2)= 0.
                   TZXMO(I,J,L,n_wms2)= 0.
                   TYYMO(I,J,L,n_wms2)= 0.
                   TYZMO(I,J,L,n_wms2)= 0.
                   TZZMO(I,J,L,n_wms2)= 0.
                endif
              end if
            end if
          end if
          ENDDO
          endif
        ENDDO
      ENDDO
      endif

      if (n_wms3 /= 0) then
      wms3_inc=1.d0
      DO L=1,LMO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            if (l.le.lmm(i,j)) then
              if (zoe(l).ge.3000.) then
                TRMO(I,J,L,n_wms3)= wms3_inc * MO(I,J,L) * oXYP(I,J)
                TXMO(I,J,L,n_wms3)= 0.
                TYMO(I,J,L,n_wms3)= 0.
                TZMO(I,J,L,n_wms3)= 0.
                if(use_qus==1) then
                   TXXMO(I,J,L,n_wms3)= 0.
                   TXYMO(I,J,L,n_wms3)= 0.
                   TZXMO(I,J,L,n_wms3)= 0.
                   TYYMO(I,J,L,n_wms3)= 0.
                   TYZMO(I,J,L,n_wms3)= 0.
                   TZZMO(I,J,L,n_wms3)= 0.
                endif
              end if
            end if
          ENDDO
        ENDDO
      ENDDO
      endif
C****
      return
      end subroutine ocn_tr_WaterMass

