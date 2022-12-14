#include "rundeck_opts.h"
      SUBROUTINE cheminit
!@sum cheminit initialize model chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls jplrts,fastj2_init,reactn

C**** GLOBAL parameters and variables:
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      USE MODEL_COM, only: Itime, ItimeI
      USE DOMAIN_DECOMP_ATM, only: getDomainBounds,grid,readt_parallel
      USE TRACER_COM, only: oh_live,no3_live
      USE TRCHEM_Shindell_COM, only: nc
      USE TRCHEM_Shindell_COM, only:
     &    prnls,prnrts,prnchg,lprn,jprn,iprn,ay,pHOx,pOx,pNOx,
     &    yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yNO3,yRXPAR,yXO2N,
#ifdef TRACERS_dCO
     &    ydC217O3,ydC218O3,yd13C2O3,
     &    yd13CXPAR,
     &    yd17OROR,yd18OROR,yd13CROR,
     &    yd17Oald,yd18Oald,yd13Cald,
     &    ydCH317O2,ydCH318O2,yd13CH3O2,
#endif  /* TRACERS_dCO */
     &    allowSomeChemReinit,pNO3,topLevelOfChemistry,nfam,ny
     &    ,pCLOx,pCLx,pOClOx,pBrOx,yCl2,yCl2O2,mostRecentNonZeroAlbedo

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var iu_data temporary unit number
!@var i,l loop dummy
      character(len=22) :: format_110
      integer :: iu_data,i,L,j,nc_10,nc_mod
      integer :: J_0,J_1,J_0S,J_1S,J_1H,J_0H,I_0,I_1
         
      call getDomainBounds(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ! Note that topLevelOfChemistry is set in 
      ! alloc_trchem_shindell_com routine

! define MOLEC reading format
 100  format(/3(50x,l1/),3(50x,i8/))
      nc_10=floor(dble(nc+1)/10.d0)
      if (nc_10 > 99) then
        ! stop the model if the format becomes invalid
        call stop_model('ERROR: Too many lines in MOLEC',255)
      endif
      nc_mod=mod(nc,10)
      write (format_110,'(a2,i1,a14,i1,a4)')
     &  '(',nc_10,'(///10a8),(///',nc_mod,'a8))'

C Read chem diagnostics parameters and molecule names
C from MOLEC file:
      call openunit('MOLEC',iu_data,.false.,.true.)
      read(iu_data,100)prnls,prnrts,prnchg,lprn,jprn,iprn
      read(iu_data,trim(format_110))ay
      call closeunit(iu_data)

! figure out first element of each family
      do i=1,ny
        select case(ay(i))
          case ('O3')    ; nfam(1)=i
          case ('NO')    ; nfam(2)=i
          case ('Cl2O2') ; nfam(3)=i
          case ('BrO')   ; nfam(4)=i
        end select
      enddo

C Read JPL chemical reactions/rates from unit JPLRX:
      call jplrts

! initialize fastj
      call fastj2_init

c Set up arrays of reaction numbers involving each molecule:
      call reactn

C Initialize a few (IM,JM,topLevelOfChemistry) arrays, first hour only:
      IF(Itime == ItimeI .and. allowSomeChemReinit == 1) THEN
        ! allowSomeChemReinit condition b/c these are in RSF files:
        pHOx(I_0:I_1,J_0:J_1,:)     =1.d0
        pOx(I_0:I_1,J_0:J_1,:)      =1.d0
        pNOx(I_0:I_1,J_0:J_1,:)     =1.d0
        pNO3(I_0:I_1,J_0:J_1,:)     =0.d0
        yCH3O2(I_0:I_1,J_0:J_1,:)   =1.d0
#ifdef TRACERS_dCO
        ydCH317O2(I_0:I_1,J_0:J_1,:)=1.d0
        ydCH318O2(I_0:I_1,J_0:J_1,:)=1.d0
        yd13CH3O2(I_0:I_1,J_0:J_1,:)=1.d0
#endif  /* TRACERS_dCO */
        yC2O3(I_0:I_1,J_0:J_1,:)    =0.d0
#ifdef TRACERS_dCO
        ydC217O3(I_0:I_1,J_0:J_1,:) =0.d0
        ydC218O3(I_0:I_1,J_0:J_1,:) =0.d0
        yd13C2O3(I_0:I_1,J_0:J_1,:) =0.d0
#endif  /* TRACERS_dCO */
        yROR(I_0:I_1,J_0:J_1,:)     =0.d0
#ifdef TRACERS_dCO
        yd17OROR(I_0:I_1,J_0:J_1,:) =0.d0
        yd18OROR(I_0:I_1,J_0:J_1,:) =0.d0
        yd13CROR(I_0:I_1,J_0:J_1,:) =0.d0
#endif  /* TRACERS_dCO */
        yXO2(I_0:I_1,J_0:J_1,:)     =0.d0
        yAldehyde(I_0:I_1,J_0:J_1,:)=0.d0
#ifdef TRACERS_dCO
        yd17Oald(I_0:I_1,J_0:J_1,:) =0.d0
        yd18Oald(I_0:I_1,J_0:J_1,:) =0.d0
        yd13Cald(I_0:I_1,J_0:J_1,:) =0.d0
#endif  /* TRACERS_dCO */
        yNO3(I_0:I_1,J_0:J_1,:)     =0.d0
        yXO2N(I_0:I_1,J_0:J_1,:)    =0.d0
        yRXPAR(I_0:I_1,J_0:J_1,:)   =0.d0
#ifdef TRACERS_dCO
        yd13CXPAR(I_0:I_1,J_0:J_1,:)=0.d0
#endif  /* TRACERS_dCO */
        oh_live(I_0:I_1,J_0:J_1,:)  =0.d0
        no3_live(I_0:I_1,J_0:J_1,:) =0.d0
        pClOx(I_0:I_1,J_0:J_1,:)    =1.d0
        pClx(I_0:I_1,J_0:J_1,:)     =0.d0
        pOClOx(I_0:I_1,J_0:J_1,:)   =0.d0
        pBrOx(I_0:I_1,J_0:J_1,:)    =1.d0
        yCl2(I_0:I_1,J_0:J_1,:)     =0.d0
        yCl2O2(I_0:I_1,J_0:J_1,:)   =0.d0
      END IF

      if(Itime == ItimeI)then
        ! First time only, read some albedo initial conditions (I,J)
        ! to be used only until rad code ALB(I,J,1) has filled in it's
        ! first non-zero values at each I,J. Array 
        ! mostRecentNonZeroAlbedo(I,J) is then saved to/read from restart
        ! files for use in rest of the run:
                       ! logicals mean: binary, old:
        call openunit('ALB_IC',iu_data,.true.,.true.)
        call readt_parallel(grid,iu_data,nameunit(iu_data),
     &  mostRecentNonZeroAlbedo,0)
        call closeunit(iu_data)
      end if

      return
      END SUBROUTINE cheminit



      SUBROUTINE jplrts
!@sum jplrts read/set up chemical reaction rates from JPL
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls lstnumc

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: pe,ea,nst,ro,
     &                               r1,sn,sb,nn,nnr,ay,rrtri
     &                              ,n_rx,n_bi,n_tri,n_nst,n_het

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
C
!@var ate temporary reactants names array
!@var i,ii,j dummy loop variable
!@var iu_data temporary unit number
      CHARACTER*8, DIMENSION(4) :: ate
      character(len=300) :: out_line
      INTEGER :: i,ii,j,iu_data,nr,nr2,nr3,nmm,nhet
      character(len=36) :: invreaction

C Read in the number of each type of reaction:
      call openunit('JPLRX',iu_data,.false.,.true.)
      read(iu_data,124)nr,nr2,nr3,nmm,nhet
      if (nr /= n_rx) then
        print*,'nr=',nr,' n_rx=',n_rx
        call stop_model('ERROR: nr (from JPLRX) /= n_rx '//
     &                  '(from TRCHEM_Shindell_COM)', 255)
      endif
      if (nr2 /= n_bi+n_nst) then
        print*,'nr2=',nr2,' n_bi=',n_bi,' n_nst=',n_nst
        call stop_model('ERROR: nr2 (from JPLRX) /= n_bi+n_nst '//
     &                  '(from TRCHEM_Shindell_COM)', 255)
      endif
      if (nr3 /= n_tri) then
        print*,'nr3=',nr3,' n_tri=',n_tri
        call stop_model('ERROR: nr3 (from JPLRX) /= n_tri '//
     &                  '(from TRCHEM_Shindell_COM)', 255)
      endif
      if (nmm /= n_nst) then
        print*,'nmm=',nmm,' n_nst=',n_nst
        call stop_model('ERROR: nmm (from JPLRX) /= n_nst '//
     &                  '(from TRCHEM_Shindell_COM)', 255)
      endif
      if (nhet /= n_het) then
        print*,'nhet=',nhet,' n_het=',n_het
        call stop_model('ERROR: nhet (from JPLRX) /= n_het '//
     &                  '(from TRCHEM_Shindell_COM)', 255)
      endif
      write(out_line,*)' '
      call write_parallel(trim(out_line))
      write(out_line,*) 'Chemical reactions used in the model: '
      call write_parallel(trim(out_line))

      do i=1,n_rx               ! >>> begin loop over total reactions <<<
        if(i <= n_rx-n_het) then !non-hetero
          if(i <= n_bi+n_nst) then   !mono or bi
            if(i == n_bi+1) read(iu_data,22)ate
            read(iu_data,16)ate,pe(i),ea(i) ! read mono and bimolecular reactions
          else                    ! read trimolecular reactions
 20         if(i == n_bi+n_nst+1) read(iu_data,22)ate
            ii=i-n_bi-n_nst
            read(iu_data,21)ate,ro(ii),sn(ii),r1(ii),sb(ii)
          end if
        else                     ! read heterogeneous reactions
          if(i == n_rx-(n_het-1)) read(iu_data,22)ate
          read(iu_data,31)ate
        end if ! (i <= n_rx-n_het)

        write(out_line,30) i,ate(1),' + ',ate(2),
     *  ' --> ',ate(3),' + ',ate(4)
        call write_parallel(trim(out_line))
c
        do j=1,2
          call lstnum(ate(j),nn(j,i))
          call lstnum(ate(j+2),nnr(j,i))
        end do

        call set_rrate_index(i, ate)
      end do                ! >>> end loop over total reactions <<<

! find the inverse reaction of a thermal decomposition
      do i=n_bi+1,n_bi+n_nst
        invreaction = trim(ay(nnr(1,i)))//'_'//
     &                trim(ay(nnr(2,i)))//'__'//
     &                trim(ay(nn(1,i)))//'_'//
     &                trim(ay(nn(2,i)))
        select case(invreaction)
          case('HO2_NO2__HO2NO2_M')
            nst(i-n_bi)=rrtri%HO2_NO2__HO2NO2_M
          case('NO3_NO2__N2O5_M')
            nst(i-n_bi)=rrtri%NO3_NO2__N2O5_M
          case('ClO_ClO__Cl2O2_M')
            nst(i-n_bi)=rrtri%ClO_ClO__Cl2O2_M
          case default
            call stop_model('ERROR: Reaction '//trim(invreaction)//
     &                      ' does not exist in the JPLRX file',255)
        end select
      enddo

 124  format(///5(/43x,i3)///)
  27  format(/(30x,i2))
  21  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f5.2,e9.2,f4.1)
  22  format(/10x,4a8/)
  25  format(//32x,2f7.1,i6)
  16  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f8.0)
  31  format(4x,a8,1x,a8,3x,a8,1x,a8)
  30  format(1x,i3,2x,a8,a3,a8,a5,a8,a3,a8)
      call closeunit(iu_data)
      return
      end SUBROUTINE jplrts



      SUBROUTINE lstnum(at,ks)
!@sum lstnum find molecule number in param list of molecules
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: nc,ay

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var at local copy of species name
!@var ks local variable to be passed back to jplrts nnr or nn array.
!@var j dummy loop variable
      INTEGER                  :: j
      INTEGER,     INTENT(OUT) :: ks
      CHARACTER*8, INTENT(IN)  :: at
      
      j=1
      do while(j <= nc)
        if(at == ay(j))then
          ks = j
          return
        else
          j = j + 1
          cycle
        endif
      enddo 
      ks = nc + 1
      if (at /= 'N2' .and. at /= 'H')
     &  call stop_model('ERROR: Tracer '//trim(at)//
     &    ' does not exist in the MOLEC file',255)

      return
      end SUBROUTINE lstnum



      subroutine fastj2_init
!@sum fastj2_init initialize fastj2 based on the currently active
!@+ chemistry scheme. It is the driver between the photolysis and
!@+ rest of chemistry code, so it can't be (in its current setup)
!@+ in the photolysis module
!@auth Kostas Tsigaridis

      use Dictionary_mod, only: sync_param
      use resolution, only: plbot, LM
      use TRCHEM_Shindell_COM, only: iprn,jprn,prnrts,n_rj
     &                              ,p_1,topLevelOfChemistry
      use photolysis, only: phtlst,inphot
     &                     ,j_iprn,j_jprn,j_prnrts,jppj,jlabel
     &                     ,jind,ks,kss,jfacta,zj 
     &                     ,aerosols_affect_photolysis
      ! also get things that are now allocatable based on top layer
      ! of chemistry and define here:
       use photolysis, only: NLGCM,ncfastj2,nbfastj,jpnl,jndlev,
     & pomegaj,fff,amf,tj2,do32,zfastj2,dmfastj2,tfastj,odcol,
     & pfastj2,o3_fastj,M__,nwfastj

      implicit none 

      NLGCM=topLevelOfChemistry
      call sync_param('aerosols_affect_photolysis',
     &                 aerosols_affect_photolysis)

      jpnl=NLGCM
      ! Stop the model if the pressure at the top of the top level of
      ! chemistry (or JPNL if someone someday sets that lower) would be
      ! 1 mb or greater (encroaching on the ozone layer). For example,
      ! see notes on the setting of colmO2 and colmO3 variables as a linear
      ! function of pressure in masterchem routine. (Exempt models with
      ! tops lower than that so that 12L models can still run):
      if(plbot(min(jpnl,topLevelOfChemistry)+1) >= 1.d0)then
        if(plbot(LM+1) < 1.d0)call stop_model
     &   ('jpnl or topLevelOfChemistry should be higher',255)
      end if
      ncfastj2=2*NLGCM+2
      nbfastj=NLGCM+1
      j_iprn=iprn
      j_jprn=jprn
      j_prnrts=prnrts
      jppj=n_rj ! needed for the photolysis routine
 
      allocate(jndlev(NLGCM))
      allocate(pomegaj(2*M__,2*NLGCM+2+1))
      allocate(fff(nwfastj,jpnl))
      allocate(amf(nbfastj,nbfastj))
      allocate(tj2(nbfastj))
      allocate(do32(nbfastj))
      allocate(zfastj2(nbfastj))
      allocate(dmfastj2(nbfastj))
      allocate(tfastj(NLGCM))
      allocate(odcol(NLGCM))
      allocate(pfastj2(NLGCM+3))
      allocate(o3_fastj(NLGCM)) ! until recently was 2*NLGCM
      allocate(jlabel(n_rj))
      allocate(jind(n_rj))
      allocate(ks(n_rj))
      allocate(kss(p_1,n_rj))
      allocate(jfacta(n_rj))
      allocate(zj(jpnl,n_rj))

C Read photolysis parameters and reactions from unit JPLPH:
      call phtlst

c fastj initialization routine:
      call inphot

      end subroutine fastj2_init


      SUBROUTINE reactn
!@sum reactn read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls guide,printls

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: nps,nds,kps,kds,nn,nnr,
     &                      npnr,ndnr,kpnr,kdnr,prnls,n_rx,n_rj
      use photolysis, only: ks,kss

      IMPLICIT NONE

c Chemical reaction lists:
      call guide(npnr,ndnr,kpnr,kdnr,nn,nnr,2,n_rx)
c Photolysis reaction lists:
      call guide( nps, nds, kps, kds,ks,kss,1,n_rj)
C Print out some diagnostics:
      if(prnls) call printls
      
      return
      end SUBROUTINE reactn



      SUBROUTINE guide(npr,ndr,kpr,kdr,nn,nnn,ns,nre)
!@sum guide read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls calcls

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: p_1,nc

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var nn   = either nn  or ks   from reactn sub
!@var nnn  = either nnr or kss  from reactn sub
!@var kpr  = either kps or kpnr from reactn sub
!@var kdr  = either kds or kdnr from reactn sub
!@var npr  = either nps or npnr from reactn sub
!@var ndr  = either nds or ndnr from reactn sub
!@var ns   = either 1   or    2 from reactn sub
!@var nre number of reactions
      INTEGER,  DIMENSION(nc)      :: kpr, kdr
      INTEGER,  DIMENSION(p_1*nre) :: npr, ndr
      INTEGER, DIMENSION(p_1,nre)  :: nn, nnn
      INTEGER                      :: ns, nre

c Chemical and photolytic destruction:
      call calcls(nn,ns,nnn,2,ndr,kdr,nre)
c Chemical and photolytic production:
      call calcls(nnn,2,nn,ns,npr,kpr,nre)
      
      return
      end SUBROUTINE guide



      SUBROUTINE calcls(nn,ns,nnn,nns,ndr,kdr,nre)
!@sum calcls Set up reaction lists for calculated gases (1 to ny)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: ny, numfam, p_1, nc, nfam,
     &                               prnls

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var kdr  = either kdr or kpr from guide sub
!@var ndr  = either ndr or npr from guide sub
!@var nre  number of reactions
!@var nns number of partic_ on opposite side of reaction
!@var ns   = either ns or   2 from guide sub
!@var nn   = either nn or nnn from guide sub
!@var nnn  = either nn or nnn from guide sub
!@var ii,k,j,i,ij,i2,newfam,ifam dummy variables
      INTEGER, DIMENSION(nc)     :: kdr
      INTEGER, DIMENSION(p_1*nre):: ndr
      INTEGER :: nre, nns, ns, k, j, i, ij, i2, newfam, ifam, ii
      INTEGER, DIMENSION(ns,nre) :: nn 
      INTEGER, DIMENSION(nns,nre):: nnn
      character(len=300) :: out_line

      k=1
      do j=1,numfam      !families, list only interfamily reactions
        newfam=0
        kdr(j)=k
        i_loop: do i=1,nre    ! 1 to # chem or phot reactions
          ij_loop: do ij=1,ns !ns # partic (prod & chem dest=2,phot dest=1)
            ! check if molecule # nn() is element of family j:
            ! if not, do nothing and move the next molecule in nn
            if(nn(ij,i) >= nfam(j).and.nn(ij,i) < nfam(j+1))then
              ! check if reaction is intrafamily:
              ! if yes, do nothing and move to next reaction
              ! this skips reactions with at least one reactant and one
              ! product of reaction i are from the same family j
              do i2=1,nns  ! nns # partic on opposite side of reac.
                if(nnn(i2,i) >= nfam(j).and.nnn(i2,i) < nfam(j+1))
     &          cycle i_loop
              enddo
              ! don't write same reaction twice:
              if(k /= 1)then
                if(ndr(k-1) == i.and.newfam /= 0) cycle ij_loop
              endif
              ndr(k)=i
              k=k+1
              newfam=1
            endif
          enddo ij_loop
        enddo i_loop
      enddo

      do j=numfam+1,nfam(1)-1     ! individual non-family molecules
        kdr(j)=k
        do i=1,nre                ! 1 to # chem or phot reactions
          do ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i) /= j) cycle ! nn is mol # of participant
            ndr(k)=i
            k=k+1
          enddo
        enddo
      enddo

      do 100 j=nfam(1),ny !indiv family mols.,list only intrafamily
        do ii=1,numfam-1
          if(j < nfam(ii+1))then
            ifam=ii
            goto 110
          endif
        enddo
        ifam=numfam
 110    kdr(j)=k
        do 100 i=1,nre          ! 1 to # chem or phot reactions
          do 100 ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i) /= j)goto100       ! nn is mol # of participant
c           check that reaction is intrafamily
            do i2=1,nns  ! nns # participants on opposite side of reac.
              if(nnn(i2,i) >= nfam(ifam).and.nnn(i2,i) < nfam(ifam+1))
     &        then
                 ndr(k)=i
                 k=k+1
                 goto 100
              endif
            enddo
 100  continue
      kdr(ny+1)=k

      if(prnls)then
        write(*,*) 'nn array size :',k
        write(out_line,*) 'nn array size :',k
        call write_parallel(trim(out_line))
      endif
      
      return
      end SUBROUTINE calcls



      SUBROUTINE printls
!@sum printls print out some chemistry diagnostics (reaction lists)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: kpnr,npnr,kdnr,ndnr,kps,nps,
     &                         ny,nn,nnr,ay,kds,nds,nc
      use photolysis, only: ks,kss

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var ireac,igas,ichange,ii dummy variables
      INTEGER :: ireac,igas,ichange,ii
      character(len=300) :: out_line

c Print reaction lists:
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ CHEMICAL PRODUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kpnr(igas+1)-kpnr(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            if (nnr(2,npnr(ireac)) > nc) then
              write(out_line,20)
     &        ' Reaction # ',npnr(ireac),' produces ',
     &        ay(nnr(1,npnr(ireac))),' and  ','X'
              call write_parallel(trim(out_line))
            else
              write(out_line,20)
     &        ' Reaction # ',npnr(ireac),' produces ',
     &        ay(nnr(1,npnr(ireac))),' and  ',ay(nnr(2,npnr(ireac)))
              call write_parallel(trim(out_line))
            end if
          enddo
        end if
      end do
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ CHEMICAL DESTRUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kdnr(igas+1)-kdnr(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,20)
     &      ' Reaction # ',ndnr(ireac),' destroys ',
     *      ay(nn(1,ndnr(ireac))),' and  ',ay(nn(2,ndnr(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
      write(out_line,*)
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ PHOTOLYTIC PRODUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kps(igas+1)-kps(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,20) ' Reaction # ',nps(ireac),' produces ',
     *      ay(kss(1,nps(ireac))),' and  ', ay(kss(2,nps(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     & '______________ PHOTOLYTIC DESTRUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kds(igas+1)-kds(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,30) ' Reaction # ',nds(ireac),' destroys ',
     *      ay(ks(nds(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
  10  format(1x,a8)
  20  format(a12,i4,a10,a8,a6,a8)
  30  format(a12,i4,a10,a8) 
       
      return
      end SUBROUTINE printls

