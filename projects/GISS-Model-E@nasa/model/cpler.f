#include "rundeck_opts.h"
#include "hycom_mpi_hacks.h"

#ifndef STANDALONE_HYCOM

#ifndef CUBED_SPHERE

      module hycom_cpler
      USE CONSTANT, only : radius,tf
      USE HYCOM_DIM_GLOB, only : iia,jja,iio,jjo,isp,ifp,ilp,ii,jj,ip
      USE HYCOM_SCALARS, only : flnma2o,flnmo2a,flnmcoso,huge,
     &                          itest,jtest
      USE HYCOM_ARRAYS, only : scp2
      USE FLUXES, only : focean_loc => focean
      USE HYCOM_DIM, only : ogrid
     &    ,aJ_0, aJ_1, aJ_0H, aJ_1H,
     &      J_0,  J_1,  J_0H,  J_1H
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : am_i_root,pack_data,unpack_data
     .    ,dist_grid
      use GEOM, only : adxyp          		! atm.cell size
      use filemanager, only : findunit
      implicit none
      integer i,j,l,n,ia,ja,jb,io,jo
      integer :: iotest=-1,jotest=-1,iatest=-1,jatest=-1
c     integer :: iotest=135,jotest=345,iatest=-1,jatest=-1
      save iatest,jatest
      logical vrbos,fail
      logical :: glbdiag=.false., chkbig=.false.
c     logical :: glbdiag=.true.,  chkbig=.true.
      integer ibig,jbig
      real scale,totwgt,peak,big
c
      private

      public fld_o2a,vec_o2a,fld_a2o,vec_a2o,tempr_o2a,cpl_wgt

      public  ilista2o, jlista2o, nlista2o, wlista2o
     .       ,ilisto2a, jlisto2a, nlisto2a, wlisto2a
     .       ,coso, sino

      public nwgta2o,nwgto2a
#ifdef ATM4x5
#ifdef HYCOM2deg
      integer, parameter:: nwgta2o=18,nwgto2a=39
#endif
#endif
#ifdef ATM2x2h
#ifdef HYCOM2deg
      integer, parameter:: nwgta2o=36,nwgto2a=19
#endif
#endif
#ifdef ATM2x2h
#ifdef HYCOM1degRefined
      integer, parameter:: nwgta2o=20,nwgto2a=80
#endif
#ifdef HYCOM1degUnrefined
      integer, parameter:: nwgta2o=20,nwgto2a=80
#endif
#endif
c
      real wlista2o(iio,jjo,nwgta2o), wlisto2a(iia,jja,nwgto2a)
     .      ,coso(iio,jjo),sino(iio,jjo),scp2_glb(iio,jjo)
      real, allocatable :: focean(:,:)
      integer ilista2o(iio,jjo,nwgta2o),jlista2o(iio,jjo,nwgta2o)
     .                                 ,nlista2o(iio,jjo)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a(iia,jja,nwgto2a)
     .                                 ,nlisto2a(iia,jja)
      integer iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8
      contains


      subroutine fld_o2a(fldo_loc,flda_loc,string5)
c --- mapping fluxes from ogcm A grid to agcm A grid
c     input: fldo (W/m^2), output: flda (W/m^2)

      character,intent(IN) :: string5*5
      real,intent(IN)    :: fldo_loc(iio,J_0H:J_1H)
      real,intent(INOUT) :: flda_loc(iia,aJ_0H:aJ_1H)
      real, allocatable  :: flda(:,:),fldo(:,:),wgto(:,:)
      real               :: glbintgl_a, glbintgl_o, epsil1
      integer iprox,jprox			! distance from iatest,jatest

      vrbos = iotest > 0 .or. iatest > 0
      if (am_i_root()) then
        allocate(flda(iia,jja),fldo(iio,jjo),wgto(iio,jjo))
      else
        allocate(flda(1,1),fldo(1,1))
      endif

      call pack_data(ogrid,fldo_loc,fldo)
c
      if (am_i_root()) then
       if (vrbos) then
        if (iotest==-1)
     .    call aij2oij(iatest,jatest,iotest,jotest)
        if (iatest==-1)
     .    call oij2aij(iotest,jotest,iatest,jatest)
        peak=1.e-8
       end if

! --- compute global integral of input field
       if (glbdiag) then
        call pack_data(ogrid,scp2,scp2_glb)
        call pack_data(agrid,focean_loc,focean)
        glbintgl_a = 0.
        glbintgl_o = 0.
        do 5 j=1,jj
        do 5 l=1,isp(j)
        do 5 i=ifp(j,l),ilp(j,l)
 5      glbintgl_o = glbintgl_o + fldo(i,j)*scp2_glb(i,j)
       end if

       big=0.
       fail=.false.
       do 16 ja=1,jja
       if (vrbos) jprox=abs(ja-jatest)			! distance from jatest
       do 16 ia=1,iia
       if (vrbos) then
        iprox=min(abs(ia-iatest),abs(ia+iia-iatest),abs(ia-iia-iatest))
        if (max(iprox,jprox).lt.5) peak=max(peak,abs(flda(ia,ja)))
       end if
       if (ia.eq.iatest .and. ja.eq.jatest)
     .  print '(a,2i5/(4(i7,i4,es9.2)))',
     .   '(fld_o2a) interpolation weights at ia,ja=',ia,ja,
     .   (ilisto2a(ia,ja,n),jlisto2a(ia,ja,n),wlisto2a(ia,ja,n),
     .   n=1,nlisto2a(ia,ja))

       if (nlisto2a(ia,ja).gt.0) then
        flda(ia,ja)=0.
        do 17 n=1,nlisto2a(ia,ja)
        io=ilisto2a(ia,ja,n)
        jo=jlisto2a(ia,ja,n)
        if (vrbos) then
         if (ip(io,jo).eq.0) then
          print '(2a,4i5)','(fld_o2a) error: reference to land point:',
     .    '  ia,ja,io,jo=',ia,ja,io,jo
          fail=.true.
         end if
        end if			! vrbos
 17     flda(ia,ja)=flda(ia,ja)+wlisto2a(ia,ja,n)*fldo(io,jo)

! --- compute global integral of output field
        if (glbdiag) glbintgl_a = glbintgl_a + flda(ia,ja)*adxyp(ia,ja)
     &                                                   *focean(ia,ja)

        if (chkbig.and.abs(flda(ia,ja)).gt.big) then
         big=abs(flda(ia,ja))
         ibig=ia
         jbig=ja
        end if
       end if			! nlisto2a > 0
 16    continue
c
       if (chkbig) print '(a,2i5,es11.3)',
     .  '(fld_o2a) largest ingtl contrib. at ia,ja=',
     .  ibig,jbig,flda(ibig,jbig)

       if (vrbos) then
        scale=10.**(2-int(log10(peak)))
        print 101,trim(string5),scale
 101    format ('(cpler) shown below:  ',a,'  scaled by',es9.1)
        call pr_9x9(fldo,iio,jjo,iotest,jotest,0.,scale,
     .    'o2a O  '//string5)
        call pr_9x9(flda,iia,jja,iatest,jatest,0.,scale,
     .    'o2a A  '//string5)

        if (string5=='s s t') then
         wgto(:,:)=0.
         i=iatest
         j=jatest
         do 7 n=1,nlisto2a(i,j)
 7       wgto(ilisto2a(i,j,n),jlisto2a(i,j,n))=wlisto2a(i,j,n)
         scale=1000.
         print 101,'coupler weights',scale
         call pr_9x9(wgto,iio,jjo,iotest,jotest,0.,scale,'o2a weights')
         print 101,'focean',scale
         call pr_9x9(focean,iia,jja,iatest,jatest,0.,scale,'focean')
        end if
       end if			! vrbos

       if (glbdiag) print '(/3a,es15.6/a,es15.6/)',
     .      "(fld_o2a) ",string5," glbintgl_o=",glbintgl_o,
     .            "                glbintgl_a=",glbintgl_a

c --- atmo model ignores -flda- values at (2:iia,jja) but expects a zonal
c --- average of those values in (1,jja).
      flda(1,jja)=sum(flda(1:iia,jja))/float(iia)
c
      end if			! am_i_root

      call unpack_data(agrid,flda,flda_loc)

      deallocate(flda,fldo)
      return
      end subroutine fld_o2a


      subroutine vec_a2o(tauxa_loc,tauya_loc,tauxo_loc,tauyo_loc)
c --- mapping vector like stress from agcm to ogcm, both on A grid
c --- input  tauxa/tauya: E_/N_ward on agcm A grid
c --- output tauxo/tauyo: +i_/+j_ward on ogcm A grid (S_/E_ward in Mercator domain)
c
      real, dimension(iia,aJ_0H:aJ_1H) :: tauxa_loc,tauya_loc
      real, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real, dimension(:,:), allocatable :: tauxa,tauya,tauxo,tauyo,
     &     sward,eward

      if (am_i_root()) then
        allocate(
     &       tauxa(iia,jja),tauya(iia,jja),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       sward(iio,jjo),eward(iio,jjo)
     &       )
      else
        allocate(
     &       tauxa(1,1),tauya(1,1),
     &       tauxo(1,1),tauyo(1,1),
     &       sward(1,1),eward(1,1)
     &       )
      endif
      call pack_data(agrid,tauxa_loc,tauxa)
      call pack_data(agrid,tauya_loc,tauya)
      call pack_data(ogrid,scp2,scp2_glb)

      tauxo(:,:)=0.
      tauyo(:,:)=0.
      if (am_i_root()) then
c
c --- mapping tauxa/tauya to ogcm grid
       do 6 j=1,jj
       do 6 l=1,isp(j)
       do 6 i=ifp(j,l),ilp(j,l)
       eward(i,j)=0.
       sward(i,j)=0.
c
       do 7 n=1,nlista2o(i,j)
       eward(i,j)=eward(i,j)+tauxa(ilista2o(i,j,n),jlista2o(i,j,n))
     .                                            *wlista2o(i,j,n)
 7     sward(i,j)=sward(i,j)-tauya(ilista2o(i,j,n),jlista2o(i,j,n))
     .                                            *wlista2o(i,j,n)
 6     continue
c
c --- rotate sward/eward to fit onto Panam grid
       do 9 j=1,jj
       do 9 l=1,isp(j)
       do 9 i=ifp(j,l),ilp(j,l)
       tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
       tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9     continue
      endif ! am_i_root

      call unpack_data(ogrid,tauxo,tauxo_loc)
      call unpack_data(ogrid,tauyo,tauyo_loc)

      deallocate(tauxa,tauya,tauxo,tauyo,sward,eward)
      return
      end subroutine vec_a2o


      subroutine fld_a2o(flda_loc,fldo_loc,string5)
c --- mapping flux fields from agcm A grid to ogcm A grid
c     input: flda (W/m^2), output: fldo (W/m^2)
c
      character,intent(IN) :: string5*5
      real,intent(IN)    :: flda_loc(iia,aJ_0H:aJ_1H)
      real,intent(INOUT) :: fldo_loc(iio,J_0H:J_1H)
      real, allocatable  :: flda(:,:),fldo(:,:),wgta(:,:)
      real               :: glbintgl_a, glbintgl_o, epsil1
      integer iprox,jprox			! distance from iotest,jotest

      vrbos = iotest > 0 .or. iatest > 0
      if (am_i_root()) then
        allocate(flda(iia,jja),fldo(iio,jjo),wgta(iia,jja))
      else
        allocate(flda(1,1),fldo(1,1))
      endif

      call pack_data(agrid,flda_loc,flda)
c
      if (am_i_root()) then
       if (vrbos) then
        if (iotest==-1)
     .    call aij2oij(iatest,jatest,iotest,jotest)
        if (iatest==-1)
     .    call oij2aij(iotest,jotest,iatest,jatest)
        fldo(:,:)=huge			! turn land points into **** in printout
        peak=1.e-8
       end if

! --- compute global integral of input field
       big=0.
       if (glbdiag) then
        call pack_data(ogrid,scp2,scp2_glb)
        call pack_data(agrid,focean_loc,focean)
        glbintgl_a = 0.
        glbintgl_o = 0.
        do 5 ja=1,jja
        do 5 ia=1,iia

        if (chkbig.and.abs(flda(ia,ja)).gt.big.and.focean(ia,ja).gt.0.) 
     .                                    then
         big=abs(flda(ia,ja))
         ibig=ia
         jbig=ja
        end if
 5      glbintgl_a = glbintgl_a
     .      + flda(ia,ja)*adxyp(ia,ja)*focean(ia,ja)
       end if

       if (chkbig) print '(2a,2i5,es11.3)',
     . '(fld_a2o) largest a field of ',
     .  trim(string5),ibig,jbig,flda(ibig,jbig)

       big=0.
       do 8 j=1,jj
       do 8 l=1,isp(j)
       do 8 i=ifp(j,l),ilp(j,l)
       if (i.eq.iotest .and. j.eq.jotest)
     .  print '(a,2i5/(4(i7,i4,f9.2)))',
     .   '(fld_a2o) interpolation weights at io,jo =',i,j,
     .   (ilista2o(i,j,n),jlista2o(i,j,n),wlista2o(i,j,n),
     .   n=1,nlista2o(i,j))

       fldo(i,j)=0.
       do 9 n=1,nlista2o(i,j)
       ia=ilista2o(i,j,n)
       ja=jlista2o(i,j,n)
 9     fldo(i,j)=fldo(i,j)+wlista2o(i,j,n)*flda(ia,ja)

! --- compute global integral of output field
       if (glbdiag) glbintgl_o = glbintgl_o + fldo(i,j)*scp2_glb(i,j)

        if (chkbig.and.abs(fldo(i,j)).gt.big) then
         big=abs(fldo(i,j))
         ibig=i
         jbig=j
        end if
 8     continue

       if (vrbos) then
         do 4 j=1,jj
         jprox=min(abs(j-jotest),abs(j+jj-jotest),abs(j-jj-jotest))
         do 4 l=1,isp(j)
         do 4 i=ifp(j,l),ilp(j,l)
           iprox=abs(i-iotest)                    ! distance from iotest
           if (max(iprox,jprox).lt.5) peak=max(peak,abs(fldo(i,j)))
 4       continue
       end if

       if (chkbig) print '(a,2i5,es11.3)',
     .  '(fld_a2o) largest ingtl contrib. at io,jo=',
     .  ibig,jbig,fldo(ibig,jbig)

       if (vrbos) then
        scale=10.**(2-int(log10(peak)))
        print 101,trim(string5),scale
 101    format ('(cpler) shown below:  ',a,'  scaled by',es9.1)
        call pr_9x9(flda,iia,jja,iatest,jatest,0.,scale,
     .   'a2o A  '//string5)
        call pr_9x9(fldo,iio,jjo,iotest,jotest,0.,scale,
     .   'a2o O  '//string5)

        if (string5=='icecv') then
         wgta(:,:)=0.
         i=iotest
         j=jotest
         do 7 n=1,nlista2o(i,j)
 7       wgta(ilista2o(i,j,n),jlista2o(i,j,n))=wlista2o(i,j,n)
         scale=1000.
         print 101,'coupler weights',scale
         call pr_9x9(wgta,iia,jja,iatest,jatest,0.,scale,'a2o  weights')
         print 101,'focean',scale
         call pr_9x9(focean,iia,jja,iatest,jatest,0.,scale,'focean')
        end if
       end if			! vrbos

       if (glbdiag) print '(/3a,es15.6/a,es15.6/)',
     .      "(fld_a2o) ",string5," glbintgl_a=",glbintgl_a,
     .            "                glbintgl_o=",glbintgl_o

      endif			! am_i_root
 100  format (' (cpler) 'a,(4(a8,'=',es10.3)))
c
      call unpack_data(ogrid,fldo,fldo_loc)

      deallocate(flda,fldo)
      return
      end subroutine fld_a2o


      subroutine vec_o2a(tauxo_loc,tauyo_loc,tauxa_loc,tauya_loc)
c --- mapping velocity/stress vector from C-grid ogcm to A-grid agcm.
c --- input  tauxo/tauyo: +i_/+j_ward (S_/E_ward in Mercador domain) on
c --- ocean C grid (@ i-1/2 & j-1/2)
c --- output tauxa/tauya: E_/N_ward on agcm A grid
c
      real,intent(in),   dimension(iio, J_0H: J_1H)::tauxo_loc,tauyo_loc
      real,intent(inout),dimension(iia,aJ_0H:aJ_1H)::tauxa_loc,tauya_loc
      real,dimension(:,:),allocatable:: tauxa,tauya,tauxo,tauyo,
     &     nward,eward
      real sine

      if (am_i_root()) then
        allocate(
     &       tauxa(iia,jja),tauya(iia,jja),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       nward(iio,jjo),eward(iio,jjo)
     &       )
      else
        allocate(
     &       tauxa(1,1),tauya(1,1),
     &       tauxo(1,1),tauyo(1,1),
     &       nward(1,1),eward(1,1)
     &       )
      endif

      call pack_data(ogrid,tauxo_loc,tauxo)
      call pack_data(ogrid,tauyo_loc,tauyo)
      if(am_i_root()) then
c
c --- average tauxo/tauyo from C to A grid & rotate to n/e orientation at A grid
c --- check velocity bounds
      do 12 j=1,jj
      jb=mod(j,jj)+1
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue
c
c --- weights are for mapping nward/eward from ogcm to agcm, both on A grid
c
      big=0.
      do 16 ja=1,jja
      do 16 ia=1,iia
      tauxa(ia,ja)=0.
      tauya(ia,ja)=0.
c
      if (nlisto2a(ia,ja).gt.0) then
        do 17 n=1,nlisto2a(ia,ja)
        tauxa(ia,ja)=tauxa(ia,ja)+eward(ilisto2a(ia,ja,n)
     .            ,jlisto2a(ia,ja,n))*wlisto2a(ia,ja,n)
 17     tauya(ia,ja)=tauya(ia,ja)+nward(ilisto2a(ia,ja,n)
     .            ,jlisto2a(ia,ja,n))*wlisto2a(ia,ja,n)

        if (chkbig.and.abs(tauxa(ia,ja)).gt.big) then
          big=abs(tauxa(ia,ja))
          ibig=ia
          jbig=ja
        end if
      end if
 16   continue

      if (vrbos) then
c       scale=10.**(2-int(log10(peak)))
        scale=1000.
        print 101,scale
 101    format ('(cpler) shown below: veco2a scaled by',es9.1)
        call pr_9x9(nward,iio,jjo,iotest,jotest,0.,scale,
     .    'o2a O  ')
        call pr_9x9(tauxa,iia,jja,iatest,jatest,0.,scale,
     .    'o2a A  ')

!       call prtfld(tauxa,iia,iia,jja,0.,1000.,'tauxa')
!       call prtfld(tauya,iia,iia,jja,0.,1000.,'tauya')
      endif ! vrbos

      if (chkbig) print '(a,2i5,es11.3)',
     . ' largest ingtl contrib. veco2a',ibig,jbig,tauxa(ibig,jbig)

      endif ! am_i_root

      call unpack_data(agrid,tauxa,tauxa_loc)
      call unpack_data(agrid,tauya,tauya_loc)

      deallocate(tauxa,tauya,tauxo,tauyo,nward,eward)
      return
      end subroutine vec_o2a


      subroutine tempr_o2a(fldo_loc,flda_loc)
c --- mapping sqrt(sqrt(temp**4)) from ogcm A grid to agcm A grid
c --- input: fldo in deg C; outout: flda in deg K
c
      real fldo_loc(iio,J_0H:J_1H),flda_loc(iia,aJ_0H:aJ_1H)
      real, allocatable :: flda(:,:),fldo(:,:)

      if (am_i_root()) then
        allocate(flda(iia,jja),fldo(iio,jjo))
      else
        allocate(flda(1,1),fldo(1,1))
      endif

      call pack_data(ogrid,fldo_loc,fldo)
      if (am_i_root()) then
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      if (nlisto2a(ia,ja) > 0) then
        do 17 n=1,nlisto2a(ia,ja)
 17     flda(ia,ja)=flda(ia,ja)+(fldo(ilisto2a(ia,ja,n),
     .    jlisto2a(ia,ja,n))+tf)**4*wlisto2a(ia,ja,n)
        flda(ia,ja)=sqrt(sqrt(flda(ia,ja)))       ! Kelvin for radiation
      end if
 16   continue
c
c --- atmo model ignores -flda- values at (2:iia,jja) but expects a zonal
c --- average of those values in (1,jja).
      flda(1,jja)=sum(flda(1:iia,jja))/float(iia)
c
      endif ! am_i_root

      call unpack_data(agrid,flda,flda_loc)
      deallocate(flda,fldo)
      return
      end subroutine tempr_o2a


      subroutine cpl_wgt
      USE HYCOM_SCALARS, only : flnmcellsz
      integer :: iz,jz
      real :: totwgt,radsqu,scp2_glb(iio,jjo),epsil1
c
#ifdef  ATM4x5
#ifdef  HYCOM2deg
       integer, parameter :: nsize1=10249200, nsize2=2079936
#endif
#endif
#ifdef ATM2x2h
#ifdef HYCOM2deg
       integer, parameter :: nsize1=20358000, nsize2=3991680
#endif
#endif
#ifdef ATM2x2h
#ifdef HYCOM1degRefined
       integer, parameter :: nsize1=45139680, nsize2=16640640  ! 387x360
#endif
#ifdef HYCOM1degUnrefined
       integer, parameter :: nsize1=41873760, nsize2=16640640  ! 359x360
#endif
#endif
      if(am_i_root()) then
        allocate(focean(iia,jja))
        call findunit(iu1)
        open(iu1,file=flnmcellsz,form='unformatted',status='old')
        read(iu1) scp2_glb
        close(iu1)
        write(*,*) 'below scp2 in wgt'
        call zebra(scp2_glb,iio,iio,jjo)
      else
        allocate(focean(1,1))
      endif
      call pack_data(agrid,focean_loc,focean)
c
      if (.not. am_i_root()) return

c --- read in all weights

      if (iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8).ne.nsize1 .or.
     .    iia*jja*((nwgto2a*2+1)*4+nwgto2a*8).ne.nsize2) then
        write(*,'(a,2i12,a,2i12)') 'wrong size in cpler1 '
     .  ,iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8)
     .  ,iia*jja*((nwgto2a*2+1)*4+nwgto2a*8)
     .  ,' should be ',nsize1,nsize2
        stop '(wrong array size in cpl_wgt)'
      endif
c
c --- read a2o interpolation weights
      call findunit(iu1)
      open(iu1,file=flnma2o,form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(iu1,rec=1) ilista2o,jlista2o,wlista2o,nlista2o
      close(iu1)
c
c --- read o2a interpolation weights
      call findunit(iu2)
      open(iu2,file=flnmo2a,form='unformatted',status='old',
     .  access='direct',recl=nsize2)
      read(iu2,rec=1) ilisto2a,jlisto2a,wlisto2a,nlisto2a
      close(iu2)
c
c --- weights represent overlap areas on unit sphere. convert to areas on
c --- actual earth (dimensioned m^2)

      radsqu=radius*radius
      print *,'(cpl_wgt) earth radius:',radius
      do ja=1,jja
       do ia=1,iia
        do n=1,nlisto2a(ia,ja)
         wlisto2a(ia,ja,n)=wlisto2a(ia,ja,n)*radsqu/(adxyp(ia,ja)
     &                                             *focean(ia,ja))
        end do
       end do
      end do
      do jo=1,jjo
       do io=1,iio
        do n=1,nlista2o(io,jo)
         wlista2o(io,jo,n)=wlista2o(io,jo,n)*radsqu/scp2_glb(io,jo)
        end do
       end do
      end do

c --- check for mismatch between sum of weights and grid cell size using epsil1

      epsil1=1.e-4
      fail=.false.
c --- check a2o weights
      do ja=1,jja
      do ia=1,iia
        if (nlisto2a(ia,ja).gt.0) then
          if (ia==iatest.and.ja==jatest) 
     .    write(*,'(a,2i4,f6.2,a,/,(2i4,f6.2))')
     .    '(cpl_wgt a2o)',ia,ja,focean(ia,ja),' => '
     .    ,((ilisto2a(ia,ja,n),jlisto2a(ia,ja,n),
     .     wlisto2a(ia,ja,n)),n=1,nlisto2a(ia,ja))
          totwgt=0.
          do n=1,nlisto2a(ia,ja)
           totwgt=totwgt+wlisto2a(ia,ja,n)
          end do
          if (totwgt.gt.0. .and. abs(totwgt-1.).gt.epsil1) then
           print '(a,2i4,a,f13.7)','ia,ja=',ia,ja,
     .      ' (cpl_wgt) sum of o2a intp.weights not 1:',totwgt
           fail=.true.
          end if
        end if
      end do
      end do

c --- check o2a weights
      do jo=1,jjo
      do io=1,iio
        if (nlista2o(io,jo).gt.0) then
          if (io==iotest.and.jo==jotest) 
     .    write(*,'(a,2i4,a,/,(2i4,2f6.2))')
     .    '(cpl_wgt o2a)',io,jo,' => '
     .    ,((ilista2o(io,jo,n),jlista2o(io,jo,n),wlista2o(io,jo,n)
     .    ,focean(ilista2o(io,jo,n),jlista2o(io,jo,n)))
     .    ,n=1,nlista2o(io,jo))
          totwgt=0.
          do n=1,nlista2o(io,jo)
           totwgt=totwgt+wlista2o(io,jo,n)
          end do
          if (abs(totwgt-1.).gt.epsil1) then
           print '(a,2i4,a,f13.7)','io,jo=',io,jo,
     .      ' (cpl_wgt) sum of a2o intp.weights not 1:',totwgt
c          fail=.true.
          end if
        end if
      end do
      end do
      if (fail) stop '(cpl_wgt.f)'

c --- in grid cells adjacent to the poles (ja=1 and 90), results are only
c --- provided on a single meridian (ia=1). delete references to cells ia>1.

      do jo=1,jjo
      do io=1,iio
        vrbos=.false.
        do n=1,nlista2o(io,jo)
          if (jlista2o(io,jo,n).eq.1 .or. jlista2o(io,jo,n).eq.jja)
     .     vrbos=.true.
        end do

        if (vrbos)
     .   print 106,'(cpl_wgt) io,jo=',io,jo,
     .     '   weight list BEFORE removing references to ia > 1',
     .     (ilista2o(io,jo,n),jlista2o(io,jo,n),wlista2o(io,jo,n),
     .     n=1,nlista2o(io,jo))
 106     format (a,2i4,a/(4(i6,i4,es10.2)))

        do n=1,nlista2o(io,jo)
         if (jlista2o(io,jo,n).eq.1 .or. jlista2o(io,jo,n).eq.jja)
     .    ilista2o(io,jo,n)=1
        end do

        if (vrbos)
     .   print 106,'(cpl_wgt) io,jo=',io,jo,
     .     '   weight list  AFTER removing references to ia > 1',
     .     (ilista2o(io,jo,n),jlista2o(io,jo,n),wlista2o(io,jo,n),
     .     n=1,nlista2o(io,jo))
      end do  
      end do  

      call findunit(iu3)
      open(iu3,file=flnmcoso,form='unformatted',status='old')
      read(iu3) iz,jz,coso,sino
      close(iu3)
      if (iz.ne.iio .or. jz.ne.jjo) then
        write(*,*) ' iz,jz=',iz,jz
        stop '(wrong iz/jz in cososino.8bin)'
      endif
      write (*,*) 'shown below: coso'
      call zebra(coso,iio,iio,jjo)
      write (*,*) 'shown below: sino'
      call zebra(sino,iio,iio,jjo)
c
      return
      end subroutine cpl_wgt

      subroutine oij2aij(io,jo,ia,ja)
c --- find 'a' point nearest prescribed 'o' point
      integer,intent(IN)  :: io,jo
      integer,intent(OUT) :: ia,ja
      real wmax
      wmax=0.
      ia=-1
      ja=-1
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      do 8 n=1,nlista2o(i,j)
      if (i.eq.io .and. j.eq.jo .and.wlista2o(i,j,n).gt.wmax) then
        wmax=wlista2o(i,j,n)
        ia=ilista2o(i,j,n)
        ja=jlista2o(i,j,n)
      end if
 8    continue
      if (ia.lt.0) then
        print '(a,2i5)','cannot find a-point near o-point',io,jo
        stop '(oij2aij error)'
      end if
      print '(2(a,2i5))','a-point nearest o-point',io,jo,'  is',ia,ja
      return
      end subroutine oij2aij


      subroutine aij2oij(ia,ja,io,jo)
c --- find 'o' point nearest prescribed 'a' point
      integer,intent(IN)  :: ia,ja
      integer,intent(OUT) :: io,jo
      real wmax
      wmax=0.
      io=-1
      jo=-1
      do 8 j=1,jja
      do 8 i=1,iia
      do 8 n=1,nlisto2a(i,j)
      if (i.eq.ia .and. j.eq.ja .and.wlisto2a(i,j,n).gt.wmax) then
        wmax=wlisto2a(i,j,n)
        io=ilisto2a(i,j,n)
        jo=jlisto2a(i,j,n)
      end if
 8    continue
      if (io.lt.0) then
        print '(a,2i5)','cannot find o-point near a-point',ia,ja
        stop '(aij2oij error)'
      end if
      print '(2(a,2i5))','o-point nearest a-point',ia,ja,'  is',io,jo
      return
      end subroutine aij2oij
      end module hycom_cpler

!--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>---<>--<>--<>--<>--<>--<>--<>

      module hycom_dynsi_cpler  ! latlon atm not using this yet
      implicit none
      contains

      subroutine init_hycom_dynsi_cpler
      end subroutine init_hycom_dynsi_cpler
      end module hycom_dynsi_cpler

#else  /* coupler with cubed-sphere atm is in cplercs.f */

      module hycom_dynsi_cpler
      USE HYCOM_DIM_GLOB, only : iio,jjo,isp,ifp,ilp,ii,jj,ip
      USE DOMAIN_DECOMP_1D, only : am_i_root,get,pack_data,unpack_data
      USE HYCOM_DIM, only : ogrid, J_0H, J_1H
      USE ICEDYN, only : igrid=>grid_ICDYN,iii=>imicdyn,jji=>jmicdyn
      USE HYCOM_CPLER, only : coso=>coso_glob, sino=>sino_glob
      USE HYCOM_ATM, only : dmui_loc,dmvi_loc
      use filemanager, only : findunit
      private
      save
      public :: init_hycom_dynsi_cpler,veci2o,veco2i,scai2o
      public :: reset_dynsi_accum,do_dynsi_accum,
     &     taux_dynsi,tauy_dynsi,ustar_dynsi

      real, dimension(:,:,:), allocatable :: wlisto2i,wlisti2o
      integer, dimension(:,:,:), allocatable ::
     &     ilisto2i,jlisto2i, ilisti2o,jlisti2o
      integer, dimension(:,:), allocatable :: nlisto2i,nlisti2o

      real, dimension(:,:), allocatable ::
     &     taux_dynsi,tauy_dynsi,ustar_dynsi

      contains

      subroutine reset_dynsi_accum
      integer :: i_j0h,i_j1h
      if(.not.igrid%have_domain) return
      if(.not.allocated(taux_dynsi)) then
        i_j0h = igrid%j_strt_halo
        i_j1h = igrid%j_stop_halo
        allocate(taux_dynsi(iii,i_j0h:i_j1h)
     &          ,tauy_dynsi(iii,i_j0h:i_j1h)
     &          ,ustar_dynsi(iii,i_j0h:i_j1h)
     &       )
      endif
      taux_dynsi = 0.
      tauy_dynsi = 0.
      ustar_dynsi = 0.
      return
      end subroutine reset_dynsi_accum


      subroutine do_dynsi_accum(dt,dtacc)
      use hycom_scalars, only : thref
      real :: dt,dtacc
      integer :: i,j
      if(.not.igrid%have_domain) return
      do j=igrid%j_strt,igrid%j_stop
        do i=1,iii
c taux,tauy are time integrals over the time interval dt
          taux_dynsi(i,j) = taux_dynsi(i,j) + dmui_loc(i,j)/dtacc
          tauy_dynsi(i,j) = tauy_dynsi(i,j) + dmvi_loc(i,j)/dtacc
          ustar_dynsi(i,j) = ustar_dynsi(i,j) + (dt/dtacc)*
     &         sqrt(sqrt(dmui_loc(i,j)**2 +dmvi_loc(i,j)**2)*thref/dt)
        enddo
      enddo
      return
      end subroutine do_dynsi_accum


      subroutine init_hycom_dynsi_cpler
      integer :: iu
      integer :: nwgti2o,nwgto2i,nsize1,nsize2

      if(.not.am_i_root()) return

      if(iio==387 .and. jjo==360 .and. iii==180 .and. jji==180) then
        nwgti2o=80; nwgto2i=50
        nsize1=178886880; nsize2=10419840
      elseif(iio==359 .and. jjo==360 .and. iii==180 .and. jji==180) then
        nwgti2o=80; nwgto2i=50
        nsize1=165944160; nsize2=10419840
      else
        call stop_model('hycom_dynsi_cpler: unknown resolution',255)
      endif

      if (iio*jjo*((nwgti2o*2+1)*4+nwgti2o*8).ne.nsize1 .or.
     .    iii*jji*((nwgto2i*2+1)*4+nwgto2i*8).ne.nsize2) then
        write(6,'(a,2i12,a,2i12)') 'wrong size in cpler '
     .  ,iio*jjo*((nwgti2o*2+1)*4+nwgti2o*8)
     .  ,iii*jji*((nwgto2i*2+1)*4+nwgto2i*8)
     .  ,' should be ',nsize1,nsize2
        call stop_model(' wrong size in hycom_dynsi_cpler',255)
      endif

      allocate(wlisto2i(iii,jji,nwgto2i)
     &        ,ilisto2i(iii,jji,nwgto2i)
     &        ,jlisto2i(iii,jji,nwgto2i)
     &        ,nlisto2i(iii,jji))
      allocate(wlisti2o(iio,jjo,nwgti2o)
     &        ,ilisti2o(iio,jjo,nwgti2o)
     &        ,jlisti2o(iio,jjo,nwgti2o)
     &        ,nlisti2o(iio,jjo))

      call findunit(iu)
      open(iu,file='taui2o',form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(iu,rec=1) ilisti2o,jlisti2o,wlisti2o,nlisti2o
      close(iu)
c
      call findunit(iu)
      open(iu,file='tauo2i',form='unformatted',status='old',
     .  access='direct',recl=nsize2)
      read(iu,rec=1) ilisto2i,jlisto2i,wlisto2i,nlisto2i
      close(iu)

      return
      end subroutine init_hycom_dynsi_cpler


      subroutine veci2o(tauxi_loc,tauyi_loc,tauxo_loc,tauyo_loc)
c --- mapping vector like stress from B grid ice model to A grid ogcm
c --- input tauxi/tauyi: E_/N_ward on ice B grid located at i+1/2 & j+1/2 corner
c --- output tauxo/tauyo: +i_/+j_ward on ogcm A grid (S_/E_ward in Mercador domain)
c
      real, dimension(iii,igrid%J_STRT_HALO:igrid%J_STOP_HALO) ::
     &     tauxi_loc,tauyi_loc
      real, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real, dimension(:,:), allocatable :: tauxi,tauyi,tauxo,tauyo,
     &     sward,eward
      integer i,j,l,n

      if(am_i_root()) then
        allocate(
     &       tauxi(iii,jji),tauyi(iii,jji),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       sward(iio,jjo),eward(iio,jjo)
     &       )
      endif
      if(igrid%have_domain) then
        call pack_data(igrid,tauxi_loc,tauxi)
        call pack_data(igrid,tauyi_loc,tauyi)
        call pack_data(ogrid,scp2,scp2_glb)
      endif
      if(am_i_root()) then
c
c --- mapping B-grid tauxi/tauyi to A-grid ogcm
      do 6 j=1,jj
      do 6 l=1,isp(j)
      do 6 i=ifp(j,l),ilp(j,l)
      eward(i,j)=0.
      sward(i,j)=0.
c
      do 7 n=1,nlisti2o(i,j)
      eward(i,j)=eward(i,j)+tauxi(ilisti2o(i,j,n),jlisti2o(i,j,n))
     .                                           *wlisti2o(i,j,n)
 7    sward(i,j)=sward(i,j)-tauyi(ilisti2o(i,j,n),jlisti2o(i,j,n))
     .                                           *wlisti2o(i,j,n)
 6    continue
c
c --- rotate sward/eward to fit onto Panam grid
      do 9 j=1,jj
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
      tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
      tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9    continue
      endif ! am_i_root
      call unpack_data(ogrid,tauxo,tauxo_loc)
      call unpack_data(ogrid,tauyo,tauyo_loc)
      if(am_i_root()) then
        deallocate(tauxi,tauyi,tauxo,tauyo,sward,eward)
      endif
      return
      end subroutine vec_i2o


      subroutine vec_o2i(tauxo_loc,tauyo_loc,tauxi_loc,tauyi_loc)
c --- mapping vector like velocity from C grid ocean model to B grid ice model
c --- input tauxo/tauyo: +i_/+j_ward (S_/E_ward in Mercador domain) on ocean C grid (@ i-1/2 & j-1/2)
c --- output tauxi/tauyi: E_/N_ward on ice B grid
c
      real, dimension(iii,igrid%J_STRT_HALO:igrid%J_STOP_HALO) ::
     &     tauxi_loc,tauyi_loc
      real, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc

      real, dimension(:,:), allocatable :: tauxi,tauyi,tauxo,tauyo,
     &     nward,eward
      real sine
      integer i,j,l,n

      if(am_i_root()) then
        allocate(
     &       tauxi(iii,jji),tauyi(iii,jji),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       nward(iio,jjo),eward(iio,jjo)
     &       )
      endif
      call pack_data(ogrid,tauxo_loc,tauxo)
      call pack_data(ogrid,tauyo_loc,tauyo)
      if(am_i_root()) then
      nward(:,:)=0.
      eward(:,:)=0.
c --- average tauxo/tauyo from C to A grid & rotate to n/e orientation at A grid
c --- check velocity bounds
      do 12 j=1,jj
      jb=mod(j,jj)+1
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue
c
c --- mapping nward/eward from A-grid ogcm to B-grid ice model
c
      do 16 j=1,jji
      do 16 i=1,iii
      tauxi(i,j)=0.
      tauyi(i,j)=0.
c
      if (nlisto2i(i,j) > 0 ) then
      do 17 n=1,nlisto2i(i,j)
      tauxi(i,j)=tauxi(i,j)+eward(ilisto2i(i,j,n)
     .          ,jlisto2i(i,j,n))*wlisto2i(i,j,n)
 17   tauyi(i,j)=tauyi(i,j)+nward(ilisto2i(i,j,n)
     .          ,jlisto2i(i,j,n))*wlisto2i(i,j,n)
      endif
 16   continue
c
      endif ! am_i_root
      if(igrid%have_domain) then
        call unpack_data(igrid,tauxi,tauxi_loc)
        call unpack_data(igrid,tauyi,tauyi_loc)
      endif
      if(am_i_root()) then
        deallocate(tauxi,tauyi,tauxo,tauyo,nward,eward)
      endif
      return
      end subroutine vec_o2i


      subroutine fld_i2o(fldi_loc,fldo_loc)
c --- mapping scalar from B-grid ice model to A-grid ogcm
c --- input fldi: quantity on ice B grid located at i+1/2 & j+1/2 corner
c --- output fldo: quantity on ogcm A grid
c
      integer i,j,l,n
      real, dimension(iii,igrid%J_STRT_HALO:igrid%J_STOP_HALO) ::
     &     fldi_loc
      real, dimension(iio, J_0H: J_1H) :: fldo_loc
      real, dimension(:,:), allocatable :: fldi,fldo

      if(am_i_root()) then
        allocate(
     &       fldi(iii,jji),
     &       fldo(iio,jjo)
     &       )
      endif
      if(igrid%have_domain) then
        call pack_data(igrid,fldi_loc,fldi)
        call pack_data(ogrid,scp2,scp2_glb)
      endif
      if(am_i_root()) then
c
c --- mapping B-grid fldi to A-grid ogcm
      do 6 j=1,jj
      do 6 l=1,isp(j)
      do 6 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
      do n=1,nlisti2o(i,j)
        fldo(i,j)=fldo(i,j)+fldi(ilisti2o(i,j,n),jlisti2o(i,j,n))
     .                                          *wlisti2o(i,j,n)
      enddo
 6    continue
      endif ! am_i_root
      call unpack_data(ogrid,fldo,fldo_loc)
      if(am_i_root()) then
        deallocate(fldi,fldo)
      endif
      return
      end subroutine fld_i2o

      end module hycom_dynsi_cpler

!--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>---<>--<>--<>--<>--<>--<>--<>

#endif /* cubed-sphere versus lat-lon atm */

#else /* code below is a trivial coupler for the standalone ocean */

      module hycom_cpler
      USE HYCOM_DIM_GLOB, only : iio,jjo,jj
      USE DOMAIN_DECOMP_1D, only : dist_grid
      USE HYCOM_DIM, only : ogrid, J_0,J_1, J_0H,J_1H
     &    ,isp,ifp,ilp,ip
c
      implicit none
      private
      type(dist_grid), pointer, public :: agrid
      public ssto2a,veca2o,flxa2o,veco2a,tempro2a,cpl_wgt
      public coso_glob, sino_glob

      real, dimension(:,:), allocatable :: coso,sino
      real, dimension(:,:), allocatable :: coso_glob,sino_glob

      contains

      subroutine ssto2a(fldo,flda)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo, output: flda
c
      real fldo(iio,J_0H:J_1H),flda(iio,J_0H:J_1H)
      flda(:,:) = fldo(:,:)
      return
      end subroutine ssto2a


      subroutine veca2o(tauxa,tauya,tauxo,tauyo)
c --- mapping vector like stress from agcm to ogcm, both on A grid
c --- input  tauxa/tauya: E_/N_ward on agcm A grid
c --- output tauxo/tauyo: +i_/+j_ward on ogcm A grid (S_/E_ward in Mercator domain)
c
      real, dimension(iio,J_0H:J_1H) :: tauxa,tauya,tauxo,tauyo
      real, dimension(iio,J_0H:J_1H) :: sward,eward
      integer i,j,l

      eward(:,:) = +tauxa(:,:)
      sward(:,:) = -tauya(:,:)
c
c --- rotate sward/eward to fit onto Panam grid
      do 9 j=j_0,j_1
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
      tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
      tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9    continue
      return
      end subroutine veca2o


      subroutine flxa2o(flda,fldo)
c --- mapping flux-like field from agcm A grid to ogcm A grid
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      real flda(iio,J_0H:J_1H),fldo(iio,J_0H:J_1H)
      fldo(:,:) = flda(:,:)
      return
      end subroutine flxa2o


      subroutine veco2a(tauxo,tauyo,tauxa,tauya)
c --- mapping vector like velocity from C grid ogcm to A grid agcm
c --- input  tauxo/tauyo: +i_/+j_ward (S_/E_ward in Mercador domain) on ocean C grid (@ i-1/2 & j-1/2)
c --- output tauxa/tauya: E_/N_ward on agcm A grid
c
      use domain_decomp_1d, only : halo_update
      integer i,j,l
      real, dimension(iio,J_0H:J_1H) :: tauxo,tauyo,tauxa,tauya
      real sine

      call halo_update(ogrid,tauyo)

      tauxa(:,:) = 0.
      tauya(:,:) = 0.

c --- average tauxo/tauyo from C to A grid & rotate to n/e orientation at A grid
c --- check velocity bounds
      do 12 j=j_0,j_1
      jb = PERIODIC_INDEX(j+1, jj) !mod(j,jj)+1
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      tauya(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      tauxa(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue

      return
      end subroutine veco2a


      subroutine tempro2a(fldo,flda)
c --- mapping sqrt(sqrt(temp**4)) from ogcm A grid to agcm A grid
c --- input: fldo in deg C; outout: flda in deg K
c
      real fldo(iio,J_0H:J_1H),flda(iio,J_0H:J_1H)
      flda(:,:) = fldo(:,:)+tf
      return
      end subroutine tempro2a


      subroutine cpl_wgt
      USE HYCOM_SCALARS, only : flnmcoso,lp
      use filemanager, only : findunit
      USE DOMAIN_DECOMP_1D, only : am_i_root,unpack_data
      integer :: iu1,iz,jz

c read rotation coeffs between hycom gridlines and geographic north
      allocate(coso(iio,j_0h:j_1h),sino(iio,j_0h:j_1h))
      if(am_i_root()) then
        allocate(coso_glob(iio,jjo),sino_glob(iio,jjo))
        call findunit(iu1)
        open(iu1,file=flnmcoso,form='unformatted',status='old')
        read(iu1) iz,jz,coso_glob,sino_glob
        close(iu1)
        if (iz.ne.iio .or. jz.ne.jjo) then
          write(6,*) ' iz,jz=',iz,jz
          stop '(wrong iz/jz in cososino.8bin)'
        endif
      endif
      call unpack_data(ogrid,coso_glob,coso)
      call unpack_data(ogrid,sino_glob,sino)

      return
      end subroutine cpl_wgt

      end module hycom_cpler

!--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>---<>--<>--<>--<>--<>--<>--<>

      module hycom_dynsi_cpler
      contains

      subroutine init_hycom_dynsi_cpler
      end subroutine init_hycom_dynsi_cpler
      end module hycom_dynsi_cpler

!--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>---<>--<>--<>--<>--<>--<>--<>

#endif /* standalone ocean */
