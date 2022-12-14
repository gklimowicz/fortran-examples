      module mdul_glbdia
! --- hycom version 0.9
      USE HYCOM_DIM,ONLY : ogrid,isp,ifp,ilp,kk,idm,jchunk
     .   ,I_0,I_1,I_0H,I_1H
     .   ,J_0,J_1,J_0H,J_1H
      USE HYCOM_SCALARS, only : baclin,thref,nstep,nstep0
     & ,diagno,area,avgbot,ocnvol,spcifh,g,onem,itest,jtest
     & ,tmean1,smean1
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, GLOBALSUM
      implicit none
      integer :: i,j,k,l,km,kn

      contains

      subroutine glbsum(mm,nn,text,h_glb_cum,s_glb_cum)
!
! --- compare area-integrated surface fluxes to volume-integrated tendencies
!
      integer  ,intent(IN) :: mm,nn
      real     ,intent(IN) :: h_glb_cum,s_glb_cum  ! global accumulation
      character,intent(IN) :: text*(*)
      real,parameter :: offset=35.
      real tmeam,tmean,smeam,smean,glbhflx,glbsflx
     .  ,tmeamj(J_0H:J_1H),tmeanj(J_0H:J_1H),heatfxj(J_0H:J_1H)
     .  ,smeamj(J_0H:J_1H),smeanj(J_0H:J_1H),saltfxj(J_0H:J_1H)

!     if (nstep.ne.nstep0+1 .and. .not.diagno) return
! --- compute global integrals of T/S, incl. changes since t=0.
! --- print time-integrated fluxes (if provided) for comparison.
      do 85 j=J_0,J_1
      tmeamj(j)=0.
      tmeanj(j)=0.
      smeamj(j)=0.
      smeanj(j)=0.
!
      do 85 k=1,kk
      km=k+mm
      kn=k+nn
      do 85 l=1,isp(j)
!
      do 85 i=ifp(j,l),ilp(j,l)
        tmeamj(j)=tmeamj(j)+ temp(i,j,km)        *dp(i,j,km)*scp2(i,j)
        tmeanj(j)=tmeanj(j)+ temp(i,j,kn)        *dp(i,j,kn)*scp2(i,j)
        smeamj(j)=smeamj(j)+(saln(i,j,km)-offset)*dp(i,j,km)*scp2(i,j)
        smeanj(j)=smeanj(j)+(saln(i,j,kn)-offset)*dp(i,j,kn)*scp2(i,j)
 85   continue
!
      call GLOBALSUM(ogrid,  tmeamj,  tmeam, all=.true.)
      call GLOBALSUM(ogrid,  tmeanj,  tmean, all=.true.)
      call GLOBALSUM(ogrid,  smeamj,  smeam, all=.true.)
      call GLOBALSUM(ogrid,  smeanj,  smean, all=.true.)
!
      if (nstep .eq. nstep0+1) then
        tmean1=tmean
        smean1=smean
      end if

      if( AM_I_ROOT() ) then
       print 100,nstep,text,
     .   tmeam/(ocnvol*onem),smeam/(ocnvol*onem),
     .   tmean/(ocnvol*onem),smean/(ocnvol*onem),
     .   nint((tmeam-tmean1)/(ocnvol*onem)*1.e10),
     .   nint((smeam-smean1)/(ocnvol*onem)*1.e10),
     .   nint((tmean-tmean1)/(ocnvol*onem)*1.e10),
     .   nint((smean-smean1)/(ocnvol*onem)*1.e10)
 100  format (i7,1x,a12,4f15.10/20x,4i15)

c     if (h_glb_cum.gt.0.) then
        print '(a20,30x,2i15,f6.1)','srf.fluxes',
     .  nint(h_glb_cum/(ocnvol*onem*spcifh)*1.e10),
     .  nint(s_glb_cum/(ocnvol*onem       )*1.e10)

      print '(i9,2(a,2f9.4),a)',nstep,' global mean temp,saln:',
     .tmean/(ocnvol*onem),smean/(ocnvol*onem),
     .'  (initl:',tmean1/(ocnvol*onem),smean1/(ocnvol*onem),')'
      print '(i9,a,2es11.3)',nstep,' global temp,saln drift 
     .(1.e6deg,1.e6psu):',
     .(tmean-tmean1)*1.e6/(ocnvol*onem),
     .(smean-smean1)*1.e6/(ocnvol*onem)
      print '(i9,a,2es11.3)',nstep,' srf.flux-induced drift 
     .(1.e6deg,1.e6psu):',
     .h_glb_cum*1.e6/(spcifh*ocnvol*onem),
     .s_glb_cum*1.e6/(       ocnvol*onem)
c       end if !h_glb_cum > 0
      end if !AM_I_ROOT
c
      return
      end subroutine glbsum

      real function glob2d(field)

! --- compute global area integral of 'field'
!
      real,intent(IN) :: field(I_0H:I_1H,J_0H:J_1H)
      real    :: fldcol(J_0H:J_1H)
      integer :: i,j,l

      do 85 j=J_0,J_1
      fldcol(j)=0.
      do 85 l=1,isp(j)
      do 85 i=ifp(j,l),ilp(j,l)
      fldcol(j)=fldcol(j)+field(i,j)*scp2(i,j)
 85   continue

      call GLOBALSUM(ogrid,fldcol,glob2d, all=.true.)
!     glob2d=glob2d/area  ! global-mean

      return
      end function glob2d

      real function ocn_surface_area_integral(field,weight)
        real,intent(IN) :: field(I_0H:I_1H,J_0H:J_1H)
        real,intent(IN) :: weight(I_0H:I_1H,J_0H:J_1H)
        real    :: fldcol(J_0H:J_1H)
        integer :: i,j,l

        do j=J_0,J_1
          fldcol(j)=0.
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            fldcol(j)=fldcol(j)+field(i,j)*weight(i,j)*scp2(i,j)
          end do
          end do
        end do

        call GLOBALSUM(ogrid,fldcol,ocn_surface_area_integral,
     &       all=.true.)

        return
      end function ocn_surface_area_integral

      end module mdul_glbdia
