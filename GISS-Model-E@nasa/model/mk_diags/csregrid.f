      module regrid_to_latlon_mod
      implicit none
      private

      integer, dimension(:), allocatable :: tile
      integer, dimension(:,:), allocatable :: ijcs,ijll
      real*8, dimension(:), allocatable :: xgrid_area,xp,yp
      real*8, dimension(:,:), allocatable :: areall

      integer :: imcub,imlon,jmlat,ncells

      public :: setup_remap,regrid1,regrid2,regrid_4d

      interface regrid1
        module procedure regrid1_2d
      end interface regrid1

      interface regrid2
        module procedure regrid2_2d
      end interface regrid2

      contains

      subroutine setup_remap(remap_fid,imlon_out,jmlat_out)
      integer :: n
      include 'netcdf.inc'
      integer :: remap_fid,imlon_out,jmlat_out
c
c read the remap file
c
      call get_dimsize(remap_fid,'imcub',imcub)
      call get_dimsize(remap_fid,'lon',imlon)
      call get_dimsize(remap_fid,'lat',jmlat)
      call get_dimsize(remap_fid,'ncells',ncells)
      allocate(tile(ncells),ijcs(2,ncells),ijll(2,ncells))
      allocate(xgrid_area(ncells),xp(ncells),yp(ncells))
      call get_var_int(remap_fid,'tile1',tile)
      call get_var_int(remap_fid,'tile1_cell',ijcs)
      call get_var_int(remap_fid,'tile2_cell',ijll)
      call get_var_r8(remap_fid,'xgrid_area',xgrid_area)
      call get_var_r8(remap_fid,'xprime',xp)
      call get_var_r8(remap_fid,'yprime',yp)

c
c calculate areas on the latlon grid
c
      allocate(areall(imlon,jmlat))
      do n=1,ncells
        areall(ijll(1,n),ijll(2,n)) = areall(ijll(1,n),ijll(2,n))
     &       + xgrid_area(n)
      enddo

      imlon_out = imlon
      jmlat_out = jmlat

      return
      end subroutine setup_remap

      subroutine get_var_r8(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      real*8 :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      endif
      status = nf_get_var_double(fid,varid,var)
      return
      end subroutine get_var_r8

      subroutine regrid1_2d(qcs,qll)
      real*4, dimension(imcub,imcub,6) :: qcs
      real*4, dimension(imlon,jmlat) :: qll
      integer :: ics,jcs,ill,jll,n
      qll(:,:) = 0.
      do n=1,ncells
        ics = ijcs(1,n)
        jcs = ijcs(2,n)
        ill = ijll(1,n)
        jll = ijll(2,n)
        qll(ill,jll) = qll(ill,jll) + xgrid_area(n)*
     &       qcs(ics,jcs,tile(n))
      enddo
      qll(:,:) = qll(:,:)/areall(:,:)
      return
      end subroutine regrid1_2d

      subroutine regrid1_3d(qcs,qll,nl)
      integer :: nl
      real*4, dimension(nl,imcub,imcub,6) :: qcs
      real*4, dimension(nl,imlon,jmlat) :: qll
      integer :: ics,jcs,ill,jll,n
      qll = 0.
      do n=1,ncells
        ics = ijcs(1,n)
        jcs = ijcs(2,n)
        ill = ijll(1,n)
        jll = ijll(2,n)
        qll(:,ill,jll) = qll(:,ill,jll) + xgrid_area(n)*
     &       qcs(:,ics,jcs,tile(n))
      enddo
      do jll=1,jmlat
        do ill=1,imlon
          qll(:,ill,jll) = qll(:,ill,jll)/areall(ill,jll)
        enddo
      enddo
      return
      end subroutine regrid1_3d

      subroutine regrid_4d(qcs,qll,nl,nk,order)
      integer :: nl,nk,order
      real*4, dimension(nl,imcub,imcub,nk,6) :: qcs
      real*4, dimension(nl,imlon,jmlat,nk) :: qll
      integer :: ics,jcs,ill,jll,n,k
      real*4, dimension(:,:,:), allocatable :: q2d
      real*4, dimension(:,:,:,:), allocatable :: q3d
      if(order.ne.1 .and. order.ne.2) then
        write(6,*) 'regrid_4d: order must be 1 or 2'
      endif
      if(nl.eq.1) then
        allocate(q2d(imcub,imcub,6))
        do k=1,nk
          q2d = qcs(1,:,:,k,:)
          if(order.eq.1) then
            call regrid1_2d(q2d,qll(1,:,:,k))
          elseif(order.eq.2) then
            call regrid2_2d(q2d,qll(1,:,:,k))
          endif
        enddo
        deallocate(q2d)
      else
        allocate(q3d(nl,imcub,imcub,6))
        do k=1,nk
          q3d = qcs(:,:,:,k,:)
          if(order.eq.1) then
            call regrid1_3d(q3d,qll(:,:,:,k),nl)
          elseif(order.eq.2) then
            call regrid2_3d(q3d,qll(:,:,:,k),nl)
          endif
        enddo
        deallocate(q3d)
      endif
      return
      end subroutine regrid_4d

      subroutine regrid2_2d(qcs,qll)
      real*4, dimension(imcub,imcub,6) :: qcs
      real*4, dimension(imlon,jmlat) :: qll
      real*4, dimension(:,:,:), allocatable :: dqi,dqj
      integer :: ics,jcs,ill,jll,n,itile
      real*4 :: dx,gradfac
      allocate(dqi(imcub,imcub,6),dqj(imcub,imcub,6))
      dx = 2./real(imcub,kind=4)
      gradfac = 1./(2.*dx)
      do itile=1,6
        do jcs=1,imcub
          dqi(1,jcs,itile) = 0.
          do ics=2,imcub-1
            dqi(ics,jcs,itile) = 
     &           (qcs(ics+1,jcs,itile)-qcs(ics-1,jcs,itile))*gradfac
          enddo
          dqi(imcub,jcs,itile) = 0.
        enddo
        dqj(:,1,itile) = 0.
        do jcs=2,imcub-1
          do ics=1,imcub
            dqj(ics,jcs,itile) =
     &           (qcs(ics,jcs+1,itile)-qcs(ics,jcs-1,itile))*gradfac
          enddo
        enddo
        dqj(:,imcub,itile) = 0.
      enddo
      qll(:,:) = 0.
      do n=1,ncells
        ics = ijcs(1,n)
        jcs = ijcs(2,n)
        ill = ijll(1,n)
        jll = ijll(2,n)
c todo: if q is non-negative, do a recon pass over cs cells
c and set dqi=dqj=0 in cells in which any qcs+dqi*xp+dqj*yp < 0
        qll(ill,jll) = qll(ill,jll) + (
     &       qcs(ics,jcs,tile(n))
     &      +dqi(ics,jcs,tile(n))*xp(n)
     &      +dqj(ics,jcs,tile(n))*yp(n)
     &       )*xgrid_area(n)
      enddo
      qll(:,:) = qll(:,:)/areall(:,:)
      deallocate(dqi,dqj)
      return
      end subroutine regrid2_2d

      subroutine regrid2_3d(qcs,qll,nl)
      integer :: nl
      real*4, dimension(nl,imcub,imcub,6) :: qcs
      real*4, dimension(nl,imlon,jmlat) :: qll
      real*4, dimension(:,:,:,:), allocatable :: dqi,dqj
      integer :: ics,jcs,ill,jll,n,itile
      real*4 :: dx,gradfac
      allocate(dqi(nl,imcub,imcub,6),dqj(nl,imcub,imcub,6))
      dx = 2./real(imcub,kind=4)
      gradfac = 1./(2.*dx)
      do itile=1,6
        do jcs=1,imcub
          dqi(:,1,jcs,itile) = 0.
          do ics=2,imcub-1
            dqi(:,ics,jcs,itile) = 
     &           (qcs(:,ics+1,jcs,itile)-qcs(:,ics-1,jcs,itile))*gradfac
          enddo
          dqi(:,imcub,jcs,itile) = 0.
        enddo
        dqj(:,:,1,itile) = 0.
        do jcs=2,imcub-1
          do ics=1,imcub
            dqj(:,ics,jcs,itile) =
     &           (qcs(:,ics,jcs+1,itile)-qcs(:,ics,jcs-1,itile))*gradfac
          enddo
        enddo
        dqj(:,:,imcub,itile) = 0.
      enddo
      qll(:,:,:) = 0.
      do n=1,ncells
        ics = ijcs(1,n)
        jcs = ijcs(2,n)
        ill = ijll(1,n)
        jll = ijll(2,n)
c todo: if q is non-negative, do a recon pass over cs cells
c and set dqi=dqj=0 in cells in which any qcs+dqi*xp+dqj*yp < 0
        qll(:,ill,jll) = qll(:,ill,jll) + (
     &       qcs(:,ics,jcs,tile(n))
     &      +dqi(:,ics,jcs,tile(n))*xp(n)
     &      +dqj(:,ics,jcs,tile(n))*yp(n)
     &       )*xgrid_area(n)
      enddo
      do jll=1,jmlat
        do ill=1,imlon
          qll(:,ill,jll) = qll(:,ill,jll)/areall(ill,jll)
        enddo
      enddo
      deallocate(dqi,dqj)
      return
      end subroutine regrid2_3d

      end module regrid_to_latlon_mod
