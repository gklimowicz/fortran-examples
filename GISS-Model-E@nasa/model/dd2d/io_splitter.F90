module io_splitter
!@sum io_splitter provides functions for I/O of a rank N>2 array
!@+   as a collection of rank N-1 arrays having distinct netcdf names.
!@+
  !use dist_grid_mod, only : dist_grid
  use dd2d_utils, only : dist_grid
  use pario, only : defvar,read_dist_data,write_dist_data,variable_exists
  implicit none
  private

  public :: defsplitvar,io_splitvar

  interface defsplitvar
     module procedure defsplitvar_2d
     module procedure defsplitvar_3d
     module procedure defsplitvar_4d
  end interface

  interface io_splitvar
     module procedure io_splitvar_2d
     module procedure io_splitvar_3d
     module procedure io_splitvar_4d
  end interface

contains

  subroutine setup_3d(jdim,splitdim,jdim_,d3)
    integer :: jdim,splitdim,jdim_,d3
!
    integer :: td
    if(splitdim < 1 .or. splitdim > 4) then
      call stop_model('bad splitdim',255)
    endif
    if(splitdim < jdim) then
      jdim_ = jdim - 1
    else
      jdim_ = jdim
    endif
    td = mod(5 + splitdim - jdim, 4) + 1 ! rotate to jdim=2
    if(td < 3) call stop_model('bad splitdim or jdim',255)
    d3 = mod(4 - td + jdim, 4) + 1
  end subroutine setup_3d

  subroutine setup_4d(jdim,splitdim,jdim_,d3,d4)
    integer :: jdim,splitdim,jdim_,d3,d4
!
    integer :: td
    if(splitdim < 1 .or. splitdim > 5) then
      call stop_model('bad splitdim',255)
    endif
    if(splitdim < jdim) then
      jdim_ = jdim - 1
    else
      jdim_ = jdim
    endif
    td = mod(6 + splitdim - jdim, 5) + 1 ! rotate to jdim=2
    if(td < 3) call stop_model('bad splitdim or jdim',255)
    d3 = mod(4 - min(td,4) + jdim, 5) + 1
    d4 = mod(6 - max(td,4) + jdim, 5) + 1
    if(d3 > d4) then ! swap
      td = d3; d3 = d4; d4 = td
    endif
  end subroutine setup_4d

  subroutine defsplitvar_2d(grid,fid,prefix,sname,dimstr, &
       splitdim,var,filter,r4_on_disk)
    type(dist_grid) :: grid
    character(len=*) :: prefix,dimstr
    character(len=*), dimension(:) :: sname
    integer :: fid,splitdim
    logical, dimension(:), intent(in), optional :: filter
    logical, optional :: r4_on_disk
    real*8, dimension(:,:,:) :: var
!
    integer :: il,iu,jl,ju,n,jdim,nsplit
    real*8, dimension(:,:), allocatable :: tr1
    character(len=128) :: varstr
    logical :: r4_on_disk_

    nsplit = size(var,splitdim)

    if(size(sname) < nsplit) then
      call stop_model('size(sname) != size(var)',255)
    endif
    if(present(filter)) then
      if(size(filter) < nsplit) then
        call stop_model('size(filter) != size(var)',255)
      endif
    endif
    if(splitdim.ne.1 .and. splitdim.ne.3) then
      call stop_model('bad splitdim',255)
    endif

    if(splitdim == 1) then
      jdim = 3
    else
      jdim = 2
    endif

    il = lbound(var,jdim-1)
    iu = ubound(var,jdim-1)
    jl = lbound(var,jdim)
    ju = ubound(var,jdim)

    allocate(tr1(il:iu,jl:ju))

    r4_on_disk_ = .false.
    if(present(r4_on_disk)) r4_on_disk_ = r4_on_disk

    do n=1,nsplit
      if(present(filter)) then
        if(.not.filter(n)) cycle
      endif
      if(trim(sname(n)).eq.'unused') cycle ! hack hack hack
      varstr = trim(prefix)//'_'//trim(sname(n))//trim(dimstr)
      call defvar(grid,fid,tr1,varstr,r4_on_disk=r4_on_disk_)
    enddo

  end subroutine defsplitvar_2d

  subroutine io_splitvar_2d(grid,fid,actstr,prefix,sname, &
       splitdim,var,filter)
    type(dist_grid) :: grid
    character(len=*) :: actstr,prefix
    character(len=*), dimension(:) :: sname
    integer :: fid,splitdim
    logical, dimension(:), intent(in), optional :: filter
    real*8, dimension(:,:,:) :: var
!
    integer :: il,iu,jl,ju,n,jdim,nsplit
    real*8, dimension(:,:), allocatable :: tr1
    character(len=128) :: vname

    nsplit = size(var,splitdim)
    if(size(sname) < nsplit) then
      call stop_model('size(sname) != size(var)',255)
    endif
    if(present(filter)) then
      if(size(filter) < nsplit) then
        call stop_model('size(filter) != size(var)',255)
      endif
    endif
    if(splitdim.ne.1 .and. splitdim.ne.3) then
      call stop_model('bad splitdim',255)
    endif

    if(splitdim == 1) then
      jdim = 3
    else
      jdim = 2
    endif

    il = lbound(var,jdim-1)
    iu = ubound(var,jdim-1)
    jl = lbound(var,jdim)
    ju = ubound(var,jdim)

    allocate(tr1(il:iu,jl:ju))

    do n=1,nsplit
      if(present(filter)) then
        if(.not.filter(n)) cycle
      endif
      if(trim(sname(n)).eq.'unused') cycle ! hack hack hack
      vname = trim(prefix)//'_'//trim(sname(n))
      if(trim(actstr).eq.'write') then
        if(splitdim == 1) then
          tr1 = var(n,:,:)
        else
          tr1 = var(:,:,n)
        endif
        call write_dist_data(grid, fid, trim(vname), tr1)
      elseif(trim(actstr).eq.'read') then
        call read_dist_data(grid, fid, trim(vname), tr1)
        if(variable_exists(grid,fid,trim(vname))) then
          if(splitdim == 1) then
            var(n,:,:) = tr1
          else
            var(:,:,n) = tr1
          endif
        endif
      endif
    enddo

  end subroutine io_splitvar_2d

  subroutine defsplitvar_3d(grid,fid,prefix,sname,dimstr, &
       jdim,splitdim,var,filter,r4_on_disk)
    type(dist_grid) :: grid
    character(len=*) :: prefix,dimstr
    character(len=*), dimension(:) :: sname
    integer :: fid,jdim,splitdim
    logical, dimension(:), intent(in), optional :: filter
    logical, optional :: r4_on_disk
    real*8, dimension(:,:,:,:) :: var
!
    integer :: il,iu,jl,ju,d3,jdim_,n,nsplit
    real*8, dimension(:,:,:), allocatable :: tr1
    character(len=128) :: varstr
    logical :: r4_on_disk_

    call setup_3d(jdim,splitdim,jdim_,d3)
    nsplit = size(var,splitdim)
    if(size(sname) < nsplit) then
      call stop_model('size(sname) != size(var)',255)
    endif
    if(present(filter)) then
      if(size(filter) < nsplit) then
        call stop_model('size(filter) != size(var)',255)
      endif
    endif
    d3 = size(var,d3)

    il = lbound(var,jdim-1)
    iu = ubound(var,jdim-1)
    jl = lbound(var,jdim)
    ju = ubound(var,jdim)

    if(jdim_ == 2) then
      allocate(tr1(il:iu,jl:ju,d3))
    else
      allocate(tr1(d3,il:iu,jl:ju))
    endif

    r4_on_disk_ = .false.
    if(present(r4_on_disk)) r4_on_disk_ = r4_on_disk

    do n=1,nsplit
      if(present(filter)) then
        if(.not.filter(n)) cycle
      endif
      if(trim(sname(n)).eq.'unused') cycle ! hack hack hack
      varstr = trim(prefix)//'_'//trim(sname(n))//trim(dimstr)
      call defvar(grid,fid,tr1,varstr,r4_on_disk=r4_on_disk_)
    enddo

  end subroutine defsplitvar_3d

  subroutine io_splitvar_3d(grid,fid,actstr,prefix,sname, &
       jdim,splitdim,var,filter)
    type(dist_grid) :: grid
    character(len=*) :: actstr,prefix
    character(len=*), dimension(:) :: sname
    integer :: fid,jdim,splitdim
    logical, dimension(:), intent(in), optional :: filter
    real*8, dimension(:,:,:,:) :: var
!
    integer :: il,iu,jl,ju,d3,d3l,d3u,jdim_,n,nsplit
    real*8, dimension(:,:,:), allocatable :: tr1
    character(len=128) :: vname


    call setup_3d(jdim,splitdim,jdim_,d3)
    nsplit = size(var,splitdim)
    if(size(sname) < nsplit) then
      call stop_model('size(sname) != size(var)',255)
    endif
    if(present(filter)) then
      if(size(filter) < nsplit) then
        call stop_model('size(filter) != size(var)',255)
      endif
    endif

    d3l = lbound(var,d3)
    d3u = ubound(var,d3)

    il = lbound(var,jdim-1)
    iu = ubound(var,jdim-1)
    jl = lbound(var,jdim)
    ju = ubound(var,jdim)

    if(jdim_ == 2) then
      allocate(tr1(il:iu,jl:ju,d3l:d3u))
    else
      allocate(tr1(d3l:d3u,il:iu,jl:ju))
    endif

    do n=1,nsplit
      if(present(filter)) then
        if(.not.filter(n)) cycle
      endif
      if(trim(sname(n)).eq.'unused') cycle ! hack hack hack
      vname = trim(prefix)//'_'//trim(sname(n))
      if(trim(actstr).eq.'write') then
        if(splitdim == 1) then
          tr1 = var(n,:,:,:)
        elseif(splitdim == 2) then
          tr1 = var(:,n,:,:)
        elseif(splitdim == 3) then
          tr1 = var(:,:,n,:)
        else
          tr1 = var(:,:,:,n)
        endif
        call write_dist_data(grid, fid, trim(vname), tr1, jdim=jdim_)
      elseif(trim(actstr).eq.'read') then
        call read_dist_data(grid, fid, trim(vname), tr1, jdim=jdim_)
        if(variable_exists(grid,fid,trim(vname))) then
          if(splitdim == 1) then
            var(n,:,:,:) = tr1
          elseif(splitdim == 2) then
            var(:,n,:,:) = tr1
          elseif(splitdim == 3) then
            var(:,:,n,:) = tr1
          else
            var(:,:,:,n) = tr1
          endif
        endif
      endif
    enddo

  end subroutine io_splitvar_3d

  subroutine defsplitvar_4d(grid,fid,prefix,sname,dimstr, &
       jdim,splitdim,var,filter,r4_on_disk)
    type(dist_grid) :: grid
    character(len=*) :: prefix,dimstr
    character(len=*), dimension(:) :: sname
    integer :: fid,jdim,splitdim
    logical, dimension(:), intent(in), optional :: filter
    logical, optional :: r4_on_disk
    real*8, dimension(:,:,:,:,:) :: var
!
    integer :: il,iu,jl,ju,d3,d4,jdim_,n,nsplit
    real*8, dimension(:,:,:,:), allocatable :: tr1
    character(len=128) :: varstr
    logical :: r4_on_disk_

    call setup_4d(jdim,splitdim,jdim_,d3,d4)
    nsplit = size(var,splitdim)
    if(size(sname) < nsplit) then
      call stop_model('size(sname) != size(var)',255)
    endif
    if(present(filter)) then
      if(size(filter) < nsplit) then
        call stop_model('size(filter) != size(var)',255)
      endif
    endif

    d3 = size(var,d3)
    d4 = size(var,d4)

    il = lbound(var,jdim-1)
    iu = ubound(var,jdim-1)
    jl = lbound(var,jdim)
    ju = ubound(var,jdim)

    if(jdim_ == 2) then
      allocate(tr1(il:iu,jl:ju,d3,d4))
    elseif(jdim_ == 3) then
      allocate(tr1(d3,il:iu,jl:ju,d4))
    else
      allocate(tr1(d3,d4,il:iu,jl:ju))
    endif

    r4_on_disk_ = .false.
    if(present(r4_on_disk)) r4_on_disk_ = r4_on_disk

    do n=1,nsplit
      if(present(filter)) then
        if(.not.filter(n)) cycle
      endif
      if(trim(sname(n)).eq.'unused') cycle ! hack hack hack
      varstr = trim(prefix)//'_'//trim(sname(n))//trim(dimstr)
      call defvar(grid,fid,tr1,varstr,r4_on_disk=r4_on_disk_)
    enddo

  end subroutine defsplitvar_4d

  subroutine io_splitvar_4d(grid,fid,actstr,prefix,sname, &
       jdim,splitdim,var,filter)
    type(dist_grid) :: grid
    character(len=*) :: actstr,prefix
    character(len=*), dimension(:) :: sname
    integer :: fid,jdim,splitdim
    logical, dimension(:), intent(in), optional :: filter
    real*8, dimension(:,:,:,:,:) :: var
!
    integer :: il,iu,jl,ju,d3,d4,d3l,d3u,d4l,d4u,jdim_,n,nsplit
    real*8, dimension(:,:,:,:), allocatable :: tr1
    character(len=128) :: vname

    call setup_4d(jdim,splitdim,jdim_,d3,d4)
    nsplit = size(var,splitdim)
    if(size(sname) < nsplit) then
      call stop_model('size(sname) != size(var)',255)
    endif
    if(present(filter)) then
      if(size(filter) < nsplit) then
        call stop_model('size(filter) != size(var)',255)
      endif
    endif

    d3l = lbound(var,d3)
    d3u = ubound(var,d3)
    d4l = lbound(var,d4)
    d4u = ubound(var,d4)
    
    il = lbound(var,jdim-1)
    iu = ubound(var,jdim-1)
    jl = lbound(var,jdim)
    ju = ubound(var,jdim)

    if(jdim_ == 2) then
      allocate(tr1(il:iu,jl:ju,d3l:d3u,d4l:d4u))
    elseif(jdim_ == 3) then
      allocate(tr1(d3l:d3u,il:iu,jl:ju,d4l:d4u))
    else
      allocate(tr1(d3l:d3u,d4l:d4u,il:iu,jl:ju))
    endif

    do n=1,nsplit
      if(present(filter)) then
        if(.not.filter(n)) cycle
      endif
      if(trim(sname(n)).eq.'unused') cycle ! hack hack hack
      vname = trim(prefix)//'_'//trim(sname(n))
      if(trim(actstr).eq.'write') then
        if(splitdim == 1) then
          tr1 = var(n,:,:,:,:)
        elseif(splitdim == 2) then
          tr1 = var(:,n,:,:,:)
        elseif(splitdim == 3) then
          tr1 = var(:,:,n,:,:)
        elseif(splitdim == 4) then
          tr1 = var(:,:,:,n,:)
        else
          tr1 = var(:,:,:,:,n)
        endif
        call write_dist_data(grid, fid, trim(vname), tr1, jdim=jdim_)
      elseif(trim(actstr).eq.'read') then
        call read_dist_data(grid, fid, trim(vname), tr1, jdim=jdim_)
        if(variable_exists(grid,fid,trim(vname))) then
          if(splitdim == 1) then
            var(n,:,:,:,:) = tr1
          elseif(splitdim == 2) then
            var(:,n,:,:,:) = tr1
          elseif(splitdim == 3) then
            var(:,:,n,:,:) = tr1
          elseif(splitdim == 4) then
            var(:,:,:,n,:) = tr1
          else
            var(:,:,:,:,n) = tr1
          endif
        endif
      endif
    enddo

  end subroutine io_splitvar_4d

end module io_splitter
