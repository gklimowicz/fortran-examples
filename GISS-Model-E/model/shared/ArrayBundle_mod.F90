
module bundle_maker
  implicit none

  private
  public :: make_bundle,make_bundle_lij

contains

  subroutine make_bundle(arr,is,ie,js,je, &
              ptr01,ptr02,ptr03,ptr04, &
        ptr05,ptr06,ptr07,ptr08,ptr09, &
        ptr10,ptr11,ptr12,ptr13,ptr14, &
        ptr15,ptr16,ptr17,ptr18,ptr19, &
        ptr20,ptr21,ptr22,ptr23,ptr24, &
        ptr25,ptr26,ptr27,ptr28,ptr29  &
       )
    integer :: is,ie,js,je
    real*8, pointer :: arr(:,:,:)
    real*8, pointer, dimension(:,:), optional :: &
                ptr01,ptr02,ptr03,ptr04, &
          ptr05,ptr06,ptr07,ptr08,ptr09, &
          ptr10,ptr11,ptr12,ptr13,ptr14, &
          ptr15,ptr16,ptr17,ptr18,ptr19, &
          ptr20,ptr21,ptr22,ptr23,ptr24, &
          ptr25,ptr26,ptr27,ptr28,ptr29
    integer :: n
!
! the size of the bundle is equal to the number of optional
! arguments present
!
    n=0
    if(present(ptr01)) n=n+1
    if(present(ptr02)) n=n+1
    if(present(ptr03)) n=n+1
    if(present(ptr04)) n=n+1
    if(present(ptr05)) n=n+1
    if(present(ptr06)) n=n+1
    if(present(ptr07)) n=n+1
    if(present(ptr08)) n=n+1
    if(present(ptr09)) n=n+1
    if(present(ptr10)) n=n+1
    if(present(ptr11)) n=n+1
    if(present(ptr12)) n=n+1
    if(present(ptr13)) n=n+1
    if(present(ptr14)) n=n+1
    if(present(ptr15)) n=n+1
    if(present(ptr16)) n=n+1
    if(present(ptr17)) n=n+1
    if(present(ptr18)) n=n+1
    if(present(ptr19)) n=n+1
    if(present(ptr20)) n=n+1
    if(present(ptr21)) n=n+1
    if(present(ptr22)) n=n+1
    if(present(ptr23)) n=n+1
    if(present(ptr24)) n=n+1
    if(present(ptr25)) n=n+1
    if(present(ptr26)) n=n+1
    if(present(ptr27)) n=n+1
    if(present(ptr28)) n=n+1
    if(present(ptr29)) n=n+1

!
! allocate the bundle
!
    if(.not. associated(arr)) then
      allocate(arr(is:ie,js:je,n))
      arr = 0.
    endif

    if(n> 0) call set_ptr(arr(:,:, 1), ptr01 ,is,js)
    if(n> 1) call set_ptr(arr(:,:, 2), ptr02 ,is,js)
    if(n> 2) call set_ptr(arr(:,:, 3), ptr03 ,is,js)
    if(n> 3) call set_ptr(arr(:,:, 4), ptr04 ,is,js)
    if(n> 4) call set_ptr(arr(:,:, 5), ptr05 ,is,js)
    if(n> 5) call set_ptr(arr(:,:, 6), ptr06 ,is,js)
    if(n> 6) call set_ptr(arr(:,:, 7), ptr07 ,is,js)
    if(n> 7) call set_ptr(arr(:,:, 8), ptr08 ,is,js)
    if(n> 8) call set_ptr(arr(:,:, 9), ptr09 ,is,js)
    if(n> 9) call set_ptr(arr(:,:,10), ptr10 ,is,js)
    if(n>10) call set_ptr(arr(:,:,11), ptr11 ,is,js)
    if(n>11) call set_ptr(arr(:,:,12), ptr12 ,is,js)
    if(n>12) call set_ptr(arr(:,:,13), ptr13 ,is,js)
    if(n>13) call set_ptr(arr(:,:,14), ptr14 ,is,js)
    if(n>14) call set_ptr(arr(:,:,15), ptr15 ,is,js)
    if(n>15) call set_ptr(arr(:,:,16), ptr16 ,is,js)
    if(n>16) call set_ptr(arr(:,:,17), ptr17 ,is,js)
    if(n>17) call set_ptr(arr(:,:,18), ptr18 ,is,js)
    if(n>18) call set_ptr(arr(:,:,19), ptr19 ,is,js)
    if(n>19) call set_ptr(arr(:,:,20), ptr20 ,is,js)
    if(n>20) call set_ptr(arr(:,:,21), ptr21 ,is,js)
    if(n>21) call set_ptr(arr(:,:,22), ptr22 ,is,js)
    if(n>22) call set_ptr(arr(:,:,23), ptr23 ,is,js)
    if(n>23) call set_ptr(arr(:,:,24), ptr24 ,is,js)
    if(n>24) call set_ptr(arr(:,:,25), ptr25 ,is,js)
    if(n>25) call set_ptr(arr(:,:,26), ptr26 ,is,js)
    if(n>26) call set_ptr(arr(:,:,27), ptr27 ,is,js)
    if(n>27) call set_ptr(arr(:,:,28), ptr28 ,is,js)
    if(n>28) call set_ptr(arr(:,:,29), ptr29 ,is,js)
    return
  end subroutine make_bundle

  subroutine make_bundle_lij(arr,l,is,ie,js,je, &
              ptr01,ptr02,ptr03,ptr04, &
        ptr05,ptr06,ptr07,ptr08,ptr09, &
        ptr10,ptr11,ptr12,ptr13,ptr14, &
        ptr15,ptr16,ptr17,ptr18,ptr19 &
       )
    integer :: l,is,ie,js,je
    real*8, pointer :: arr(:,:,:,:)
    real*8, pointer, dimension(:,:,:), optional :: &
                ptr01,ptr02,ptr03,ptr04, &
          ptr05,ptr06,ptr07,ptr08,ptr09, &
          ptr10,ptr11,ptr12,ptr13,ptr14, &
          ptr15,ptr16,ptr17,ptr18,ptr19
    integer :: n
!
! the size of the bundle is equal to the number of optional
! arguments present
!
    n=0
    if(present(ptr01)) n=n+1
    if(present(ptr02)) n=n+1
    if(present(ptr03)) n=n+1
    if(present(ptr04)) n=n+1
    if(present(ptr05)) n=n+1
    if(present(ptr06)) n=n+1
    if(present(ptr07)) n=n+1
    if(present(ptr08)) n=n+1
    if(present(ptr09)) n=n+1
    if(present(ptr10)) n=n+1
    if(present(ptr11)) n=n+1
    if(present(ptr12)) n=n+1
    if(present(ptr13)) n=n+1
    if(present(ptr14)) n=n+1
    if(present(ptr15)) n=n+1
    if(present(ptr16)) n=n+1
    if(present(ptr17)) n=n+1
    if(present(ptr18)) n=n+1
    if(present(ptr19)) n=n+1

!
! allocate the bundle
!
    if(.not. associated(arr)) then
      allocate(arr(l,is:ie,js:je,n))
      arr = 0.
    endif

    if(n> 0) call set_ptr_lij(arr(:,:,:, 1), ptr01 ,is,js)
    if(n> 1) call set_ptr_lij(arr(:,:,:, 2), ptr02 ,is,js)
    if(n> 2) call set_ptr_lij(arr(:,:,:, 3), ptr03 ,is,js)
    if(n> 3) call set_ptr_lij(arr(:,:,:, 4), ptr04 ,is,js)
    if(n> 4) call set_ptr_lij(arr(:,:,:, 5), ptr05 ,is,js)
    if(n> 5) call set_ptr_lij(arr(:,:,:, 6), ptr06 ,is,js)
    if(n> 6) call set_ptr_lij(arr(:,:,:, 7), ptr07 ,is,js)
    if(n> 7) call set_ptr_lij(arr(:,:,:, 8), ptr08 ,is,js)
    if(n> 8) call set_ptr_lij(arr(:,:,:, 9), ptr09 ,is,js)
    if(n> 9) call set_ptr_lij(arr(:,:,:,10), ptr10 ,is,js)
    if(n>10) call set_ptr_lij(arr(:,:,:,11), ptr11 ,is,js)
    if(n>11) call set_ptr_lij(arr(:,:,:,12), ptr12 ,is,js)
    if(n>12) call set_ptr_lij(arr(:,:,:,13), ptr13 ,is,js)
    if(n>13) call set_ptr_lij(arr(:,:,:,14), ptr14 ,is,js)
    if(n>14) call set_ptr_lij(arr(:,:,:,15), ptr15 ,is,js)
    if(n>15) call set_ptr_lij(arr(:,:,:,16), ptr16 ,is,js)
    if(n>16) call set_ptr_lij(arr(:,:,:,17), ptr17 ,is,js)
    if(n>17) call set_ptr_lij(arr(:,:,:,18), ptr18 ,is,js)
    if(n>18) call set_ptr_lij(arr(:,:,:,19), ptr19 ,is,js)

    return
  end subroutine make_bundle_lij

  subroutine set_ptr(arr,ptr,is,js)
    integer :: is,js
    real*8, target :: arr(is:,js:)
    real*8, pointer :: ptr(:,:)
#ifdef INTEL18_WITH_TRAPS
    ! Handle pointer bounds not being inherited from target.
    ! This syntax was not used initially because not all
    ! compilers accepted it when this routine was written.
    ! Compiler acceptance to be re-assessed.
    ptr(is:,js:) => arr
#else
    ptr => arr
#endif
    return
  end subroutine set_ptr

  subroutine set_ptr_lij(arr,ptr,is,js)
    integer :: is,js
    real*8, target :: arr(:,is:,js:)
    real*8, pointer :: ptr(:,:,:)
#ifdef INTEL18_WITH_TRAPS
    ! Handle pointer bounds not being inherited from target.
    ! This syntax was not used initially because not all
    ! compilers accepted it when this routine was written.
    ! Compiler acceptance to be re-assessed.
    ptr(1:,is:,js:) => arr
#else
    ptr => arr
#endif
    return
  end subroutine set_ptr_lij

end module bundle_maker

module ArrayBundle_mod
!@sum subroutines for packing/unpacking a selection of arrays
!@+   to/from a single bundle
!@auth I. Aleinov, D. Gueyffier
  implicit none
  private

  public lookup_str
  public ab_init, ab_add, ab_bundle, ab_unbundle, ab_copy, get_bounds

  integer, parameter :: N_LOOKUP_MAX=256

  type lookup_record
    integer :: lm,km
    real*8, pointer :: src(:), dest(:)
    real*8, pointer :: src_w(:,:), dest_w(:,:)
  end type lookup_record

  type lookup_str
    private
    integer :: si_0, si_1, sj_0, sj_1  ! source bounds
    integer :: di_0, di_1, dj_0, dj_1  ! destination bounds
    integer  :: n_lookup=0
    type (lookup_record) :: lr(N_LOOKUP_MAX)
  end type lookup_str

contains

  subroutine get_bounds(lstr, si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1)
!@sum accessor, returns domain bouds for both source and dest. grids
    implicit none
    type (lookup_str) :: lstr
    integer :: si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1

    si_0 = lstr%si_0
    si_1 = lstr%si_1
    sj_0 = lstr%sj_0
    sj_1 = lstr%sj_1
    di_0 = lstr%di_0
    di_1 = lstr%di_1
    dj_0 = lstr%dj_0
    dj_1 = lstr%dj_1

  end subroutine get_bounds

  subroutine ab_init( lstr, si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1 )
    implicit none
    type (lookup_str) :: lstr
    integer :: si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1 

    lstr%si_0 = si_0
    lstr%si_1 = si_1
    lstr%sj_0 = sj_0
    lstr%sj_1 = sj_1
    lstr%di_0 = di_0
    lstr%di_1 = di_1
    lstr%dj_0 = dj_0
    lstr%dj_1 = dj_1

    lstr%n_lookup = 0
  end subroutine ab_init


  subroutine ab_add( lstr, src, dest, shp, flag, src_w, dest_w )
    implicit none
    type (lookup_str) :: lstr
    real*8, dimension(*), target :: src, dest
    integer :: shp(:)
    character*(*) :: flag
    real*8, dimension(:,:), target, optional :: src_w, dest_w

    integer :: si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1 
    integer :: sim, sjm, dim, djm, lm, km

    lstr%n_lookup = lstr%n_lookup + 1
    if ( lstr%n_lookup > N_LOOKUP_MAX ) call stop_model("ab_add: increase N_LOOKUP_MAX", 255)

    si_0 = lstr%si_0
    si_1 = lstr%si_1
    sj_0 = lstr%sj_0
    sj_1 = lstr%sj_1
    di_0 = lstr%di_0
    di_1 = lstr%di_1
    dj_0 = lstr%dj_0
    dj_1 = lstr%dj_1

    sim = si_1 - si_0 + 1
    sjm = sj_1 - sj_0 + 1
    dim = di_1 - di_0 + 1
    djm = dj_1 - dj_0 + 1

    select case( flag )
    case('ij')
      lm = 1
      km = 1
    case('lij')
      lm = shp(1)
      km = 1
    case('ijk')
      lm = 1
      km = shp(3)
    case('lijk')
      lm = shp(1)
      km = shp(4)
    case default
      call stop_model("ab_add: unexpected flag", 255)
    end select

    lstr%lr(lstr%n_lookup)%lm = lm
    lstr%lr(lstr%n_lookup)%km = km

    lstr%lr(lstr%n_lookup)%src => src(1:lm*sim*sjm*km)
    lstr%lr(lstr%n_lookup)%dest => dest(1:lm*dim*djm*km)

    if ( present(src_w) .and. present(dest_w) ) then
      lstr%lr(lstr%n_lookup)%src_w => src_w(:,:)
      lstr%lr(lstr%n_lookup)%dest_w => dest_w(:,:)
    else
      if ( present(src_w) ) call stop_model("ab_add: use both weights or none", 255)
      nullify( lstr%lr(lstr%n_lookup)%src_w )
      nullify( lstr%lr(lstr%n_lookup)%dest_w )
    endif

  end subroutine ab_add


  subroutine ab_bundle( lstr, buf_s, buf_d )
    !USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
    implicit none
    type (lookup_str) :: lstr
    real*8, dimension(:,:,:), pointer :: buf_s, buf_d

    integer :: si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1 
    integer im,jm,km,lm
    integer i,j,k,l,m,n,ind


    si_0 = lstr%si_0
    si_1 = lstr%si_1
    sj_0 = lstr%sj_0
    sj_1 = lstr%sj_1
    di_0 = lstr%di_0
    di_1 = lstr%di_1
    dj_0 = lstr%dj_0
    dj_1 = lstr%dj_1

    im = si_1 - si_0 + 1
    jm = sj_1 - sj_0 + 1

    n = 0
    do k=1,lstr%n_lookup
      n = n+lstr%lr(k)%km*lstr%lr(k)%lm
    enddo

    allocate( buf_s(n, si_0:si_1, sj_0:sj_1) )
    allocate( buf_d(n, di_0:di_1, dj_0:dj_1) )

    n = 0
    do m = 1,lstr%n_lookup
      km = lstr%lr(m)%km
      lm = lstr%lr(m)%lm
      do k=0,km-1
        do l=0,lm-1
          n = n+1

          do j=0,jm-1
            do i=0,im-1
              ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
              buf_s(n,i+si_0,j+sj_0) = lstr%lr(m)%src(ind)
            enddo
          enddo

          if ( associated(lstr%lr(m)%src_w) ) then
            do j=0,jm-1
              do i=0,im-1
                ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                buf_s(n,i+si_0,j+sj_0) = buf_s(n,i+si_0,j+sj_0) &
                     & * lstr%lr(m)%src_w(i+1,j+1)
              enddo
            enddo
          endif

        enddo
      enddo
    enddo

  end subroutine ab_bundle


  subroutine ab_unbundle( lstr, buf_s, buf_d )
    implicit none
    type (lookup_str) :: lstr
    real*8, dimension(:,:,:), pointer :: buf_s, buf_d

    integer :: di_0, di_1, dj_0, dj_1
    integer im,jm,km,lm
    integer i,j,k,l,m,n,ind

    di_0 = lstr%di_0
    di_1 = lstr%di_1
    dj_0 = lstr%dj_0
    dj_1 = lstr%dj_1

    im = di_1 - di_0 + 1
    jm = dj_1 - dj_0 + 1

    n = 0
    do m = 1,lstr%n_lookup
      km = lstr%lr(m)%km
      lm = lstr%lr(m)%lm
      do k=0,km-1
        do l=0,lm-1
          n = n+1

          do j=0,jm-1
            do i=0,im-1
              ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
              lstr%lr(m)%dest(ind) = buf_d(n,i+di_0,j+dj_0)
            enddo
          enddo

        enddo
      enddo
    enddo

    n = 0
    do m = 1,lstr%n_lookup
      km = lstr%lr(m)%km
      lm = lstr%lr(m)%lm
      do k=0,km-1
        do l=0,lm-1
          n = n+1

          if ( associated(lstr%lr(m)%dest_w ) ) then
            do j=0,jm-1
              do i=0,im-1
                ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                if ( lstr%lr(m)%dest_w(i+1,j+1) .ne. 0.d0 ) then
                  lstr%lr(m)%dest(ind) = lstr%lr(m)%dest(ind) &
                       &                   / lstr%lr(m)%dest_w(i+1,j+1)
                else
                  lstr%lr(m)%dest(ind) = 0.d0
                endif
              enddo
            enddo
          endif

        enddo
      enddo
    enddo

    deallocate( buf_s )
    deallocate( buf_d )     

  end subroutine ab_unbundle


  subroutine ab_copy( lstr )
    implicit none
    type (lookup_str) :: lstr

    integer :: si_0, si_1, sj_0, sj_1, di_0, di_1, dj_0, dj_1
    integer im,jm,km,lm,nm
    integer m

    si_0 = lstr%si_0
    si_1 = lstr%si_1
    sj_0 = lstr%sj_0
    sj_1 = lstr%sj_1
    di_0 = lstr%di_0
    di_1 = lstr%di_1
    dj_0 = lstr%dj_0
    dj_1 = lstr%dj_1

    im = si_1 - si_0 + 1
    jm = sj_1 - sj_0 + 1

    ! make sure we copy arrays of the same size
    if ( im .ne. di_1 - di_0 + 1 .or. jm .ne. dj_1 - dj_0 + 1 )  &
         & call stop_model("ab_copy: dims of arrays differ", 255)

    do m = 1,lstr%n_lookup
      km = lstr%lr(m)%km
      lm = lstr%lr(m)%lm
      nm = lm*im*jm*km
      lstr%lr(m)%dest(1:nm) = lstr%lr(m)%src(1:nm)
    enddo

  end subroutine ab_copy

end module ArrayBundle_mod

