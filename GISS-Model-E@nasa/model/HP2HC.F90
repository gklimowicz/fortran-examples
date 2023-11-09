module hp2hc_types
	type double_2_ptr
		real*8, dimension(:,:), pointer :: p
	end type double_2_ptr
end module hp2hc_types



module hp2hc
	use glint2_modele

	implicit none

	! The matrix used to convert height points to height classes
	type(glint2_modele_matrix) :: hp_to_hc

contains

! -----------------------------------------------------------
subroutine multiply_by_bundle(mat, ival, oval)
use hp2hc_types
implicit none

	type(glint2_modele_matrix), intent(in) :: mat
	type(double_2_ptr), dimension(:,:), intent(in) :: ival ! (nvars, nhc)
	type(double_2_ptr), dimension(:,:) :: oval	! (nvars, nhc)
	! --------------------------------

	integer :: i,j,k
	integer :: nvars, nhc, nmat
	integer :: ri, rj, rk
	integer :: ci, cj, ck

print *,'BEGIN multiply_by_bundle'
	nvars = size(oval,1)
	nhc = size(oval,2)
	nmat = size(mat%vals)	! Number of non-zero elements in the matrix

! Dummy: just copy over
do j=1,nvars
	do k=1,nhc
		oval(j,k)%p = ival(j,k)%p
	end do
end do
return


print *,'multiply_by_bundle: 1'
	! Clear the arrays
	do j=1,nvars
		do k=1,nhc
			oval(j,k)%p = 0
		end do
	end do

!do i=1,5
!	print *,'hp2hc',i,mat%rows_i(i),mat%rows_j(i),mat%rows_k(i),mat%cols_i(i),mat%cols_j(i),mat%cols_k(i)
!end do


!print *,'multiply_by_bundle: 2'
	! Do the sparse matrix multiply on our bundled values
	do i=1,nmat
		! Index for height points (row in matrix)
		ri = mat%rows_i(i)		! 1..im
		rj = mat%rows_j(i)		! 1..jm
		rk = mat%rows_k(i)		! Height class

		! Index for height classes (column in matrix)
		ci = mat%cols_i(i)
		cj = mat%cols_j(i)
		ck = mat%cols_k(i)		! Height class

!print *,'multiply_by_bundle',i
!print *,'row',ri,rj,rk
!print *,'col',ci,cj,ck
		do j=1,nvars
!print *,'oval 1:', lbound(oval,1), ubound(oval,1)
!print *,'oval 2:', lbound(oval,2), ubound(oval,2)
!print *,'ival 1:', lbound(ival,1), ubound(ival,1)
!print *,'ival 2:', lbound(ival,2), ubound(ival,2)
!print *,'multiply_by_bundle',i,j,ci,cj,ri,rj,oval(j,ck)%p,ival(j,rk)%p
			oval(j,ck)%p(ci,cj) = oval(j,ck)%p(ci,cj) + &
				mat%vals(i) * ival(j,rk)%p(ri,rj)
		end do
	end do

	! Copy HP=1 to HC=1
	! (this is not explicitly encoded in the matrix)
	do j=1,nvars
		oval(j,1)%p = ival(j,1)%p
	end do

print *,'END multiply_by_bundle'
end subroutine multiply_by_bundle

! ----------------------------------------------------
subroutine bundle_init_li(atmglas, vars)
use hp2hc_types
USE EXCHANGE_TYPES, only : atmgla_xchng_vars
	type(atmgla_xchng_vars), dimension(:) :: atmglas
	type(double_2_ptr), dimension(:,:), allocatable :: vars

	integer :: nhc, iv, ih
	integer :: nvars = 4	! Must match below
#ifdef TRACERS_ON
	ntm = size(atmglas%gtracer,1)	! Or any variable with ntm dimension
	nvars = nvars + 1 * ntm		! Must match below
#endif
	print *,'bundle_init_li'
	nhc = size(atmglas,1)

	allocate(vars(nvars,nhc))

	do ih=1,nhc
		iv = 0
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%snow
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp2
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtempr

		! ntm = atmglas(ih)%ntm
		! gtracer(ntm, im, jm)
		! truno(ntm, im, jm)
#ifdef TRACERS_ON
		do j=1,ntm
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)%gtracer(j,:,:)
		end do
#endif
	end do

	! Now iv should equal nvars
	if (iv /= nvars) then
		write (6,*) 'bundle_init_li: nvars = ',nvars,'must be set to',iv
		call stop_model( 'bundle_init_li::nvars too big or small', 255 )
	end if
end subroutine bundle_init_li
! ---------------------------------------------------------
subroutine bundle_precip_li(atmglas, vars)
use hp2hc_types
USE EXCHANGE_TYPES, only : atmgla_xchng_vars
	type(atmgla_xchng_vars), dimension(:) :: atmglas
	type(double_2_ptr), dimension(:,:), allocatable :: vars

	integer :: nhc, iv, ih
	integer :: nvars = 6	! Must match below
#ifdef TRACERS_ON
	ntm = size(atmglas%gtracer,1)	! Or any variable with ntm dimension
	nvars = nvars + 2 * ntm		! Must match below
#ifdef TRACERS_WATER
	nvars = nvars + ntm
#endif
#endif
	nhc = size(atmglas,1)
	print *,'bundle_precip_li nhc = ',nhc
	allocate(vars(nvars,nhc))

	iv = 0
	do ih=1,nhc
		iv = 0
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%runo
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%implm
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%implh
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp2
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtempr

		! ntm = atmglas(ih)%ntm
		! gtracer(ntm, im, jm)
		! truno(ntm, im, jm)
#ifdef TRACERS_ON
		do j=1,ntm
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)%gtracer(j,:,:)
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)%truno(j,:,:)
#ifdef TRACERS_WATER
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)%implt(j,:,:)
#endif
		end do
#endif
	end do

	! Now iv should equal nvars
	if (iv /= nvars) then
		write (6,*) 'bundle_precip_li: nvars = ',nvars,'must be set to',iv
		call stop_model( 'bundle_precip_li::nvars too big or small', 255 )
	end if
end subroutine bundle_precip_li
! ---------------------------------------------------------
subroutine bundle_ground_li(atmglas, vars)
use hp2hc_types
USE EXCHANGE_TYPES, only : atmgla_xchng_vars
	type(atmgla_xchng_vars), dimension(:) :: atmglas
	type(double_2_ptr), dimension(:,:), allocatable :: vars

	integer :: nhc, iv, ih
	integer :: nvars = 11	! Must match below
#ifdef TRACERS_ON
	ntm = size(atmglas%gtracer,1)	! Or any variable with ntm dimension
	nvars = nvars + 2 * ntm		! Must match below
#ifdef TRACERS_WATER
	nvars = nvars + 1 * ntm
#endif
#endif
	print *,'bundle_ground_li'
	nhc = size(atmglas,1)

	allocate(vars(nvars,nhc))

	iv = 0
	do ih=1,nhc
		iv = 0
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%e1
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%snow
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%snowfr
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%snowdp
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%implm
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%implh
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%runo
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%eruno
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp2
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtempr

		! ntm = atmglas(ih)%ntm
		! gtracer(ntm, im, jm)
		! truno(ntm, im, jm)
#ifdef TRACERS_ON
		do j=1,ntm
#ifdef TRACERS_WATER
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)implt(j,:,:)
#endif
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)truno(j,:,:)
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)gtracer(j,:,:)
		end do
#endif
	end do

	! Now iv should equal nvars
	if (iv /= nvars) then
		write (6,*) 'bundle_ground_li: nvars = ',nvars,'must be set to',iv
		call stop_model( 'bundle_ground_li::nvars too big or small', 255 )
	end if
end subroutine bundle_ground_li
! ---------------------------------------------------------
subroutine bundle_surface_landice(atmglas, vars)
use hp2hc_types
USE EXCHANGE_TYPES, only : atmgla_xchng_vars
	type(atmgla_xchng_vars), dimension(:) :: atmglas
	type(double_2_ptr), dimension(:,:), allocatable :: vars

	integer :: nhc, iv, ih
	integer :: nvars = 17	! Must match below
#ifdef TRACERS_WATER
	ntm = size(atmglas%gtracer,1)	! Or any variable with ntm dimension
	nvars = nvars + 2 * ntm		! Must match below
#endif
	print *,'bundle_surface_landice'
	nhc = size(atmglas,1)

	allocate(vars(nvars,nhc))

	iv = 0
	do ih=1,nhc
		iv = 0
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%dth1
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%dq1
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%e0
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%e1
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%evapor
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemp
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%gtemps
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%latht
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%sensht
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%solar
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%trheat
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%tgrnd
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%tgr4
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%dmua
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%dmva
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%uflux1
		iv = iv + 1
		vars(iv,ih)%p => atmglas(ih)%vflux1

		! ntm = atmglas(ih)%ntm
		! gtracer(ntm, im, jm)
		! truno(ntm, im, jm)
#ifdef TRACERS_WATER
		do j=1,ntm
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)trsrfflx(j,:,:)
			iv = iv + 1
			vars(iv,ih)%p => atmglas(ih)trevapor(j,:,:)
		end do
#endif
	end do

	! Now iv should equal nvars
	if (iv /= nvars) then
		write (6,*) 'bundle_surface_landice: nvars = ',nvars,'must be set to',iv
		call stop_model( 'bundle_surface_landice::nvars too big or small', 255 )
	end if
end subroutine bundle_surface_landice
! ---------------------------------------------------------
subroutine bundle_hp_to_hc(bundler, atmglas_hp, atmglas)
	USE EXCHANGE_TYPES, only : atmgla_xchng_vars
	use hp2hc_types
	implicit none
	type(atmgla_xchng_vars), dimension(:) :: atmglas_hp
	type(atmgla_xchng_vars), dimension(:) :: atmglas
	interface
		subroutine bundler(atmglas, vars)
			USE EXCHANGE_TYPES, only : atmgla_xchng_vars
			use hp2hc_types
			type(atmgla_xchng_vars), dimension(:) :: atmglas
			type(double_2_ptr), dimension(:,:), allocatable :: vars
		end subroutine bundler
	end interface
	! ----------------------------------------------
	type(double_2_ptr), dimension(:,:), allocatable :: hp, hc

!	! Dummy: Just copy over
!	atmglas(:) = atmglas_hp(:)
!	return

print *,'BEGIN bundle_hp_to_hc'
	call bundler(atmglas_hp, hp)
	call bundler(atmglas, hc)
	
	call multiply_by_bundle(hp_to_hc, hp, hc)

print *,'END bundle_hp_to_hc'
end subroutine bundle_hp_to_hc
! ----------------------------------------------------





#if 0
! ----------------------------------------------------
#define HP_TO_HC_X(varbundle) \
subroutine hp_to_hc_/**/varbundle(atmglas_hp, atmglas) \
	USE EXCHANGE_TYPES, only : atmgla_xchng_vars \
	implicit none \
	! type(glint2_modele_matrix), intent(in) :: hp_to_hc \
	!(1-min(nhc,2)/2:nhc) ! glacial ice \
	type(atmgla_xchng_vars), allocatable, dimension(:) :: atmglas_hp \
	type(atmgla_xchng_vars), allocatable, dimension(:) :: atmglas \
 \
	call bundle_arrays_/**/varbundle(atmglas_hp, hp) \
	call bundle_arrays_/**/varbundle(atmglas, hc) \
 \
	call multiply_by_bundle(hp_to_hc, hp, hc) \
end subroutine hp_to_hc_/**/varbundle

HP_TO_HC_X(init_li)
HP_TO_HC_X(precip_li)
HP_TO_HC_X(ground_li)
HP_TO_HC_X(surface_landice)



! ----------------------------------------------------
#endif



end module hp2hc
