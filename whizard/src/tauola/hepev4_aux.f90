module hepev4_aux

  use hep_common

  implicit none
  private

  integer :: nhep_original
  integer :: nhepev4_call=0

contains

  subroutine hepev4_fill
    integer :: i
    integer :: j
    double precision  :: sumdiff
    nhepev4_call = nhepev4_call + 1
    spinlh(:,1:nhep_original) = 0.0D0
    icolorflowlh(:,1:nhep_original) = 0
    loop_i: do i = 1,nhep_original
       loop_j: do j = 1,nup
          check_consist: if(idup(j) == idhep(i)) then
             sumdiff = sum((pup(1:3,j) - phep(1:3,i))**2)
             check_sumdiff: if(sumdiff < 1.d-6) then
                spinlh(3,i) = spinup(j)
                icolorflowlh(:,i) = icolup(:,j)
               exit loop_j
             end if check_sumdiff
          end if check_consist
       end do loop_j
    end do loop_i

  end subroutine hepev4_fill

  subroutine hepev4_update(tauspin_pyjets)
    double precision, dimension(:), intent(in) :: tauspin_pyjets
    spinlh(1:2,nhep_original+1:nhep) = 0.0D0
    spinlh(3,nhep_original+1:nhep) = tauspin_pyjets(nhep_original+1:nhep)
    icolorflowlh(:,nhep_original+1:nhep) = 0
  end subroutine hepev4_update

end module hepev4_aux
