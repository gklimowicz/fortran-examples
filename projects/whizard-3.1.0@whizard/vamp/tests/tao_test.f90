program tao_test
  use tao_random_numbers, only: test30 => tao_random_test
  use tao52_random_numbers, only: test52 => tao_random_test
  implicit none
  call test30 ("tmp.tao")
  call test52 ("tmp.tao")
  stop 0
end program tao_test
