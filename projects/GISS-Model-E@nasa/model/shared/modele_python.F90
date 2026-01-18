#include "rundeck_opts.h"

module modele_python_mod

implicit none

contains

! ------------------------------------------
! Configures ModelE in ways good when running from Python
subroutine init_for_python()
	use stop_model_mod

#ifdef USE_FEXCEPTION
	call set_stop_model_ptr(stop_model_fexception)
#endif
end subroutine init_for_python
! ------------------------------------------

end module modele_python_mod
