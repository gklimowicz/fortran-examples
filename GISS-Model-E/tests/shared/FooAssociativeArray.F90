module Foo_mod
   implicit none
   private
   
   public :: Foo

   type Foo
     integer :: value
   end type Foo

end module Foo_mod

#define TYPE Foo
#define PRINT_IT
#define DO_THIS=1

#include <AssociativeArrayTemplate.h>


