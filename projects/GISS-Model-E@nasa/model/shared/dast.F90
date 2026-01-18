module vector_integer_mod
#define _entry integer
#define EQUAL_DEFINED
#include "containers/vector.fh"
end module vector_integer_mod

module vector_real_mod
#define _entry real
#define EQUAL_DEFINED
#include "containers/vector.fh"
end module vector_real_mod

module vector_real8_mod
#define _entry real(kind=kind(1.d0))
#define EQUAL_DEFINED
#include "containers/vector.fh"
end module vector_real8_mod

module set_integer_mod
#define _entry integer
#include "containers/set.fh"
end module set_integer_mod

module map_integer_integer_mod
#define _key integer
#define _value integer
#include "containers/map.fh"
end module map_integer_integer_mod
