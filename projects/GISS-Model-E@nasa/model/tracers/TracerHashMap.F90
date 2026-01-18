#define HAS_PRINT
#define TYPE Tracer
#define TYPE_NAME Tracer
#include "../shared/AssociativeArrayTemplate.h"

#define VALUE_TYPE Tracer
#undef ITERATOR_TYPE
#define ITERATOR_TYPE TracerIterator
#include "../shared/HashMapTemplate.h"

