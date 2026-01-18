#define TYPE Attribute
#define USE_MODULE AbstractAttribute_mod
#define TYPE_NAME AbstractAttribute
#define HAS_PRINT
!#define WRAPPED_TYPE

#include "AssociativeArrayTemplate.h"

#define VALUE_TYPE AbstractAttribute
#define ASSOCIATIVE_ARRAY_TYPE AttributeAssociativeArray
#undef ITERATOR_TYPE
#define HASH_TYPE AttributeHashMap

#include "HashMapTemplate.h"

