! #define ASSERT_ALWAYS(expr,str...) If (.not.(expr))Call assert_(__LINE__,__FILE__, ##str)
#define ASSERT_ALWAYS(expr,str) If (.not.(expr))Call assert_(__LINE__,__FILE__, str)


#ifndef NDEBUG
#define ASSERT(expr,str) ASSERT_ALWAYS(expr, str)
#else
#define ASSERT(expr,str)
#endif

#ifndef NDEBUG
#define WARN(expr,str) If (.not. (expr)) Call warn_(__LINE__,__FILE__, str)
#else
#define WARN(expr,str)
#endif

