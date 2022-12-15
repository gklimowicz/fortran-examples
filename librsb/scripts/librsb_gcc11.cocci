// usage: spatch --local-includes --sp-file librsb_gcc11.cocci rsb_krnl_bcss_spmv_u.c > librsb_gcc11.patch && patch  < librsb_gcc11.patch
// description: apply coccinelle rule to insert corrective GCC pragmas
// Patch for a https://gcc.gnu.org/bugzilla/show_bug.cgi?id=103995 workaround.
@pragma_inject@
identifier i =~ "rsb__BCSR_spmv_sasa_double_complex_[CH]__t[NTC]_r1_c1_uu_s[HS]_dE_uG";
type T;
@@
+ #pragma GCC push_options
+ #pragma GCC optimize "-O3", "-fno-tree-loop-vectorize"
T i(...)
{
 ...
}
+ #pragma GCC pop_options
