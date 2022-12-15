#!/usr/bin/awk -f
# Patch for a https://gcc.gnu.org/bugzilla/show_bug.cgi?id=103995 workaround.
# $ diff rsb_krnl_bcss_spmv_u.c  <( awk -f librsb_gcc11.awk < rsb_krnl_bcss_spmv_u.c )
BEGIN { wp=1; of=""; }
/^rsb_err_t rsb__BCSR_spmv_sasa_double_complex_[CH]__t[NTC]_r1_c1_uu_s[HS]_dE_uG/ { wp=0; system("echo '#pragma GCC push_options';echo '#pragma GCC optimize \"-O3\", \"-fno-tree-loop-vectorize\"'"); }
/^}/ {  if ( wp == 0 ) { system("echo '#pragma GCC pop_options'"); } ;wp=1;}
/.*/ { print; }
