#! /bin/sh
ver=2
rel=0
mode=fast
root=`pwd`
indir=${root}/input
tmpdir=${root}/tmp
outroot=${root}/output
outdir=${outroot}/v${ver}/r${rel}
# acc="sband tesla nlc_a nlc_b nlc_c"
acc="sband tesla nlc_a"
xmkdir () {
  for d in "$@"; do
    mkdir $d 2>/dev/null || true
  done
}
rm -fr ${tmpdir}
xmkdir ${outroot} ${outroot}/v${ver} ${outdir} ${tmpdir}
cd ${tmpdir}
cat /dev/null >${outdir}/Params.f
for a in $acc; do
  case "$a" in
    *1000*) energy=TEV1;;
         *) energy=GEV500;;
  esac
  sed -n -e 's/;//' -e 's/lumi_[eg][eg]=//p' \
    ${indir}/${a}_${mode}/${a}_${mode} >Totals
  echo 5.0D0 >Powers
  cp ${indir}/${a}_${mode}/lumidiff-??.dat .
  ${root}/circe1_int.bin
  rm -fr ${outdir}/${a}_${mode}
  xmkdir ${outdir}/${a}_${mode}
  sed -e "s/@ENERGY@/$energy/g" \
      -e "s/@ACC@/`echo $a | tr a-z A-Z`/g" Parameters \
    >>${outdir}/Params.f
done
cd ${root}
rm -fr ${tmpdir}
exit 0
