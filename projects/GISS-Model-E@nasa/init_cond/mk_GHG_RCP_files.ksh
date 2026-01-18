#!/bin/ksh
echo "Before executing this program"
echo "- Copy the original files into a subdirectory called: original_files"
echo "- compile GHG_IPCC_to_GISS.f into GHG_IPCC_to_GISS.exe"
echo "Copies of the orig. files are e.g. on /archive/u/rruedy/AR5_GHG_forcings/original_files"
echo
echo "ready to go ? (y=yes)"
read reply
if [[ $reply != [yY]* ]]
then echo "no action" ; exit ; fi

# Create header files
cat <<EOF > head_RCP85
Global-Mean Greenhouse Gas Mixing Ratios Used for IPCC_AR5 runs
RCP 8.5 -------------------------------------------------------
       CO2    N2O    CH4   CFC-11 CFC-12 others
Year   ppm    ppm    ppm     ppb    ppb    ppb
-------------------------------------------------
EOF

cat <<EOF > head_RCP6
Global-Mean Greenhouse Gas Mixing Ratios Used for IPCC_AR5 runs
RCP 6.0 -------------------------------------------------------
       CO2    N2O    CH4   CFC-11 CFC-12 others
Year   ppm    ppm    ppm     ppb    ppb    ppb
-------------------------------------------------
EOF

cat <<EOF > head_RCP6SCP6TO45
Global-Mean Greenhouse Gas Mixing Ratios Used for IPCC_AR5 runs
RCP 6.0 -> 4.5 ------------------------------------------------
       CO2    N2O    CH4   CFC-11 CFC-12 others
Year   ppm    ppm    ppm     ppb    ppb    ppb
-------------------------------------------------
EOF

cat <<EOF > head_RCP45
Global-Mean Greenhouse Gas Mixing Ratios Used for IPCC_AR5 runs
RCP 4.5 -------------------------------------------------------
       CO2    N2O    CH4   CFC-11 CFC-12 others
Year   ppm    ppm    ppm     ppb    ppb    ppb
-------------------------------------------------
EOF

cat <<EOF > head_RCP3PD
Global-Mean Greenhouse Gas Mixing Ratios Used for IPCC_AR5 runs
RCP3 PD (2.6) -------------------------------------------------
       CO2    N2O    CH4   CFC-11 CFC-12 others
Year   ppm    ppm    ppm     ppb    ppb    ppb
-------------------------------------------------
EOF

cd original_files ; gunzip *gz
for x in *DAT
do a=$( grep THISFILE_FIRSTDATAROW $x ) ; nhead=${a##* }
   nh=${nhead%%,*} ; echo $x ${nh}  ; skip=$(( $nh - 1 ))
   ../GHG_IPCC_to_GISS.exe $x $skip
   cat ../head_${x%_MIDYR_CONC.DAT} fort.2 > ../GHG_${x%_MIDYR_CONC.DAT}.txt ; rm -f fort.2
done
gzip * ; rm -f ../head_RCP*
