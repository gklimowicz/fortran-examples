#!/usr/bin/ksh
#
FC=ifort

NetCDF=/usr/local/other/netcdf/4.0.1_gnu
UDUNITS=/usr/local/other/udunits/2.0.3_intel-10.1.021
HDF=/usr/local/other/hdf5/1.8.3_serialGNU
UUID=/usr/local/other/uuid/1.6.2

#CMOR=/usr/local/other/cmor/2.0rc5_GNU              # GNU   compiler
CMOR=/usr/local/other/cmor/2.0rc5_intel-10.1.017   # Intel compiler

#Dont forget include next lines at your profile:
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/other/netcdf/4.0.1_gnu/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/other/udunits/2.0.3_intel-10.1.021/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/other/hdf5/1.8.3_serialGNU/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/other/uuid/1.6.2/lib

fileFort=$1
fileExe=${fileFort%.*}.exe 
rm -rf local_subs.mod ${fileExe} 
mkdir ./Test

rm -rf ${fileExe}
${FC} ${fileFort} -convert big_endian -O2 -traceback   \
     -I${CMOR} -I${CMOR}/include -L${CMOR}/lib -lcmor  \
     -I${NetCDF}/include -L${NetCDF}/lib -lnetcdf      \
     -I${UDUNITS}/include -L${UDUNITS}/lib -ludunits2  \
     -I${HDF}/include -L${HDF}/lib -lhdf5 -lhdf5_hl -lm -lz \
     -I${UUID}/include -L${UUID}/lib -luuid -o ${fileExe}






