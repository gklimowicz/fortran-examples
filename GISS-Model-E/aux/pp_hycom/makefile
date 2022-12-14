SHELL=/bin/sh

#FC=    ifort -assume byterecl -convert big_endian
#LD=    ifort -assume byterecl -convert big_endian

NC=/usr/local/other/netcdf/3.6.2_intel-10.1.013

LEVEL_OPT = -g  -traceback -C
#LEVEL_OPT = -O2 -traceback

#CPPDEFS = -DHYCOM_RES_387x360x26
CPPDEFS = -DHYCOM_RES_387x360x32
#CPPDEFS = -DHYCOM_RES_359x360x26
#CPPDEFS = -DHYCOM_RES_359x360x32
CPPDEFS = -D$(RESOLUTION)

FC=     ifort $(LEVEL_OPT) $(CPPDEFS) -convert big_endian -openmp
FC=     ifort $(LEVEL_OPT) $(CPPDEFS) -convert big_endian

LD=     ifort $(LEVEL_OPT) $(CPPDEFS)  -convert big_endian -I${NC}/include -openmp
LD=     ifort $(LEVEL_OPT) $(CPPDEFS)  -convert big_endian -I${NC}/include

LDFLAGS =
CPPFLAGS =
LIBS = -L${NC}/lib -lnetcdf

avg1_obj=hycom_dimen.o const_proc.o hycom_arrays.o avg.o getdat_nc.o gtdpth.o meshsz.o
avg2_obj=coagflx.o extrct.o prtmsk.o bigrid.o restep.o findmx.o psmoo.o thruflow.o zebra.o

avg_obj = $(avg1_obj) $(avg2_obj)
avg_src = $(avg1_obj:.o=.F90) $(avg2_obj: .o=.f)

ov2_obj=const_proc.o hycom_dimen.o hycom_arrays.o pltov_mon.o getdat_nc.o \
	extrct.o prtmsk.o bigrid.o gtdpth.o meshsz.o restep.o \
	findmx.o thruflow.o zebra.o

latlonz1_obj=hycom_dimen.o const_proc.o hycom_arrays.o cpler.o \
	latlonz.o getdat_nc.o gtdpth.o meshsz.o
latlonz2_obj= bigrid.o extrct.o restep.o grdtrns.o findmx.o zebra.o

latlonz_obj = $(latlonz1_obj) $(latlonz2_obj)
latlonz_src = $(latlonz1_obj:.o=.F90) $(latlonz2_obj:.o=.f)

get2d_obj=const_proc.o hycom_dimen.o get2d.o
get2d_src = $(get2d_obj:.o=.F90)

avg:    $(avg_obj)
	$(LD) -o avg $(avg_obj) $(LIBS) $(LDFLAGS)

pltov_mon:  $(ov2_obj)
	$(LD) -o pltov_mon $(ov2_obj) $(LIBS) $(LDFLAGS)

latlonz:    $(latlonz_obj)
	$(LD) -o latlonz $(latlonz_obj) $(LIBS) $(LDFLAGS)

get2d:  $(get2d_obj)
	$(LD) -o get2d $(get2d_obj) $(LIBS) $(LDFLAGS)

avg.o: avg.F90 const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c avg.F90
bigrid.o: bigrid.f hycom_dimen.F90
	$(LD) $(LDFLAGS) -c bigrid.f
coagflx.o: coagflx.f const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c coagflx.f
const_proc.o: const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c const_proc.F90
cpler.o: cpler.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c cpler.F90
extrct.o: extrct.f const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c extrct.f
get2d.o: get2d.F90 const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c get2d.F90
getdat_nc.o: getdat_nc.F90 const_proc.F90 hycom_dimen.F90 hycom_arrays.F90
	$(LD) $(LDFLAGS) -c getdat_nc.F90
gtdpth.o: gtdpth.F90 const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c gtdpth.F90
grdtrns.o: grdtrns.f
	$(LD) $(LDFLAGS) -c grdtrns.f
hycom_dimen.o: hycom_dimen.F90
	$(LD) $(LDFLAGS) -c hycom_dimen.F90
hycom_arrays.o: hycom_arrays.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c hycom_arrays.F90
latlonz.o: latlonz.F90 const_proc.F90 hycom_dimen.F90 hycom_arrays.F90
	$(LD) $(LDFLAGS) -c latlonz.F90
meshsz.o: meshsz.F90 hycom_arrays.F90 hycom_dimen.F90 const_proc.F90
	$(LD) $(LDFLAGS) -c meshsz.F90
prtmsk.o: prtmsk.f
	$(LD) $(LDFLAGS) -c prtmsk.f
psmoo.o: psmoo.f hycom_dimen.F90
	$(LD) $(LDFLAGS) -c psmoo.f
restep.o: restep.f const_proc.F90 hycom_dimen.F90
	$(LD) $(LDFLAGS) -c restep.f
zebra.o: zebra.f
	$(LD) $(LDFLAGS) -c zebra.f


tar: makefile input readme $(latlonz_src) $(get2d_src) $(avg_src)
	tar -c.F90 proc.tar mak.F90ile input readme $(latlonz_src) $(get2d_src) $(avg_src)

clean:
	\rm *.o *.mod
