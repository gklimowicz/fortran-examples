#!/usr/bin/env python
#
# This script add height classes to exiting input files.

import sys
import argparse
import io
import struct
import numpy
import re
import giss.io.giss
from giss.util import *
import netCDF4
import numpy as np
from odict import odict
import os.path
import snowdrift

ifname = "TOPO"


# -------------------------------------------------------------
# Reads info about a variable, no matther whether that variable is
# sitting in a netCDF file, or was read in from a GISS-format file
# @return (val, sdims, dtype)
#   val = The np.array values of the variable
#   sdim = Tuple of names of the variable's dimensions
#   dtype = Type of the variable
def gread(handle, var_name) :
	if isinstance(handle, odict) :	# Just fetch the variable from topo
		return handle[var_name]
	else :		# We have a netcdf handle
		val = read_ncvar(handle, var_name)
		var = handle.variables[var_name]
		return (val, var.dimensions, var.dtype)

# Reads info about a variable, no matther whether that variable is
# sitting in a netCDF file, or was read in from a GISS-format file
# @return (sdims, shape, dtype)
#   shape = Dimensions of the np.array that would hold the variable
#   sdim = Tuple of names of the variable's dimensions
#   dtype = Type of the variable
def gread_dims(handle, name) :
	if isinstance(handle, odict) :	# Just fetch the variable from topo
		var = handle[name]		# (value, sdims) pairs
		sdims = var[1]
		shape = var[0].shape
		dtype = var[2]
	else :		# We have a netcdf handle
		var = handle.variables[name]
		sdims = var.dimensions
		shape = var.shape
		dtype = var.dtype
	return (sdims, shape, dtype)
# -------------------------------------------------------------
# Creates a new variable with an extra dimension prepended
# @param varpair = (val, sdims, dtype) input
# @return (ovar, odims).  Space is allocated, user must fill it in
def prepend_dim(varpair, dimname, dimval) :
	ivar = varpair[0]
	isdims = varpair[1]
	idtype = varpair[2]
	odims = (dimname,) + isdims
	ovar = np.zeros((dimval,) + ivar.shape)
	return (ovar, odims, idtype)
# -------------------------------------------------------------
# Given a bunch of variables that have been determined, writes
# them to a netCDF file.  Defines the dimensions, etc. as needed
# @param ofname Name of netCDF file to write
# @param wvars = list of variables to write
#        wvars[i] = (handle, name)
#        handle = Way to look up the variable (TOPO array, netCDF handle, etc)
#        name = Name of variable
def write_netcdf(ofname, wvars) :
	nc = netCDF4.Dataset(ofname, 'w')

	# Collect together the dimension names and lengths for each variable
	wdims = odict()
	for vv in wvars :		# (handle, name) pairs
		handle = vv[0]
		name = vv[1]
		sdims, shape, dtype = gread_dims(handle, name)

		for i in range(0, len(sdims)) :
			sdim = sdims[i]
			shp = shape[i]
			if sdim in wdims :
				if shp != wdims[sdim] :
					raise Exception('Inconsistent dimension for %s (%d vs %d)' % (sdim, shp, wdims[sdim]))
			else :
				wdims[sdim] = shp
	print wdims

	# Define the dimensions
	for name, val in wdims.iteritems() :
		print name, val
		nc.createDimension(name, val)

	# Define the variables
	for vv in wvars :		# (handle, name) pairs
		handle = vv[0]
		name = vv[1]
		sdims, shape, dtype = gread_dims(handle, name)
		nc.createVariable(name, dtype, sdims)

	# Copy in the data
	for vv in wvars :		# (handle, name) pairs
		handle = vv[0]
		name = vv[1]
		val, sdims, dtype = gread(handle, name)
		nc.variables[name][:] = val[:]

	nc.close()
# -------------------------------------------------------------

# Height-classify a TOPO and GIC file
# @param nhc Number of height classes to use
def hcgic(TOPO_iname, GIC_iname, TOPO_oname, GIC_oname, hcfunc, nhc) :

	ivars = {}

	# ================ Open the GIC file
	# Open the GIC file and read the things we need to change
	ncgic = netCDF4.Dataset(GIC_iname, 'r')

#	snowli_dims = ncgic.variables['snowli'].dimensions
#	print ncgic.dimensions[snowli_dims[0]]
#	print ncgic.dimensions[ncgic.variables['snowli'].dimensions[1]]

	# ================= Read the TOPO file
	topo = odict()
	for rec in giss.io.giss.reader(TOPO_iname) :
		val = np.zeros(rec.data.shape)	# Promote to double
		val[:] = rec.data[:]
		item = (val, (u'jm', u'im'), 'f8')
		topo[rec.var.lower()] = item

	# ================= Read out dimensions
#	for key in topo.iterkeys() :
#		print key
#
#	shape = topo['snowli'][0].shape
#	jm = shape[0]
#	im = shape[1]
#	n1 = im*jm

	# ================ Height-classify variables
	ovars = hcfunc(topo, ncgic, nhc)

	# ============= Write out new TOPO

	# Collect variables to write
	ovars_remain = odict(ovars)
	wvars = []	# Holds tuples (source, name)
	for tl in topo.iteritems() :
		name = tl[0]
		# Fetch variable from correct place
		if name in ovars :
			wvars.append((ovars, name))
			del ovars_remain[name]
		else :
			wvars.append((topo, name))

	# Define and write the variables
	write_netcdf(TOPO_oname, wvars)


	# ============= Write out new GIC
	# Collect variables to write (from original GIC file)
	wvars = []	# Holds tuples (source, name)
	for name in ncgic.variables :
		# Fetch variable from correct place
		if name in ovars :
			wvars.append((ovars, name))
			del ovars_remain[name]
		else :
			wvars.append((ncgic, name))

	# Add on new variable(s) we've created
	for name in ovars_remain.iterkeys() :
		wvars.append((ovars, name))

	# Define and write the variables
	write_netcdf(GIC_oname, wvars)

# ================================================================
# ----------------------------------------------------------------
# Height-classify the variables...
# @param topo TOPO input file.  The following variables are needed:
#        fgice, fgrnd, zatmo
# @param ncgic GIC input file  The following variables are needed:
#        tlandi, snowli
# @param nhc Number of height classes to add
# @return The varialbes created (in a dict):
#         fgice, fgrnd, elevhc, tlandi, snowli, fhc
def dummy_hc(topo, ncgic, nhc) :
	# Copy over fgice and fgrnd
	# The values of these might have to change a bit, but
	# not the structure.
	ovars = odict()
	ovars['fgice'] = gread(topo, 'fgice')
	ovars['fgrnd'] = gread(topo, 'fgrnd')

	# Allocate but only read things we need now
	zatmo = gread(topo, 'zatmo')
	tlandi = gread(ncgic, 'tlandi')
	snowli = gread(ncgic, 'snowli')
	ovars['elevhc'] = prepend_dim(zatmo, 'nhc', nhc)
	ovars['tlandi'] = prepend_dim(tlandi, 'nhc', nhc)
	ovars['snowli'] = prepend_dim(snowli, 'nhc', nhc)
	ovars['fhc'] = (np.zeros(ovars['snowli'][0].shape), ovars['snowli'][1], 'f8')

	# Make height-classified versions of vars by copying
	for i in range(0,nhc) :
		ovars['elevhc'] [0] [i,:] = zatmo[0][:]
		ovars['tlandi'] [0] [i,:] = tlandi[0][:]
		ovars['snowli'] [0] [i,:] = snowli[0][:]
		ovars['fhc'] [0] [:] = 1.0 / float(nhc)

	return ovars
# -------------------------------------------------------------
# Height-classify the variables...
# @param topo TOPO input file.  The following variables are needed:
#        fgice, fgrnd, zatmo
# @param ncgic GIC input file  The following variables are needed:
#        tlandi, snowli
# @param nhc Number of height classes to add
# @return The varialbes created (in a dict):
#         fgice, fgrnd, elevhc, tlandi, snowli, fhc
def test3_hc(topo, ncgic, nhc) :
	# Copy over fgice and fgrnd
	# The values of these might have to change a bit, but
	# not the structure.
	ovars = odict()
	ovars['fgice'] = gread(topo, 'fgice')
	ovars['fgrnd'] = gread(topo, 'fgrnd')

	fgice = ovars['fgice'][0]
	fgrnd = ovars['fgrnd'][0]

	# Allocate but only read things we need now
	zatmo_t = gread(topo, 'zatmo')
	tlandi_t = gread(ncgic, 'tlandi')
	snowli_t = gread(ncgic, 'snowli')
	ovars['elevhc'] = prepend_dim(zatmo_t, 'nhc', nhc)
	ovars['tlandi'] = prepend_dim(tlandi_t, 'nhc', nhc)
	ovars['snowli'] = prepend_dim(snowli_t, 'nhc', nhc)
	ovars['fhc'] = (np.zeros(ovars['snowli'][0].shape), ovars['snowli'][1], 'f8')

	elevhc = ovars['elevhc'][0]
	tlandi = ovars['tlandi'][0]
	snowli = ovars['snowli'][0]
	fhc = ovars['fhc'][0]


	# Make height-classified versions of vars by copying
	for ihc in range(0,nhc) :
		elevhc[ihc,:] = zatmo_t[0][:]
		tlandi[ihc,:] = tlandi_t[0][:]
		snowli[ihc,:] = snowli_t[0][:]
		fhc[:] = 1.0 / float(nhc)

	# ---------------------------------------------
	jm = snowli.shape[1]
	im = snowli.shape[2]
	print fgrnd.shape, jm, im
	for j in range(0,jm) :
		for i in range(0, im) :
			if fgice[j,i] == 0 : continue
			elevhc[0,j,i] = max(0, zatmo_t[0][j,i] - 300)
			elevhc[1,j,i] = zatmo_t[0][j,i]
			elevhc[2,j,i] = zatmo_t[0][j,i] + 300
#			elevhc[0,j,i] = 0
#			elevhc[1,j,i] = 1000
#			elevhc[2,j,i] = 2000


	fhc[0,:] = .25
	fhc[1,:] = .5
	fhc[2,:] = .25

	return ovars
# -------------------------------------------------------------
# This is moved to giss.snowdrift package
# # @return .mask2, .elevation2, .sd, .n1
# def read_searise_ice(overlap_fname, searise_fname, n1) :
# 	ret = {}
# 
# 	# =============== Read stuff from ice grid (mask2, elevation2)
# 	# fname = os.path.
# # oin(data_root, 'searise/Greenland_5km_v1.1.nc')
# 	print 'Opening ice data file %s' % searise_fname
# 	searise_nc = netCDF4.Dataset(searise_fname)
# 	mask2 = get_landmask(searise_nc)
# 	ret['mask2'] = mask2
# 	topg = np.array(searise_nc.variables['topg'], dtype='d').flatten('C')
# 	thk = np.array(searise_nc.variables['thk'], dtype='d').flatten('C')
# 	elevation2 = topg + thk
# 	ret['elevation2'] = elevation2
# 	searise_nc.close()
# 
# 	# ========= Set up height class categories
# 	tops = np.array([200,400,700,1000,1300,1600,2000,2500,3000,10000], dtype='d')
# #	tops = np.array([10000], dtype='d')
# 	nhc = tops.shape[0]
# 	height_max1 = np.tile(tops, (n1,1))		# Produces an n1 x nhc array
# 
# 	# ================ Read overlap matrix
# 	# See p. 14 of "The CESM Land Ice Model: Documentation and User's Guide"
# 	# by William Lipscomb (June 2010)
# 	print 'overlap_fname = %s' % overlap_fname
# 	sd = snowdrift.Snowdrift(overlap_fname)
# 	sd.init(elevation2, mask2, height_max1)
# 	grid1 = sd.grid1()
# 	print 'Loaded grid1, n1=%d' % grid1.n
# 
# 	ret['nhc'] = nhc
# 	ret['sd'] = sd
# 	ret['grid1'] = grid1
# #	ret['n1'] = grid1.n1
# 
# 	return giss.util.Struct(ret)

# -------------------------------------------------------------------
def get_landmask(searise_nc) :
# landcover:ice_sheet = 4 ;
# landcover:land = 2 ;
# landcover:local_ice_caps_not_connected_to_the_ice_sheet = 3 ;
# landcover:long_name = "Land Cover" ;
# landcover:no_data = 0 ;
# landcover:ocean = 1 ;
# landcover:standard_name = "land_cover" ;
	mask2 = np.array(searise_nc.variables['landcover'], dtype=np.int32).flatten('C')
	mask2 = np.where(mask2==4,np.int32(1),np.int32(0))
	return mask2
# ---------------------------
# Height-classify the variables...
# @param topo TOPO input file.  The following variables are needed:
#        fgice, fgrnd, zatmo
# @param ncgic GIC input file  The following variables are needed:
#        tlandi, snowli
# @param nhc Number of height classes to add
# @return The varialbes created (in a dict):
#         fgice, fgrnd, elevhc, tlandi, snowli, fhc
# overlap_fname = 'searise_ll_overlap-4x5-5.nc'
def overlap_hc(overlap_fname, topo, ncgic, nhc) :

	data_root = os.environ['SNOWDRIFT_FIG_DATA']
	jm = topo['fgice'][0].shape[0]
	im = topo['fgice'][0].shape[1]
	n1 = jm * im


	# ============ Start off with a "default" version of the
	# height-classified vars
	# Copy over fgice and fgrnd
	# The values of these might have to change a bit, but
	# not the structure.
	ovars = odict()
#	ovars['fgice'] = gread(topo, 'fgice')
#	ovars['fgrnd'] = gread(topo, 'fgrnd')

#	fgice = ovars['fgice'][0]
#	fgrnd_orig = ovars['fgrnd'][0]

	# Allocate but only read things we need now
	zatmo_t = gread(topo, 'zatmo')
	tlandi_t = gread(ncgic, 'tlandi')
	snowli_t = gread(ncgic, 'snowli')
	ovars['tlandi'] = prepend_dim(tlandi_t, 'nhc', nhc)
	ovars['snowli'] = prepend_dim(snowli_t, 'nhc', nhc)
	snowli = ovars['snowli'][0]
	tlandi = ovars['tlandi'][0]

	# Make height-classified versions of vars by copying
	for ihc in range(0,nhc) :
		tlandi[ihc,:] = tlandi_t[0][:]
		snowli[ihc,:] = snowli_t[0][:]

	# ========== Initialize FHC and related arrays
	# These will have nan off-ice-sheet
	fhc1h = np.zeros((nhc, n1))
#	fhc1h[:] = np.nan
	fhc1h[:] = 0
	fhc1h[0,:] = 1

	elev1hxy = np.zeros((nhc, jm, im))
	for ihc in range(0,nhc) :
		elev1hxy[ihc,:,:] = zatmo_t[0][:,:]
	elev1h = elev1hxy.reshape((nhc,n1))
#	elev1h = np.zeros((nhc, n1))
#	elev1h[:] = np.nan

	# These will have a value everywhere, defaulting to value
	# read from ModelE data files.
	fgice_t = gread(topo, 'fgice')		# "original" tuple
	# Convert to double
	fgice_od = np.zeros(fgice_t[0].shape)
	fgice_od[:] = fgice_t[0][:]
	fgice1 = fgice_od.reshape((n1,))
#	print '*******************',fgice_od.dtype, fgice1.dtype

	# ========== Add in FHC, etc. for Greenland grid
	greenland = giss.snowdrift.read_searise_ice(overlap_fname, os.path.join(data_root, 'searise/Greenland_5km_v1.1.nc'), n1)
	if (n1 != greenland.grid1.n) :
		raise Exception('Greenland file has different GCM grid size than GCM files: %d vs %d' % (greenland.grid1.n, n1))
	greenland.sd.compute_fhc(fhc1h, elev1h, fgice1)

	# ========== Add in FHC, etc. for Antarctica grid

	

	# ============ Convert to (im, jm) format and store away
#	if np.not_equal(np.isnan(fhc1h), np.isnan(elev1h)).any() :
#		raise Exception('Sanity Check: Nan extents differ on fhc1h vs. elev1h')
	fhc1hxy = fhc1h.reshape((nhc, jm,im))
	elev1hxy = elev1h.reshape((nhc, jm,im))
	fgice1xy = fgice1.reshape((jm,im))
#	fgice1xy[:,:] = fgice_t[0][:,:]		# TEMP
	ovars['fhc'] = (fhc1hxy, ovars['snowli'][1], 'f8')
	ovars['elevhc'] = (elev1hxy, ovars['snowli'][1], 'f8')
	ovars['fgice'] = (fgice1xy, fgice_t[1], 'f8')

	# ======= Adjust fgrnd accordingly, to keep (FGICE + FGRND) constant
	fgrnd_t = gread(topo, 'fgrnd')		# "original" tuple
	flake_t = gread(topo, 'flake')		# "original" tuple
	focean_t = gread(topo, 'focean')		# "original" tuple
	fgrnd = np.ones(fgrnd_t[0].shape)
#	fgrnd[:,:] = fgice_t[0][:,:] + fgrnd_t[0][:,:] - fgice1xy[:,:]
	# Re-compute invariant, since we promoted from single to double precision
	# (as well as changing fgice)
	fgrnd[:,:] = np.ones(fgrnd_t[0].shape) - focean_t[0] - flake_t[0] - fgice1xy

#	fgrnd[:,:] = fgrnd_t[0][:,:]	#TEMP
	ovars['fgrnd'] = (fgrnd, fgrnd_t[1], 'f8')

	# ======= Compute zatmo and merge into existing
	zatmo1xy = np.nansum(fhc1hxy * elev1hxy, 0)
	mask = np.logical_not(np.isnan(zatmo1xy))
	zatmo = np.zeros(zatmo_t[0].shape)
	zatmo[:,:] = zatmo_t[0][:,:]
	zatmo[mask] = zatmo1xy[mask]
	ovars['zatmo'] = (zatmo, zatmo_t[1], 'f8')

	#np.set_printoptions(threshold='nan')
	#print np.nansum(fhc1hxy,0)		# SHOULD BE 1


#	# =========== Store it away (temporary)
#	nc = netCDF4.Dataset('fhc.nc', 'w')
#	nc.createDimension('nhc', nhc)
#	nc.createDimension('jm', jm)
#	nc.createDimension('im', im)
#	nc.createVariable('fhc', 'f8', ('nhc', 'jm', 'im'))
#	nc.variables['fhc'][:] = fhc1hxy[:]
#	nc.close()
#
#	print 'fhc = ', fhc1hxy[:,42,22]
#
#
#	# Check that dimensions in overlap matrix match with dimensions
#	# in TOPO and GIC files
#	if im*jm != grid1.n :
#		raise Exception('n1=%d from ovlerap file not compatible with n1=%d from TOPO file' % (grid1.n,im*jm))

	return ovars
# =============================================================

dir = '/Users/rpfische/cmrun'
TOPO_in  = os.path.join(dir, 'Z72X46N.cor4_nocasp')
TOPO_out = os.path.join(dir, 'Z72X46N.cor4_nocasp_hc.nc')
GIC_in  = os.path.join(dir, 'GIC.E046D3M20A.1DEC1955.ext.nc')
GIC_out = os.path.join(dir, 'GIC.E046D3M20A.1DEC1955.ext_hc.nc')

# Curry overlap_hc()
class MyOverlapHC() :
	def __init__(self) :
		self.overlap_fname = 'searise_ll_overlap-4x5-5.nc'
	def __call__(self, topo, ncgic, nhc) :
		return overlap_hc(self.overlap_fname, topo, ncgic, nhc)

hcgic(TOPO_in, GIC_in, TOPO_out, GIC_out, MyOverlapHC(), 10)
#hcgic(TOPO_in, GIC_in, TOPO_out, GIC_out, test3_hc, 3)
#hcgic(TOPO_in, GIC_in, TOPO_out, GIC_out,  dummy_hc, 10)
