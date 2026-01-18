import netCDF4
from giss.util import *
import giss.io.giss



# Does some sanity checks of height-classified input files to ModelE



def compare_old_new_relative(v0, v1, vname) :
	zero = (v1-v0)/(v0+1e-20)	# Avoid 0/0
	errors = []

	for x in range(0,20) :
		epsilon = np.exp(-x*np.log(10))
		# Count grid cells where invariant is violated by epsiolon
		errors.append(np.sum(np.abs(zero) > epsilon))

	print vname, errors
# -------------------------------------
def compare_old_new_absolute(v0, v1, vname) :
	zero = v1-v0			# We're all scaled 0-1, so this works
	errors = []

	for x in range(0,20) :
		epsilon = np.exp(-x*np.log(10))
		# Count grid cells where invariant is violated by epsiolon
		errors.append(np.sum(np.abs(zero) > epsilon))

	print vname, errors
# -------------------------------------
def test_one(one) :
	errors = []

	for x in range(0,20) :
		epsilon = np.exp(-x*np.log(10))
		# Count grid cells where invariant is violated by epsiolon
		errors.append(np.sum(np.abs(one-1.0) > epsilon))
	return errors
# -------------------------------------
topo0_fname = 'TOPO'
topo_fname = 'TOPO1.nc'
gic_fname = 'GIC1.nc'



topo_nc = netCDF4.Dataset(topo_fname, 'r')
gic_nc = netCDF4.Dataset(gic_fname, 'r')

focean = read_ncvar(topo_nc, 'focean')
flake = read_ncvar(topo_nc, 'flake')
fgrnd = read_ncvar(topo_nc, 'fgrnd')
fgice = read_ncvar(topo_nc, 'fgice')

# Read old topo file
topo0 = dict()
for rec in giss.io.giss.reader(topo0_fname) :
	val = np.zeros(rec.data.shape)	# Promote to double
	val[:] = rec.data[:]
	topo0[rec.var.lower()] = val

focean0 = topo0['focean']
flake0 = topo0['flake']
fgrnd0 = topo0['fgrnd']
fgice0 = topo0['fgice']



# ========================================================
# Check invariant, and count # of grid cells changed from original files

xone0 = focean0 + flake0 + fgrnd0 + fgice0
print 'xone0',test_one(xone0)

print focean0.dtype, focean.dtype


compare_old_new_absolute(focean0, focean, 'focean')
compare_old_new_absolute(flake0, flake, 'flake')
compare_old_new_absolute(fgrnd0, fgrnd, 'fgrnd')
compare_old_new_absolute(fgice0, fgice, 'fgice')

xone = focean + flake + fgrnd + fgice
print 'xone',test_one(xone)

# ========================================================
# Check zatmo consistent with elevhc and fhc

zatmo0 = topo0['zatmo']
zatmo = read_ncvar(topo_nc, 'zatmo')
print 'zatmo0 nans =', np.sum(np.isnan(zatmo0))
print 'zatmo nans =', np.sum(np.isnan(zatmo))
compare_old_new_relative(zatmo0, zatmo, 'zatmo')

elevhc = read_ncvar(gic_nc, 'elevhc')
fhc = read_ncvar(gic_nc, 'fhc')

print 'fhc non-nans =', np.sum(np.logical_not(np.isnan(fhc)))
fhc_sum = np.sum(fhc,0)
print 'fhc_sum',test_one(fhc_sum)

zatmo_new = np.sum(fhc * elevhc,0)
mask = np.isnan(zatmo_new)
zatmo_new[mask] = zatmo[mask]
print 'zatmo_new nans =', np.sum(np.isnan(zatmo_new))
compare_old_new_relative(zatmo, zatmo_new, 'zatmo_new')

# ==================================================
# Check for nan in TLANDI
tlandi = read_ncvar(gic_nc, 'tlandi')
print 'tlandi has %d nans' % np.sum(np.isnan(tlandi))
