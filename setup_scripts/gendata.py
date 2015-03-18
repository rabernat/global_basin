import numpy as np
from netCDF4 import Dataset

# the lowest resolution
#power_of_two = 5
#use_lowres = True

power_of_two = 8
use_lowres = False


tilesize = 64

Nx = 2**power_of_two
Ny = 4*Nx

# width at equator in degrees
basin_width = 45.
fact = Nx / basin_width

lon0 = 0.
rSphere = 6370e3
meters_per_degree_eq = 2*np.pi*rSphere/360.

dlon = 1/fact
XG = lon0 + dlon * np.arange(Nx+1)
XC = lon0 + (0.5*dlon) + dlon * np.arange(Nx)

# asymmetry between north and south
Nassym = Ny/16

lat = [0]
dlat = []
for j in xrange(Ny/2 + Nassym):
    dlat.append(dlon * np.cos(np.radians(lat[-1])))    
    lat.append(lat[-1] + dlat[-1])
lat = np.array(lat)
lat_full = np.hstack([-lat[::-1], lat[1:(-2*Nassym)]])
dlat_full = np.diff(lat_full)



latmin = lat_full[0]
YG = np.hstack([latmin, latmin+np.cumsum(dlat_full)])
YC = np.hstack([latmin + dlat_full[0]/2, latmin+np.cumsum(dlat_full[:-1])])

dely = np.diff(YG)
delx = np.diff(XG)

xc, yc = np.meshgrid(XC, YC)
xg, yg = np.meshgrid(XG, YG)

#Ny = len(YG)

# vertical grid
DRF_surf = 10.
DRF_deep = 200.
DRF_top = DRF_surf*np.ones(20)
scale_fac = 1.15
DRF_tc = np.logspace(np.log10(DRF_surf), np.log10(DRF_deep), 30)
DRF_abyss = DRF_deep*np.ones(14)
DRF = np.hstack([DRF_top, DRF_tc, DRF_abyss])
RF = np.hstack([0, np.cumsum(DRF)])

# a half-resoution version
DRF_lowres = DRF[::2] + DRF[1::2]
RF_lowres = np.hstack([0, np.cumsum(DRF_lowres)])

if use_lowres:
    delR = DRF_lowres
    Nr = 32
    RC = 0.5*(RF_lowres[1:] + RF_lowres[:-1])
else:
    delR = DRF
    Nr = 64
    RC = 0.5*(RF[1:] + RF[:-1])

# bathymetry
# start with full depth
bathy = RF[-1]*np.ones_like(xc)

# shelves
shelfwidth_deg = 2.
slopewidth_deg = 2.
shelfdepth = 140. # meters
shelfmask = np.zeros_like(xc, 'bool')
slopemask = np.zeros_like(xc, 'bool')
# southern shelf
shelfmask[yc <= (latmin + shelfwidth_deg)] = 1
slopemask[yc <= (latmin + shelfwidth_deg + slopewidth_deg)] = 1
# western shelf
shelfmask[xc <= (lon0 + shelfwidth_deg)] = 1
slopemask[xc <= (lon0 + shelfwidth_deg + slopewidth_deg)] = 1
# eastern shelf
shelfmask[xc >= (basin_width - shelfwidth_deg)] = 1
slopemask[xc >= (basin_width - shelfwidth_deg - slopewidth_deg)] = 1

# mid ocean ridge
ridgewidth_deg = 3.
ridgedepth = RF[50]
ridgemask = np.zeros_like(xc, 'bool')
rslopemask = np.zeros_like(xc, 'bool')
ridgemask[:, Nx/2] = 1
ridgemask[yc <= (latmin + shelfwidth_deg + slopewidth_deg/4)] = 0
rw = ridgewidth_deg / dlon
rslopemask[:, (Nx/2-rw):(Nx/2+rw+2)] = 1
# special mask for the ridge / shelf transition
#rstransmask = ridgemask & slopemask & ~shelfmask
#Nrs = rstransmask.sum()

# drake passage
dplat = [-61, -56]
dpdepth = RF[50]
dpmask = shelfmask & (yc > dplat[0]) & (yc < dplat[1])

# fill in known values
bathy[ridgemask] = ridgedepth
#bathy[rstransmask] = np.linspace(shelfdepth, ridgedepth, Nrs)
bathy[shelfmask] = shelfdepth
# now add side boundary
bathy[:,-1] = 0
# open drake passage
bathy[dpmask] = dpdepth


# where we need to interpolate
imask = (slopemask | rslopemask) & ~shelfmask & ~ridgemask
#imask[-1] = 0

## try to use diffusion
def diffuse_2D(f):
    fe = np.roll(f, -1, axis=-1)
    fw = np.roll(f, 1, axis=-1)
    fn = np.roll(f, -1, axis=-2)
    fs = np.roll(f, 1, axis=-2)
    # get rid of periodic bc in latitude
    fn[-1] = f[-1]
    fs[0] = f[0]
    return (-4*f + fe + fw + fn + fs)

bnew = bathy.copy()
for n in xrange(10000):
    bprev = bnew.copy()
    bnew += 0.2 * diffuse_2D(bnew)
    bnew[~imask] = bathy[~imask]
    delta = ((bnew - bprev)**2).mean()
    print n, delta
    if delta < 1e-10:
        break

# close top boundary
bnew[-1] = 0
# bathy needs to be negative
bnew *= -1.

# bathy
output_fname = 'bathy_basin_%g.f4' % Nx
tp = '>f4'
bnew.astype(tp).tofile('../bin_files/%s' % output_fname)

# dely
dely.astype(tp).tofile('../bin_files/dely_basin_%g.f4' % Nx)

# stuff for data
print "xgOrigin = %g," % XG[0]
print "delX = %g*%g," % (Nx, dlon)
print "ygOrigin = %g," % YG[0]
print "delYfile = 'dely_basin_%g.f4'," % Nx
print "delR = ", np.array2string(
                        delR,
                        max_line_width=50,
                        separator=', ')[2:-1]


# initial conditions from WOA
WOA_dir = '/Users/rpa/RND/Data/NODC_WOA98'

nct = Dataset(WOA_dir + '/otemp.anal1deg.nc')
ncs = Dataset(WOA_dir + '/salt.anal1deg.nc')
tvar = nct.variables['otemp']
svar = ncs.variables['salt']
lat = nct.variables['lat'][:]
lev = nct.variables['level'][:]
# just use the atlantic
irange = np.r_[300:350]

tbar = tvar[0,:,:,irange].mean(axis=-1)
sbar = svar[0,:,:,irange].mean(axis=-1)
# interpolate first in lat
t1 = -99*np.ones((len(lev),Ny))
s1 = -99*np.ones((len(lev),Ny))

for k in xrange(len(lev)):
    idx = ~tbar[k,:].mask
    if any(idx):
        t1[k,:] = np.interp(YC, lat[idx][::-1], tbar[k,idx][::-1] )
        s1[k,:] = np.interp(YC, lat[idx][::-1], sbar[k,idx][::-1] )
t1 = np.ma.masked_array(t1, t1==-99)
s1 = np.ma.masked_array(s1, s1==-99)
t2 = np.zeros((Nr,Ny))
s2 = np.zeros((Nr,Ny))
for j in xrange(Ny):
    idx = ~t1[:,j].mask
    if any(idx):
        t2[:,j] = np.interp(RC, lev[idx], t1[idx,j])
        s2[:,j] = np.interp(RC, lev[idx], s1[idx,j])
        
# save data
np.tile(t2[:,:,np.newaxis], [1,1,Nx]).astype(tp).tofile(
        '../bin_files/theta_woa_%gx%gx%g.f4' % (Nx,Ny,Nr))
np.tile(s2[:,:,np.newaxis], [1,1,Nx]).astype(tp).tofile(
        '../bin_files/salt_woa_%gx%gx%g.f4' % (Nx,Ny,Nr))       
