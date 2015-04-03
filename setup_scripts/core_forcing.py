from pylab import *
#from scipy.io import netcdf
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from mpl_toolkits.basemap import maskoceans
from astropy.convolution import Gaussian1DKernel, convolve
import os

outdir = os.path.join('../bin_files')

### CORE ###
core_dir = os.path.join(os.environ['D'],'COREv2')
core_suff = '15JUNE2009'
core_outsuff = 'corev2_zon_avg'

#### SCOW ####
scow_dir = os.path.join(os.environ['D'],'scow/monthly_fields')
scow_suff = 'monthly_maps'
scow_outsuff = 'scow_zon_avg'

# binary format for output files
tp = '>f4'

Nt = 12
startdate1 = '00000101'
startdate2 = '000000'
repeat_period = 24*60*60*360 # our definition of a year
forc_period = repeat_period/Nt # our definition of a month

lon0 = 0.
dlon = 45.

# big list of core variables for exf 
# don't do stress, instead get from scow
fields_core =[['t_10', 'T_10', 1., 'atemp',],
         ['q_10', 'Q_10', 1., 'aqh'],
#         ['u_10', 'U_10', 1., 'uwind'],
#         ['v_10', 'V_10', 1., 'vwind'],
         ['ncar_precip', 'PRC_MOD', 1e-3, 'precip'],
         #['ncar_precip', 'SNOW', 1e-3, 'precipsnow'],
         ['ncar_rad', 'SWDN', 1., 'swdown'],
         ['ncar_rad', 'LWDN', 1., 'lwdown']]

fields_scow = [['wind_stress_zonal', 1., 'ustress'],
               ['wind_stress_meridional', 1., 'vstress']]
         
def average_data_core(ncfilename, varname, lonname='LON',latname='LAT', fac=1.):
    fname = os.path.join(core_dir, ncfilename)
    ncf = Dataset(fname)
    ncv = ncf.variables[varname]
    (nt,ny,nx) = ncv.shape
    print (nt, ny, nx)
    # need to interpolate
    lon, lat = ncf.variables[lonname][:], ncf.variables[latname][:]
    lon[lon>=180] -=360.
    lons, lats = np.meshgrid(lon,lat)
    
    dt = nt / Nt
    
    # indices for atlantic and pacific
    regions = {
        'globe': np.r_[0:nx],
        'atl': np.r_[155:192],
        'pac': np.r_[76:140]
    }
    
    data = {}
    for r in regions:
        data[r] = np.zeros((Nt, ny))
    
    # 1.5 degree smoother 
    ker = Gaussian1DKernel(1)
    
    for n in xrange(Nt):
        # zonal avg of month
        d = maskoceans(lons, lats, ncv[(n*dt):((n+1)*dt)].mean(axis=0), inlands=False)
        d.mask = ~d.mask
        for r, iidx in regions.iteritems():  
            dz = fac*d[:,iidx].mean(axis=1)
            dz_sm = convolve(dz.filled(np.nan), ker, boundary='extend')
            data[r][n] = dz_sm                      
    return data['globe'], data['atl'], data['pac'], lat

def average_data_scow(ncfilename, lonname='longitude',latname='latitude', fac=1.,
                        stress_min=-0.03):
    ncf = Dataset(os.path.join(scow_dir, ncfilename))
    months = ['january', 'february', 'march', 'april', 'may', 'june',
              'july', 'august', 'september', 'october', 'november', 'december']

    nt = 12

    lon, lat = ncf.variables[lonname][:], ncf.variables[latname][:]
    #lon[lon>=180] -=360.
    #lons, lats = np.meshgrid(lon,lat)
    dlat = np.diff(lat)[0]
    lat0 = lat[0]
    
    # pad lat to extend farther south (need to get to 71S)
    padlen = 5
    lat = np.pad(lat, (padlen,0), mode='constant')
    lat[:padlen] = -dlat*np.arange(1,6)[::-1] + lat0
    
    ny, nx = len(lat), len(lon)
    print (nt, ny, nx)
    
    # indices for atlantic and pacific
    regions = {
        'globe': np.r_[0:nx],
        'atl': np.r_[1160:1440],
        'pac': np.r_[554:1054]
    }
    
    data = {}
    for r in regions:
        data[r] = np.zeros((Nt, ny)) 
    
    # 1.5 degree smoother 
    ker = Gaussian1DKernel(6)

    for n in range(len(months)):
        ncv = ncf.variables[months[n]]
        d = np.pad(ncv[:], ((padlen,0),(0,0)), mode='constant',
                        constant_values=-999)
        # ocean mask is already built in
        #d = maskoceans(lons, lats, ncv[:], inlands=False)
        #d.mask = ~d.mask
        d = np.ma.masked_less(d, -99)
        
        for reg, iidx in regions.iteritems():
        
            dz = fac*d[:,iidx].mean(axis=1)

            # hack to produce easterlies
            dz_mask = (dz==0.) # will pick up the sea ice margin
            jmax = np.where(dz_mask[:100])[0].max()+1
            jidx = np.arange(0,jmax)
            dz[jidx] = np.interp(jidx, [0, jmax], [stress_min, dz[jmax]])

            # smooth data with 1 degree gaussian kernel
            dz_sm = convolve(dz, ker, boundary='extend')

            data[reg][n] = dz_sm
        
    return data['globe'], data['atl'], data['pac'], lat


glob, atl, pac, outsuff = {}, {}, {}, {}
lat = {}

for ncname, varname, fac, mname in fields_core:
    glob[mname], atl[mname], pac[mname], lat[mname] = average_data_core(
                '%s.%s.nc' % (ncname, core_suff), varname, fac=fac)
    outsuff[mname] = core_outsuff

Nlon = 2
# write data.exf

header = """ &EXF_NML_01
  repeatPeriod = 31104000.,
  ocean_emissivity = 0.97,
  ice_emissivity = 0.95,
  snow_emissivity = 0.95,
  exf_iprec = 32,
  hu = 10.,
  ht = 10.,
 &
 &EXF_NML_02
"""

myflds = sort(glob.keys())

with open('../input/data.exf', 'w') as f:
    f.write(header)
    
    for k in myflds:
        f.write("  %sfile = '%s.%s'\n" % (k, k, outsuff[k]))
        f.write("  %sstartdate1 = %s\n" % (k, startdate1))
        f.write("  %sstartdate2 = %s\n" % (k, startdate2))
        f.write("  %speriod = %g\n" % (k, forc_period))
        f.write("\n")
    f.write(' &\n')
    f.write(' &EXF_NML_03\n')
    f.write(' &\n')
    f.write(' &EXF_NML_04\n')
        
    for k in myflds:
        Nlat = len(lat[k])      
        dlat = np.diff(lat[k])
        if np.any(dlat != dlat[0]):
            lat_str = np.array2string(
                            dlat,
                            max_line_width=50,
                            separator=', ')[2:-1]
        else:
            lat_str = '%g*%g' % (len(dlat), dlat[0])    
        f.write("  %s_nlon = %g\n" % (k, Nlon))
        f.write("  %s_nlat = %g\n" % (k, Nlat))
        f.write("  %s_lon0 = %g\n" % (k, lon0))
        f.write("  %s_lon_inc = %g\n" % (k, dlon))
        f.write("  %s_lat0 = %g\n" % (k, lat[k][0]))
        f.write("  %s_lat_inc = %s\n" % (k, lat_str))
        f.write("\n")
    f.write(' &')

# create a rescaled precipation 
for scalefac in [0.,0.25,0.5,0.75, 0.9]:
    lat0 = 40.
    dl = 8.
    fac = scalefac/2 + 0.5
    scalefun = -np.tanh((lat['precip'] - lat0)/dl) * (1-fac) + fac
    pscaled = glob['precip'] * scalefun
    newkey = 'precip_scaled_%03d' % (100*scalefac)
    glob[newkey] = pscaled 
    outsuff[newkey] = outsuff['precip']
    
# create a rescaled N. Atl. temp
for tfac in [2., 5., 10., 20.]:
    dl = 8.
    lat0 = 60.
    scalefun = np.exp((lat['atemp'] - lat0)/dl)
    atemp_scaled = glob['atemp'] - tfac*scalefun
    newkey = 'atemp_cold_north_%02ddeg' % tfac
    glob[newkey] = atemp_scaled
    outsuff[newkey] = outsuff['atemp']


for ncname, fac, mname in fields_scow:
    glob[mname], atl[mname], pac[mname], lat[mname] = average_data_scow(
                '%s_%s.nc' % (ncname, scow_suff), fac=fac)
    outsuff[mname] = scow_outsuff

# write output files
# duplicate once in longitude to make sure interpolation works
for k, data in glob.iteritems():
    d2d = np.tile(data[:,:,np.newaxis],[1,1,Nlon])
    d2d.astype(tp).tofile('%s/%s.%s' % (outdir, k, outsuff[k]))
