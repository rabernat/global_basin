from pylab import *
#from scipy.io import netcdf
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import os

ddir = os.path.join(os.environ['D'],'COREv2')
outdir = os.path.join('../bin_files')

suff = '15JUNE2009'
outsuff = 'corev2_zon_avg'
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
fields =[['t_10', 'T_10', 1., 'atemp',],
         ['q_10', 'Q_10', 1., 'aqh'],
         ['u_10', 'U_10', 1., 'uwind'],
         ['v_10', 'V_10', 1., 'vwind'],
         ['ncar_precip', 'PRC_MOD', 1e-3, 'precip'],
         #['ncar_precip', 'SNOW', 1e-3, 'precipsnow'],
         ['ncar_rad', 'SWDN', 1., 'swdown'],
         ['ncar_rad', 'LWDN', 1., 'lwdown']]
         
def average_data(ncfilename, varname, lonname='LON',latname='LAT', fac=1.):
    ncf = Dataset(os.path.join(ddir,ncfilename))
    ncv = ncf.variables[varname]
    (nt,ny,nx) = ncv.shape
    print (nt, ny, nx)
    # need to interpolate
    lon, lat = ncf.variables[lonname][:], ncf.variables[latname][:]
    dt = nt / Nt
    data_out = zeros((Nt, ny))
    for n in xrange(Nt):
        # zonal avg of month
        data_out[n] = ncv[(n*dt):((n+1)*dt)].mean(axis=0).mean(axis=1)
    return fac*data_out, lat

res = {}
lat = {}
for ncname, varname, fac, mitgcmname in fields:
    res[mitgcmname], lat[mitgcmname] = average_data(
                '%s.%s.nc' % (ncname, suff), varname, fac=fac)

# write output files
# duplicate once in longitude to make sure interpolation works
Nlon = 2
for k, data in res.iteritems():
    d2d = np.tile(data[:,:,np.newaxis],[1,1,Nlon])
    d2d.astype(tp).tofile('%s/%s.%s' % (outdir, k, outsuff))


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

myflds = sort(res.keys())

with open('../input/data.exf', 'w') as f:
    f.write(header)
    
    for k in myflds:
        Nlat = len(lat[k])
        lat_str = np.array2string(
                        np.diff(lat[k]),
                        max_line_width=50,
                        separator=', ')[2:-1]
        f.write("  %sfile = '%s.%s'\n" % (k, k, outsuff))
        f.write("  %sstartdate1 = %s\n" % (k, startdate1))
        f.write("  %sstartdate2 = %s\n" % (k, startdate2))
        f.write("  %speriod = %g\n" % (k, forc_period))
        f.write("\n")
    f.write(' &\n')
    f.write(' &EXF_NML_03\n')
    f.write(' &\n')
    f.write(' &EXF_NML_04\n')
        
    for k in myflds:    
        f.write("  %s_nlon = %g\n" % (k, Nlon))
        f.write("  %s_nlat = %g\n" % (k, Nlat))
        f.write("  %s_lon0 = %g\n" % (k, lon0))
        f.write("  %s_lon_inc = %g\n" % (k, dlon))
        f.write("  %s_lat0 = %g\n" % (k, lat[k][0]))
        f.write("  %s_lat_inc = %s\n" % (k, lat_str))
        f.write("\n")
    f.write(' &')

