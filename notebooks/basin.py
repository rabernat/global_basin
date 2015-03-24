import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from mitgcmdata import MITgcmmodel
from MITgcmutils import jmd95
import os


class BasinModel(MITgcmmodel.ModelInstance):
    """Subclass of MITgcmdata.ModelSetup specified for basin model."""

    def __init__(self, run_name, base_dir='../', fig_dir='../figures/', *args, **kwargs):
        
        self.run_name = run_name
        self.fig_dir = fig_dir
        output_dir = os.path.join(base_dir, 'run.%s' % run_name)
        
        MITgcmmodel.ModelInstance.__init__(self, output_dir,
                *args, **kwargs)
        
        # some variables to make plotting and pcolors easier  
        self.X = self.xc[0,0,:]
        self.Y = self.yc[0,:,0]
        self.Z = self.rc[:,0,0]
        self.YG = np.hstack([self.yg[0,:,0], 2*self.yg[0,-1,0]-self.yg[0,-2,0]])
        self.XG = np.hstack([self.xg[0,0], 2*self.xg[0,0,-1]-self.xg[0,0,-2]])
        self.RF = self.rf[:,0,0]
        
        self.set_plot_defaults()
    
    ### Analysis Function ###
    
    def load_tave_fields(self):
        """Load standard time-average output data from default_iter."""   
        
        # 3D fields
        self.T = self.mask_field(self.rdmds('Ttave'))
        self.S = self.mask_field(self.rdmds('Stave'))
        self.PH = self.mask_field(self.rdmds('PhHytave'))
        self.U = self.rdmds('uVeltave')
        self.V = self.rdmds('vVeltave')
        self.W = self.rdmds('wVeltave')    
        
        # surface fluxes
        self.empmr = self.rdmds('EmPmRtave')* self.hFacC[0]
        self.qnet = self.rdmds('QNETtave')* self.hFacC[0]
        self.fu = self.rdmds('FUtave') * self.hFacC[0]
        self.fv = self.rdmds('FVtave') * self.hFacC[0]
        
        # Gent-McWilliams streamfunction
        self.GM_PsiY = None
        self.GM_PsiX = None
        try:
            self.GM_PsiY = self.rdmds('GM_PsiYtave')
            self.GM_PsiX = self.rdmds('GM_PsiXtave')
        except IOError:
            pass
        
        self.Kr = None
        try:
            self.Kr = self.rdmds('GGL90diffKr-T')
        except IOError:
            try:
                self.Kr = self.rdmds('KPPdiffKzT-T')
            except IOError:
                pass
        
        # ice
        self.heff = self.rdmds('HEFFtave')
            
        self.calc_moc()
        self.calc_barotropic_streamfunction()
        self.calc_sigma()
        
    def calc_sigma(self):        
        self.sig0 = self.mask_field(jmd95.densjmd95(self.S, self.T, 0.))-1000.
        self.sig2 = self.mask_field(jmd95.densjmd95(self.S, self.T, 2000.))-1000.
    
    def calc_barotropic_streamfunction(self):
        Ubt = self.integrate(self.U, axes=(-3,))
        self.PsiZ_bt = -np.cumsum(Ubt * self.dyg, axis=-2)
        
    def calc_moc(self, norm=1e6):
        """Calculate Eulerian and GM-induced overturning streamfunctions."""
        Vbar = self.integrate(self.V, axes=(-1,))
        self.PsiY_bar = np.cumsum(self.drf*Vbar, axis=-3)/norm

        if self.GM_PsiY is not None:
            self.PsiY_eddy = self.integrate(self.GM_PsiY, axes=(-1,))/norm
            self.PsiY_res = self.PsiY_bar + self.PsiY_eddy
            
    # utility functions
    def find_i_lon(self, lon):
        return find_closest(self.X, lon)
    def find_j_lat(self, lat):
        return find_closest(self.Y, lat)
    def find_k_depth(self, depth):
        return find_closest(self.Z, depth)
            
    ### Plotting Function ###
    def set_plot_defaults(self,
        Slevs = np.arange(33,37,0.2),
        Tlevs = np.arange(-2,25),
        sig2levs = np.hstack([np.arange(32,36), np.arange(36,38,0.2)]),
        psilevs = np.arange(-30,31)+0.5,
        Krlevs = 10.**(np.arange(-5,-0.9,0.2)),
        icefac=200.):
        
        self.Slevs = Slevs
        self.Tlevs = Tlevs
        self.sig2levs = sig2levs
        self.psilevs = psilevs
        self.icefac = 200.
    
    def generate_standard_plots(self, save=False, itersuff=False, dpi=200., fmt='png'):
        
        if save:
            # for saving
            plt.rcParams['font.size'] = 8
            hfigsize = (6.5,3.)
            vfigsize = (3.,6.5)
        else:
            # just for on screen
            plt.rcParams['font.size'] = 12
            hfigsize = (12,5.)
            vfigsize = (5.,12.)

        
        myfigs = {}

        myfigs['zm_theta'] = plt.figure(figsize=hfigsize)
        ax = plt.subplot(111)
        self.plot_zm_theta(ax)
        plt.tight_layout()

        myfigs['zm_salt'] = plt.figure(figsize=hfigsize)
        ax = plt.subplot(111)
        self.plot_zm_salt(ax)
        plt.tight_layout()
        
        myfigs['zm_Kr'] = plt.figure(figsize=hfigsize)
        ax = plt.subplot(111)
        self.plot_zm_Kr(ax)
        plt.tight_layout()
        
        myfigs['surface_flux'] = plt.figure(figsize=hfigsize)
        ax1 = plt.subplot(121)
        self.plot_zm_qnet(ax1)
        ax2 = plt.subplot(122)
        self.plot_zm_empmr(ax2)

        myfigs['v_sections'] = plt.figure(figsize=hfigsize)
        j1 = self.find_j_lat(-30)
        j2 = self.find_j_lat(30)
        ax1 = plt.subplot(121)
        self.plot_v_section(j1, ax1)
        ax2 = plt.subplot(122)
        self.plot_v_section(j2, ax2)
        plt.tight_layout()
        
        myfigs['barotropic_sf'] = plt.figure(figsize=vfigsize)
        ax = plt.subplot(111)
        self.plot_barotropic_streamfunction(ax)
        plt.tight_layout()
        
        myfigs['w100'] = plt.figure(figsize=vfigsize)
        k = self.find_k_depth(-100)
        ax = plt.subplot(111)
        self.plot_w_map(k, ax)
        plt.tight_layout()
        
        myfigs['w500'] = plt.figure(figsize=vfigsize)
        k = self.find_k_depth(-500)
        ax = plt.subplot(111)
        self.plot_w_map(k, ax)
        plt.tight_layout()
        
                        
        myfigs['w2000'] = plt.figure(figsize=vfigsize)
        k = self.find_k_depth(-2000)
        ax = plt.subplot(111)
        self.plot_w_map(k, ax)
        plt.tight_layout()
        
        if save:
            out_dir = os.path.join(self.fig_dir, self.run_name)
            if itersuff:
                isuff = '.%010d' % self.default_iter
            else:
                isuff = ''
            try:
                os.makedirs(out_dir)
            except OSError:
                pass
            for name, fig in myfigs.iteritems():
                fname = os.path.join(out_dir, '%s%s.%s' % (name, isuff, fmt))
                fig.savefig(fname, dpi=dpi)
    
    def plot_zm_qnet(self, ax, lims=[-100,100]):
        p = ax.plot(self.Y, self.qnet.mean(axis=-1), 'k-')
        ax.set_ylim(lims)
        ax.set_xlim([self.YG[0], self.YG[-1]])
        ax.grid()
        ax.set_xlabel('latitude')
        ax.set_ylabel(r'flux (W m$^{-2}$)')
        ax.set_title('Zonal Mean Surface Heat Flux')
        
    def plot_zm_empmr(self, ax, scale=(24*60*60*360/1000.), lims=[-2,2]):
        p = ax.plot(self.Y, scale*self.empmr.mean(axis=-1), 'k-')
        ax.set_ylim(lims)
        ax.set_xlim([self.YG[0], self.YG[-1]])
        ax.grid()
        ax.set_xlabel('latitude')
        ax.set_ylabel('flux (m / year)')
        ax.set_title('Zonal Mean Surface Freshwater Flux')     
        
    def plot_zm_theta(self, ax):
        c = self.contour_meridional_section(self.T.mean(axis=-1), ax,
                            self.Tlevs, cmap='RdBu_r', extend='both')
        plt.colorbar(c, ax=ax)
        self.decorate_meridional_section(ax)
        ax.set_title('Zonal Mean Potential Temperature (deg. C)')    
    
    def plot_zm_salt(self, ax):
        c = self.contour_meridional_section(self.S.mean(axis=-1), ax,
                        self.Slevs, cmap='PuOr_r', extend='both')
        plt.colorbar(c, ax=ax)
        self.decorate_meridional_section(ax)
        ax.set_title('Zonal Mean Salinity (PSU)')   
        
    def plot_zm_Kr(self, ax):
        c = self.contour_meridional_section(self.Kr.mean(axis=-1), ax, 
                                 10.**(np.arange(-5,-0.9,0.25)),
                                 cmap='YlOrRd', #extend='both',
                                 locator=ticker.LogLocator())
        plt.colorbar(c, ax=ax)
        self.decorate_meridional_section(ax)
        ax.set_title(r'$log_{10}(K_r)$ (m$^2$ s$^{-2}$)')
        
    def plot_v_section(self, j, ax, vscale=5, dv=0.5):
        vlevs=np.arange(-vscale,vscale,dv)+dv/2
        vticks=np.arange(-vscale,vscale+dv,2*dv)
        cf = self.contour_zonal_section(100*self.V[:,j], ax, vlevs,
                                    cmap='RdBu_r', extend='both')
        plt.colorbar(cf, ax=ax, orientation='horizontal', ticks=vticks)
        self.decorate_zonal_section(ax, j)
        ax.set_title('Meridional Velocity (cm/s) | %2.0f Latitude' % self.YG[j])
                
    def plot_barotropic_streamfunction(self, ax, norm=1e6, psilims=[-80,10], psilevs=5):
        levs = np.arange(psilims[0],psilims[-1]+psilevs,psilevs)
        c = self.latlon_map_contourf(self.PsiZ_bt/1e6, ax, levs,
                                     cmap='RdYlBu', extend='both')
        plt.colorbar(c, ax=ax, orientation='horizontal', pad=0.1,
                     ticks=ticker.MaxNLocator(nbins=7), shrink=0.7)
        self.decorate_latlon_map(ax)
        ax.set_title('Barotropic Streamfunction (Sv)')
        
    def plot_w_map(self, k, ax, norm=24*60*60*360., lim=200.):
        pc = self.latlon_map_pcolormesh(norm * self.W[k], ax,
                              cmap='RdBu_r')
        plt.colorbar(pc, ax=ax, orientation='horizontal', pad=0.1,
                     ticks=ticker.MaxNLocator(nbins=7), shrink=0.7)
        self.decorate_latlon_map(ax)
        pc.set_clim([-lim, lim])
        ax.set_title('w (m/year) | %4.0f m' % -self.RF[k+1])
        
    def contour_meridional_section(self, f, ax, *args, **kwargs):
        cf = ax.contourf(self.Y, self.Z, f, *args, **kwargs)
        return cf

    def decorate_meridional_section(self, ax, bathy=True, sig2=True, moc=True, ice=True):
        """Add bathymetrey, zonal mean density, and MOC streamfunction to section plot."""
       
        if sig2:
            c = ax.contour(self.Y, self.Z, self.sig2.mean(axis=-1), self.sig2levs, colors='0.5')
            plt.clabel(c, fmt='%3.1f')
        
        if moc:
            ax.contour( self.Y, self.Z, self.PsiY_res.squeeze(), self.psilevs,
                     extend='both', colors='k' )
        if bathy:
            ax.set_frame_on(False)
            # drake passage and surface
            ax.plot(self.Y, -self.depth[:,-1], 'k-', linewidth=2)
            # ridge and shelf
            #ax.plot(self.Y, -self.depth[:,self.Nx/2], 'k-', linewidth=2)
            ax.fill_between(self.Y, -self.depth[:, self.Nx/4], -self.depth.max(),
                            facecolor='0.5')
            ax.fill_between(self.Y, -self.depth[:, self.Nx/2], -self.depth.max(),
                            alpha=0.2, facecolor='0.5')
        
        # ice
        if ice:
            ax.fill_between(self.Y, self.icefac*self.heff.mean(axis=-1), facecolor='1.0')
        
        ax.set_xlabel('latitude')
        ax.set_ylabel('depth (m)')

    def decorate_zonal_section(self, ax, j=None, sig2=True, bathy=True):
        """Add bathymetrey, zonal mean density, and MOC streamfunction to section plot."""
       
        if (j is not None) and sig2:
            c = ax.contour(self.X, self.Z, self.sig2[:,j], self.sig2levs, colors='0.5')
            plt.clabel(c, fmt='%3.1f')
            
        if (j is not None) and bathy:
            #ax.set_frame_on(False)
            ax.fill_between(self.X, -self.depth[j], -self.depth.max(),
                            facecolor='0.5')
        
        ax.set_xlabel('longitude')
        ax.set_ylabel('depth (m)')        
        
    def contour_zonal_section(self, f, ax, *args, **kwargs):
        cf = ax.contourf(self.X, self.Z, f, *args, **kwargs)
        return cf
    
    def decorate_latlon_map(self, ax):
        ax.set_aspect(1)
        ax.set_ylim((self.YG[0], self.YG[-1]))
        ax.set_xlim((self.XG[0], self.XG[-1]))
        # land
        pc = self.latlon_map_pcolormesh(
                np.ma.masked_greater(self.hFacC[0],0.), ax, cmap='Greys')
        pc.set_clim([-1,0])
        ax.contour(self.X, self.Y, self.depth, [140,4000], colors='0.5')
        
        ax.set_ylabel('latitude')
        ax.set_xlabel('longitude')
        ax.set_yticks(range(-70,61,10))
        ax.set_xticks(range(0,45,10))
        
    def latlon_map_pcolormesh(self, f, ax, *args, **kwargs):
        pc = ax.pcolormesh(self.XG, self.YG, f.squeeze(), *args, **kwargs)
        return pc
    
    def latlon_map_contourf(self, f, ax, *args, **kwargs):
        cf = ax.contourf(self.X, self.Y, f.squeeze(), *args, **kwargs)
        return cf
 
        
    
def find_closest(arr, val):
    """Utitlity for finding closest integer value."""
    return np.argmin((arr-val)**2)
        

