import  os
import  sys
import  glob
import  fitsio
import  matplotlib
import  pylab                                  as       pl
import  pandas                                 as       pd
import  numpy                                  as       np
import  astropy.io.fits                        as       fits
import  matplotlib.pyplot                      as       plt
import  numpy.lib.recfunctions                 as       rfn
import  healpy                                 as       hp

from    mpl_toolkits.axes_grid1                import   make_axes_locatable
from    fast_scatter                           import   fast_scatter
from    matplotlib                             import   rc
from    astropy.table                          import   Table, vstack
from    desitarget.targets                     import   encode_targetid
from    desitarget.geomask                     import   is_in_box
from    desitarget.targetmask                  import   desi_mask
from    mpl_toolkits.axes_grid1                import   make_axes_locatable
from    mpl_toolkits.axes_grid1.inset_locator  import   inset_axes
from    scipy.stats                            import   pearsonr
from    phot_sys                               import   set_photsys


rc('font', **{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

##
nside           = np.int(sys.argv[1])
parea           = hp.nside2pixarea(nside, degrees = True)

def read_elgs(nrandom=20000):
  cols          = ['NOBS_G', 'NOBS_R', 'NOBS_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2']
  
  randoms       = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-2.fits')
  randoms       = randoms[1][cols + ['RA', 'DEC', 'PHOTSYS']][:nrandom]
    
  return  randoms, cols
  
    
if __name__ == '__main__':    
  nside         = 512

  randoms, cols = read_elgs(nrandom=90000)

  '''
  set_photsys(randoms)

  sys           = np.unique(randoms['PHOTSYS'])
  split         = []

  for _ in sys:
    split.append(randoms[randoms['PHOTSYS'] == _])
  
  print(sys)
  '''
  
  ##  randoms       = Table(data=randoms, names=['RA', 'DEC'] + cols)  
  ##  randoms.pprint()

  binary       = np.load('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity_{}.npy'.format(nside))
  rhpind       = binary[:,0]
  rhpra        = binary[:,1]
  rhpdec       = binary[:,2]
  rtdensity    = binary[:,3]

  ##
  fig, axarr = plt.subplots(nrows=7, ncols=2, figsize=(12, 25))

  plt.subplots_adjust(left = 0.05, right = 0.95, hspace=0.4, wspace=0.2, top = 0.955, bottom = 0.025)

  ##
  hppix        = hp.ang2pix(nside, (90. - randoms['DEC']) * np.pi / 180., randoms['RA'] * np.pi / 180., nest=False)
  hpind, cnts  = np.unique(hppix, return_counts=True)

  theta, phi   = hp.pix2ang(nside, hpind, nest=False)
  hpra, hpdec  = 180. / np.pi * phi, 90. -180. / np.pi * theta

  hpra[hpra > 300.] -= 360.
  hpra  += 60.

  common, x_ind, y_ind = np.intersect1d(rhpind, hpind, return_indices=True) 
  
  ##
  for count, toplot in enumerate(cols):
    row          = count % 7
    col          = np.int(count >= 7) 

    print(row, col, toplot)
    
    values       = np.array([np.mean(randoms[toplot][hppix == x]) for x in hpind])

    rr           = pearsonr(rtdensity[x_ind], values[y_ind])
    
    vmin         = np.quantile(values, 0.05)
    vmax         = np.quantile(values, 0.95)
    step         = (vmax - vmin) / 50.
      
    fast_scatter(axarr[row][col], hpra, hpdec, values, vmin, vmax, 50, cmap='jet', printit=False)

    toplot       = toplot.split('_')

    if toplot[-1][0] != 'W':
      toplot[-1]      = toplot[-1].lower()

    ##
    axarr[row][col].set_title(r'${}$-{}'.format(toplot[-1], toplot[0].upper()))
    axarr[row][col].set_xlim(360., 0.)
    
    label        = r'$r \simeq$ {:.2e}'.format(rr[0]) 
    
    axarr[row][col].text(75., 70., label)
    
  ##   
  pl.savefig('../../genfiles/elgs/skydepths/depths.png')
