import  os
import  sys
import  json
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
from    astropy.table                          import   Table, vstack, Column
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
  toload        = ['RA', 'DEC'] + cols
  
  randoms       = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-2.fits')
  randoms       = randoms[1][toload][:nrandom]
    
  return  randoms, cols
  
    
if __name__ == '__main__':    
  nside         = 512

  randoms, cols = read_elgs(nrandom=1000)

  names                  = ['RA', 'DEC'] + cols
  randoms                = Table(data=randoms, names=names)

  randoms['PHOTSYS']     = Column(data=np.array(['N'] * len(randoms)), name='PHOTSYS', dtype='S32')
  randoms['PHOTSYS'][randoms['DEC'] < 32.375] = 'S'
    
  randoms            = set_photsys(randoms)

  ##  randoms.sort('PSFDEPTH_W2')
  
  randoms.pprint()
  
  sys           = np.unique(randoms['PHOTSYS'][:])
  split         = []

  for _ in sys:
    split.append(randoms[randoms['PHOTSYS'] == _])
  
  ##  print(sys)

  ##
  binary       = np.load('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity_{}.npy'.format(nside))
  rhpind       = binary[:,0]
  rhpra        = binary[:,1]
  rhpdec       = binary[:,2]
  rtdensity    = binary[:,3]

  for ss, randoms in zip(sys, split):
    ##
    hppix        = hp.ang2pix(nside, (90. - randoms['DEC']) * np.pi / 180., randoms['RA'] * np.pi / 180., nest=False)
    hpind, cnts  = np.unique(hppix, return_counts=True)

    theta, phi   = hp.pix2ang(nside, hpind, nest=False)
    hpra, hpdec  = 180. / np.pi * phi, 90. -180. / np.pi * theta
  
    hpra[hpra > 300.] -= 360.
    hpra  += 60.

    common, x_ind, y_ind = np.intersect1d(rhpind, hpind, return_indices=True) 

    output         = {}
    
    ##
    for count, toplot in enumerate(cols):
      row          = count % 7
      col          = np.int(count >= 7) 

      ##  print(row, col, toplot)

      values       = np.array([np.mean(randoms[toplot].quantity[hppix == x]) for x in hpind])
      rr           = pearsonr(rtdensity[x_ind], values[y_ind])

      output[toplot] = rr[0]
      
      print('{: <10} \t {: <10} \t {:.2f}'.format(ss, toplot, rr[0]))

    with open('dat/{}_pearson.json'.format(ss), 'w') as fp:
      json.dump(output, fp, sort_keys=True, indent=4)
            
  print('\n\nDone.\n\n')
