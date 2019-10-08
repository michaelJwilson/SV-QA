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


rc('font', **{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

##
nside          = np.int(sys.argv[1])
parea          = hp.nside2pixarea(nside, degrees = True)

def read_elgs(nrandom=20000, _all=False):
  cols         = ['NOBS_G', 'NOBS_R', 'NOBS_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2', 'MASKBITS']

  if _all:
    nrandom    = -1
    randoms    = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits')

  else:
    randoms    = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-2.fits')

  ##
  randoms      = randoms[1][cols + ['RA', 'DEC']][:nrandom]
    
  return  randoms, cols
  
    
if __name__ == '__main__':
  nside         = 512
  _all          = True
  save          = False
  
  ##
  randoms, cols = read_elgs(nrandom=5000, _all=_all)

  isin          = np.ones(len(randoms), dtype=np.float)

  for bit in [1, 5, 6, 7, 11, 12, 13]:
    isin        = isin * (1.0 - np.clip(np.bitwise_and(randoms['MASKBITS'], np.ones_like(randoms['MASKBITS']) * 2 ** bit), a_min=0.0, a_max=1.0))
    ##  isin    = isin * np.array([(x & 2 ** bit) == 0 for x in randoms['MASKBITS']]).astype(np.int)

  isin          = isin.astype(np.bool)

  remaining     = randoms[isin]

  print('Input:  {};  Remaining:  {}  ({})'.format(len(randoms), len(remaining), np.count_nonzero(isin)))
  
  ##  randoms   = Table(data=randoms, names=['RA', 'DEC'] + cols)  
  ##  randoms.pprint()
  '''
  binary       = np.load('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity_{}.npy'.format(nside))
  rhpind       = binary[:,0]
  rhpra        = binary[:,1]
  rhpdec       = binary[:,2]
  rtdensity    = binary[:,3]
  '''

  ##
  npix           = hp.pixelfunc.nside2npix(nside)
  indices        = np.arange(npix)
  result         = np.zeros_like(indices)
  denom          = np.zeros_like(indices)
  
  ##
  hppix          = hp.ang2pix(nside, (90. - remaining['DEC']) * np.pi / 180., remaining['RA'] * np.pi / 180., nest=False)
  hpind, cnts    = np.unique(hppix, return_counts=True)
  
  for i, ind in enumerate(hpind):
    result[ind] += cnts[i]  

  hppix          = hp.ang2pix(nside, (90. - randoms['DEC']) * np.pi / 180., randoms['RA'] * np.pi / 180., nest=False)
  hpind, cnts    = np.unique(hppix, return_counts=True)

  for i, ind in enumerate(hpind):
    denom[ind]  += cnts[i]

  mask           = np.array(denom > 0.0).astype(np.int)

  result         = mask * result / denom

  occupied       = result[result > 0.0]

  print('\n\nDeviation level: {:.3} per cent.'.format(100. * np.std(occupied)))
  
  ##
  theta, phi     = hp.pix2ang(nside, range(npix), nest=False)
  hpra, hpdec    = 180. / np.pi * phi, 90. -180. / np.pi * theta

  if save:
    np.save(os.environ['CSCRATCH'] + '/BGS/SV-ASSIGN/elgs/pix_area.npy', np.c_[np.arange(npix), hpra, hpdec, result])

  print('\n\nDone.\n\n')
    
  if not _all:
    fig, axarr   = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
      
    plt.subplots_adjust(left = 0.05, right = 0.95, hspace=0.4, wspace=0.2, top = 0.955, bottom = 0.025)
      
    hpra[hpra > 300.] -= 360.
    hpra        += 60.

    fast_scatter(axarr, hpra, hpdec, mask, -0.1, 1.1, 50, cmap='BuPu', printit=False)
      
    axarr.set_xlim(365., -5.)
    
    pl.savefig('pix_area_test.png')
