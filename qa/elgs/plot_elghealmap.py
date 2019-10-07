import  os
import  sys
import  glob
import  fitsio
import  matplotlib
import  pylab                                 as       pl
import  pandas                                as       pd
import  numpy                                 as       np
import  astropy.io.fits                       as       fits
import  matplotlib.pyplot                     as       plt
import  numpy.lib.recfunctions                as       rfn
import  healpy                                as       hp

from    mpl_toolkits.axes_grid1               import   make_axes_locatable
from    fast_scatter                          import   fast_scatter
from    matplotlib                            import   rc
from    astropy.table                         import   Table, vstack
from    desitarget.targets                    import   encode_targetid
from    desitarget.geomask                    import   is_in_box
from    desitarget.targetmask                 import   desi_mask
from    mpl_toolkits.axes_grid1               import   make_axes_locatable
from    mpl_toolkits.axes_grid1.inset_locator import   inset_axes


rc('font', **{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

##
nside          = np.int(sys.argv[1])
parea          = hp.nside2pixarea(nside, degrees = True)

def read_elgs(mtype, snr, ext):
  ##  [hpind, hpra, hpdec, tdensity].
  return  np.load('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity_{}_{:d}{}_{:d}.npy'.format(mtype, snr, ext, nside))
  
    
if __name__ == '__main__':    
  ext        = '_r_z'

  snrs       = [3, 4, 6, 9, 12, 20]
  
  mtypes     = ['DEV', 'COMP', 'EXP', 'REX', 'PSF']

  fig, axarr = plt.subplots(nrows=5, ncols=6, figsize=(20, 20))

  plt.subplots_adjust(left = 0.05, right = 0.95, hspace=0.6, wspace=0.4, top = 0.925, bottom = 0.05)
  
  ##
  for i, mtype in enumerate(mtypes):
    for j, snr in enumerate(snrs):
      hmap                         = read_elgs(mtype, snr, ext)
      hpind, hpra, hpdec, tdensity = np.hsplit(hmap, 4)

      hpra, hpdec, tdensity        = hpra.flatten(), hpdec.flatten(), tdensity.flatten()
              
      print(mtype, snr, np.unique(tdensity, return_counts=True))
      
      hpra[hpra > 300.] -= 360.                                                                                                                                                                               
      hpra  += 60.
      
      vmin   = np.quantile(tdensity, 0.001)
      vmax   = np.quantile(tdensity, 0.999)

      print(mtype, snr, vmin, vmax)
      
      fast_scatter(axarr[i][j], hpra, hpdec, tdensity, vmin, vmax, 50, cmap='jet', printit=False)

      axarr[i][j].set_title(r'{} ELG $r,z$ SNR $\geq {})$'.format(mtype, snr))
      axarr[i][j].set_xlim(360., 0.)

  pl.savefig('../../genfiles/elgs/skydepths/elgs{}.png'.format(ext))
