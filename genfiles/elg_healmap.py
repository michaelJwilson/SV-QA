import  os
import  glob
import  fitsio
import  pylab                  as      pl
import  pandas                 as      pd
import  numpy                  as      np
import  astropy.io.fits        as      fits
import  matplotlib.pyplot      as      plt
import  numpy.lib.recfunctions as      rfn
import  healpy                 as      hp

from    astropy.table         import   Table, vstack
from    desitarget.targets    import   encode_targetid
from    desitarget.geomask    import   is_in_box
from    desitarget.targetmask import   desi_mask


nside        = 512
parea        = hp.nside2pixarea(nside, degrees = True)

cols         = ['RA', 'DEC', 'MORPHTYPE', 'DESI_TARGET']

def write_elgs(elgs, write_elgs=False):
  print('Writing ELGs.')

  hppix        = hp.ang2pix(nside, (90. - elgs['DEC']) * np.pi / 180., elgs['RA'] * np.pi / 180., nest=False)
  hpind, cnts  = np.unique(hppix, return_counts=True)

  theta,phi    = hp.pix2ang(nside, hpind, nest=False)
  hpra, hpdec  = 180. / np.pi * phi, 90. -180. / np.pi * theta

  tdensity     = cnts / parea  ##  per sq. deg.                                                                                                                                                                   

  np.save('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity_{}.npy'.format(nside), np.c_[hpind, hpra, hpdec, tdensity])
  
  if write_elgs:
    elgs.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/elgs/elgs.fits', format='fits', overwrite=True)
  
targets      = glob.glob('/project/projectdirs/desi/target/catalogs/dr8/0.32.0/targets/main/resolve/dark/*.fits')

elgs         = Table([fits.open(targets[0])[1].data[x] for x in cols], names=cols)
elgs         = elgs[(elgs['DESI_TARGET'] & desi_mask.mask('ELG')) != 0]

for i, _targets in enumerate(targets[1:]):
  print('Joining {} of {} ELG file.  Currently {} ELGs.'.format(i, len(targets) - 1, len(elgs)))

  _in        = Table([fits.open(_targets)[1].data[x] for x in cols], names=cols)
  _in        = _in[(_in['DESI_TARGET'] & desi_mask.mask('ELG')) != 0]

  _in.pprint(max_width = -1)

  elgs       = vstack([elgs, _in])
  
  ##  Now checking indiv. obj. in the tbc pixels
  
  if i % 5 == 0:
    write_elgs(elgs)
    
##
write_elgs(elgs, write_elgs=True)
