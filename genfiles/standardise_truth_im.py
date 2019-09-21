import os
import glob
import numpy            as      np
import pylab            as      pl
import astropy.io.fits  as      fits
import astropy.units    as      u

from   astropy.table    import  Table
from   astropy          import  constants as const


scratch   = os.environ['CSCRATCH']
_all      = ['hsc_pdr1_deep.forced.reduced-match.fits', 'hsc_pdr1_wide.forced.reduced-match.fits', 'hsc_pdr1_udeep.forced.reduced-match.fits']

for _file in _all:
  __file    = scratch + '/BGS/SV-ASSIGN/truth/' + _file

  for gcap in ['north', 'south']:
    ##  If north exists, replace with current south. 
    _file   = __file.replace('-north', '-' + gcap)

    ##  North not present, must be match. 
    _file   = _file.replace('-match', '-' + gcap)
    
    dat     = Table(fits.open(_file)[1].data)
    survey  = _file.split('/')[-1].split('.')[0]
  
    ##
    dat['RA']    = dat['ra']
    dat['DEC']   = dat['dec']
    dat['Z']     = dat['frankenz_photoz_best']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'].quantity)
    dat['ZWARN'] = -99  * np.ones_like(dat['Z'].quantity, dtype=np.int)
 
    print(dat)
    
    '''
    a_g          = 1.155636 * (1 - MW_TRANSMISSION_G) - 0.001767
    a_r          = 0.818918 * (1 - MW_TRANSMISSION_G) - 0.001252
    a_i          = 0.584431 * (1 - MW_TRANSMISSION_G) - 0.000893
    a_z          = 0.450745 * (1 - MW_TRANSMISSION_G) - 0.000689
    a_y          = 0.384616 * (1 - MW_TRANSMISSION_G) - 0.000588
    '''
  
    del  dat['ra']
    del  dat['dec']

    fout = _file.split('.fits')
    fout = fout[0] + '-standard.fits'
    
    fout = fout.split('/')
    fout = '/'.join(x for x in fout[:-1]) + '/standard/' + fout[-1]
    
    print('Writing ... {}'.format(fout))

    try:
      dat.write(fout, format='fits', overwrite=True)
    
    except:
      print('\n\nFailed to write: {}'.format(survey))
  
print('\n\nDone.\n\n')

