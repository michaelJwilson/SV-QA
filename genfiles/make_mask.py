import os
import glob
import numpy           as     np
import astropy.io.fits as     fits

from   astropy.table   import Table


def get_gaia():
  files  = glob.glob('/project/projectdirs/cosmo/work/gaia/chunks-gaia-dr2-astrom/*.fits')
  ##  print(files)

  result = np.zeros([3]).reshape(1, 3)

  ##  https://github.com/legacysurvey/legacypipe/blob/master/py/legacypipe/reference.py#L48
  for i, _file in enumerate(files):
    hdr    = fits.open(_file)[1].header
    dat    = Table(fits.open(_file)[1].data)['ra', 'dec', 'phot_g_mean_mag']
    dat    = dat[dat['phot_g_mean_mag'] < 16.]

    arr    = np.c_[np.array(dat['ra']), np.array(dat['dec']), np.array(dat['phot_g_mean_mag'])]
    
    result = np.vstack([result, arr])

    print('Retrieved GAIA {} of {}.'.format(i, len(files)))
    
  return  result[1:]


def get_tycho():
  _file    = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/tycho2.kd.fits'
  dat      = Table(fits.open(_file)[1].data)
  dat      = dat[dat['MAG_VT'] < 13.]

  dat['MRADIUS'] = np.minimum(1800., 150. * 2.5**((11. - dat['MAG_VT']) / 3.)) * 0.262
  
  dat      = dat['RA', 'DEC', 'MAG_VT', 'MRADIUS']

  return  dat


if __name__ == '__main__':
  '''
  ##   Ignores proper motions and Tycho overlap. 
  cat = Table(get_gaia(), names=['RA', 'DEC', 'G'])

  cat['ISBRIGHT'] = cat['G'] < 13.
  cat['ISMEDIUM'] = cat['G'] < 16.

  ##   https://github.com/legacysurvey/legacypipe/blob/master/py/legacypipe/reference.py#L196                                                                             
  ##   Removed 1./3600. for deg.
  cat['MRADIUS']  = np.minimum(1800., 150. * 2.5**((11. - cat['G'])/3.)) * 0.262
  
  print(cat)
  
  scratch = os.environ['CSCRATCH']
  cat.write(scratch + '/BGS/SV-ASSIGN/masks/gaia_stellar_mask.fits', format='fits', overwrite=True)
  '''

  tycho = get_tycho()
  tycho.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/tycho_stellar_mask.fits')
  
  print('\n\nDone.\n\n')
