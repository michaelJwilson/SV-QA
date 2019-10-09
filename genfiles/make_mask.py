import os
import glob
import ephem
import fitsio
import astropy
import astropy.units        as     u
import pylab                as     pl 
import healpy               as     hp
import numpy                as     np
import astropy.io.fits      as     fits
import matplotlib.pyplot    as     plt

from   astropy.table        import Table
from   fast_scatter         import fast_scatter
from   astropy.coordinates  import SkyCoord


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
  
  dat            = dat['RA', 'DEC', 'MAG_VT', 'MRADIUS']

  return  dat

def get_HI(nside=256, plotit=False, printit=False):
  ''' 
  Reads Lenz et. al. (1610.06175) HI column density. 
  
  Native resolution:  16.2 arcmin -- nside=256 = 13.74 arcmin. 

  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/AFJNWJ#
  https://nbviewer.jupyter.org/github/DanielLenz/ebv_tools/blob/master/examples.ipynb
  '''

  _ebv               = hp.read_map(os.environ['SCRATCH'] + '/BGS/SV-ASSIGN/ebv_lhd.hpx.fits', verbose=False)

  _nside             = hp.get_nside(ebv)
  _npix              = hp.nside2npix(_nside)
  _ordering          = 'ring'

  if plotit:  
    hp.mollview(np.log10(_ebv), title='Lenz++', unit='log(E(B-V) [mag])')
    pl.show()

  ##  Cast to new resolution.
  ebv                = hp.pixelfunc.ud_grade(_ebv, nside, order_in=_ordering) 

  ##  np.clip(ebv, a_min=ebv[ebv > 0.0].min(), a_max = None)

  assert  np.all(ebv > 0.0)

  if printit:
    print(nside, np.count_nonzero(ebv < 0.0), len(ebv))
    print(np.sort(ebv))

  ## 
  theta, phi  = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=False)
  hpra, hpdec = 180. / np.pi * phi, 90. - 180. / np.pi * theta

  s           = SkyCoord(ra = hpra * u.degree, dec = hpdec * u.degree, frame='icrs')
  s           = s.galactic

  glon        = s.l.value
  glat        = s.l.b.value
        
  pix         = hp.ang2pix(nside, glon, glat, lonlat=True)

  return  hpra, hpdec, np.log10(ebv_map[pix])
  

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
  '''
  tycho = get_tycho()
  tycho.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/tycho_stellar_mask.fits')
  '''

  hpra, hpdec, loghi  = get_HI(1024)
  
  print(hpra.min(), hpra.max())
  
  ##  Wrap randoms                                                                                                                                                                                                             
  hpra[hpra > 300.]  -= 360.
  hpra               += 60.
  
  vmin                = loghi.min() ##  np.percentile(loghi, 0.001)
  vmax                = loghi.max() ##  np.percentile(loghi, 0.999)
  
  fig, ax             = plt.subplots(nrows=1, ncols=1)

  ##
  fast_scatter(ax, hpra, hpdec, loghi, vmin, vmax, 50, markersize=1., cmap='viridis', printit=False, alpha=0.25)

  ##
  ls                  = np.arange(0.1, 359.5, 1.)
  bs                  = np.zeros_like(ls)
  
  gc                  = SkyCoord(l=ls*u.degree, b=bs*u.degree, frame='galactic', equinox='J2000')
  gc                  = gc.transform_to('icrs') 

  ind                 = np.argsort(gc.ra.deg)

  ra                  = gc.ra.deg[ind]
  dec                 =	gc.dec.deg[ind]
  
  ax.plot(ra, dec, 'k-')

  ax.set_xlim(360., 0.)
  
  pl.savefig('hone.png')
  
  print('\n\nDone.\n\n')
