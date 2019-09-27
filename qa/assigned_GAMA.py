from   __future__        import  with_statement

import sys
import glob
import random
import contextlib
import matplotlib;  matplotlib.use('pdf')
import numpy              as       np
import pylab              as       pl
import matplotlib.pyplot  as       plt
import astropy.io.fits    as       fits
import astropy.units      as       u

from   astropy.table      import  Table, join, hstack
from   astropy            import  constants as const
from   desitarget.geomask import  circles

try:
    from urllib.parse import urlencode

except ImportError:
    from urllib import urlencode

try:
    from urllib.request import urlopen

except ImportError:
    from urllib2 import urlopen

plt.style.use('dark_background')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{hyperref}')

##  plt.rcParams['pgf.preamble'] = [r'\usepackage{hyperref}']

def make_tiny(url):
  request_url = ('http://tinyurl.com/api-create.php?' + 
  urlencode({'url':url}))

  with contextlib.closing(urlopen(request_url)) as response:
    return response.read().decode('utf-8').split('//')[1]

##
nrow       = 2
ncol       = 5

rasterized = False

##  
fig, axarr = plt.subplots(nrows=nrow, ncols=ncol, figsize=(20, 10))

tiles      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)

##  Legacy matched truth. 
legacys    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/legacy/ls-GAMA-south.fits')[1].data)

##  Assigned truth with legacy columns. 
files      = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/GAMA*.fits')

for i, _file in enumerate(files):
  tile     = _file.split('/')[-1].split('_')[-1].split('.fits')[0]

  row      = i % nrow
  col      = i % 5
  
  dat      = Table(fits.open(_file)[1].data)
  
  if len(dat) > 0:
    axarr[row][col].scatter(dat['RA'].quantity, dat['DEC'].quantity, s=1, rasterized=rasterized, alpha=0.0, marker='.')

    xlim = axarr[row][col].get_xlim()
    ylim = axarr[row][col].get_ylim()
  
    axarr[row][col].scatter(legacys['RA'], legacys['DEC'], s=3, rasterized=rasterized, alpha=0.3, marker='.', c='white') 

    ##
    axarr[row][col].scatter(dat['RA'], dat['DEC'], s=7, rasterized=rasterized, alpha=1., marker='.', c='gold')

    axarr[row][col].set_xlim(xlim)
    axarr[row][col].set_ylim(ylim)
  
    ##  Hyperlink title.
    _t    = tiles[tiles['TILEID'] == np.int(tile)] 
    ra    = _t['RA'].quantity[0]
    dec   = _t['DEC'].quantity[0]

    hlink = r'http://legacysurvey.org/viewer?ra={:.4f}&dec={:.4f}&layer=dr8&zoom=12&desifiber={:.4f},{:.4f}'.format(ra, dec, ra, dec)

    ##  wlink = r'http://www.astro.utah.edu/~u6022465/SV/tiles/SV_BGS/fits_files/tile-{}.fits'.format(tile)
    alink = r'https://portal.nersc.gov/project/desi/users/mjwilson/SV-ASSIGN/tile-{}.fits'.format(tile)
    tlink = r'https://portal.nersc.gov/project/desi/users/mjwilson/SV-TARGETS/tile-targets-{}.fits'.format(tile)

    axarr[row][col].set_title('(GAMA) Tile:  {}\n'.format(tile) + r'\url{%s}' % make_tiny(hlink) + '\n' + r'\url{%s}' % make_tiny(alink) + '\n' + r'\url{%s}' % make_tiny(tlink))

##  
plt.tight_layout()

pl.savefig('assigned_GAMA.png')

print('\n\nDone.\n\n')
