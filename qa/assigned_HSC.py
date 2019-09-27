from    __future__          import   with_statement

import  sys
import  glob
import  random
import  contextlib
import  matplotlib
import  numpy               as       np
import  pylab               as       pl
import  matplotlib.pyplot   as       plt
import  astropy.io.fits     as       fits
import  astropy.units       as       u

from    astropy.table       import  Table, join, hstack
from    astropy             import  constants as const
from    desitarget.geomask  import  circles

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


nrow       = 2
ncol       = 6

rasterize  = False

##  
fig, axarr = plt.subplots(nrows=nrow, ncols=ncol, figsize=(16, 7.5))


tiles      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)

legacyn    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/legacy/ls-hsc_pdr1_wide.forced.reduced-north.fits')[1].data)
legacys    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/legacy/ls-hsc_pdr1_wide.forced.reduced-south.fits')[1].data)

files      = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/hsc_pdr1_wide*.fits')

row        = 0

for i, _file in enumerate(files):
  tile     = _file.split('/')[-1].split('_')[-1].split('.fits')[0]

  if i > 5:
    row    = 1

  ##
  col      = i % ncol
  dat      = Table(fits.open(_file)[1].data)

  print(i, row, col, _file, len(dat))
  
  if len(dat) > 0:
    axarr[row][col].set_title('Tile {}'.format(tile))

    ## 
    axarr[row][col].scatter(dat['RA'], dat['DEC'], s=1, rasterized=rasterize, alpha=0.1, marker='.', c='k')
  
    axarr[row][col].set_xlim(axarr[row][col].get_xlim())
    axarr[row][col].set_ylim(axarr[row][col].get_ylim())
  
    axarr[row][col].scatter(legacyn['RA'], legacyn['DEC'], s=1, rasterized=rasterize, alpha=0.5, marker='x', c='k')
    axarr[row][col].scatter(legacys['RA'], legacys['DEC'], s=1, rasterized=rasterize, alpha=0.5, marker='x', c='k') 

    ##
    axarr[row][col].scatter(dat['RA'], dat['DEC'], s=2, rasterized=True, alpha=1., marker='.', c='gold')

    ##  Hyperlink title.                                                                                                                                                                                                                
    _t    = tiles[tiles['TILEID'] == np.int(tile)]
    ra    = _t['RA'].quantity[0]
    dec   = _t['DEC'].quantity[0]
    
    hlink = r'http://legacysurvey.org/viewer?ra={:.4f}&dec={:.4f}&layer=dr8&zoom=12&desifiber={:.4f},{:.4f}'.format(ra, dec, ra, dec)

    ##  wlink = r'http://www.astro.utah.edu/~u6022465/SV/tiles/SV_BGS/fits_files/tile-{}.fits'.format(tile)                                                                                                                             
    wlink = r'https://portal.nersc.gov/project/desi/users/mjwilson/SV-ASSIGN/tile-{}.fits'.format(tile)

    print(hlink)

    axarr[row][col].set_title('(HSC) Tile:  {}\n'.format(tile) + r'\url{%s}' % make_tiny(hlink) + '\n' + r'\url{%s}' % make_tiny(wlink))
  
plt.tight_layout()

##  pl.show()
pl.savefig('assigned_HSC.png')

print('\n\nDone.\n\n')

