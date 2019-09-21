import glob
import random
import numpy             as       np
import pylab             as       pl
import matplotlib.pyplot as       plt
import astropy.io.fits   as       fits
import astropy.units     as       u

from   astropy.table      import  Table, join, hstack
from   astropy            import  constants as const
from   desitarget.geomask import  circles


nrow       = 2
ncol       = 5

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
  
  axarr[row][col].set_title('Tile {}'.format(tile))

  ## 
  axarr[row][col].scatter(dat['RA'], dat['DEC'], s=1, rasterized=True, alpha=0.5, marker='.', c='gold')

  axarr[row][col].set_xlim(axarr[row][col].get_xlim())
  axarr[row][col].set_ylim(axarr[row][col].get_ylim())

  axarr[row][col].scatter(legacys['RA'], legacys['DEC'], s=1, rasterized=True, alpha=0.5, marker='x', c='k') 

  ##
  axarr[row][col].scatter(dat['RA'], dat['DEC'], s=1, rasterized=True, alpha=0.5, marker='.', c='gold')
  
##  
plt.tight_layout()

pl.show()
##  pl.savefig('spec_truth.png')

print('\n\nDone.\n\n')

