import os
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


scratch    = os.environ['CSCRATCH']

tiles      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)

sv_mtl     = Table(fits.open(scratch + '/BGS/SV-ASSIGN/mtls/svmtl_074036.fits')[1].data)

##  Limit to only BGS.                                                                                                                                                                                                       
##  sv_mtl = sv_mtl[sv_mtl['SV1_BGS_TARGET'].quantity.value > 0]

_assigned  = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/assigned_targets.fits')[1].data)
_assigned.sort('BRICK_OBJID')

##  Those that were assigned in this tile.                                                                                                                                                                                   
assigned   = join(_assigned, sv_mtl, keys=['BRICKID', 'BRICK_OBJID'], join_type='inner')
atypes, _  = np.unique(assigned['FA_TYPE'].quantity, return_counts=True)

print(atypes)
print(_)

exit(1)

gaia       = Table(fits.open(scratch + '/BGS/SV-ASSIGN/masks/gaia_stellar_mask.fits')[1].data, names=['RA', 'DEC', 'G'])

legacys    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/GAMA-south.fits')[1].data)

_file      = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/standard/bytile/GAMA-south-standard_074036.fits'  
dat        = Table(fits.open(_file)[1].data)

pl.title('Tile {:06d}'.format(74056))
  
## 
pl.scatter(dat['RA'], dat['DEC'], s=1, rasterized=True, alpha=0.5, marker='.', c='gold')

ax = pl.gca()

pl.xlim(ax.get_xlim())
pl.ylim(ax.get_ylim())

##  axarr[row][col].scatter(legacys['RA'], legacys['DEC'], s=1, rasterized=True, alpha=0.5, marker='x', c='k') 
  
pl.scatter(gaia['RA'], gaia['DEC'], s=1, rasterized=True, alpha=0.5, marker='.', c='cyan')
    
plt.tight_layout()

pl.show()
##  pl.savefig('spec_truth.png')

print('\n\nDone.\n\n')

