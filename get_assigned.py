import os
import sys
import pylab           as     pl
import numpy           as     np
import astropy.io.fits as     fits

from   astropy.table   import Table, vstack


scratch  = os.environ['CSCRATCH']

tiles    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)
utiles   = np.unique(tiles['TILEID'].quantity)

result   = Table(names=['TARGETID', 'FA_TYPE', 'BRICK_OBJID', 'TARGET_RA', 'TARGET_DEC', 'BRICKID'], dtype=('i8', 'i1', 'i4', 'f8', 'f8', 'i4'))

for tile in utiles:
  ##  Scrape from the tile picker site. 
  ##  cmd = 'wget http://www.astro.utah.edu/~u6022465/SV/tiles/SV_BGS/fits_files/tile-{:06}.fits -O /global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile, tile)
  ##  os.system(cmd)

  _fits  = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile))
  
  dat    = Table(_fits[1].data)['TARGETID', 'FA_TYPE', 'BRICK_OBJID', 'TARGET_RA', 'TARGET_DEC', 'BRICKID']
  result = vstack([result, dat])

  print(dat['BRICKID'])
  print(result)
  
  ##  exit(1)
  
result.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/assigned_targets.fits', format='fits', overwrite=True)

