import  os
import  glob
import  random
import  numpy                  as      np
import  pylab                  as      pl
import  matplotlib.pyplot      as      plt
import  astropy.io.fits        as      fits
import  astropy.units          as      u
import  healpy                 as      hp

from    astropy.table          import  Table, join, hstack, vstack
from    astropy                import  constants as const
from    desitarget.geomask     import  circles
from    desiutil.bitmask       import  BitMask
from    desitarget.targetmask  import  load_mask_bits
from    astropy.coordinates    import  SkyCoord
from    astropy                import  units     as u
from    desimodel.footprint    import  pix2tiles, find_tiles_over_point, is_point_in_desi


##  https://faun.rc.fas.harvard.edu/eschlafly/desi/tiling/dr8/note.pdf
##  AIRMASS:        Airmass if observed 15. deg. from transit.
##  STAR_DENSITY:   Median # of GAIA stars with G < 19.5 per sq. deg. in tile. 
##  IMAGEFRAC R:    Fraction of this tile within 1.605 deg. of r imaging.     
##  IMAGEFRAC GRZ:  Fraction of this tile within 1.605 deg. of grz imaging. 
_file      = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')
tiles      = Table(_file[1].data)['AIRMASS', 'STAR_DENSITY', 'IMAGEFRAC_R', 'IMAGEFRAC_GRZ', 'TILEID', 'RA', 'DEC']

##                                                                                                                                                                                                                                      
scratch    = os.environ['CSCRATCH']
_sweeps    = glob.glob(scratch + '/BGS/SV-ASSIGN/sweeps/*.fits')

cols       = ['PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'ANYMASK_R', 'ANYMASK_G', 'ANYMASK_Z']
load       = cols + ['TARGETID']

sweep      = Table(fits.open(_sweeps[0])[1].data)[load]

for i, _sweep in enumerate(_sweeps[1:]):
  sweep    = vstack(sweep, Table(fits.open(_sweep)[1].data)[load])

  print('{} of {}.'.format(i, len(_sweeps)))

  if i > 4:
    break
  
##
##  files  = glob.glob(scratch + '/BGS/SV-ASSIGN/mtls/svmtl_*.fits')
files      = glob.glob(scratch + '/BGS/SV-ASSIGN/fiberassign/tile-*.fits')
                       
for fname in files:
  tile     = Table(fits.open(fname)[1].data)

  for x in cols:
    try:
      del tile[x]
      
    except:
      print('Column {} not available.'.format(x))
      
  tile     = join(tile, sweep, keys='TARGETID', join_type='left')

  print(tile)

  ##  tile.write(fname, format='fits', overwrite=True)
  
print('\n\nDone.\n\n')
