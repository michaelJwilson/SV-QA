import  os
import  glob
import  random
import  numpy                  as      np
import  pylab                  as      pl
import  matplotlib.pyplot      as      plt
import  astropy.io.fits        as      fits
import  astropy.units          as      u
import  healpy                 as      hp

from    astropy.table          import  Table, join, hstack
from    astropy                import  constants as const
from    desitarget.geomask     import  circles
from    desiutil.bitmask       import  BitMask
from    desitarget.targetmask  import  load_mask_bits
from    astropy.coordinates    import  SkyCoord
from    astropy                import  units     as u
from    desimodel.footprint    import  pix2tiles, find_tiles_over_point, is_point_in_desi


scratch    = os.environ['CSCRATCH']

## 
gaia       = Table(fits.open(scratch + '/BGS/SV-ASSIGN/masks/gaia_stellar_mask.fits')[1].data, names=['ra', 'dec', 'G'])
print('GAIA (G < 13) mask stars: {}'.format(len(gaia)))

##  http://legacysurvey.org/dr8/external/
tycho      = Table(fits.open('/global/project/projectdirs/cosmo/staging/tycho2/tycho2.kd.fits')[1].data)
tycho      = tycho[tycho['MAG_VT'] < 13.]
print('Tycho (VT < 13) mask stars: {}'.format(len(tycho)))

##  Large galaxy cat.
lgal       = Table(fits.open('/global/project/projectdirs/cosmo/staging/largegalaxies/v2.0/LSLGA-v2.0.fits')[1].data)
print('LSLGA mask galaxies: {}'.format(len(lgal)))

##
ngcsc      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/NGC-star-clusters.fits', ignore_missing_end=True)[1].data)
print('NGCSC mask stars: {}'.format(len(ngcsc)))

##
_assigned  = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/assigned_targets.fits')[1].data)
_assigned.sort('BRICK_OBJID')

##  Retain only the science targets; Dropping standards etc. 
_assigned  = _assigned[_assigned['FA_TYPE'].quantity == 1]

##  Load one tile, currently. 
sv_mtl     = Table(fits.open(scratch + '/BGS/SV-ASSIGN/mtls/svmtl_074057.fits')[1].data)

##  Limit to only BGS.
sv_mtl     = sv_mtl[sv_mtl['SV1_BGS_TARGET'].quantity.value > 0]

##  
_bitdefs   = load_mask_bits("sv1")
bgs_mask   = BitMask('sv1_bgs_mask', _bitdefs)

types      =  bgs_mask.names()
bits       = [bgs_mask.bitnum(x) for x in types]

##  
assigned   = join(_assigned, sv_mtl, keys=['BRICKID', 'BRICK_OBJID'], join_type='inner')
assert  len(_assigned) == len(assigned)

splits     = {}
title      = ''

for _type in ['BGS_FAINT', 'BGS_BRIGHT', 'BGS_FAINT_EXT', 'BGS_LOWQ', 'BGS_FIBMAG']:
  sample        = _type.split('_')
  sample.pop(0)
  sample        = ''.join(x for x in sample)

  splits[_type] = assigned[(assigned['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum(_type)) != 0]
  title        += sample + ': {}'.format(len(splits[_type])) + ';  '
  
  print(sample, len(splits[_type]))

##
sep_limit                  = .5

##  Separation from bright objects. 
cassigned                  = SkyCoord(assigned['RA'].quantity * u.deg, assigned['DEC'].quantity * u.deg, frame='icrs')    

##
print('Solving for Tycho separation.')
ctycho                     = SkyCoord(   tycho['RA'].quantity * u.deg,    tycho['DEC'].quantity * u.deg, frame='icrs')
idxc, idxcatalog, d2d, d3d = cassigned.search_around_sky(ctycho, sep_limit * u.deg)
d2d                        = np.array([x.degree for x in d2d])

pl.hist(60. * d2d, bins=50, alpha=0.5)
pl.xlabel('Tycho (VT < 13.) separation [arcmin.]')
pl.xscale('log')
pl.yscale('log')
pl.savefig('plots/tycho.png')
pl.clf()

##
print('Solving for LSLGA separation.')
clgal                      = SkyCoord(    lgal['RA'].quantity * u.deg,     lgal['DEC'].quantity * u.deg, frame='icrs')
idxc, idxcatalog, d2d, d3d = cassigned.search_around_sky(clgal, sep_limit * u.deg)
d2d                        = np.array([x.degree for x in d2d])

pl.hist(60. * d2d, bins=50, alpha=0.5)
pl.xlabel('LSLGA separation [arcmin]')
pl.xscale('log')
pl.yscale('log')
pl.savefig('plots/LSLGA.png')
pl.clf()

##
print('Solving for NGCSC separation.')
cngcsc                     = SkyCoord(   ngcsc['ra'].quantity * u.deg,    ngcsc['dec'].quantity * u.deg, frame='icrs')
idxc, idxcatalog, d2d, d3d = cassigned.search_around_sky(cngcsc, sep_limit * u.deg)                                                                                                     
d2d                        = np.array([x.degree for x in d2d])
pl.hist(60. * d2d, bins=50, alpha=0.5)
pl.xlabel('NGCSC separation [arcmin]')
pl.xscale('log')
pl.yscale('log')
pl.savefig('plots/NGCSC.png')    
pl.clf()

##
print('Solving for GAIA separation.')                
cgaia                      = SkyCoord(    gaia['ra'].quantity * u.deg,     gaia['dec'].quantity * u.deg, frame='icrs')
idxc, idxcatalog, d2d, d3d = cassigned.search_around_sky(cgaia, sep_limit * u.deg)
d2d                        = np.array([x.degree for x in d2d])
pl.hist(60. * d2d, bins=50, alpha=0.5)
pl.xlabel('GAIA (G < 13) separation [arcmin]')
pl.xscale('log')
pl.yscale('log')
pl.savefig('plots/GAIA.png')
pl.clf()

