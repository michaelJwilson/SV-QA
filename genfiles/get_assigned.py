import os
import sys
import pylab                  as      pl
import numpy                  as      np
import astropy.io.fits        as      fits

from   astropy.table          import  Table, vstack, Column
from   desitarget.targetmask  import  load_mask_bits
from   desiutil.bitmask       import  BitMask
from   desimodel.io           import  load_desiparams, load_platescale


scratch      = os.environ['CSCRATCH']

tiles        = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)
utiles       = np.unique(tiles['TILEID'].quantity)

result       = Table(names=['TARGETID', 'FA_TYPE', 'BRICK_OBJID', 'TARGET_RA', 'TARGET_DEC', 'BRICKID'], dtype=('i8', 'i1', 'i4', 'f8', 'f8', 'i4'))

##
_bitdefs     = load_mask_bits("sv1")

desi_mask    = BitMask('sv1_desi_mask', _bitdefs)
bgs_mask     = BitMask('sv1_bgs_mask',  _bitdefs)

types        =  bgs_mask.names()
bits         = [bgs_mask.bitnum(x) for x in types]

##  L428 of https://github.com/desihub/desimodel/blob/master/py/desimodel/focalplane/geometry.py
params       =  load_desiparams()
fiber_dia    = params['fibers']['diameter_um']

#- Platescales in um/arcsec
ps           = load_platescale()

for tile in utiles:
  print('Solving for Tile {}.'.format(tile))
  
  ##  Scrape from the tile picker site. 
  ##  cmd = 'wget http://www.astro.utah.edu/~u6022465/SV/tiles/SV_BGS/fits_files/tile-{:06}.fits -O /global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile, tile)
  ##  os.system(cmd)

  _fits   = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile))
  
  dat     = Table(_fits[1].data)['TARGETID', 'FA_TYPE', 'BRICK_OBJID', 'TARGET_RA', 'TARGET_DEC', 'BRICKID', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'BGS_TARGET', 'DESIGN_X','DESIGN_Y', 'DESI_TARGET']
  
  result  = vstack([result, dat])

  isbgs   = ((dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('BGS_ANY')) != 0)

  ##  BRIGHT_OBJECT, IN_BRIGHT_OBJECT, NEAR_BRIGHT_OBJECT. 
  ##  exclude = (dat['SV1_DESI_TARGET'] & desi_mask.mask('BRIGHT_OBJECT') != 0) | (dat['SV1_DESI_TARGET'] & desi_mask.mask('IN_BRIGHT_OBJECT') != 0) | (dat['SV1_DESI_TARGET'] & desi_mask.mask('NEAR_BRIGHT_OBJECT') != 0)
  exclude = (dat['SV1_DESI_TARGET'] & desi_mask.mask('BRIGHT_OBJECT|IN_BRIGHT_OBJECT|NEAR_BRIGHT_OBJECT') != 0)
    
  ##  print(dat.columns)
  ##  print(result)

  ##  Calculate fiber size.
  x            = np.asarray(dat['DESIGN_X'])
  y            = np.asarray(dat['DESIGN_Y'])
  r            = np.sqrt(x**2. + y**2.)

  radial_scale = np.interp(r, ps['radius'], ps['radial_platescale'])
  az_scale     = np.interp(r, ps['radius'], ps['az_platescale'])

  ##  radial and azimuthal fiber radii in arcsec
  rr           = 0.5 * fiber_dia / radial_scale
  raz          = 0.5 * fiber_dia / az_scale

  fiber_area   = (np.pi * rr * raz)

  geomean      = np.sqrt(rr * raz)
  
  ##  Light weight viewer friendly version.
  viewer                 = dat['TARGET_RA', 'TARGET_DEC']

  viewer['RA']           = viewer['TARGET_RA']
  viewer['DEC']          = viewer['TARGET_DEC']

  del viewer['TARGET_RA']
  del viewer['TARGET_DEC']
  
  ##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L1080
  ##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L593
  viewer['color']        = Column(data = ['black']  * len(dat), dtype='S16', length=len(dat))
  ##  viewer['border-style'] = Column(data = ['dashed'] * len(dat), dtype='S16', length=len(dat))
  ##  viewer['opacity']      = 0.3
  viewer['radius']       = geomean
  viewer['lineWidth']    = 2.
  viewer['name']         = ''
  
  ##  print(viewer)

  ##  FA_TYPE definition:  https://github.com/desihub/fiberassign/blob/master/src/targets.h#L28
  ##
  ##  --  Skies and Standards --                                                                                                                                                                            
  viewer['color'][(dat['DESI_TARGET'] & 2 ** desi_mask.bitnum('BAD_SKY'))         != 0 ]  = 'teal'
  viewer['color'][(dat['DESI_TARGET'] & 2 ** desi_mask.bitnum('SKY'))             != 0 ]  = 'blue'
  viewer['color'][(dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('MWS_ANY')  )   != 0 ]  = 'acqua'
  viewer['color'][(dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('STD_WD')  )    != 0 ]  = 'white'
  viewer['color'][(dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('STD_BRIGHT'))  != 0 ]  = 'lime'

  ##  Safe.
  viewer['color'][(dat['FA_TYPE'] & 2 ** 3) != 0]                                         = 'Fuchsia'
  
  ##
  ##  --  Science --
  ##
  ##  Note:  reverse order to expected redshift success (previous colors overwritten for overlap.)
  ##
  ##  BGS_LOWQ
  viewer['color'][~exclude & ((dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_LOWQ'))      != 0)]  = 'purple'

  ##  BGS_FAINT_EXT                                                                                                                                                                                               
  viewer['color'][~exclude & ((dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_FAINT_EXT')) != 0)]  = 'maroon'

  ##  BGS_FIBMAG                                                                                                                                                                                                  
  viewer['color'][~exclude & ((dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_FIBMAG'))    != 0)]  = 'green'

  ##  BGS_FAINT                                                                                                                                                                                                   
  viewer['color'][~exclude & ((dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_FAINT'))     != 0)]  = 'silver'
  
  ##  BGS_BRIGHT      
  viewer['color'][~exclude & ((dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_BRIGHT'))    != 0)]  = 'yellow'

  if len(dat[viewer['color'] == 'black']) > 0:
    print( dat[viewer['color'] == 'black'])
  
  try:
    ##  Append mask sources.
    _fits   = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/tycho_{:06}.fits'.format(tile))
    mask    = Table(_fits[1].data)

    mviewer              = mask['RA', 'DEC', 'MRADIUS']
    mviewer['radius']    = mviewer['MRADIUS']
    mviewer['color']     = Column(data = ['red']  * len(mask), dtype='S16', length=len(mask))    
    mviewer['lineWidth'] = 2.
    mviewer['name']      = ''
    
    del mviewer['MRADIUS']
  
    viewer = vstack([viewer, mviewer])
    
  except:
    print('Unable to retrieve Tycho mask for Tile {}'.format(tile))

    
  try:
    ##  Append mask sources.                                                                                                                                                                                    
    _fits                = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/gaia_{:06}.fits'.format(tile))
    mask                 = Table(_fits[1].data)
    mask                 = mask[mask['ISBRIGHT'] == True]
    
    mviewer              = mask['RA', 'DEC', 'MRADIUS']
    mviewer['radius']    = mviewer['MRADIUS']
    mviewer['color']     = Column(data = ['red']  * len(mask), dtype='S16', length=len(mask))
    mviewer['lineWidth'] = 2.
    mviewer['name']      = ''
    
    del mviewer['MRADIUS']

    viewer = vstack([viewer, mviewer])

  except:
    print('Unable to retrieve Tycho mask for Tile {}'.format(tile))

  try:
    randoms = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/randoms/randoms_{:06}.fits'.format(tile))[1].data)
    isin    = (randoms['MASKBITS'] & 2**1) != 0
    randoms = randoms[isin]

    mviewer = randoms['RA', 'DEC']

    mviewer['color'] = Column(data = ['red']  * len(randoms), dtype='S16', length=len(randoms))
    mviewer['name']  = ''

    viewer  = vstack([viewer, mviewer])

  except:
    print('Unable to retrieve randoms for Tile {}'.format(tile))
    
  try:
    clouds           = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/clouds/clouds_{:06}.fits'.format(tile))[1].data)

    mviewer          = clouds['RA', 'DEC']
    mviewer['color'] = Column(data = ['lime']  * len(clouds), dtype='S16', length=len(clouds))
    mviewer['name']  = Column(data = clouds['NAME'], dtype='S16', length=len(clouds))

    viewer           = vstack([viewer, mviewer])

  except:
    print('Unable to retrieve molecular clouds for Tile {}'.format(tile))
  
  viewer.write('/project/projectdirs/desi/www/users/mjwilson/SV-ASSIGN/tile-{:06}.fits'.format(tile), format='fits', overwrite=True)

  
##  Complete catalogue. 
result.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/assigned_targets.fits', format='fits', overwrite=True)

##  Set permissions on viewer files.
os.system('chmod --reference=/project/projectdirs/desi/www/users/mjwilson/plots/visibility-nofullmoon-26-0.pdf /project/projectdirs/desi/www/users/mjwilson/SV-ASSIGN/*')

print('Number of safe: {}'.format(len(result[result['FA_TYPE'] == 4])))
print('\n\nDone.\n\n')
