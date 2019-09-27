import os
import sys
import pylab                  as      pl
import numpy                  as      np
import astropy.io.fits        as      fits

from   astropy.table          import  Table, vstack, Column
from   desitarget.targetmask  import  load_mask_bits
from   desiutil.bitmask       import  BitMask
from   desimodel.io           import  load_desiparams, load_platescale


scratch   = os.environ['CSCRATCH']

tiles     = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)
utiles    = np.unique(tiles['TILEID'].quantity)

##
_bitdefs  = load_mask_bits("sv1")

desi_mask = BitMask('sv1_desi_mask', _bitdefs)
bgs_mask  = BitMask('sv1_bgs_mask',  _bitdefs)

types     =  bgs_mask.names()
bits      = [bgs_mask.bitnum(x) for x in types]

##  L428 of https://github.com/desihub/desimodel/blob/master/py/desimodel/focalplane/geometry.py
params    =  load_desiparams()
fiber_dia =  params['fibers']['diameter_um']

#- Platescales in um/arcsec
ps        = load_platescale()

##  Add in GAMA labels.                                                                                                                                                                                                        
_fits     = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/legacy/ls-GAMA-south.fits')
gama      = Table(_fits[1].data)

for tile in utiles:  
  ##  Scrape from the tile picker site. 
  ##  cmd = 'wget http://www.astro.utah.edu/~u6022465/SV/tiles/SV_BGS/fits_files/tile-{:06}.fits -O /global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile, tile)
  ##  os.system(cmd)
  try:
    print('Solving for Tile {}.'.format(tile))

    fname  = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/mtls/svmtl_{:06}.fits'.format(tile)
    _fits  = fits.open(fname)  
    dat    = Table(_fits[1].data)['TARGETID', 'BRICK_OBJID', 'RA', 'DEC', 'BRICKID', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET']
    
  except:
    print('Unable to retrieve {}.'.format(fname))
    
    continue
  
  isbgs    = (dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('BGS_ANY')) != 0

  exclude  = (dat['SV1_DESI_TARGET'] & desi_mask.mask('BRIGHT_OBJECT|IN_BRIGHT_OBJECT|NEAR_BRIGHT_OBJECT') != 0)

  print('Tile {} excludes {} for bright mask.'.format(tile, len(dat[exclude])))
  
  ##  print(dat.columns)
  ##  print(result)
  
  ##  Light weight viewer friendly version.
  viewer                  = dat['RA', 'DEC']
  
  ##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L1080
  ##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L593
  viewer['color']        = Column(data = ['black']  * len(dat), dtype='S16', length=len(dat))
  viewer['radius']       = 0.0   ## [arcsecond] ~ 2x DESI projected fiber size. 
  viewer['fillOpacity']  = 0.5
  viewer['lineWidth']    = 2.
  viewer['name']         = Column(data = ['']  * len(dat), dtype='S18', length=len(dat))

  ##  BGS_LOWQ                                                                                                                                                                                                                          
  viewer['color'][~exclude & ((dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_LOWQ'))      != 0)]  = 'purple'


    ##  Append mask sources.                                                                                                                                                                                     
    _fits   = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/gaia_{:06}.fits'.format(tile))
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
    randoms = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/randoms/randoms_{:06}.fits'.format(tile))[1].data)
    isin    = (randoms['MASKBITS'] & 2**1) != 0  
    randoms = randoms[isin]

    mviewer          = randoms['RA', 'DEC']

    mviewer['RA']    = mviewer['RA']
    mviewer['DEC']   = mviewer['DEC']
  
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

    print(mviewer)
    
    viewer           = vstack([viewer, mviewer])

  except:
    print('Unable to retrieve molecular clouds for Tile {}'.format(tile))    
    
  viewer.write('/project/projectdirs/desi/www/users/mjwilson/SV-TARGETS/tile-targets-{:06}.fits'.format(tile), format='fits', overwrite=True)
  
##  Set permissions on viewer files.
os.system('chmod --reference=/project/projectdirs/desi/www/users/mjwilson/plots/visibility-nofullmoon-26-0.pdf /project/projectdirs/desi/www/users/mjwilson/SV-TARGETS/*')

