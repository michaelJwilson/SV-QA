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

##  L428 of https://github.com/desihub/desimodel/blob/master/py/desimodel/focalplane/geometry.py
params    =  load_desiparams()
fiber_dia =  params['fibers']['diameter_um']

#- Platescales in um/arcsec
ps        = load_platescale()

##  Add in GAMA labels.                                                                                                                                                                                                        
_fits     = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/elgs/hsc.fits')
dat       = Table(_fits[1].data)

##  dat   = dat[dat['MORPHTYPE'] == 'REX']
dat       = dat[dat['SNR']       <=    12]

des       = (0. < dat['RA']) & (dat['RA'] < 45.) & (-30. < dat['DEC']) & (dat['DEC'] < 5.)
isin      = des  ##  isin | des

##                                                                                                                                                                                                                                      
dat       = dat[isin]

dat.sort('SNR')

dat.pprint()

##  Light weight viewer friendly version.
viewer                = dat['RA', 'DEC']
  
##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L1080
##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L593
viewer['color']       = Column(data = ['yellow']  * len(dat), dtype='S16', length=len(dat))
viewer['radius']      = 10.   ##  [arcsecond] ~ 2x DESI projected fiber size. 
viewer['fillOpacity'] = 0.5
viewer['lineWidth']   = 1.

names                 = ['{:.1f}'.format(x) for x in dat['FRANKENZ']]  ##  [''] * len(names)

viewer['name']        = Column(data = names, dtype='S18', length=len(dat))

##                                                                                                                                                                                                  
_fits     = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/elgs/elgs.fits')
dat       = Table(_fits[1].data)

##  dat   = dat[dat['MORPHTYPE'] == 'REX']
dat       = dat[dat['SNR']       <=    12]

##  [248.5, 242.6, 42., 44.4]
isin      = (242.6 < dat['RA']) & (dat['RA'] < 248.5) & (42. < dat['DEC']) & (dat['DEC'] < 44.4)

##  [215.5, 213., 51.7, 53.0]
isin      = isin | ((dat['RA'] < 215.5) & (dat['RA'] > 213.0) & (51.7 < dat['DEC']) & (dat['DEC'] < 53.0))

des       = (0. < dat['RA']) & (dat['RA'] < 45.) & (-30. < dat['DEC']) & (dat['DEC'] < 5.)
isin      = des  ##  isin | des

##
dat       = dat[isin]

##
dat.pprint()

##  Light weight viewer friendly version.                                                                                                                                                                                            
mviewer               =  dat['RA', 'DEC']

##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L1080                                                                                                                                                 
##  https://github.com/legacysurvey/decals-web/blob/master/templates/index.html#L593                                                                                                                                                   
mviewer['color']       = Column(data = ['green']  * len(dat), dtype='S16', length=len(dat))
mviewer['radius']      = 5.   ##  [arcsecond] ~ 1x DESI projected fiber size.                                                                                                                                                          
mviewer['fillOpacity'] = 0.5
mviewer['lineWidth']   = 1.
mviewer['name']        = Column(data = ['']  * len(dat), dtype='S18', length=len(dat))

mviewer['color'][dat['MORPHTYPE'] == 'REX'] = 'red'
mviewer['color'][dat['MORPHTYPE'] == 'DEV'] = 'blue'
mviewer['color'][dat['MORPHTYPE'] == 'EXP'] = 'green'
mviewer['color'][dat['MORPHTYPE'] == 'PSF'] = 'white'
mviewer['color'][dat['MORPHTYPE'] == 'DUP'] = 'cyan'

viewer                 = vstack([viewer, mviewer])

print(viewer)

'''
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
'''    

##
viewer.write('/project/projectdirs/desi/www/users/mjwilson/SV-ELGTARGETS/tile-targets.fits', format='fits', overwrite=True)

##  Set permissions on viewer files.
os.system('chmod --reference=/project/projectdirs/desi/www/users/mjwilson/plots/visibility-nofullmoon-26-0.pdf /project/projectdirs/desi/www/users/mjwilson/SV-ELGTARGETS/*')

