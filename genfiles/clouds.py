import os
import numpy               as      np
import pylab               as      pl
import astropy.units       as      u
import astropy.io.fits     as      fits
import matplotlib.pyplot   as      plt

from   astropy.table       import  Table, vstack, unique, join
from   astropy.coordinates import  SkyCoord
from   desimodel.footprint import  pix2tiles, find_tiles_over_point, is_point_in_desi


##  https://ui.adsabs.harvard.edu/abs/2014ApJ...786...29S/abstract
def get_clouds(fname, plotit=False):
  root  = os.environ['CSCRATCH'] + '/BGS/SV-ASSIGN/clouds/'

  dat   = Table(fits.open(root + fname)[1].data) 
  names = dat['name']
  stars = dat['nstar_tot']

  gc    = SkyCoord(l=dat['l'] * u.degree, b=dat['b'] * u.degree, frame='galactic')
  gc    = gc.transform_to('icrs') 
      
  result = Table()

  result['l']     = dat['l']
  result['b']     = dat['b']
  result['RA']    = gc.ra.deg
  result['DEC']   = gc.dec.deg
  result['NAME']  = names

  return result


if __name__ == '__main__':
  scratch    = os.environ['CSCRATCH']
  _file      = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')
  tiles      = Table(_file[1].data)['AIRMASS', 'STAR_DENSITY', 'IMAGEFRAC_R', 'IMAGEFRAC_GRZ', 'TILEID', 'RA', 'DEC']
  
  bclouds    = get_clouds('bigcloud.fits')
  mclouds    = get_clouds('mbmcloud.fits')

  clouds     = vstack([bclouds, mclouds])   
  clouds.write(os.environ['CSCRATCH'] + '/BGS/SV-ASSIGN/clouds/clouds.fits', format='fits', overwrite=True)

  ##                                                                                                                                                                                                                                    
  nside      = 64

  (indesi, tindex)           = is_point_in_desi(tiles, clouds['RA'].quantity.value, clouds['DEC'].quantity.value, return_tile_index=True)

  ##  Prevent 'duplication' of e.g. RA cols after join.                                                                                                                                                                                
  del tiles['RA']
  del tiles['DEC']

  ##  Some objects in the .mtl are not in the tiles;  rough cut on creation I guess.                                                                                                                                                   
  clouds['TILEID']          = -99 * np.ones_like(clouds['RA'].quantity.value, dtype=np.int)
  clouds['TILEID'][indesi]  = tiles['TILEID'].quantity[tindex[indesi]]
  clouds.sort('TILEID')

  clouds     = clouds[clouds['TILEID'] > 0]
  clouds     = join(clouds, tiles, keys=['TILEID'], join_type='left')
  
  ##  Unique TILEIDS.                                                                                                                                                                                                                  
  utiles     = np.unique(clouds['TILEID'].quantity)

  for _tile in utiles:
    fname    = scratch + '/BGS/SV-ASSIGN/clouds/clouds_{:06d}.fits'.format(_tile)
    tile_cut = clouds[clouds['TILEID'] == _tile]
    
    print(tile_cut)
    
    tile_cut.write(fname, format='fits', overwrite=True)

print('\n\nDone.\n\n')

  
