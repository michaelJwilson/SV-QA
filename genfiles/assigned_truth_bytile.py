import  os
import  glob
import  numpy                 as      np
import  astropy.io.fits       as      fits

from    astropy.table         import  Table
from    desiutil.bitmask      import  BitMask
from    desitarget.targetmask import  load_mask_bits
from    desimodel.footprint   import  is_point_in_desi


##
nside   = 64

_file   = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')
tiles   = Table(_file[1].data)['AIRMASS', 'STAR_DENSITY', 'IMAGEFRAC_R', 'IMAGEFRAC_GRZ', 'TILEID', 'RA', 'DEC']

##  https://desi.lbl.gov/trac/wiki/TargetSelectionWG/TargetingTruthTables/MatchedTruthCatalogs
scratch = os.environ['CSCRATCH']

files   = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/*.fits')

for _file in files:    
  ##  Truth spectroscopic catalogue.  Look for a given survey in both north and south. 
  truth      = Table(fits.open(_file)[1].data)
  
  (indesi, tindex)        = is_point_in_desi(tiles, truth['RA'].quantity.value, truth['DEC'].quantity.value, return_tile_index=True)

  ##  Prevent 'duplication' of RA cols after join.                                                                                                                                                                                      
  ##  del tiles['RA']
  ##  del tiles['DEC']

  ##  Some objects in the .mtl are not in the tiles;  rough cut on creation I guess.                                                                                                                                                    
  truth['TILEID']         = -99 * np.ones_like(truth['RA'].quantity.value, dtype=np.int)
  truth['TILEID'][indesi] = tiles['TILEID'].quantity[tindex[indesi]]
  truth = truth[truth['TILEID'].quantity > 0]
  
  truth.sort('TILEID')

  if len(truth) > 0:
    print(truth)
    
    splits = _file.split('/')
    fname  = '/'.join(x for x in splits[:-1])
    tail   = '/' + splits[-1].split('.fits')[0]
  
    ##  Unique TILEIDS.
    utiles     = np.unique(truth['TILEID'].quantity)

    for _tile in utiles:
      out  = truth[truth['TILEID'] == _tile]

      print('Writing {} ...'.format(fname + '/bytile/' + tail + '_{:06d}.fits'.format(_tile)))

      if len(out) > 0:
        out.write(fname + '/bytile/' + tail + '_{:06d}.fits'.format(_tile), format='fits', overwrite=True)
