import  os
import  glob
import  random
import  fitsio
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


##  https://faun.rc.fas.harvard.edu/eschlafly/desi/tiling/dr8/note.pdf
##  AIRMASS:        Airmass if observed 15. deg. from transit.
##  STAR_DENSITY:   Median # of GAIA stars with G < 19.5 per sq. deg. in tile. 
##  IMAGEFRAC R:    Fraction of this tile within 1.605 deg. of r imaging.     
##  IMAGEFRAC GRZ:  Fraction of this tile within 1.605 deg. of grz imaging. 
scratch    = os.environ['CSCRATCH']
_file      = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')
tiles      = Table(_file[1].data)['AIRMASS', 'STAR_DENSITY', 'IMAGEFRAC_R', 'IMAGEFRAC_GRZ', 'TILEID', 'RA', 'DEC']

print(is_point_in_desi(tiles, 347.2099, 20.7548, return_tile_index=True))
print(tiles[31])

exit(1)

dat        = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits')
randoms    = Table(dat[1][:])
  
##
nside      = 64
  
(indesi, tindex)           = is_point_in_desi(tiles, randoms['RA'].quantity.value, randoms['DEC'].quantity.value, return_tile_index=True)

##  Prevent 'duplication' of e.g. RA cols after join. 
del tiles['RA']
del tiles['DEC']

##  Some objects in the .mtl are not in the tiles;  rough cut on creation I guess. 
randoms['TILEID']          = -99 * np.ones_like(randoms['RA'].quantity.value, dtype=np.int)
randoms['TILEID'][indesi]  = tiles['TILEID'].quantity[tindex[indesi]]
randoms.sort('TILEID')
  
##  Not associated to a tile?
##  notin = sv_mtl[sv_mtl['TILEID'] < 0]
##  notin.sort('DEC')
##  print(notin)

##  pl.plot(notin['RA'], notin['DEC'], c='k', marker='.', lw=0, markersize=1)
##  pl.show()

##  Only keep the targets in a given tile. 
randoms    = randoms[randoms['TILEID'] > 0]
randoms    = join(randoms, tiles, keys=['TILEID'], join_type='left')

##  Unique TILEIDS.
utiles     = np.unique(randoms['TILEID'].quantity)

for _tile in utiles:
  fname    = scratch + '/BGS/SV-ASSIGN/randoms/randoms_{:06d}.fits'.format(_tile)
  tile_cut = randoms[randoms['TILEID'] == _tile]

  print(tile_cut)
  
  tile_cut.write(fname, format='fits', overwrite=True)


  print('\n\nDone.\n\n')
