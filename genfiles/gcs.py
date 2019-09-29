import numpy               as     np
import astropy.io.fits     as     fits
import astropy.units       as     u

from astropy.table         import Table, vstack
from desimodel.footprint   import is_point_in_desi, tiles2pix
from desimodel.io          import load_tiles 
from desitarget.targetmask import desi_mask


gcs            = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/masks/NGC-star-clusters.fits')[1].data)

tiles          = np.array([list(x) for x in load_tiles(onlydesi=True)])
tiles          = Table([tiles[:, 0].astype(np.int), tiles[:, 1].astype(np.float32), tiles[:,2].astype(np.float32)], names=['TILEID', 'RA', 'DEC'])

isin, index    = is_point_in_desi(tiles, gcs['ra'], gcs['dec'], return_tile_index=True)

gcs['TILEID']       = -99
gcs['TILEID'][isin] = tiles['TILEID'][index[isin]]

gcs            = gcs[gcs['name'] == 'NGC5904']
gcs['RA']      = gcs['ra']
gcs['DEC']     = gcs['dec']
gcs            = gcs['RA', 'DEC', 'name', 'TILEID']

tiles          = tiles[tiles['TILEID'] == gcs['TILEID']]

gcs.pprint()
tiles.pprint()

pix            = tiles2pix(nside=2, tiles=tiles)

for p in pix:
  targets                 = Table(fits.open('/project/projectdirs/desi/target/catalogs/dr8/0.32.0/targets/main/resolve/bright/targets-dr8-hp-{}.fits'.format(p))[1].data)
  targets                 = targets[(targets['DESI_TARGET'] & 2 ** desi_mask.bitnum('BGS_ANY')) != 0]
  isin, index             = is_point_in_desi(tiles, targets['RA'], targets['DEC'], return_tile_index=True)
  targets['TILEID']       = -99
  targets['TILEID'][isin] = tiles['TILEID'][index[isin]]

  targets         = vstack([gcs, targets])
  targets['name'] = ''
  targets         = targets['RA', 'DEC', 'TILEID', 'name']
  targets         = targets[targets['TILEID'] > 0]
  
  print(targets)

  targets.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/main/targets-dr8-hp-{}.fits'.format(p), format='fits', overwrite=True)
  
