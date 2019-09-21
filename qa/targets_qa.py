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

##  Load one tile, currently. 
Tile       = 74057
sv_mtl     = Table(fits.open(scratch + '/BGS/SV-ASSIGN/mtls/svmtl_{:06d}.fits'.format(Tile))[1].data)

##  Limit to only BGS.
sv_mtl     = sv_mtl[sv_mtl['SV1_BGS_TARGET'].quantity.value > 0]
print(len(sv_mtl))

##  
_bitdefs   = load_mask_bits("sv1")
bgs_mask   = BitMask('sv1_bgs_mask', _bitdefs)

types      =  bgs_mask.names()
bits       = [bgs_mask.bitnum(x) for x in types]

splits     = {}
title      = ''

##
for _type in ['BGS_FAINT', 'BGS_BRIGHT', 'BGS_FAINT_EXT', 'BGS_LOWQ', 'BGS_FIBMAG']:
  sample        = _type.split('_')
  sample.pop(0)
  sample        = ''.join(x for x in sample)

  splits[_type] = sv_mtl[(sv_mtl['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum(_type)) != 0]
  title        += sample + ': {}'.format(len(splits[_type])) + ';  '
  
  print(sample, len(splits[_type]))
  
##
colors     = plt.rcParams['axes.prop_cycle'].by_key()['color']
  
##  Source types.
morphtypes = set(sv_mtl['MORPHTYPE'])

## 
rmag       = 22.5 - 2.5 * np.log10(sv_mtl['FLUX_R'].quantity / sv_mtl['MW_TRANSMISSION_R'].quantity)
rfibmag    = 22.5 - 2.5 * np.log10(sv_mtl['FIBERFLUX_R'].quantity / sv_mtl['MW_TRANSMISSION_R'].quantity)

for mtype, color in zip(morphtypes, colors):
  isin     = sv_mtl['MORPHTYPE'] == mtype
  
  pl.plot(rmag[isin], rfibmag[isin], '.', markersize=1, c=color, label=mtype+'-{}'.format(len(sv_mtl[sv_mtl['MORPHTYPE'] == mtype])))   
  
  print(mtype, len(sv_mtl[sv_mtl['MORPHTYPE'] == mtype]))

for lim in [19.5, 20.1, 20.5]:
  pl.axvline(lim, ymin=0., ymax=1., c='k', alpha=0.7)
  
rfiblim    =  21.0511
pl.axhline(rfiblim, xmin=0., xmax=1., c='k', alpha=0.7)

pl.xlim(18.5, 21.0)
pl.ylim(18.5, 22.5)
pl.xlabel(r'$r$ mag.')
pl.ylabel(r'$r$ fiber mag.')
pl.title('Tile {}'.format(Tile))
pl.legend(frameon=False)

pl.show()
##  pl.savefig('plots/fibmag.png')
pl.clf()

##  (z - W1) - (3. / 2.5) * (g - r) + 1.2;
gmag  = 22.5 - 2.5 * np.log10(sv_mtl['FLUX_G'].quantity  / sv_mtl['MW_TRANSMISSION_G'].quantity)
zmag  = 22.5 - 2.5 * np.log10(sv_mtl['FLUX_Z'].quantity  / sv_mtl['MW_TRANSMISSION_Z'].quantity)
w1mag = 22.5 - 2.5 * np.log10(sv_mtl['FLUX_W1'].quantity / sv_mtl['MW_TRANSMISSION_W1'].quantity)

color = (zmag - w1mag) - (3. / 2.5) * (gmag - rmag) + 1.2
gd_zs = (color > (rfibmag - 19.6) / 2.5) & (rfibmag < 21.0)

pl.xlim(18., 23.)
pl.ylim(-1.,  3.)
pl.xlabel(r'$r$ fiber mag.')
pl.ylabel(r'($z-W1) - 3/2.5(g-r) + 1.2$')
pl.plot(rfibmag, color, '.', markersize=1, c='k')
pl.plot(rfibmag[gd_zs], color[gd_zs], 'rx', markersize=3)
pl.title(title, fontsize=9)
pl.show()

exit(1)

##
bins = np.arange(0., 10., 0.25)
    
for band in ['G', 'R', 'Z']:
  fracflux = assigned['FRACFLUX_{}'.format(band)].quantity
  masked   = len(fracflux[fracflux > 5.])

  pl.hist(fracflux, bins=bins, label=band + '-{}'.format(masked), alpha=0.5)

pl.axvline(5., c='k')
pl.xlim(0., 10)
pl.xlabel('FRACFLUX')
pl.yscale('log')
pl.legend(frameon=False)
pl.title(title, fontsize=9)
pl.savefig('plots/fracflux.png')
pl.clf()

##
bins = np.arange(0., 2., 0.025)

for band in ['G', 'R', 'Z']:
  fracmasked = assigned['FRACMASKED_{}'.format(band)].quantity
  masked     = len(fracmasked[fracmasked > .4])

  pl.hist(fracmasked, bins=bins, label=band + '-{}'.format(masked), alpha=0.5)

pl.axvline(.4, c='k')
pl.xlim(0., 1.1)
pl.xlabel('FRACMASKED')
pl.yscale('log')
pl.legend(frameon=False)
pl.title(title, fontsize=9) 
pl.savefig('plots/fracmasked.png')
pl.clf()

##
bins = np.arange(0., 2., 0.025)

for band in ['G', 'R', 'Z']:                                                                                                                                                                                        
  fracin = assigned['FRACIN_{}'.format(band)].quantity                                                                                                                                                        
  masked = len(fracin[fracin < .3])
  
  pl.hist(fracin, bins=bins, label=band + '-{}'.format(masked), alpha=0.5)                                                                                                                                         

pl.axvline(.4, c='k')
pl.xlim(0., 1.1)
pl.xlabel('FRACIN') 
pl.yscale('log')                                                                                                                                                                                                    
pl.legend(frameon=False)                                                                                                                                                                                            
pl.savefig('plots/fracin.png')
pl.clf()
