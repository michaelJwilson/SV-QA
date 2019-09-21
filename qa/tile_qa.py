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
from    desimodel.io           import  load_tiles

##
scratch    = os.environ['CSCRATCH']
root       = scratch + '/svdc2019c2/survey/'

##                                                                                                                                                                                                                                      
atfile     = '/global/cscratch1/sd/mjwilson/svdc2019c2/survey/basetiles/original/schlafly-tiles.fits'
alltiles   = Table(fits.open(atfile)[1].data) 
print(alltiles.columns)

tiles      = Table(fits.open(scratch + '/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)

for x, label in zip(['AIRMASS', 'STAR_DENSITY', 'EBV_MED'], ['Airmass at 15. deg. from zenith', 'Median # of G < 19.5 GAIA stars / sq. deg.', 'EBV']):
  sample = tiles[x].quantity                                                                                                                                                                               
  pl.hist(sample, bins=50, alpha=0.5)  

  pl.xlabel(label)
  pl.yscale('linear')                                                                                                                                                                                                                       
  pl.legend(frameon=False)

  mn   = np.median(alltiles[x])
  std  = alltiles[x].std()
  
  ax   = pl.gca()
  uax  = ax.twiny()  
  axt  = ax.get_xticks()

  uaxt = [(np.float(x) - mn) / std for x in axt]
  uaxt = ['{:.1f}'.format(x) for x in uaxt]
    
  uax.set_xticks(axt)
  uax.set_xbound(ax.get_xbound())
  uax.set_xticklabels(uaxt)
  uax.set_xlabel(r'$\Delta$ {} / (IN DESI) STD. DEV.'.format(x))
  
  ##  pl.show()
  pl.savefig('plots/' + x + '.png')
  pl.clf()

for x in ['R', 'GRZ']:
  sample = tiles['IMAGEFRAC_{}'.format(x)].quantity
  sample = 100. * (1. - sample)
  
  pl.hist(sample, bins=50, alpha=0.5)
  pl.xlabel('Fraction of tile within 1.605 deg without {} imaging [%]'.format(x))  
  pl.legend(frameon=False)  
  '''
  mn   = np.median(alltiles['IMAGEFRAC_{}'.format(x)])
  std  = alltiles[x].std()

  ax   = pl.gca()
  uax  = ax.twiny()
  axt  = ax.get_xticks()

  uaxt = [(np.float(x) - mn) / std for x in axt]
  uaxt = ['{:.1f}'.format(x) for x in uaxt]

  uax.set_xticks(axt)
  uax.set_xbound(ax.get_xbound())
  uax.set_xticklabels(uaxt)
  uax.set_xlabel(r'$\Delta$ {} / (IN DESI) STD. DEV.'.format(x))
  '''
  ##  pl.show()
  pl.savefig('plots/' + x + '.png')
  pl.clf()
