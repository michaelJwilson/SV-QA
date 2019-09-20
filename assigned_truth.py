import glob
import random
import numpy             as      np
import pylab             as      pl
import matplotlib.pyplot as      plt
import astropy.io.fits   as      fits
import astropy.units     as      u

from   astropy.table      import  Table, join, hstack
from   astropy            import  constants as const
from   desitarget.geomask import circles


plt.figure(figsize=(7, 10))

compute    = True

ax1        = plt.subplot(4, 1, 1)
ax2        = plt.subplot(4, 1, 2)
ax3        = plt.subplot(4, 1, 3)
ax4        = plt.subplot(4, 1, 4)

tiles      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)

_assigned  = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/assigned_targets.fits')[1].data)
_assigned.sort('BRICK_OBJID')

if compute:
  files    = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/standard/*.fits')

else:
  files    = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/matched/*.fits')
    
for _file in files:
  survey       = _file.split('/')[-1].split('.')[0].split('-standard')[0]

  print('Solving for {}'.format(survey))
  
  if compute:
    ##  Spectroscopic.                                                                                                                                                                                                               
    dat          = Table(fits.open(_file)[1].data, meta={'name': survey})['Z', 'ZERR', 'ZWARN']  ##  ['RA', 'DEC']

    legacy       = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/legacy/ls-' + survey + '.fits'
    legacy       = Table(fits.open(legacy)[1].data)

    ##  Line matched.
    legacy       = hstack([legacy, dat])
    legacy['BRICK_OBJID'] = legacy['OBJID']

    ##  
    legacy.sort('BRICK_OBJID')

    ##  print(legacy.columns)
    
    matched      = join(_assigned, legacy, keys=['BRICKID', 'BRICK_OBJID'])
  
    ##  print(len(targetids), np.sum(isin))
     
    ##  print('\n\n{}\n'.format(survey))
    ##  print(legacy)
    print(matched)

    if len(matched) == 0:
      continue

    matched.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/' + survey + '.fits', format='fits', overwrite=True)

  else:
    matched    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/' + survey +'.fits')[1].data)
    
  label        = survey + '-{:d}'.format(len(matched))
  
  ##  Randomise the order (to prevent correlation with e.g. ra).                                                                                                                                                                     
  sample       = random.sample(range(0, len(matched)), len(matched))
  matched      = matched[sample]

  ##  Subsample.                                                                                                                                                                                                                     
  matched      = matched[::2]
  
  if survey in ['sdss-north', 'sdss-south']:
    ax1.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label, alpha=0.5)

  elif survey in ['2dFGRS-north', '2dFGRS-south', '6dFGS-north', '6dFGS-south', 'ages_reduced-north', 'ages_reduced-south']:
    ax2.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label, alpha=0.5)

  elif survey in ['GAMA-north', 'GAMA-south']:
    ax3.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label, alpha=0.5)
    
  else:
    ax3.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label, alpha=0.5)

##  Last plot show SV tiles.                                                                                                                                                                                                         
##  kwargs = {'ec': 'k', 'lw': 5, 'alpha': 0.5}
##  circles(tiles['RA'].quantity, tiles['DEC'].quantity, s=1.6, c='k', **kwargs)
ax4.plot(_assigned['TARGET_RA'], _assigned['TARGET_DEC'], c='gold', marker='x', lw=0, markersize=1, label='DESI', alpha=0.3)
    
##
for ax in [ax1, ax2, ax3, ax4]:  
  ax.axhline(y=0.,     c='k', alpha=0.5)
  ax.axhline(y=32.375, c='k', alpha=0.5)
  
  ax.set_xlim(360.,   0.)
  ax.set_ylim(-40.,  80.)

  ax.legend(frameon=False, fontsize=10, ncol=2, loc=1)
  
plt.tight_layout()
  
pl.savefig('spec_truth.png')

print('\n\nDone.\n\n')

