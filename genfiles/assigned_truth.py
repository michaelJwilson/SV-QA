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

compute    = False

ax1        = plt.subplot(4, 1, 1)
ax2        = plt.subplot(4, 1, 2)
ax3        = plt.subplot(4, 1, 3)
ax4        = plt.subplot(4, 1, 4)

tiles      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)

_assigned  = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/assigned_targets.fits')[1].data)
_assigned.sort('BRICK_OBJID')

ax4.plot(_assigned['TARGET_RA'][::10], _assigned['TARGET_DEC'][::10], c='gold', marker='.', lw=0, markersize=1, label='DESI', alpha=0.1, rasterized=False, zorder=2)

if compute:
  files    = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/standard/*.fits')

else:
  files    = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/*.fits')
  
for _file in files:
  survey   = _file.split('/')[-1].split('.')[0].split('-standard')[0]

  if survey in ['hsc_pdr1_wide', 'hsc_pdr1_deep', 'hsc_pdr1_udeep']:
    survey = _file.split('/')[-1]
    survey = survey.split('.fits')[0]
        
  print('Solving for {}'.format(survey))
  
  if compute:
    ##  Spectroscopic.
    dat          = Table(fits.open(_file)[1].data, meta={'name': survey})['Z', 'ZERR', 'ZWARN']  ##  ['RA', 'DEC']

    legacy       = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/legacy/ls-' + survey.replace('-standard', '') + '.fits'
    legacy       = Table(fits.open(legacy)[1].data)

    ##  Line matched.
    legacy                = hstack([legacy, dat])
    legacy['BRICK_OBJID'] = legacy['OBJID']

    ##  
    ##  legacy.sort('BRICK_OBJID')

    ##  print(legacy.columns)
    
    matched      = join(_assigned, legacy, keys=['BRICKID', 'BRICK_OBJID'], join_type='inner')
  
    ##  print(len(targetids), np.sum(isin))
     
    if len(matched) == 0:
      continue

    print(matched)
    
    matched.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/' + survey.replace('-standard', '') + '.fits', format='fits', overwrite=True)

  else:
    matched    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/' + survey +'.fits')[1].data)
    
  label        = survey + '-{:d}'.format(len(matched))
  label        = label.replace('north', 'N').replace('south', 'S')
  
  ##  Randomise the order (to prevent correlation with e.g. ra).                                                                                                                                                                     
  sample       = random.sample(range(0, len(matched)), len(matched))
  matched      = matched[sample]

  ##  Subsample.                                                                                                                                                                                                                     
  matched      = matched[::2]
  
  if survey in ['sdss-north', 'sdss-south']:
    ax1.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label.upper(), alpha=0.5)

  elif survey in ['2dFGRS-north', '2dFGRS-south', '6dFGS-north', '6dFGS-south', 'ages_reduced-north', 'ages_reduced-south']:
    ax2.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label.upper(), alpha=0.5)

  elif survey in ['GAMA-north', 'GAMA-south', 'primus-north', 'primus-south']:
    ax1.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label.upper(), alpha=0.5)

  elif survey in ['hsc_pdr1_udeep.forced.reduced-north', 'hsc_pdr1_udeep.forced.reduced-south', 'hsc_pdr1_deep.forced.reduced-north', 'hsc_pdr1_deep.forced.reduced-south',\
                  'hsc_pdr1_wide.forced.reduced-north',  'hsc_pdr1_wide.forced.reduced-south']:

    label = label.replace('hsc_pdr1', 'HSC1').replace('.forced', '').replace('.reduced', '').replace('_', '-').replace('udeep', 'UD').replace('deep', 'D').replace('wide', 'W')
    
    ax4.scatter(matched['RA'], matched['DEC'], s=10, label=label.upper(), alpha=1., zorder=1)
    
  else:    
    ax3.scatter(matched['RA'], matched['DEC'], s=5, rasterized=True, label=label.upper(), alpha=0.5)
    
##
for ax, ncol in zip([ax1, ax2, ax3, ax4], [3, 2, 2, 3]):  
  ax.axhline(y=0.,     c='k', alpha=0.5)
  ax.axhline(y=32.375, c='k', alpha=0.5)
  
  ax.set_xlim(360.,   0.)
  ax.set_ylim(-40.,  80.)

  ax.legend(frameon=False, fontsize=9, ncol=ncol, loc=3, columnspacing=0, labelspacing=0, handletextpad=0)
  
plt.tight_layout()
  
pl.savefig('../plots/spec_truth.png')

print('\n\nDone.\n\n')

