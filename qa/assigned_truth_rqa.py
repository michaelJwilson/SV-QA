import glob
import random
import numpy             as      np
import pylab             as      pl
import matplotlib.pyplot as      plt
import astropy.io.fits   as      fits
import astropy.units     as      u

from   astropy.table      import  Table, join
from   astropy            import  constants as const
from   desitarget.geomask import circles
from   matplotlib.ticker  import FormatStrFormatter


files       = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/*.fits')
    
nsurvey     = len(files)
nrow        = np.int(np.ceil(nsurvey / 2))

fig, axarr  = plt.subplots(nrow, 2, sharex=True)
fig.set_size_inches((7, 10))

colors      = plt.rcParams['axes.prop_cycle'].by_key()['color']
uspec       = []
    
for i, _file in enumerate(files):
  survey    = _file.split('/')[-1].split('.fits')[0].split('-standard')[0]
  matched   = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/' + survey +'.fits')[1].data)
  
  targetids = matched['TARGETID'].quantity.value.tolist()
  uspec    += targetids

  rs        = matched['FLUX_R'].quantity / matched['MW_TRANSMISSION_R'].quantity
  rs        = 22.5 - 2.5 * np.log10(rs)
  rs        = np.sort(rs)
  rs        = rs[np.isfinite(rs)]

  ##  print(survey, len(rs))
  
  row       = i % nrow
  col       = i % 2

  label     = survey
  label     = label.replace('north', 'N').replace('south', 'S')
  label     = label.replace('hsc_pdr1', 'HSC1').replace('.forced', '').replace('.reduced', '').replace('_', '-').replace('udeep', 'UD').replace('deep', 'D').replace('wide', 'W')
  label    += ' ({})'.format(len(rs))
  
  axarr[row][col].hist(rs, bins=100, label=label.upper(), alpha=0.3, color=colors[i % len(colors)])

  ##  Note problem with AGES. 
  axarr[row][col].set_xlim(19.,  21.)
  axarr[row][col].set_ylim(10., 5.e4)

  axarr[row][col].set_yscale('log')

  axarr[row][col].legend(frameon=False, loc=1, ncol=2, fontsize=9)

##
uspec = set(uspec)

print('Unique spectra: {}'.format(len(uspec)))

axarr[0][0].set_title('Unique  spectra:  {}'.format(len(uspec)))

axarr[-1][0].set_xlabel(r'$r_{AB}$')
axarr[-1][1].set_xlabel(r'$r_{AB}$')

axarr[-1][0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axarr[-1][1].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

plt.tight_layout()
pl.show()
##  pl.savefig('plots/assigned_truth_rhist.png')
