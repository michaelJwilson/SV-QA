import os
import sys
import copy
import glob
import random
import contextlib
import matplotlib;  matplotlib.use('pdf')

import numpy                         as       np
import pylab                         as       pl
import healpy                        as       hp 
import matplotlib.pyplot             as       plt
import astropy.io.fits               as       fits
import astropy.units                 as       u

from   astropy.table                 import   Table, join, hstack, Column
from   astropy                       import   constants as const
from   desitarget.geomask            import   circles
from   in_des                        import   in_des
from   phot_sys                      import   set_photsys
from   desitarget.sv1.sv1_targetmask import   desi_mask, bgs_mask


try:
    from urllib.parse import urlencode

except ImportError:
    from urllib import urlencode

try:
    from urllib.request import urlopen

except ImportError:
    from urllib2 import urlopen

plt.style.use('dark_background')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{hyperref}')

##  plt.rcParams['pgf.preamble'] = [r'\usepackage{hyperref}']

def make_tiny(url):
  request_url = ('http://tinyurl.com/api-create.php?' + 
  urlencode({'url':url}))

  with contextlib.closing(urlopen(request_url)) as response:
    return response.read().decode('utf-8').split('//')[1].replace('tinyurl.com/', '')

    
##
nrow       = 4
ncol       = 5
nside      = 1028

rasterized = False

print('\n\nWelcome.')

##  
fig, axarr = plt.subplots(nrows=nrow, ncols=ncol, figsize=(35, 35))

scratch    = os.environ['CSCRATCH']
_file      = fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')
tiles      = Table(_file[1].data)['RA', 'DEC', 'TILEID', 'AIRMASS', 'EBV_MED', 'STAR_DENSITY', 'IMAGEFRAC_R', 'IMAGEFRAC_GRZ', 'BRIGHTVTMAG', 'BRIGHTRA', 'BRIGHTDEC']

##  
rand_files = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/randoms/*.fits')

##  
names      = ['TILEID', 'PHOTSYS', 'NOBS_G2R', 'BRIGHT', 'SATUR_R', 'ALLMASK_R', 'BAILOUT', 'MEDIUM', 'GALAXY', 'CLUSTER', 'GALDEPTH_R',\
              'PSFSIZE_R', 'TARDEN_MAX', 'TARDEN_STD', 'TARDEN_MAXTURL', 'LOQDEN_MAX', 'LOQDEN_STD', 'LOQDEN_MAXTURL']

data       = {}

for name in names:
  data[name] = []

print('\n\nResolving target fluctuations on an angular scale of {} arcmin.\n\n'.format(hp.pixelfunc.nside2resol(nside, arcmin=True)))
  
for i, _file in enumerate(rand_files):
  tileid     = _file.split('_')[-1].split('.fits')[0]

  randoms    = Table(fits.open(r'/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/randoms/randoms_{}.fits'.format(tileid))[1].data)
  set_photsys(randoms, verbose=False)

  randoms                =  randoms['RA', 'DEC', 'PHOTSYS', 'NOBS_R', 'MASKBITS', 'GALDEPTH_R', 'PSFSIZE_R']
  randoms['NOBS_G2R']    =  randoms['NOBS_R'] >= 2

  del randoms['NOBS_R']

  randoms['BRIGHT']      = (randoms['MASKBITS'] & 2** 1) > 0  
  randoms['SATUR_R']     = (randoms['MASKBITS'] & 2** 3) > 0
  randoms['ALLMASK_R']   = (randoms['MASKBITS'] & 2** 6) > 0
  randoms['BAILOUT']     = (randoms['MASKBITS'] & 2**10) > 0
  randoms['MEDIUM']      = (randoms['MASKBITS'] & 2**11) > 0
  randoms['GALAXY']      = (randoms['MASKBITS'] & 2**12) > 0
  randoms['CLUSTER']     = (randoms['MASKBITS'] & 2**13) > 0

  assert  not np.any((randoms['MASKBITS'] & 2**10) > 0)
  assert  not np.any((randoms['MASKBITS'] & 2**13) > 0)
  
  data['TILEID'].append(np.int(tileid))
  data['PHOTSYS'].append('/'.join(np.unique(randoms['PHOTSYS'])))

  data['GALDEPTH_R'].append(-2.5 * (np.log10(5. / np.sqrt(np.median(np.array(randoms['GALDEPTH_R'])))) -9.))
  data['PSFSIZE_R'].append(np.median(np.array(randoms['PSFSIZE_R'])))

  for name in names:
    if not name in ['TILEID', 'PHOTSYS', 'TARDEN_MAX', 'TARDEN_STD', 'TARDEN_MAXTURL', 'LOQDEN_MAX', 'LOQDEN_STD', 'LOQDEN_MAXTURL', 'GALDEPTH_R', 'PSFSIZE_R', 'ASSIGN_TURL', 'TARGETS_TURL']:
      data[name].append(100.0 * np.count_nonzero(randoms[name] == True) / len(randoms))

  ##  Assigned truth with legacy columns.                                                                                                                                                                                        
  _file       = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/mtls/svmtl_{:06d}.fits'.format(np.int(tileid))
  dat         = Table(fits.open(_file)[1].data)  
  dat         = dat[(dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('BGS_ANY')) != 0]

  ##  Do not include low quality. 
  dat         = dat[dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_LOWQ') == 0]
  
  ## 
  pix         = hp.ang2pix(nside, (90. - dat['DEC'].quantity.value) * np.pi / 180., dat['RA'].quantity.value * np.pi / 180., nest=False)
  upix, cnts  = np.unique(pix, return_counts=True)

  maxpix      = upix[np.argmax(cnts)]
  
  theta, phi  = hp.pix2ang(nside, maxpix, nest=False)
  mra, mdec   = 180. / np.pi * phi, 90. -180. / np.pi * theta
  
  link        = 'http://legacysurvey.org/viewer/viewer/?ra={:.4f}&dec={:.4f}&layer=dr8&zoom=12'.format(mra, mdec)
  
  data['TARDEN_MAX'].append(np.sort(cnts)[-3:])
  data['TARDEN_STD'].append(np.std(cnts))
  data['TARDEN_MAXTURL'].append(make_tiny(link))

  ##  And again for low quality. 
  _file       = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/mtls/svmtl_{:06d}.fits'.format(np.int(tileid))
  dat         = Table(fits.open(_file)[1].data)
  dat         = dat[(dat['SV1_DESI_TARGET'] & 2 ** desi_mask.bitnum('BGS_ANY')) != 0]
  dat         = dat[dat['SV1_BGS_TARGET'] & 2 ** bgs_mask.bitnum('BGS_LOWQ') != 0]

  pix         = hp.ang2pix(nside, (90. - dat['DEC'].quantity.value) * np.pi / 180., dat['RA'].quantity.value * np.pi / 180., nest=False)
  upix, cnts  = np.unique(pix, return_counts=True)

  maxpix      = upix[np.argmax(cnts)]

  theta, phi  = hp.pix2ang(nside, maxpix, nest=False)
  mra, mdec   = 180. / np.pi * phi, 90. -180. / np.pi * theta

  link        = 'http://legacysurvey.org/viewer/viewer/?ra={:.4f}&dec={:.4f}&layer=dr8&zoom=12'.format(mra, mdec)

  data['LOQDEN_MAX'].append(np.sort(cnts)[-3:])
  data['LOQDEN_STD'].append(np.std(cnts))
  data['LOQDEN_MAXTURL'].append(make_tiny(link))

##  
tile_fails = Table([data[x] for x in data.keys()], names=names)
tile_fails = join(tile_fails, tiles, join_type='left', keys='TILEID')

tile_link    = [make_tiny('http://legacysurvey.org/viewer/viewer/?ra={:.4f}&dec={:.4f}&layer=dr8&zoom=12'.format(x['RA'], x['DEC'])) for x in tile_fails]
assign_link  = [make_tiny('https://portal.nersc.gov/project/desi/users/mjwilson/SV-ASSIGN/tile-{:06d}.fits'.format(np.int(x['TILEID']))) for x in tile_fails]
targets_link = [make_tiny('https://portal.nersc.gov/project/desi/users/mjwilson/SV-TARGETS/tile-targets-{:06d}.fits'.format(np.int(x['TILEID']))) for x in tile_fails]

tile_fails['TILE_TURL']    = Column(data=tile_link,    name='TILE_TURL',    dtype='S32', length=len(tile_fails))
tile_fails['ASSIGN_TURL']  = Column(data=assign_link,  name='ASSIGN_TURL',  dtype='S32', length=len(tile_fails))
tile_fails['TARGETS_TURL'] = Column(data=targets_link, name='TARGETS_TURL', dtype='S32', length=len(tile_fails))

##  Lower the precision for pprint.
for col in ['NOBS_G2R', 'GALDEPTH_R', 'BRIGHT', 'SATUR_R', 'ALLMASK_R', 'BAILOUT', 'MEDIUM', 'GALAXY', 'CLUSTER',\
            'TARDEN_STD', 'LOQDEN_STD', 'EBV_MED', 'STAR_DENSITY', 'IMAGEFRAC_R', 'IMAGEFRAC_GRZ', 'BRIGHTVTMAG',\
            'PSFSIZE_R']:

  tile_fails[col]   = Column(data=np.array(tile_fails[col]), name=col, dtype=np.float16, length=len(tile_fails[col]))

tile_fails['NGAMA'] = Column(data=['--'] * len(tile_fails), name='NGAMA', dtype='S16', length=len(tile_fails))
tile_fails['NAGES'] = Column(data=['--'] * len(tile_fails), name='NAGES', dtype='S16', length=len(tile_fails))
tile_fails['NSDSS'] = Column(data=['--'] * len(tile_fails), name='NSDSS', dtype='S16', length=len(tile_fails))

##  Append.
for i, _file in enumerate(rand_files):
  tileid   = np.int(_file.split('_')[-1].split('.fits')[0])

  ##  GAMA
  try:
    dat    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/GAMA-south_{:06d}.fits'.format(tileid))[1].data)
    tile_fails['NGAMA'][tile_fails['TILEID'] == tileid] = '{:d}'.format(len(dat))
    
  except:
    tile_fails['NGAMA'][tile_fails['TILEID'] == tileid] = '{:d}'.format(0)
      
  ##  AGES
  try:
    dat   = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/ages_reduced-north_{:06d}.fits'.format(tileid))[1].data)
    ngal  = len(dat)

  except:
    ngal  = 0  

  try:
    dat   = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/ages_reduced-south_{:06d}.fits'.format(tileid))[1].data)
    ngal += len(dat)

  except:
    pass
    
  tile_fails['NAGES'][tile_fails['TILEID'] == tileid] = '{:d}'.format(ngal)

  ##  SDSS                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  try:
    dat    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/sdss-north_{:06d}.fits'.format(tileid))[1].data)
    ngal   = len(dat)

  except:
    ngal   = 0

  try:
    dat    = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/assigned/bytile/sdss-south_{:06d}.fits'.format(tileid))[1].data)
    ngal  += len(dat)

  except:
    pass

  tile_fails['NSDSS'][tile_fails['TILEID'] == tileid] = '{:d}'.format(ngal)

tile_fails.sort('TILEID')
  
tile_fails.pprint(max_lines=-1, max_width=-1)

tile_fails.write(os.environ['CSCRATCH'] + '/BGS/SV-ASSIGN/tiles/tile-quality.fits', format='fits', overwrite=True)

print('\n\nDone.\n\n')
