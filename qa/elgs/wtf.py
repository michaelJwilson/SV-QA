import  os
import  glob
import  numpy               as      np 
import  astropy.io.fits     as      fits
 
from    astropy.table       import  Table, join
from    desitarget.targets  import  encode_targetid
from    phot_sys            import  set_photsys


print('\n\nWelcome.\n\n')

scratch = os.environ['CSCRATCH']

##
dat     = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_north_lite.fits')[1].data)
dat     = set_photsys(dat)
##  dat     = dat[dat['PHOTSYS'] == 'DES']

dat.sort('SNR')

dat.pprint(max_lines=50, max_width=-1)

print('\n\n')

## 
##  scratch = os.environ['CSCRATCH']
##  elgs    = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_north_lite.fits')[1].data)
##  dat.pprint()
##  print('\n\n')

##
radii    = ['_0p5', '_0p75', '_1p0', '_1p5', '_2p0', '_3p5', '_5p0', '_7p0']

apf_namg = ['apflux_resid_g'.upper() + x for x in radii]
apf_namr = ['apflux_resid_r'.upper() + x for x in radii]
apf_namz = ['apflux_resid_z'.upper() + x for x in radii]

names    = ['TARGETID', 'BRICKID', 'OBJID'] + apf_namg + apf_namr + apf_namz + ['PSF_FLUXG', 'PSF_FLUXG_IVAR', 'PSF_FLUXR', 'PSF_FLUXR_IVAR',\
            'PSF_FLUXZ', 'PSF_FLUXZ_IVAR', 'REX_FLUXG', 'REX_FLUXG_IVAR', 'REX_FLUXR', 'REX_FLUXR_IVAR', 'REX_FLUXZ', 'REX_FLUXZ_IVAR']

tractors                      =   glob.glob(scratch + '/BGS/SV-ASSIGN/new_tractors/north/249//*.fits')

for tractor in tractors:
  tractor                     =  Table(fits.open(tractor)[0].data, names=names)
  columns                     =  [x.upper() for x in tractor.columns]
  
  ##  arcsec.
  tractor['APFLUX_RADII']       =  np.array([[0.5, 0.75, 1.0, 1.5, 2.0, 3.5, 5.0, 7.0]] * len(tractor))

  ##  [NANOMAGGIES PER SQ. ARCSEC].
  tractor['APFLUX_RESID_7M5_G']      = (tractor['APFLUX_RESID_G_7p0'] - tractor['APFLUX_RESID_G_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  tractor['APFLUX_RESID_7M5_R']      = (tractor['APFLUX_RESID_R_7p0'] - tractor['APFLUX_RESID_R_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  tractor['APFLUX_RESID_7M5_Z']      = (tractor['APFLUX_RESID_Z_7p0'] - tractor['APFLUX_RESID_Z_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  '''
  tractor['APFLUX_RESID_7M5_G_IVAR'] =  tractor['APFLUX_RESID_G_7p0_IVAR'] + tractor['APFLUX_RESID_G_5p0_IVAR'] 
  tractor['APFLUX_RESID_7M5_R_IVAR'] =  tractor['APFLUX_RESID_R_7p0_IVAR'] + tractor['APFLUX_RESID_R_5p0_IVAR']
  tractor['APFLUX_RESID_7M5_Z_IVAR'] =  tractor['APFLUX_RESID_Z_7p0_IVAR'] + tractor['APFLUX_RESID_Z_5p0_IVAR']
  '''
  tractor.sort('PSF_FLUXG')

  del  tractor['APFLUX_RADII']
  
  tractor.pprint(max_width=-1)
  '''
  matched                       = join(dat, tractor, keys=['TARGETID'], join_type='inner')

  if len(matched) > 0:
    matched.pprint(max_width=-1)

    print('Success!')

    exit(0)

  else:
    print('No matches.')
  '''

  exit(1)

print('\n\nDone.\n\n')
