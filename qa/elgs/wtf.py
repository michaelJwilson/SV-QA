import  os
import  glob
import  numpy               as      np 
import  astropy.io.fits     as      fits
 
from    astropy.table       import  Table, join
from    desitarget.targets  import  encode_targetid


print('Welcome.')

scratch = os.environ['CSCRATCH']

##
dat     = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_north_lite.fits')[1].data)
##  dat     = dat[dat['PHOTSYS'] == 'DES']

##  dat.sort('GSNR')

dat.pprint(max_lines=250)

print('\n\n')

## 
##  scratch = os.environ['CSCRATCH']
##  elgs    = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_north_lite.fits')[1].data)
##  dat.pprint()
##  print('\n\n')

##  /north/
radii    = ['_0p5', '_0p75', '_1p0', '_1p5', '_2p0', '_3p5', '_5p0', '_7p0']
apf_namg = ['apflux_resid_g'.upper() + x for x in radii]
apf_namr = ['apflux_resid_r'.upper() + x for x in radii]
apf_namz = ['apflux_resid_z'.upper() + x for x in radii]

names    = ['TARGETID', 'BRICKID', 'OBJID'] + apf_namg + apf_namr + apf_namz + ['PSF_FLUXG', 'PSF_FLUXG_IVAR', 'PSF_FLUXR', 'PSF_FLUXR_IVAR',\
            'PSF_FLUXZ', 'PSF_FLUXZ_IVAR', 'REX_FLUXG', 'REX_FLUXG_IVAR', 'REX_FLUXR', 'REX_FLUXR_IVAR', 'REX_FLUXZ', 'REX_FLUXZ_IVAR']

tractors                        =  glob.glob(scratch + '/BGS/SV-ASSIGN/new_tractors/*/*.fits')

for tractor in tractors:
  try:                                                                                                                                                                                                                                 
    tractor                     =  Table(fits.open(tractor)[0].data, names=names)
    columns                     =  [x.upper() for x in tractor.columns]
  
  except:
    print('Ignoring failure on {}.'.format(tractor))
    continue

  ##  arcsec.
  tractor['APFLUX_RADII']       =  np.array([[0.5, 0.75, 1.0, 1.5, 2.0, 3.5, 5.0, 7.0]] * len(tractor))

  ##  [NANOMAGGIES PER SQ. ARCSEC].
  tractor['APFLUX_RESID_7M5_G'] = (tractor['APFLUX_RESID_G_7p0'] - tractor['APFLUX_RESID_G_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  tractor['APFLUX_RESID_7M5_R'] = (tractor['APFLUX_RESID_R_7p0'] - tractor['APFLUX_RESID_R_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  tractor['APFLUX_RESID_7M5_Z'] = (tractor['APFLUX_RESID_Z_7p0'] - tractor['APFLUX_RESID_Z_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))

  tractor.sort('APFLUX_RESID_7M5_G')

  del  tractor['APFLUX_RADII']
  
  tractor.pprint()

  matched                       = join(dat, tractor, keys=['TARGETID'], join_type='inner')

  if len(matched) > 0:
    matched.pprint()

    print('Success!')

    exit(0)

print('\n\nDone.\n\n')
