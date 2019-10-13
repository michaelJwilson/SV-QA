import  os
import  glob
import  pylab               as      pl
import  numpy               as      np 
import  astropy.io.fits     as      fits
 
from    astropy.table       import  Table, join, Column 
from    desitarget.targets  import  encode_targetid
from    phot_sys            import  set_photsys
from    make_tiny           import  make_tiny


print('\n\nWelcome.\n\n')

scratch       = os.environ['CSCRATCH']

##
hsc           = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_north.fits')[1].data)
hsc           = hsc[hsc['DEC'] > -30.]
hsc           = set_photsys(hsc)
hsc           = hsc[hsc['PHOTSYS'] == 'BMZLS']

cols          = hsc.columns
cols          = ['TARGETID', 'RA','DEC', 'DCHISQ', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'NOBS_G','NOBS_R','NOBS_Z', 'SHAPEEXP_R', 'TYPE','mizuki_photoz_best','FRANKENZ','PHOTSYS', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']

hsc           = hsc[cols]
hsc           = hsc[hsc['TYPE'] == 'REX']
hsc['MIZUKI'] = hsc['mizuki_photoz_best']
hsc           = hsc[hsc['SHAPEEXP_R'] > 2.5]

urls          = ['http://viewer.legacysurvey.org/?ra={:.4f}&dec={:.4f}&zoom=16&layer=dr8'.format(i['RA'], i['DEC']) for i in hsc]
urls          = [make_tiny(url) for url in urls]
hsc['TURL']   =  urls

del  hsc['mizuki_photoz_best']

hsc.sort('SHAPEEXP_R')
hsc.pprint(max_lines=50, max_width=-1)

pl.hist(hsc['SHAPEEXP_R'], bins=100)
pl.yscale('log')
pl.show()

'''
hsc              = hsc[hsc['SHAPEEXP_R'] > 9.0]

toview           = hsc['RA', 'DEC']
toview['radius'] = 1.6

toview.write('./toview.fits', format='fits')
'''
print('\n\n')
'''
## 
elgs                      = fits.open(scratch + '/BGS/SV-ASSIGN/elgs/elgs.fits')[1].data
elgs                      = Table(elgs) 
elgs['PHOTSYS']           = Column(data=['N'] * len(elgs), name='PHOTSYS')

##
in_south                  = elgs['DEC'] < 32.375
elgs['PHOTSYS'][in_south] = 'S'
elgs                      = set_photsys(hsc)

elgs.pprint()

##
ecols                     = list(elgs.columns)
hcols                     = list(hsc.columns)
tojoin                    = set(ecols) ^ set(hcols)
tojoin                    = list(tojoin)
'''
##  elgs                  = join(elgs, hsc[tojoin], keys='TARGETID', join_type='inner')
##  elgs.pprint()

print('\n\n')

##
radii    = ['_0p5', '_0p75', '_1p0', '_1p5', '_2p0', '_3p5', '_5p0', '_7p0']

apf_namg = ['apflux_resid_g'.upper() + x for x in radii]
apf_namr = ['apflux_resid_r'.upper() + x for x in radii]
apf_namz = ['apflux_resid_z'.upper() + x for x in radii]

apv_namg = ['apflux_ivar_g'.upper() + x for x in radii]
apv_namr = ['apflux_ivar_r'.upper() + x for x in radii]
apv_namz = ['apflux_ivar_z'.upper() + x for x in radii]

names    = ['TARGETID', 'BRICKID', 'OBJID'] + apf_namg + apv_namg + apf_namr + apv_namr + apf_namz + apv_namz + ['PSF_FLUXG', 'PSF_FLUXG_IVAR', 'PSF_FLUXR', 'PSF_FLUXR_IVAR',\
            'PSF_FLUXZ', 'PSF_FLUXZ_IVAR', 'REX_FLUXG', 'REX_FLUXG_IVAR', 'REX_FLUXR', 'REX_FLUXR_IVAR', 'REX_FLUXZ', 'REX_FLUXZ_IVAR']

tractors                             =  glob.glob(scratch + '/BGS/SV-ASSIGN/new_tractors/north/*/*.fits')

##
for tractor in tractors:
  tractor                            =  Table(fits.open(tractor)[0].data, names=names)
  columns                            =  [x.upper() for x in tractor.columns]
  
  ##  arcsec.
  tractor['APFLUX_RADII']            =  np.array([[0.5, 0.75, 1.0, 1.5, 2.0, 3.5, 5.0, 7.0]] * len(tractor))

  ##  [NANOMAGGIES PER SQ. ARCSEC].
  tractor['APFLUX_RESID_7M5_G']      = (tractor['APFLUX_RESID_G_7p0'] - tractor['APFLUX_RESID_G_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  tractor['APFLUX_RESID_7M5_R']      = (tractor['APFLUX_RESID_R_7p0'] - tractor['APFLUX_RESID_R_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  tractor['APFLUX_RESID_7M5_Z']      = (tractor['APFLUX_RESID_Z_7p0'] - tractor['APFLUX_RESID_Z_5p0']) / (np.pi * (tractor['APFLUX_RADII'][:,7]**2. - tractor['APFLUX_RADII'][:,6]**2.))
  '''
  tractor['APFLUX_RESID_7M5_G_IVAR'] =  tractor['APFLUX_RESID_G_7p0_IVAR'] + tractor['APFLUX_RESID_G_5p0_IVAR'] 
  tractor['APFLUX_RESID_7M5_R_IVAR'] =  tractor['APFLUX_RESID_R_7p0_IVAR'] + tractor['APFLUX_RESID_R_5p0_IVAR']
  tractor['APFLUX_RESID_7M5_Z_IVAR'] =  tractor['APFLUX_RESID_Z_7p0_IVAR'] + tractor['APFLUX_RESID_Z_5p0_IVAR']
  '''

  tractor['PSF_GSNR']                =  tractor['PSF_FLUXG'] * np.sqrt(tractor['PSF_FLUXG_IVAR'])
  tractor['PSF_RSNR']                =  tractor['PSF_FLUXR'] * np.sqrt(tractor['PSF_FLUXR_IVAR'])
  tractor['PSF_ZSNR']                =  tractor['PSF_FLUXZ'] * np.sqrt(tractor['PSF_FLUXZ_IVAR'])

  tractor['REX_GSNR']                =  tractor['REX_FLUXG'] * np.sqrt(tractor['REX_FLUXG_IVAR'])
  tractor['REX_RSNR']                =  tractor['REX_FLUXR'] * np.sqrt(tractor['REX_FLUXR_IVAR'])
  tractor['REX_ZSNR']                =  tractor['REX_FLUXZ'] * np.sqrt(tractor['REX_FLUXZ_IVAR'])

  del  tractor['APFLUX_RADII']

  tractor = tractor['TARGETID', 'PSF_GSNR', 'PSF_RSNR', 'PSF_ZSNR', 'REX_GSNR', 'REX_RSNR', 'REX_ZSNR', 'APFLUX_RESID_7M5_G', 'APFLUX_RESID_7M5_R', 'APFLUX_RESID_7M5_Z']
  tractor.sort('PSF_GSNR')
  
  tractor.pprint()

  exit(1)

  matches                            =  np.in1d(tractor['TARGETID'], hsc['TARGETID'])
  
  if np.any(matches):
    print('\n\nSuccess!\n\n')

    matched                          =  join(elgs, tractor, keys=['TARGETID'], join_type='inner')
    matched.sort('SNR')
    matched.pprint(max_width=-1)

    print('\n\nDone.\n\n')
    
    exit(0)

  else:
    print('No matches.')

print('\n\nDone.\n\n')
