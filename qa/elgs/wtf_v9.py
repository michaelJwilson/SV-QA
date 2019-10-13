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

scratch                              =  os.environ['CSCRATCH']

cols                                 = ['RA', 'DEC', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'SHAPEEXP_R',\
                                        'TYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'BRICK_PRIMARY', 'BRIGHTBLOB',\
                                        'APFLUX_RESID_G', 'APFLUX_RESID_R', 'APFLUX_RESID_Z', 'OBJID', 'BRICKID',\
                                        'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z']

##  /Users/MJWilson/Work/desi/scratch/BGS/SV-ASSIGN/skybackground/tractor/
tractors                             =  glob.glob(scratch + '/BGS/SV-ASSIGN/skybackground/tractor/*.fits')

##
for i, _tractor in enumerate(tractors): 
  tractor                            =  Table(fits.open(_tractor)[1].data)
  columns                            =  [x.upper().strip() for x in tractor.columns]
  tractor                            =  Table(tractor, names=columns)
  tractor                            =  tractor[cols]
  
  _len                               =  len(tractor)

  del  tractor['BRICKID']

  ## 
  _metrics                           =  os.environ['CSCRATCH'] + '/BGS/SV-ASSIGN/skybackground/metrics/all-models-' + _tractor.split('/')[-1].split('.fits')[0].split('-')[1] + '.fits'
  metrics                            =  Table(fits.open(_metrics)[1].data)
  columns                            =  [x.upper().strip() for x in metrics.columns]
  metrics                            =  Table(metrics, names=columns)
  metrics                            =  metrics['PSF_FLUX', 'REX_FLUX', 'OBJID', 'BRICKID']

  print('\n\nLoading {} and {}.\n\n'.format(_tractor, _metrics))

  tractor.pprint()

  print('\n\n')
  
  metrics.pprint()
  
  ##
  tractor                            =  join(tractor, metrics, keys='OBJID')

  tractor                            =  tractor[tractor['TYPE'] != 'PSF ']
  tractor['PSF_MIN']                 =  np.c_[tractor['PSFSIZE_G'], tractor['PSFSIZE_R'], tractor['PSFSIZE_Z']].min(axis=-1)

  tractor['SNR2']                    =  np.sqrt(tractor['FLUX_IVAR_G']) * tractor['FLUX_G'] **2. + np.sqrt(tractor['FLUX_IVAR_R']) * tractor['FLUX_R'] **2. + np.sqrt(tractor['FLUX_IVAR_Z']) * tractor['FLUX_Z'] **2.

  ##                                                                                                                                                                                                                                                                                     
  tractor                            =  tractor[tractor['SHAPEEXP_R'] < 0.4]

  urls                               =  [make_tiny('http://legacysurvey.org/viewer?ra={:.4f}&dec={:.4f}&layer=dr8&zoom=16'.format(x['RA'], x['DEC'])) for x in tractor]
  tractor['TURL']                    =  urls

  tractor.sort('SHAPEEXP_R')

  print('\n\n')

  tractor.pprint(max_width=-1, max_lines=50)

  print('\n{}    {}'.format(_len, len(tractor)))

  tractor['name']                    = tractor['TYPE']
  tractor                            = tractor['RA', 'DEC', 'name', 'PSF_MIN'] 
  tractor['radius']                  = Column(data=tractor['PSF_MIN'], name='radius')
  
  tractor.write('toview.fits', format='fits', overwrite=True)
  
  break

print('\nDone.\n\n')
