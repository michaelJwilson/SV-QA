import glob
import numpy            as      np
import pylab            as      pl
import astropy.io.fits  as      fits
import astropy.units    as      u

from   astropy.table    import  Table
from   astropy          import  constants as const


##  Only standardise the spectroscopic aspect, LS half is fine.  
files = glob.glob('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/truth/primus*.fits')

for _file in files:
  dat          = Table(fits.open(_file)[1].data)
  survey       = _file.split('/')[-1].split('.fits')[0]

  if survey in ['6dFGS-north', '6dFGS-south']:
    dat['RA']    = dat['RAJ2000']
    dat['DEC']   = dat['DEJ2000']
    dat['Z']     = np.array([x / const.c.to('km/s').value for x in dat['cz']])
    dat['ZERR']  = np.array([x / const.c.to('km/s').value for x in dat['e_cz']])
    dat['ZWARN'] = dat['q_cz']
    
    del  dat['RAJ2000']
    del  dat['DEJ2000']
    del  dat['cz']
    del  dat['e_cz']
    del  dat['q_cz']
    
    ##  print(dat)

  elif survey in ['2dFGRS-north', '2dFGRS-south']:
    dat['RA']    = dat['RAJ2000']
    dat['DEC']   = dat['DEJ2000']
    dat['Z']     = dat['z']
    dat['ZERR']  = -99. * np.ones_like(dat['z'])
    dat['ZWARN'] = dat['q_z']
    
    del  dat['RAJ2000']
    del  dat['DEJ2000']
    del  dat['z']

    ##  print(dat)
    
  elif survey in ['ages_reduced-north', 'ages_reduced-south']:
    dat['RA']    = dat['RAJ2000']
    dat['DEC']   = dat['DEJ2000']
    dat['Z']     = dat['z1']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'])
    dat['ZWARN'] = -99. * np.ones_like(dat['Z'])

    del  dat['RAJ2000']
    del  dat['DEJ2000']
    del  dat['z1']

    ##  print(dat)
    
  elif survey in ['2dflens-north', '2dflens-south']:
    dat['RA']  = dat['R.A.']
    dat['DEC'] = dat['Dec.']
    dat['Z']   = dat['z']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'])
    dat['ZWARN'] = -99. * np.ones_like(dat['Z'])

    del  dat['R.A.']
    del  dat['Dec.']
    del  dat['z']

    ##  print(dat) 

  elif survey in ['VIPERS_W1_SPECTRO_PDR2-north', 'VIPERS_W1_SPECTRO_PDR2-south', 'VIPERS_W4_SPECTRO_PDR2-north', 'VIPERS_W4_SPECTRO_PDR2-south']:
    dat['RA']  = dat['alpha']
    dat['DEC'] = dat['delta']
    dat['Z']   = dat['zspec']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'])
    dat['ZWARN'] = -99. * np.ones_like(dat['Z'])

    del  dat['alpha']
    del  dat['delta']
    del  dat['zspec']

    ##  print(dat)

  elif survey in ['sdss-north', 'sdss-south']:    
    dat['RA']  = dat['PLUG_RA']
    dat['DEC'] = dat['PLUG_DEC']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'])
    dat['ZWARN'] = -99. * np.ones_like(dat['Z'])

    del  dat['PLUG_RA']
    del  dat['PLUG_DEC']
    
    ##  print(dat)

  elif survey in ['OzDES-north', 'OzDES-south']:
    dat['Z'] = dat['z']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'])
    dat['ZWARN'] = -99. * np.ones_like(dat['Z'])
    
    del  dat['z']

    ##  print(dat)

  elif survey in ['GAMA-north', 'GAMA-south']:
    ##  print(survey)
    ##  print(dat)
    dat['ZERR']  = -99. * np.ones_like(dat['Z'])
    dat['ZWARN'] = -99. * np.ones_like(dat['Z'])

  elif survey in ['hsc_pdr1_deep.forced.reduced-north', 'hsc_pdr1_deep.forced.reduced-south', 'hsc_pdr1_udeep.forced.reduced-north', 'hsc_pdr1_udeep.forced.reduced-south',\
                  'hsc_pdr1_wide.forced.reduced-north', 'hsc_pdr1_wide.forced.reduced-south']:    
    ##                                                                                                                                                                                                                                
    dat['RA']    = dat['ra']
    dat['DEC']   = dat['dec']
    dat['Z']     = dat['frankenz_photoz_best']
    dat['ZERR']  = -99. * np.ones_like(dat['Z'].quantity)
    dat['ZWARN'] = -99  * np.ones_like(dat['Z'].quantity, dtype=np.int)

    '''                                                                                                                                                                                                  
    a_g          = 1.155636 * (1 - MW_TRANSMISSION_G) - 0.001767                                                                                                                                        
    a_r          = 0.818918 * (1 - MW_TRANSMISSION_G) - 0.001252                                                                                                                                          
    a_i          = 0.584431 * (1 - MW_TRANSMISSION_G) - 0.000893                                                                                                                                                
    a_z          = 0.450745 * (1 - MW_TRANSMISSION_G) - 0.000689                                                                                                                                       
    a_y          = 0.384616 * (1 - MW_TRANSMISSION_G) - 0.000588                                                                                                                                      
    '''

    del  dat['ra']
    del  dat['dec']

  elif survey in ['primus-north', 'primus-south']:
    ##  print(survey)                                                                                                                                                                                                                  
    ##  print(dat)
    dat['Z']     = dat['ZPRIMUS']
    dat['ZERR']  = dat['ZPRIMUS_ERR']
    dat['ZWARN'] = dat['ZPRIMUS_ZWARNING']
    
  else:
    print('Error.')
    
    print(survey)
    print(dat.columns)
    print(dat)
    
    exit(1)


  fout  = _file.split('.fits')
  fout  = fout[0] + '-standard.fits'

  fout  = fout.split('/')
  fout  = '/'.join(x for x in fout[:-1]) + '/standard/' + fout[-1]
    
  print('Writing ... {}'.format(fout))

  try:
    dat.write(fout, format='fits', overwrite=True)

  except:
    print('\n\nFailed to write: {}'.format(survey))
      
print('\n\nDone.\n\n')

