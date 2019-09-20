import  os
import  numpy                 as      np
import  astropy.io.fits       as      fits

from    astropy.table         import  Table
from    desiutil.bitmask      import  BitMask
from    desitarget.targetmask import  load_mask_bits


##  https://desi.lbl.gov/trac/wiki/TargetSelectionWG/TargetingTruthTables/MatchedTruthCatalogs
scratch   = os.environ['CSCRATCH']
truthdir  = '/project/projectdirs/desi/target/analysis/truth/dr8.0/'

_all      = ['2dFGRS-match.fits', '6dFGS-match.fits', '2dflens-match.fits', 'GAMA-DR3-SpecObj-match.fits', 'OzDES-DR1-match.fits',\
             'VIPERS_W1_SPECTRO_PDR2-match.fits', 'VIPERS_W4_SPECTRO_PDR2-match.fits', 'ages_reduced-match.fits', 'sdss-specObj-dr14-unique-trimmed-match.fits']

##  Writes desired truth catalogues to scratch.
for i, survey in enumerate(_all):
  for gcap in ['north', 'south']:
    infile     = truthdir + '/{}/matched/'.format(gcap) + survey

    if not os.path.exists(infile):
      ##  print('Cannot find: {}'.format(infile))
      continue
    
    ##  Truth spectroscopic catalogue.  Look for a given survey in both north and south. 
    truth      = Table(fits.open(truthdir + '/{}/matched/'.format(gcap) + survey)[1].data)

    ##  Line matched LS output.
    infile     = truthdir + '/{}/matched/'.format(gcap) + 'ls-dr8.0-' + survey

    if not os.path.exists(infile):
      ##  print('Cannot find: {}'.format(infile))
      continue

    lsm        = Table(fits.open(truthdir + '/{}/matched/'.format(gcap) + 'ls-dr8.0-' + survey)[1].data)    
    nsurvey    = survey.split('-')[0]
    
    print('\n\n')
    print(nsurvey.upper(), gcap.upper())
    print('\n')
    print(lsm)
    
    print('Writing.. {}'.format(scratch + '/BGS/SV-ASSIGN/truth/legacy/' + 'ls-' + nsurvey + '-' + gcap + '.fits'))
    print('Writing.. {}'.format(scratch + '/BGS/SV-ASSIGN/truth/' + nsurvey + '-' + gcap + '.fits'))
    
    lsm.write(scratch + '/BGS/SV-ASSIGN/truth/legacy/' + 'ls-' + nsurvey + '-' + gcap + '.fits', format='fits', overwrite=True)
    truth.write(scratch + '/BGS/SV-ASSIGN/truth/' + nsurvey + '-' + gcap + '.fits', format='fits', overwrite=True)
