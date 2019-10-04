import glob
import fitsio
import numpy               as      np
import astropy.io.fits     as      fits

from   astropy.table       import  Table, vstack
from   desitarget.targets  import  encode_targetid


def return_tractors(sfile, hsphere, verbose=False):
  _sname = sfile.split('.fits')[0].split('-')[1:]
  lims   = []

  for sname in _sname:
    if 'p' in sname:
      lims.append([np.float(sname.split('p')[0]), np.float(sname.split('p')[1])])

    else:
      lims.append([np.float(sname.split('m')[0]), -1.0 * np.float(sname.split('m')[1])])
      
  ##  Tractor brick is ~ 0.25 x 0.25 deg. 
  dirs     = np.arange(np.int(np.floor(lims[0][0])) - 1, np.int(np.ceil(lims[1][0])) + 1)
  isin     = []
  
  for _dir in dirs:
    ##  
    tractors = glob.glob('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/tractor/{}/*.fits'.format(hsphere, _dir))

    ##  print(_dir)
    
    ##  for x in tractors:
    ##    print(x)
  
    for tfile in tractors:
      tname  = tfile.split('.fits')[0].split('-')[1]

      if 'p' in tname:
        ##  Brick center.                                                                                                                                                                                                              
        ra, dec = np.float(tname.split('p')[0]) / 10., np.float(tname.split('p')[1]) / 10.

      else:
        ##  Brick center.                                                                                                                                                                                                             
        ra, dec = np.float(tname.split('m')[0]) / 10., -1.0 * np.float(tname.split('m')[1]) / 10.

      _isin = (ra + 0.15 >= lims[0][0]) & (dec + 0.15 >= lims[0][1]) & (ra - 0.15 <= lims[1][0]) & (dec - 0.15 <= lims[1][1])

      if _isin:
        isin.append(tfile)
   
  if verbose:
    print('Solving for {}.'.format(sfile))
    print('Sweep lower bound: {} {}'.format(lims[0][0], lims[0][1]))
    print('Sweep upper bound: {} {}'.format(lims[1][0], lims[1][1]))

    for x in tractors:
      print(x)
        
  return  isin
