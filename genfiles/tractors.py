import  os
import  glob
import  fitsio
import  pandas              as      pd
import  numpy               as      np
import  astropy.io.fits     as      fits

from    astropy.table       import  Table, vstack
from    desitarget.targets  import  encode_targetid

##
cols         = ['release', 'brick_primary', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'brickid', 'objid']

for hsphere in ['north', 'south']:  ##  ['north', 'south']
  tractors   = glob.glob('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/tractor/*/*.fits'.format(hsphere))
  
  for i, _tractor in enumerate(tractors):
    name     = _tractor.split('/')[-1].split('.fits')[0]
    
    _output  = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tractors/{}/{}.fits'.format(hsphere, name)

    done     = os.path.exists(_output)

    if not done:
      print('Solving for {} of {}.'.format(i, len(tractors)))

      name     = _tractor.split('/')[-1].split('.fits')[0]

      tractor  = fitsio.FITS(_tractor)

      rows     = tractor[1].get_nrows() 

      dtypes   = tractor[1].get_rec_dtype()[0]

      data     = np.array(tractor[1][cols][:])

      ##  encode_targetid(objid=data['objid'], brickid=data['brickid'], release=data['release'], sky=0, mock=0)
     
      output   = fitsio.FITS(_output, 'rw')
      output.write(data, names = cols)
    
      tractor.close()

    else:
      print('File {} exists.'.format(_output))
