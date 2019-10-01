import glob
import astropy.io.fits     as      fits

from   astropy.table       import  Table, vstack
from   desitarget.targets  import  encode_targetid


cols    = ['BRICKID', 'OBJID', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'ALLMASK_G', 'ALLMASK_R', 'ALLMASK_Z', 'ANYMASK_R', 'ANYMASK_G', 'ANYMASK_Z']
dtypes  = ['int32', 'int32', 'float32', 'float32', 'float32', 'int16', 'int16' , 'int16', 'int16', 'int16', 'int16']

hsphere = 'north'

files   = glob.glob('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/sweep/8.0/*.fits'.format(hsphere))

for i, _file in enumerate(files):
  name   = _file.split('.fits')[0].split('-')
  name   = name[1] + '-' + name[2]

  print('Solving for {}, {} of {}.'.format(name, i, len(files)))
  
  dat             = Table(fits.open(_file)[1].data)[cols]
  dat['TARGETID'] = encode_targetid(objid=dat['OBJID'], brickid=dat['BRICKID'], release=8, mock=0, sky=0)
  
  dat.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/sweeps/{}/sweep-{}-vadd.fits'.format(hsphere, name), format='fits', overwrite=True)
