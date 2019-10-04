import  os
import  glob
import  fitsio
import  pylab               as      pl
import  pandas              as      pd
import  numpy               as      np
import  astropy.io.fits     as      fits
import  matplotlib.pyplot   as      plt 
          
from    astropy.table       import  Table, vstack
from    desitarget.targets  import  encode_targetid
from    desitarget.geomask  import  is_in_box


camera      = b'mosaic'   ##  ['90prime', 'mosaic', 'decam']
band        = b'r'        ##  [b'g', b'r', b'z']

if camera  == b'mosaic':
  band      = b'z'

nrandom     = np.int(20000)

cols        = ['expnum', 'filter', 'ra_center', 'dec_center', 'ra0', 'dec0', 'ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3']
skies       = ['exptime', 'meansky', 'stdsky', 'ccdskysb', 'minsky', 'maxsky']

ccd         = fitsio.FITS('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/ccds/ccds-annotated-{}-dr8.fits'.format(camera.decode('UTF-8')))
dtype       = ccd[1].get_rec_dtype()[0]

assert  np.all(ccd[1]['camera'][:] == camera)

##
randoms     = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-2.fits')
randoms     = randoms[1]['RA', 'DEC'][:nrandom]

_file       = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/skies/skies_{}_{}.txt'.format(camera.decode('UTF-8'), band.decode('UTF-8'))

nfail       = 0

if not os.path.exists(_file):
  ccd        = ccd[1].read(vstorage='object')
  ccd        = np.array(ccd, dtype=dtype)[cols + skies]

  exps       = ccd[ccd['filter'] == band]

  result     = np.zeros(len(randoms) * len(skies)).reshape(len(randoms), len(skies))
  count      = np.zeros(len(randoms), dtype=np.int32)

  for i, x in enumerate(exps):
    try:
      ##  http://legacysurvey.org/ccdordering/  
      if camera   == b'decam':
        inccd  = is_in_box(randoms, [x['ra3'], x['ra1'], x['dec3'], x['dec1']])

      elif camera == b'90prime':
        ##  BASS
        inccd  = is_in_box(randoms, [x['ra2'], x['ra0'], x['dec2'], x['dec0']])

      elif camera == b'mosaic':
        ##  MzLS
        inccd  = is_in_box(randoms, [x['ra0'], x['ra2'], x['dec0'], x['dec2']])
      
      else:
        raise  ValueError('Invalid inut for camera.')

    except:
      nfail   += 1

      print('Failed for {}'.format([x['ra0'], x['ra2'], x['dec0'], x['dec2']]))

      continue
      
    ##  
    if len(inccd) > 0:
      print('Solving for {} of {}.'.format(i, len(exps)))

      toadd  = np.array(list(x[skies]))

      result[inccd] += toadd
      count[inccd]  +=     1
      
  ##
  result[count > 0] /= count[count > 0].astype(np.float)[:,None]

  np.savetxt(_file, result, fmt='%.6le')
  
else:
  result    = np.loadtxt(_file)

##
ncol        = 2
nrow        = np.int(np.ceil(len(skies) / 2))

##  
fig, axarr  = plt.subplots(nrows=np.int(np.ceil(len(skies) / 2)), ncols=2, figsize=(30, 10))

##  Wrap randoms
randoms['RA'][randoms['RA'] > 300.] -= 360.
randoms['RA'] += 60.

##  Cut to non-DES.
result      =  result[(-30. < randoms['DEC'])]
randoms     = randoms[(-30. < randoms['DEC'])] 

if camera   == 'decam':
  result, randoms = [result[(randoms['DEC'] < 30.)], randoms[(randoms['DEC'] < 30.)]]

elif camera == '90prime':
  result, randoms = [result[(randoms['DEC'] > 35.)], randoms[(randoms['DEC'] > 35.)]]

else:
  result, randoms = [result[(randoms['DEC'] > 35.)], randoms[(randoms['DEC'] > 35.)]]

##
vmaxs       = np.quantile(result, 0.55, axis=1)

for i, _ in enumerate(skies):
  row = i % nrow
  col = i % 2
  
  print(row, col)

  nresult  = result[:,i]

  isin     = np.isfinite(nresult)

  nresult  = nresult[isin]
  nresult  = nresult - np.median(nresult)
  nresult /= np.std(nresult)

  vmax     = np.quantile(nresult, 0.95)
  
  sc       = axarr[row][col].scatter(randoms['RA'][isin], randoms['DEC'][isin], c=nresult, s=1, rasterized=True, vmin=nresult.min(), vmax=vmax)

  plt.colorbar(sc, ax=axarr[row][col])

  axarr[row][col].set_title(skies[i].upper())
  axarr[row][col].set_xlim(360., 0.)

print('Number of failures: {}'.format(nfail))
  
plt.subplots_adjust(hspace=0.4)
##  pl.show()

print('\n\nDone.\n\n')
