import  os
import  sys
import  glob
import  fitsio
import  matplotlib
import  pylab                                 as       pl
import  pandas                                as       pd
import  numpy                                 as       np
import  astropy.io.fits                       as       fits
import  matplotlib.pyplot                     as       plt 
import  numpy.lib.recfunctions                as       rfn
import  healpy                                as       hp

from    mpl_toolkits.axes_grid1               import   make_axes_locatable
from    fast_scatter                          import   fast_scatter
from    matplotlib                            import   rc
from    astropy.table                         import   Table, vstack
from    desitarget.targets                    import   encode_targetid
from    desitarget.geomask                    import   is_in_box
from    desitarget.targetmask                 import   desi_mask
from    mpl_toolkits.axes_grid1               import   make_axes_locatable
from    mpl_toolkits.axes_grid1.inset_locator import   inset_axes


##  plt.style.use(['dark_background'])

rc('font', **{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

nside       = np.int(sys.argv[1])

camera      = str.encode(sys.argv[2])  ##  [b'90prime', b'mosaic', b'decam']
band        = str.encode(sys.argv[3])  ##  [b'g', b'r', b'z']

recompute   = False
noplot      = False
plot_elgs   = True

if camera  == b'mosaic':
  band      = b'z'

nrandom     = np.int(20000)

##  Sky rms for the entire image (in counts).
##  Our pipeline (not the CP) estimate of the sky level, average over the image, in ADU.
##  Standard deviation of our sky level.
##  Sky surface brightness (in AB mag/arcsec2).
##  Min. of our sky level.
##  Max. of our sky level.
##  FWHM (in pixels) measured by the CP.
##  Community pipeline number.    

cols        = ['expnum', 'camera', 'filter', 'ra_center', 'dec_center', 'ra0', 'dec0', 'ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3']
skies       = ['exptime',\
               'skyrms',\
               'meansky',\
               'stdsky',\
               'ccdskysb',\
               'minsky',\
               'maxsky',\
               'fwhm',\
               'plver']

def remap(x, printit=False):  
  uentries, cnts = np.unique(x, return_counts = True)

  result      = np.zeros(len(x), dtype=[('plverf', np.float32)])
  
  for u in uentries:
    entry     = u.decode('UTF-8')[1:].replace('.', '').strip()
    lentry    = len(entry)
    
    rentry    = np.float(entry) / 10. ** (lentry - 1)

    print('Remapping {}, {} ({}) {}.'.format(u, entry, lentry, rentry))
    
    result[x == u] = rentry

  ##  print(result['plverf'])
    
  ##
  print(np.unique(result['plverf'], return_counts=True))
  
  return  result

##
randoms     = fitsio.FITS('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-2.fits')
randoms     = randoms[1]['RA', 'DEC'][:nrandom]

_file       = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/skies/skies_{}_{}.txt'.format(camera.decode('UTF-8'), band.decode('UTF-8'))

nfail       = 0

if (not os.path.exists(_file)) | recompute:
  print('{} not found, recalculating.'.format(_file))

  ccd          = fitsio.FITS('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/ccds/ccds-annotated-{}-dr8.fits'.format(camera.decode('UTF-8')))
  dtype        = ccd[1].get_rec_dtype()[0]

  ccd          = ccd[1].read(vstorage='object')
  ccd          = np.array(ccd, dtype=dtype)[cols + skies]

  plverf        = remap(ccd['plver']) ##  np.array(remap(ccd['plver']), dtype=[('plverf', np.float32)])                                                                                                                                   
  ##                                                                                                                                                                                                                                  
  ccd           = rfn.merge_arrays([ccd, plverf], flatten = True, usemask = False)

  skies.remove('plver')

  skies         = skies + ['plverf']
  ccd           = ccd[cols + skies]
  
  exps          = ccd[ccd['filter'] == band]
 
  result        = np.zeros(len(randoms) * len(skies)).reshape(len(randoms), len(skies))
  count         = np.zeros(len(randoms), dtype=np.int32)

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

  if noplot:
    exit(0)
  
else:
  result    = np.loadtxt(_file)

  skies.remove('plver')

  skies     = skies + ['plverf']
  
##
ncol        = 2
nrow        = np.int(np.ceil(len(skies) / 2))

##  
fig, axarr  = plt.subplots(nrows=np.int(np.ceil(len(skies) / 2)), ncols=2, figsize=(10, 10))

##
plt.subplots_adjust(left = 0.05, right = 0.95, hspace=0.6, wspace=0.4, top = 0.925, bottom = 0.05)

##  Wrap randoms
randoms['RA'][randoms['RA'] > 300.] -= 360.
randoms['RA'] += 60.

##  Cut to non-DES.
result      =  result[(-30. < randoms['DEC'])]
randoms     = randoms[(-30. < randoms['DEC'])] 

if camera   == b'decam':
  result, randoms = [result[(randoms['DEC'] < 30.)], randoms[(randoms['DEC'] < 30.)]]

elif camera == b'90prime':
  result, randoms = [result[(randoms['DEC'] > 35.)], randoms[(randoms['DEC'] > 35.)]]

else:
  result, randoms = [result[(randoms['DEC'] > 35.)], randoms[(randoms['DEC'] > 35.)]]

##
for i, _ in enumerate(skies):
  row = i % nrow
  col = i % 2
  
  print(row, col)

  nresult  = result[:,i]

  isin     = np.isfinite(nresult)

  nresult  = nresult[isin]
  ##  nresult  = nresult - np.median(nresult)
  ##  nresult /= np.std(nresult)

  parea        = hp.nside2pixarea(nside, degrees = True)
  hppix        = hp.ang2pix(nside, (90. - randoms['DEC'][isin]) * np.pi / 180., randoms['RA'][isin] * np.pi / 180., nest=False)
  hpind, cnts  = np.unique(hppix, return_counts=True)
  
  theta,phi    = hp.pix2ang(nside, hpind, nest=False)
  hpra, hpdec  = 180. / np.pi * phi, 90. -180. / np.pi * theta

  values       = np.array([np.mean(nresult[hppix == x]) for x in hpind])

  vmin         = np.quantile(values, 0.05)
  vmax         = np.quantile(values, 0.95)
  step         = (vmax - vmin) / 50.

  fast_scatter(axarr[row][col], hpra, hpdec, values, vmin, vmax, step, markersize=0.7, cmap='jet')
  
  if i == 0:
    ylims      = axarr[row][col].get_ylim()
  
  axarr[row][col].set_title(skies[i].upper())
  axarr[row][col].set_xlim(360., 0.)

##  
if plot_elgs:
  nside        = 512

  binary       = np.load('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity_{}.npy'.format(nside))
  hpind        = binary[:,0]
  hpra         = binary[:,1]
  hpdec        = binary[:,2]
  tdensity     = binary[:,3]

  ##  Cut to non-DES.                                                                                                                                                                                                                 
  hpra         = hpra[hpdec > -30.]
  tdensity     = tdensity[hpdec > -30.]
  hpdec        = hpdec[hpdec > -30.]

  if camera   == b'decam':
    hpra       = hpra[hpdec < 30.]
    tdensity   = tdensity[hpdec < 30.]
    hpdec      = hpdec[hpdec < 30.]

  else:
    hpra       = hpra[hpdec > 35.]
    tdensity   = tdensity[hpdec > 35.]
    hpdec      = hpdec[hpdec > 35.]
    
  ##  Wrap randoms.                                                                                                                                                                                                                 
  hpra[hpra > 300.] -= 360.                                                                                                                                                                                                            
  hpra        += 60.  
    
  ##  Digitize. 
  vmin         =  np.quantile(tdensity, 0.001) 
  vmax         =  np.quantile(tdensity, 0.999)    
  step         = (vmax - vmin) / 50.
  
  fast_scatter(axarr[-1][-1], hpra, hpdec, tdensity, vmin, vmax, step, cmap='jet')
  
  axarr[-1][-1].set_title('ELG DENSITY')

  axarr[-1][-1].set_xlim(360., 0.)
  axarr[-1][-1].set_ylim(ylims)

##
fig.suptitle(r'{}      ${}$-band'.format(camera.decode('UTF-8').upper(), band.decode('UTF-8')), fontsize=14)
  
print('Number of failures: {}'.format(nfail))
  
pl.savefig('skydepths/skydepth_{}_{}.pdf'.format(camera.decode('UTF-8'), band.decode('UTF-8')))

print('\n\nDone.\n\n')
