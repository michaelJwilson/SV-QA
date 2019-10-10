import  os
import  glob
import  fitsio
import  pandas              as      pd
import  numpy               as      np
import  astropy.io.fits     as      fits

from    astropy.table       import  Table, vstack
from    desitarget.targets  import  encode_targetid
from    fitsio              import  FITS, FITSHDR


##
cols         = ['release', 'brick_primary', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'brickid', 'objid']
mcols        = ['psf_flux', 'psf_flux_ivar', 'rex_flux', 'rex_flux_ivar', 'brickid', 'objid']

verbose      =  False

for hsphere in ['north']:  ##  ['north', 'south']
  tractors   = glob.glob('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/tractor/*/*.fits'.format(hsphere))[::-1]
  
  for i, _tractor in enumerate(tractors):
    _dir     = _tractor.split('/')[-2]
    name     = _tractor.split('/')[-1].split('.fits')[0]
    brick    = name.split('-')[-1]
    
    metrics  = '/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/metrics/{}/all-models-{}.fits'.format(hsphere, _dir, brick)
    
    print(_tractor, _dir, name, metrics)
    
    os.system('mkdir -p {}'.format('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tractors/{}/{}'.format(hsphere, _dir)))

    _output  = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tractors/{}/{}/{}.fits'.format(hsphere, _dir, name)

    done     = os.path.exists(_output)
    
    if not done:
      print('Solving for {} of {}.'.format(i, len(tractors)))

      name     = _tractor.split('/')[-1].split('.fits')[0]

      tractor  = fitsio.FITS(_tractor)

      rows     = tractor[1].get_nrows() 

      dtypes   = tractor[1].get_rec_dtype()[0]

      data     = np.array(tractor[1][cols][:])

      tractor.close()

      ##
      mheader  = fitsio.read_header(metrics)
      metrics  = fitsio.FITS(metrics)
      
      grz      = np.array([mheader['BRICK_{}'.format(x)] for x in ['G', 'R', 'Z']])
      nbands   = np.count_nonzero(grz)

      rows     = metrics[1].get_nrows()
      
      mdata    = metrics[1][mcols][:]

      if len(data) != len(mdata):
        print(np.sort(data['objid']))
        print(np.sort(mdata['objid']))
        
        print('Missing objects:  {}'.format(set(data['objid']) ^ set(mdata['objid'])))

        data   = data[:len(mdata)] 

      ##
      tid      = encode_targetid(objid=data['objid'], brickid=data['brickid'], release=data['release'], sky=0, mock=0)

      ##
      fluxs    = np.zeros((rows, 6))
      ivflux   = np.zeros((rows, 6))

      for j, mtype in enumerate(['psf', 'rex']):
        counter          = 3 * j
        bcounter         = 0
        
        for i, isin in enumerate(grz):
          if isin: 
            if len(mdata['{}_flux'.format(mtype)].shape) == 1:
              fluxs[:,counter + i]  = mdata['{}_flux'.format(mtype)][bcounter]
              ivflux[:,counter + i] = mdata['{}_flux_ivar'.format(mtype)][bcounter]

            else:
              fluxs[:,counter + i]  = mdata['{}_flux'.format(mtype)][:,bcounter]
              ivflux[:,counter + i] = mdata['{}_flux_ivar'.format(mtype)][:,bcounter]
          
            bcounter    += 1

      ##
      
      radii    = ['_0p5', '_0p75', '_1p0', '_1p5', '_2p0', '_3p5', '_5p0', '_7p0']
      apf_namg = ['apflux_resid_g'.upper() + x for x in radii]
      apf_namr = ['apflux_resid_r'.upper() + x for x in radii]
      apf_namz = ['apflux_resid_z'.upper() + x for x in radii]
            
      names    = ['TARGETID', 'BRICKID', 'OBJID'] + apf_namg + apf_namr + apf_namz + ['PSF_FLUXG', 'PSF_FLUXG_IVAR', 'PSF_FLUXR', 'PSF_FLUXR_IVAR',\
                  'PSF_FLUXZ', 'PSF_FLUXZ_IVAR', 'REX_FLUXG', 'REX_FLUXG_IVAR', 'REX_FLUXR', 'REX_FLUXR_IVAR', 'REX_FLUXZ', 'REX_FLUXZ_IVAR']

      '''
      output   = Table(data=[tid, data['apflux_resid_g'], data['apflux_resid_r'], data['apflux_resid_z'], fluxs[:,0], ivflux[:,0],\
                             fluxs[:,1], ivflux[:,1], fluxs[:,2], ivflux[:,2], fluxs[:,3], ivflux[:,3], fluxs[:,4], ivflux[:,4], fluxs[:,5], ivflux[:,5]],\
                             names=names)
      
      output.pprint()

      output.write(_output, format='fits', overwrite=True)
      '''
      
      outdata   =  np.c_[tid, data['brickid'], data['objid'], data['apflux_resid_g'], data['apflux_resid_r'], data['apflux_resid_z'],\
                         fluxs[:,0], ivflux[:,0], fluxs[:,1], ivflux[:,1], fluxs[:,2], ivflux[:,2],\
                         fluxs[:,3], ivflux[:,3], fluxs[:,4], ivflux[:,4], fluxs[:,5], ivflux[:,5]]

      ##  print(len(tid), len(data), len(fluxs), len(ivflux))

      if verbose:
        toprint = Table(data=outdata, names=names)
        toprint.pprint(max_width=-1)
        
      output   = FITS(_output, 'rw')
      output.write(outdata, names=names)
      
      metrics.close()

      print('{} written.'.format(_output))
      
    else:
      continue
      
      ##  print('File {} exists.'.format(_output))
    
