import  os
import  sys
import  glob
import  fitsio
import  pylab                   as       pl
import  pandas                  as       pd
import  numpy                   as       np
import  astropy.io.fits         as       fits
import  matplotlib.pyplot       as       plt
import  numpy.lib.recfunctions  as       rfn
import  healpy                  as       hp

from    astropy.table           import   Table, vstack
from    desitarget.targets      import   encode_targetid
from    desitarget.geomask      import   is_in_box
from    desitarget.targetmask   import   desi_mask


nside          = np.int(sys.argv[1])
parea          = hp.nside2pixarea(nside, degrees = True)

cols           = ['RA', 'DEC', 'MORPHTYPE', 'DESI_TARGET', 'FLUX_G', 'FLUX_R', 'FLUX_Z',\
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',\
                  'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'TARGETID', 'DCHISQ']

def write_elgs(elgs, mtype=None, cutlevel=None, cutcol=None, write_elgs=False, ext=''):  
  if write_elgs:
    elgs.write('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/elgs/elgs.fits', format='fits', overwrite=True)
    exit(0)
    
  opath        = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/healmaps/elg_tdensity'
  
  if mtype is not None:
    assert  np.all(elgs['MORPHTYPE'] == mtype)
    opath     += '_{}'.format(mtype.strip())
    
  if cutlevel is not None:
    if cutcol is None:
      raise ValueError()
    
    assert  np.all(elgs[cutcol] >= cutlevel)
    opath     += '_{:.1f}'.format(cutlevel)

  opath        = opath + ext + '_{}.npy'.format(nside)
    
  print('Writing {}.'.format(opath))

  elgs.sort(cutcol)
  elgs.pprint()
  
  ##
  hppix        = hp.ang2pix(nside, (90. - elgs['DEC']) * np.pi / 180., elgs['RA'] * np.pi / 180., nest=False)
  hpind, cnts  = np.unique(hppix, return_counts=True)

  theta,phi    = hp.pix2ang(nside, hpind, nest=False)
  hpra, hpdec  = 180. / np.pi * phi, 90. -180. / np.pi * theta

  tdensity     = cnts / parea  ##  per sq. deg.                                                                                                                                                    
  
  np.save(opath, np.c_[hpind, hpra, hpdec, tdensity])
  
    
if __name__ == '__main__':
 compute        = False

 ext            = '_chi2'  ##  '_r_z'

 dchi2          = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0]
 snrs           = [3., 4., 6., 9., 12., 20.]

 if compute:
  targets       = glob.glob('/project/projectdirs/desi/target/catalogs/dr8/0.32.0/targets/main/resolve/dark/*.fits')

  ##  elgs      = Table(fits.open(targets[0])[1].data)
  elgs          = Table([fits.open(targets[0])[1].data[x] for x in cols], names=cols)
  elgs          = elgs[(elgs['DESI_TARGET'] & desi_mask.mask('ELG')) != 0]

  ##  S/N should be prior to extinction correction. 
  elgs['GSNR']  = np.sqrt(elgs['FLUX_IVAR_G']) * elgs['FLUX_G'] ##  elgs['MW_TRANSMISSION_G']
  elgs['RSNR']  = np.sqrt(elgs['FLUX_IVAR_R']) * elgs['FLUX_R'] ##  elgs['MW_TRANSMISSION_R']
  elgs['ZSNR']  = np.sqrt(elgs['FLUX_IVAR_Z']) * elgs['FLUX_Z'] ##  elgs['MW_TRANSMISSION_Z']

  ##
  elgs['RZSNR'] = np.sqrt(elgs['RSNR'] ** 2. + elgs['ZSNR'] ** 2.)
  
  elgs['SNR']   = np.sqrt(elgs['GSNR'] ** 2. + elgs['RSNR'] ** 2. + elgs['ZSNR'] ** 2.)

  elgs['BFIT']     = np.array([np.array(x).max()               for x in elgs['DCHISQ']])
  elgs['BFIT2']    = np.array([np.sort(np.array(x))[:-1].max() for x in elgs['DCHISQ']])

  elgs['CHI2DIFF'] = np.clip(np.log10(elgs['BFIT'] - elgs['BFIT2']), a_min=-6., a_max=None)

  print(list(elgs.columns))

  del  elgs['DCHISQ']
  
  elgs.sort('SNR')
  
  elgs.pprint()
  
  ##  
  mtypes       = np.unique(elgs['MORPHTYPE'])
  
  for i, _targets in enumerate(targets[1:]):
    print('Joining {} of {} ELG file.  Currently {} ELGs.'.format(i, len(targets) - 1, len(elgs)))

    _in          = Table([fits.open(_targets)[1].data[x] for x in cols], names=cols)
    _in          = _in[(_in['DESI_TARGET'] & desi_mask.mask('ELG')) != 0]

    _in['GSNR']  = np.sqrt(_in['FLUX_IVAR_G']) * _in['FLUX_G'] ##  / _in['MW_TRANSMISSION_G']
    _in['RSNR']  = np.sqrt(_in['FLUX_IVAR_R']) * _in['FLUX_R'] ##  / _in['MW_TRANSMISSION_R']
    _in['ZSNR']  = np.sqrt(_in['FLUX_IVAR_Z']) * _in['FLUX_Z'] ##  / _in['MW_TRANSMISSION_Z']

    ##
    _in['RZSNR'] = np.sqrt(_in['RSNR'] ** 2. + _in['ZSNR'] ** 2.)

    ##
    _in['SNR']   = np.sqrt(_in['GSNR'] ** 2. + _in['RSNR'] ** 2. + _in['ZSNR'] ** 2.)

    _in['BFIT']     = np.array([np.array(x).max()               for x in _in['DCHISQ']])
    _in['BFIT2']    = np.array([np.sort(np.array(x))[:-1].max() for x in _in['DCHISQ']])

    _in['CHI2DIFF'] = np.clip(np.log10(_in['BFIT'] - _in['BFIT2']), a_min=-6., a_max=None)

    del  _in['DCHISQ']

    _in.sort('SNR')
    
    _in.pprint(max_width = -1)

    elgs         = vstack([elgs, _in])
  
    ##  Now checking indiv. obj. in the tbc pixels.
    if i % 5 == 0:
      ##  continue
      write_elgs(elgs, write_elgs=False, mtype=None, cutcol=None, cutlevel=None)

  ##  Forces exit. 
  write_elgs(elgs, write_elgs=True, mtype=None, cutcol=None, cutlevel=None) 
      
 else:
   elgs   = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/elgs/elgs.fits')[1].data)  
   mtypes = np.unique(elgs['MORPHTYPE'])

##
if ext == '':
  write_elgs(elgs, write_elgs=False, mtype=None, cutcol=None, cutlevel=None)
  
##
for mtype in mtypes:
  ##  for snr in snrs:
  for _dchi2 in dchi2:
    toprint = elgs[(elgs['MORPHTYPE'] == mtype) & (elgs['CHI2DIFF'] >= _dchi2) & (elgs['CHI2DIFF'] >= _dchi2)]
    write_elgs(toprint, write_elgs=False, mtype=mtype, cutlevel=_dchi2, cutcol='CHI2DIFF', ext=ext)
