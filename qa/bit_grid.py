import os
import sys
import pandas                        as      pd
import pylab                         as      pl
import numpy                         as      np
import astropy.io.fits               as      fits

from   astropy.table                 import  Table, vstack, Column
from   desimodel.io                  import  load_desiparams, load_platescale
from   desitarget.cuts               import  _psflike, _check_BGS_targtype_sv
from   desitarget.targetmask         import  desi_mask, bgs_mask
from   desitarget.sv1.sv1_targetmask import  desi_mask                        as svdesi_mask
from   desitarget.sv1.sv1_targetmask import  bgs_mask                         as  svbgs_mask
from   collections                   import  OrderedDict


BGS_MASKBITS = OrderedDict()

BGS_MASKBITS['BMB_NOBS']        =   0x1,    # (gnobs >= 1) & (rnobs >= 1) & (znobs >= 1)
BGS_MASKBITS['BMB_FRACMASK']    =   0x2,    # (gfracmasked < 0.4) & (rfracmasked < 0.4) & (zfracmasked < 0.4)
BGS_MASKBITS['BMB_FRACFLUX']    =   0x4,    # (gfracflux < 5.0) & (rfracflux < 5.0) & (zfracflux < 5.0)
BGS_MASKBITS['BMB_FRACIN']      =   0x8,    # (gfracin > 0.3) & (rfracin > 0.3) & (zfracin > 0.3)
BGS_MASKBITS['BMB_FLUXIVAR']    =   0x10,   # (gfluxivar > 0) & (rfluxivar > 0) & (zfluxivar > 0)
BGS_MASKBITS['BMB_BRIGHTMASK']  =   0x20,   # (maskbits & 2**1) == 0
BGS_MASKBITS['BMB_RMGLO']       =   0x40,   # (rflux > gflux * 10**(-1.0/2.5))
BGS_MASKBITS['BMB_RMGHI']       =   0x80,   # (rflux < gflux * 10**(4.0/2.5))
BGS_MASKBITS['BMB_ZMRLO']       =   0x100,  # (zflux > rflux * 10**(-1.0/2.5))
BGS_MASKBITS['BMB_ZMRHI']       =   0x200,  # (zflux < rflux * 10**(4.0/2.5))
BGS_MASKBITS['BMB_GRR']         =   0x400,  # (Grr > 0.6)
BGS_MASKBITS['BMB_ZEROGAIAMAG'] =   0x800,  # (gaiagmag == 0)
BGS_MASKBITS['BMB_PSF']         =   0x1000, # (_psflike(objtype))
BGS_MASKBITS['BMB_AEN']         =   0x2000, # (Non-PSF by Gaia AEN)   

def isBGS_colors(rflux=None, rfiberflux=None, south=True, targtype=None, primary=None):
    """
    Standard set of masking cuts used by all BGS target selection classes
    (see, e.g., :func:`~desitarget.cuts.isBGS` for parameters).
    """

    if primary is None:
        primary = np.ones_like(rflux, dtype='?')

    bgs = primary.copy()

    if targtype == 'lowq':
        bgs &= rflux > 10**((22.5-20.1)/2.5)

    elif targtype == 'bright':
        bgs &= rflux > 10**((22.5-19.5)/2.5)

    elif targtype == 'faint':
        bgs &= rflux > 10**((22.5-20.1)/2.5)
        bgs &= rflux <= 10**((22.5-19.5)/2.5)

    elif targtype == 'faint_ext':
        bgs &= rflux > 10**((22.5-20.5)/2.5)
        bgs &= rflux <= 10**((22.5-20.1)/2.5)

    elif targtype == 'fibmag':
        bgs &= rflux <= 10**((22.5-20.1)/2.5)
        bgs &= rfiberflux > 10**((22.5-21.0511)/2.5)

    else:
        _check_BGS_targtype_sv(targtype)

    return  bgs

def notinBGS_mask(gflux=None, rflux=None, zflux=None, gnobs=None, rnobs=None, znobs=None, primary=None,
                  gfracmasked=None, rfracmasked=None, zfracmasked=None,
                  gfracflux=None, rfracflux=None, zfracflux=None,
                  gfracin=None, rfracin=None, zfracin=None, w1snr=None,
                  gfluxivar=None, rfluxivar=None, zfluxivar=None, Grr=None,
                  gaiagmag=None, maskbits=None, objtype=None, targtype=None):
    """
    Standard set of masking cuts used by all BGS target selection classes
    (see, e.g., :func:`~desitarget.cuts.isBGS_faint` for parameters).
    """
    _check_BGS_targtype_sv(targtype)

    if primary is None:
        primary = np.ones_like(gnobs, dtype='?')

    bgs_qcs   = primary.copy()
    bgs       = primary.copy()
    
    # quality cuts definitions
    bgs_qcs  &= (gnobs >= 1) & (rnobs >= 1) & (znobs >= 1)
    bgs_qcs  &= (gfracmasked < 0.4) & (rfracmasked < 0.4) & (zfracmasked < 0.4)
    bgs_qcs  &= (gfracflux < 5.0) & (rfracflux < 5.0) & (zfracflux < 5.0)
    bgs_qcs  &= (gfracin > 0.3) & (rfracin > 0.3) & (zfracin > 0.3)
    bgs_qcs  &= (gfluxivar > 0) & (rfluxivar > 0) & (zfluxivar > 0)
    bgs_qcs  &= (maskbits & 2**1) == 0

    # color box
    bgs_qcs  &= rflux > gflux * 10**(-1.0/2.5)
    bgs_qcs  &= rflux < gflux * 10**(4.0/2.5)
    bgs_qcs  &= zflux > rflux * 10**(-1.0/2.5)
    bgs_qcs  &= zflux < rflux * 10**(4.0/2.5)

    if targtype == 'lowq':
        bgs &= Grr > 0.6
        bgs |= gaiagmag == 0
        bgs |= (Grr < 0.6) & (~_psflike(objtype)) & (gaiagmag != 0)
        bgs &= ~bgs_qcs

    else:
        bgs &= Grr > 0.6
        bgs |= gaiagmag == 0
        bgs |= (Grr < 0.6) & (~_psflike(objtype)) & (gaiagmag != 0)
        bgs &= bgs_qcs

    return  bgs

def set_BGSMASKBITS(gflux=None, rflux=None, zflux=None, gnobs=None, rnobs=None, znobs=None, primary=None,
                    gfracmasked=None, rfracmasked=None, zfracmasked=None,
                    gfracflux=None, rfracflux=None, zfracflux=None,
                    gfracin=None, rfracin=None, zfracin=None,
                    gfluxivar=None, rfluxivar=None, zfluxivar=None, Grr=None,
                    gaiagmag=None, maskbits=None, objtype=None, gaia_aen=None):
    """                                                                                                                                                                                                                                 
    MASKBITS equivalents for quality cuts applied. 
    """
    result = np.zeros_like(gflux, dtype=np.int16)

    # quality cuts definitions                                                                                                                                                                                                          
    result[~ ((gnobs >= 1) & (rnobs >= 1) & (znobs >= 1))]                      += BGS_MASKBITS['BMB_NOBS']
    result[~ ((gfracmasked < 0.4) & (rfracmasked < 0.4) & (zfracmasked < 0.4))] += BGS_MASKBITS['BMB_FRACMASK'] 
    result[~ ((gfracflux < 5.0) & (rfracflux < 5.0) & (zfracflux < 5.0))]       += BGS_MASKBITS['BMB_FRACFLUX']
    result[~ ((gfracin > 0.3) & (rfracin > 0.3) & (zfracin > 0.3))]             += BGS_MASKBITS['BMB_FRACIN']
    result[~ ((gfluxivar > 0) & (rfluxivar > 0) & (zfluxivar > 0))]             += BGS_MASKBITS['BMB_FLUXIVAR']
    result[~ ((maskbits & 2**1) == 0)]                                          += BGS_MASKBITS['BMB_BRIGHTMASK']
    
    # color box                                                                                                                                                                                                                         
    result[~ (rflux > gflux * 10**(-1.0/2.5))]                                  += BGS_MASKBITS['BMB_RMGLO']
    result[~ (rflux < gflux * 10**(4.0/2.5))]                                   += BGS_MASKBITS['BMB_RMGHI']
    result[~ (zflux > rflux * 10**(-1.0/2.5))]                                  += BGS_MASKBITS['BMB_ZMRLO']
    result[~ (zflux < rflux * 10**(4.0/2.5))]                                   += BGS_MASKBITS['BMB_ZMRHI']

    result[  ((gaiagmag > 0.0) & (Grr   <= 0.6))]                               += BGS_MASKBITS['BMB_GRR']
    result[~ (gaiagmag == 0.0)]                                                 += BGS_MASKBITS['BMB_ZEROGAIAMAG']
    result[~ (~_psflike(objtype))]                                              += BGS_MASKBITS['BMB_PSF']


    ispsf_byaen    = (gaiagmag < 18) & (gaia_aen < 10. ** 0.5)
    ispsf_byaen    = ispsf_byaen | ((gaiagmag >= 18) & (gaia_aen < 10. ** (0.5 + 0.2 * (gaiagmag - 18.))))    

    ##  Is is GAIA to start. 
    ispsf_byaen    = ispsf_byaen & (gaiagmag > 0.0)

    result[ispsf_byaen]                                                         += BGS_MASKBITS['BMB_AEN']
    
    return  result

##  
applyto   = 'assign'

scratch   = os.environ['CSCRATCH']

tiles     = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)
utiles    = np.unique(tiles['TILEID'].quantity)

ds_dtypes = ['SKY', 'BAD_SKY', 'SUPP_SKY']

##  ['BRIGHT_OBJECT', 'IN_BRIGHT_OBJECT', 'NEAR_BRIGHT_OBJECT']
sv_dtypes = ['STD_WD', 'STD_FAINT', 'STD_BRIGHT', 'MWS_ANY', 'BGS_ANY', 'LRG', 'LRG_INIT_4PASS', 'LRG_SUPER_4PASS', 'LOWZ_FILLER', 'ELG', 'QSO', 'SCND_ANY']
sv_btypes = ['BGS_BRIGHT', 'BGS_FAINT', 'BGS_FAINT_EXT', 'BGS_FIBMAG', 'BGS_LOWQ']

names     = []
##  names    += ['SKY', 'BAD', 'SUP']
names    += ['WD', 'STD DRK', 'STD BRT', 'MWS', 'BGS', 'LRG', 'LRG INIT', 'LRG SUPER', 'LRG LOZ', 'ELG', 'QSO', 'SCD']
names    += ['BRT', 'FNT', 'FEXT', 'FMAG', 'LOQ']
names    += ['MEDIUM', 'LSLGA', 'SCLUSTER']
names    += ['NOBS', 'FMASK', 'FFLUX', 'FIN', 'IVAR', 'BMASK', 'RMGLO', 'RMGHI', 'ZMRLO', 'ZMRHI', 'GRR', 'GG', 'PSF', 'AEN PSF']
names    += ['BPIX', 'STR', 'BLEED', 'EDGE', 'OLIER']

print('\n\nWelcome.\n\n')

for ii, tile in enumerate(utiles):  
  skyrow  = []
  row     = []
  
  try:
    print('Solving for Tile {} ({} of {}).'.format(tile, ii, len(utiles)))

    if applyto == 'assign':
      fname = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile)

    elif applyto == 'targets':
      fname = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/mtls/svmtl_{:06}.fits'.format(tile)

    else:
      raise  ValueError('Must apply to either assign or targets.') 
      
    _fits = fits.open(fname)  
        
  except:
    print('Unable to retrieve {}.'.format(fname))
    continue

  ##  
  dat   = Table(_fits[1].data)  ##  ['FIBER', 'DESI_TARGET', 'BGS_TARGET', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET'] 

  ##  print(dat.columns)                                                                                                                                                                                                          

  if applyto == 'assign':
    ##  Remove skies.                                                                                                                                                                                                               
    isin  = np.ones(len(dat), dtype=bool)

    for x in ds_dtypes:
      isin = isin & ((dat['DESI_TARGET'] & desi_mask.mask(x)) == 0)

    skies = dat[~isin]
    dat   = dat[ isin]

    ##                                                                                                                                                                                                                              
    dat.sort('FIBER')

  ##
  ttypes        = ['faint', 'bright', 'faint_ext', 'lowq', 'fibmag']
  
  dat['GFLUX']  = dat['FLUX_G'] / dat['MW_TRANSMISSION_G']
  dat['RFLUX']  = dat['FLUX_R'] / dat['MW_TRANSMISSION_R']
  dat['ZFLUX']  = dat['FLUX_Z'] / dat['MW_TRANSMISSION_Z']

  dat['RFFLUX'] = dat['FIBERFLUX_R'] / dat['MW_TRANSMISSION_R']
  
  ##  Weird there's no MW transmission.
  dat['Grr']    = dat['GAIA_PHOT_G_MEAN_MAG'] - 22.5 + 2.5 * np.log10(dat['FLUX_R'])
  dat['SNR_W1'] = dat['FLUX_W1'] * np.sqrt(dat['FLUX_IVAR_W1'])

  ## 
  dat['BGS_MASKBITS'] = set_BGSMASKBITS(gnobs=dat['NOBS_G'], rnobs=dat['NOBS_R'], znobs=dat['NOBS_Z'], primary=None, gfracmasked=dat['FRACMASKED_G'], rfracmasked=dat['FRACMASKED_R'], zfracmasked=dat['FRACMASKED_Z'],\
                                        gfracflux=dat['FRACFLUX_G'], rfracflux=dat['FRACFLUX_R'], zfracflux=dat['FRACFLUX_Z'], gfracin=dat['FRACIN_G'], rfracin=dat['FRACIN_R'], zfracin=dat['FRACIN_Z'],\
                                        gfluxivar=dat['FLUX_IVAR_G'], rfluxivar=dat['FLUX_IVAR_R'], zfluxivar=dat['FLUX_IVAR_Z'], Grr=dat['Grr'], gaiagmag=dat['GAIA_PHOT_G_MEAN_MAG'], maskbits=dat['MASKBITS'],\
                                        objtype=dat['MORPHTYPE'], gflux=dat['GFLUX'], rflux=dat['RFLUX'], zflux=dat['ZFLUX'], gaia_aen=dat['GAIA_ASTROMETRIC_EXCESS_NOISE'])
  
  ##
  color_cut           = np.zeros_like(dat['RFLUX'], dtype=bool)
  in_mask             = np.zeros_like(dat['RFLUX'], dtype=bool)
  
  for target in ['lowq']:
    color_cut   = color_cut | isBGS_colors(rflux=dat['RFLUX'], rfiberflux=dat['RFFLUX'], south=True, targtype=target, primary=None)
    
  for target in ['lowq']:
    ##  Low-quality overlap with STD. BRIGHT, STD. FAINT and MWS.
    in_mask     = in_mask | notinBGS_mask(gnobs=dat['NOBS_G'], rnobs=dat['NOBS_R'], znobs=dat['NOBS_Z'], primary=None, gfracmasked=dat['FRACMASKED_G'], rfracmasked=dat['FRACMASKED_R'], zfracmasked=dat['FRACMASKED_Z'],\
                                          gfracflux=dat['FRACFLUX_G'], rfracflux=dat['FRACFLUX_R'], zfracflux=dat['FRACFLUX_Z'], gfracin=dat['FRACIN_G'], rfracin=dat['FRACIN_R'], zfracin=dat['FRACIN_Z'],\
                                          gfluxivar=dat['FLUX_IVAR_G'], rfluxivar=dat['FLUX_IVAR_R'], zfluxivar=dat['FLUX_IVAR_Z'], Grr=dat['Grr'], gaiagmag=dat['GAIA_PHOT_G_MEAN_MAG'], maskbits=dat['MASKBITS'],\
                                          targtype=target, w1snr=dat['SNR_W1'], objtype=dat['MORPHTYPE'], gflux=dat['GFLUX'], rflux=dat['RFLUX'], zflux=dat['ZFLUX'])
    
  for i, x in enumerate(color_cut):
    lowq = (dat['SV1_BGS_TARGET'][i]  & svbgs_mask.mask('BGS_LOWQ')) != 0

    if lowq == (x & in_mask[i]):
        continue

    else:
        ##  Low-quality overlap with STD. BRIGHT, STD. FAINT and MWS.
        print('\n\nMismatches between LOWQ def. and that in catalogue:  {} {}'.format(lowq, x & in_mask[i]))

        if applyto == 'assign':        
          for x in ds_dtypes:
            print(x, (dat['DESI_TARGET'][i]     &   desi_mask.mask(x)) != 0)

        for x in sv_dtypes:
            print(x, (dat['SV1_DESI_TARGET'][i]     &   svdesi_mask.mask(x)) != 0)

  if applyto == 'assign':            
    for x in ds_dtypes:
      skyrow.append((skies['DESI_TARGET']  & desi_mask.mask(x))   != 0)

  for x	in sv_dtypes:
    row.append((dat['SV1_DESI_TARGET'] & svdesi_mask.mask(x)) != 0)

  for x in sv_btypes:
    row.append((dat['SV1_BGS_TARGET']  & svbgs_mask.mask(x))  != 0)

  for bit, x in zip([11, 12, 13], ['MEDIUM', 'LSLGA', 'SC-NGC']):
    row.append((dat['MASKBITS']        & 2 ** bit)            != 0)
    
  for x in BGS_MASKBITS.keys():
    row.append((dat['BGS_MASKBITS']    & BGS_MASKBITS[x])     != 0)

  for bit, _type in zip([0, 1, 6, 8, 11], ['BPIX', 'STR', 'BLEED', 'EDGE', 'OUTLIER']):
    is_set = ((dat['ALLMASK_G'] & 2 ** bit) != 0) | ((dat['ALLMASK_R'] & 2 ** bit) != 0) | ((dat['ALLMASK_Z'] & 2 ** bit) != 0)
    row.append(is_set)
    
  ##
  skyrow = np.array(skyrow).astype(np.int).T
  row    = np.array(row).astype(np.int).T

  ##                                                                                                                                                                                                                                   
  np.savetxt('bitgrid/{}/txt/sky_tile_{:06d}.txt'.format(applyto, tile), skyrow, fmt='%d\t', header='\t'.join(ds_dtypes))
  np.savetxt('bitgrid/{}/txt/tile_{:06d}.txt'.format(applyto, tile), row, fmt='%d\t', header='\t'.join(names))
  '''
  ##
  batch = 25

  for i in range(200):  
    df  = pd.DataFrame(data=row[i * batch : (i + 1) * batch], index=np.arange(i * batch, (i + 1) * batch), columns=names)  # 1st row as the column names
    df.index.name  = '{:06d}'.format(tile)

    column_format  = '|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|'
    column_format += '|c|c|c|c|c|c|c|c|c|c|c|c|c|'
               
    table   = df.to_latex(column_format=column_format, bold_rows=False, multirow=True)
    lines   = table.split('\n')
    
    with open('bitgrid/{}/tex/tile_{:06d}_{}.tex'.format(applyto, tile, i * batch), 'w') as tf:
      tf.write(lines[0])
      tf.write(lines[1])
      tf.write(lines[2])
      tf.write(lines[3])
      tf.write(lines[4])
      
      for line in lines[5:-4]:
        tf.write(line)
        tf.write('\n')
        tf.write('\midrule')
        tf.write('\n')

      tf.write(lines[-4])
      tf.write(lines[-3])
      tf.write(lines[-2])
      tf.write('\n')
      tf.write(lines[-1])
  '''    
  ##  break

print('\n\nDone.\n\n')
  
