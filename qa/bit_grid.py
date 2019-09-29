import os
import sys
import pandas                        as      pd
import pylab                         as      pl
import numpy                         as      np
import astropy.io.fits               as      fits

from   astropy.table                 import  Table, vstack, Column
from   desimodel.io                  import  load_desiparams, load_platescale
from   desitarget.sv1.sv1_targetmask import  desi_mask, bgs_mask


scratch   = os.environ['CSCRATCH']

tiles     = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/tiles/BGS_SV_30_3x_superset60_Sep2019.fits')[1].data)
utiles    = np.unique(tiles['TILEID'].quantity)

ds_dtypes = ['SKY', 'BAD_SKY', 'SUPP_SKY']

##  ['BRIGHT_OBJECT', 'IN_BRIGHT_OBJECT', 'NEAR_BRIGHT_OBJECT']
sv_dtypes = ['STD_WD', 'STD_FAINT', 'STD_BRIGHT', 'MWS_ANY', 'BGS_ANY', 'LRG', 'ELG', 'QSO', 'SCND_ANY']
sv_btypes = ['BGS_BRIGHT', 'BGS_FAINT', 'BGS_FAINT_EXT', 'BGS_LOWQ', 'BGS_FIBMAG']

names     = ['SKY', 'BAD', 'SUP'] + ['WD', 'STD FNT', 'STD BRT', 'MWS', 'BGS', 'LRG', 'ELG', 'QSO', 'SCD'] + ['BRT', 'FNT', 'FEXT', 'LOQ', 'FMAG']

print('\n\nWelcome.\n\n')

for tile in utiles:  
  row     = []

  try:
    print('Solving for Tile {}.'.format(tile))

    fname = '/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/fiberassign/tile-{:06}.fits'.format(tile)
    _fits = fits.open(fname)  
    dat   = Table(_fits[1].data)['FIBER', 'DESI_TARGET', 'BGS_TARGET', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET']

    ##
    dat.sort('FIBER')
    
    ##  print(dat)
    
  except:
    print('Unable to retrieve {}.'.format(fname))
    continue

  for x in ds_dtypes:
    row.append((dat['DESI_TARGET']     & desi_mask.mask(x)) != 0)

  for x	in sv_dtypes:
    row.append((dat['SV1_DESI_TARGET'] & desi_mask.mask(x)) != 0)

  for x in sv_btypes:
    row.append((dat['SV1_BGS_TARGET']  & bgs_mask.mask(x)) != 0)

  ##  
  row = np.array(row).astype(np.int).T

  batch = 25

  for i in range(200):  
    df  = pd.DataFrame(data=row[i * batch : (i + 1) * batch], index=np.arange(i * batch, (i + 1) * batch), columns=names)  # 1st row as the column names
    df.index.name = '{:06d}'.format(tile)

    table = df.to_latex(column_format='|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|', bold_rows=False, multirow=True)
    lines = table.split('\n')
    
    with open('bitgrid/tile_{:06d}_{}.tex'.format(tile, i * batch), 'w') as tf:
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
          
  break

print('\n\nDone.\n\n')
  
