import glob
import numpy  as np
import pandas as pd


pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)

def process_df(i, _file):
  df            = pd.read_csv(_file, sep='\t', header=0, index_col=False)
  df            = df.rename(columns={"# SKY": "SKY"})
  df            = df.rename(columns={"# WD": "WD"})
  df.index.name = 'FIBER'
  df            = df.set_index(i*5000 + df.index.values, drop=False)

  return df

##
print('\n\nWelcome.\n\n')

##
applyto          = 'assign'

files            = glob.glob('bitgrid/{}/txt/tile_*.txt'.format(applyto))
skies            = glob.glob('bitgrid/{}/txt/sky_tile_*.txt'.format(applyto))

sky              = process_df(0, skies[0])
df               = process_df(0, files[0])

for i, _sky  in enumerate(skies[1:]):
  sky            = sky.append(process_df(i, _sky))

for i, _file in enumerate(files[1:]):
  df             = df.append(process_df(i, _file))

##
sky.index.name   = 'FIBER'
df.index.name    = 'FIBER'

sky.columns      = [x.replace('_', ' ') for x in sky.columns]

skies            = pd.DataFrame(sky.sum(axis = 0, skipna = True), dtype=np.int32)
skies.columns    = ['N' + applyto.upper()]

skies['Comment'] = [''] * len(skies)
skies['Comment'] = skies['Comment'].astype('object')

skies['Comment'][0] = 'Sky.'
skies['Comment'][1] = 'Bad sky.'
skies['Comment'][2] = 'Supplementary sky.'

##  print(skies)

##  print(skies)

##
cols            = list(df.columns)
ncol            = len(cols)

##  print(ncol, cols)


summary         = pd.DataFrame(np.zeros((ncol, ncol), dtype=np.int64), columns=cols)
summary         = summary.set_index(np.array(cols).T, drop=False)

for x in cols:
  for y in cols:
    ##  summary.at[x, x]  = len(df[df[x] + df[y] == 2])
    summary[x][y] = len(df[df[x] + df[y] == 2])

##
legend          = summary.copy()[['WD']]
legend.columns  = ['N' + applyto.upper()]

legend['N' + applyto.upper()] = np.diag(summary)

##
legend['Comment']     =  ''
legend['Comment'][0]  =  'Standard - White Dwarf.'
legend['Comment'][1]  =  'Standard - Dark.'
legend['Comment'][2]  =  'Standard - Bright.'
legend['Comment'][3]  =  'Milky Way Survey.'
legend['Comment'][4]  =  'Bright Galaxy Survey.'
legend['Comment'][5]  =  'Luminous Red Galaxy.'
legend['Comment'][6]  =  'Luminous Red Galaxy - Init.'
legend['Comment'][7]  =  'Luminous Red Galaxy - Super.'
legend['Comment'][8]  =  'Luminous Red Galaxy - Low-redshift.'
legend['Comment'][9]  =  'Emission Line Galaxy.'
legend['Comment'][10] =  'Quasar.'
legend['Comment'][11] =  'Secondary.'
legend['Comment'][12] =  'Bright.'
legend['Comment'][13] =  'Faint.'
legend['Comment'][14] =  'Faint extended.'
legend['Comment'][15] =  'Fiber magnitude limited.'
legend['Comment'][16] =  'Low quality.'
legend['Comment'][17] =  'Medium star mask.'
legend['Comment'][18] =  'Large galaxy mask.'
legend['Comment'][19] =  'Globular cluster mask.'
legend['Comment'][20] =  'NOBS $<$ 1 in any of grz.'
legend['Comment'][21] =  'Profile-weighted fraction of pixels masked.'
legend['Comment'][22] =  'Profile-weighted fraction of flux from others.'
legend['Comment'][23] =  'Fraction of a sources flux within the blob.'
legend['Comment'][24] =  'Inverse variance of flux.'
legend['Comment'][25] =  'Bright star mask.'
legend['Comment'][26] = r'$(r-g) < -1$.'
legend['Comment'][27] = r'$(r-g) > 4$.'
legend['Comment'][28] = r'$(z-r) < -1$.'
legend['Comment'][29] =	r'$(z-r) >  4$.'
legend['Comment'][30] = r'Stellar by GAIA color.'
legend['Comment'][31] = r'GAIA detection in G.'
legend['Comment'][32] = r'PSF Tractor morphology.'
legend['Comment'][33] = r'PSF by GAIA AEN.'
legend['Comment'][34] = r'Bad column, bright pixel etc (ALL).'
legend['Comment'][35] = r'Saturated (ALL).'
legend['Comment'][36] = r'Bleed trail (ALL).'
legend['Comment'][37] = r'Edge pixel (ALL).'
legend['Comment'][38] = r'{\tt legacypipe} outlier pixel (ALL).'

legend                = skies.append(legend, ignore_index=True)

print('\n\n')
print(legend)

for x in cols:
  summary.at[x, x]    = -99

print('\n\n')
print(summary)

##
column_format         = '|c|c|c|'
legend                = legend.to_latex(column_format=column_format, bold_rows=False, multirow=True, escape=False)  
lines                 = legend.split('\n')

with open('bitgrid/{}/tex/{}_legend.tex'.format(applyto, applyto), 'w') as tt:
  tt.write(lines[0])
  tt.write('\n')
  tt.write(lines[1])
  tt.write('\n')
  tt.write(lines[2])
  tt.write('\n')
  tt.write(lines[3])
  tt.write('\n')
  tt.write(lines[4])
  tt.write('\n')
  tt.write(lines[5])
  tt.write('\n')
  tt.write(lines[6])
  tt.write('\n')
  tt.write('\midrule')

  for line in lines[7:24]:
    tt.write(line)         
    tt.write('\n')

  tt.write('\midrule')

  for line in lines[24:-1]:
    tt.write(line)
    tt.write('\n')

##
column_format = '|c' * (1 + ncol) + '|'

summary       = summary.to_latex(column_format=column_format, bold_rows=False, multirow=True)
summary       = summary.replace('-99', '--')
lines         = summary.split('\n')  

with open('bitgrid/{}/tex/{}_summary.tex'.format(applyto, applyto), 'w') as tf:
  tf.write(lines[0])
  tf.write('\n')
  tf.write(lines[1])
  tf.write('\n')
  tf.write(lines[2])
  tf.write('\n')
  tf.write(lines[3])
  tf.write('\n')
  
  for line in lines[4:-4]:
    tf.write(line)                                                                                                                                                                                                                                                                  
    tf.write('\n')
    tf.write('\midrule')
    tf.write('\n')

  ##
  tf.write(lines[-4])
  tf.write('\n')
  tf.write(lines[-3])
  tf.write('\n')
  tf.write(lines[-2])                                                                                                                                                                                                                                                            
  tf.write('\n')                                                                                                                                                                                                                                                               
  tf.write(lines[-1])           

print('\n\nDone.\n\n')
