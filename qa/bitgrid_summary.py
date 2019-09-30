import glob
import numpy  as np
import pandas as pd


pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)

def process_df(i, _file):
  df            = pd.read_csv(_file, sep='\t', header=0, index_col=False)
  df            = df.rename(columns={"# SKY": "SKY"})
  df.index.name = 'FIBER'
  df            = df.set_index(i*5000 + df.index.values, drop=False)

  return df

##
files           = glob.glob('bitgrid/txt/*.txt')

df              = process_df(0, files[0])

for i, _file in enumerate(files[1:]):
  df            = df.append(process_df(i, _file))

##  
print(df)
  
##
cols          = list(df.columns)
ncol          = len(cols)

summary       = pd.DataFrame(np.zeros((ncol, ncol), dtype=np.int64), columns=cols)
summary       = summary.set_index(np.array(cols).T, drop=False)

for x in cols:
  for y in cols:
    summary[x][y] = len(df[df[x] + df[y] == 2])

print(summary)

column_format = '|c' * 34 + '|'

summary       = summary.to_latex(column_format=column_format, bold_rows=False, multirow=True)
lines         = summary.split('\n')  

with open('bitgrid/tex/summary.tex', 'w') as tf:                                                                                                                                                                                                              
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
                                                                                                                                                                                                                                                                                                
  tf.write(lines[-4])                                                                                                                                                                                                                                                                        
  tf.write('\n')
  tf.write(lines[-3])
  tf.write('\n')
  tf.write(lines[-2])                                                                                                                                                                                                                                                                        
  tf.write('\n')                                                                                                                                                                                                                                                                             
  tf.write(lines[-1])           
  
print('\n\nDone.\n\n')
