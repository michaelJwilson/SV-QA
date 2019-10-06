import os

##
nside = 128

runs  = []

for camera in ['decam', '90prime', 'mosaic']:
    for band in ['g', 'r', 'z']:
        runs.append([camera, band])

        if camera == 'mosaic':
          break
        
runs  = runs[::-1]

for run in runs:
  cmd = 'python legacysky_randoms.py {} {} {}'.format(nside, run[0], run[1])
    
  print(cmd)
  
  os.system(cmd)