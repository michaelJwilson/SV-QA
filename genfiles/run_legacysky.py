import os

##
nside = 128

for camera in ['decam', '90prime', 'mosaic']:
    for band in ['g', 'r', 'z']:
        cmd = 'python legacysky_randoms.py {} {} {}'.format(nside, camera, band)
        
        print(cmd)

        os.system(cmd)
