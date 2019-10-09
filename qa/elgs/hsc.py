import  os
import  sys
import  glob
import  corner
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
from    astropy.table                         import   Table, vstack, hstack, join, Column
from    desitarget.targets                    import   encode_targetid
from    desitarget.geomask                    import   is_in_box
from    desitarget.targetmask                 import   desi_mask
from    mpl_toolkits.axes_grid1               import   make_axes_locatable
from    mpl_toolkits.axes_grid1.inset_locator import   inset_axes
from    desitarget.targets                    import   encode_targetid
from    phot_sys                              import   set_photsys
from    in_des                                import   in_des


rc('font', **{'family':'serif', 'serif':['Times']})
rc('text', usetex=True)

##
##  nside    = np.int(sys.argv[1])
##  parea    = hp.nside2pixarea(nside, degrees = True)  
    
if __name__ == '__main__':    
 scratch   = os.environ['CSCRATCH']
 truthdir  = '/project/projectdirs/desi/target/analysis/truth/dr8.0/'

 gcap      = 'north'
 snr       =      9
 compute   =  False
 
 if compute:
  ##
  elgs      = Table(fits.open('/global/cscratch1/sd/mjwilson/BGS/SV-ASSIGN/elgs/elgs.fits')[1].data)                                                                                                                                   
  mtypes    = np.unique(elgs['MORPHTYPE'])    
  
  ## 
  ##  _all  = ['hsc_pdr1_wide.forced.reduced-match.fits', 'hsc_pdr1_deep.forced.reduced-match.fits', 'hsc_pdr1_udeep.forced.reduced-match.fits']
  _all      = ['hsc_pdr1_wide.forced.reduced-match.fits']

  cols      = ['RA', 'DEC', 'frankenz_photoz_best', 'SURVEY', 'TARGETID']
  
  ##  Writes desired truth catalogues to scratch.                                                                                                                                                                                      
  for i, survey in enumerate(_all):
      infile     = truthdir + '/{}/matched/'.format(gcap) + survey

      if not os.path.exists(infile):
        print('Cannot find: {}'.format(infile))                                                                                                                                                                                     
        continue

      ##  Truth spectroscopic catalogue.  Look for a given survey in both north and south.                                                                                                                                             
      hsc        = Table(fits.open(truthdir + '/{}/matched/'.format(gcap) + survey)[1].data)
      hsc        = hsc['mizuki_photoz_best', 'frankenz_photoz_best']
      
      ##  Line matched LS output.                                                                                                                                                                                                       
      infile     = truthdir + '/{}/matched/'.format(gcap) + 'ls-dr8.0-' + survey

      if not os.path.exists(infile):
        print('Cannot find: {}'.format(infile))                                                                                                                                                                                     
        continue

      lsm        = Table(fits.open(truthdir + '/{}/matched/'.format(gcap) + 'ls-dr8.0-' + survey)[1].data)

      ##  hsc_pdr1_wide.forced.reduced
      depth      = survey.split('-')[0].replace('.forced.reduced', '').split('_')[-1].upper()
      nsurvey    = 'HSC1-' + depth
            
      ##
      lsm['TARGETID'] = encode_targetid(objid=lsm['OBJID'], brickid=lsm['BRICKID'], release=lsm['RELEASE'], sky=0, mock=0)
      lsm['SURVEY']   = depth

      hsc        = hstack([lsm, hsc])
      
      hsc['FRANKENZ'] = hsc['frankenz_photoz_best']

      del  hsc['frankenz_photoz_best']
            
      print('Solving for {}.'.format(nsurvey))

      ##  hsc.pprint()

  cols  =  hsc.columns
  ecols = elgs.columns
  
  diff  = list(set(cols)  ^ set(ecols))
  ecols = list(set(ecols) & set(diff))
  ecols = ecols + ['TARGETID']

  ##
  ##  hsco = join(hsc, elgs[ecols], join_type='left', keys='TARGETID')
  hsco  = hsc
  
  hsco.pprint(max_width=-1)
  hsco.write(scratch + '/BGS/SV-ASSIGN/elgs/hsc_{}.fits'.format(gcap), format='fits', overwrite=True)
  
  exit(1)
  
  ##
  elgs = join(elgs, hsc, join_type='inner', keys='TARGETID')
  elgs.write(scratch + '/BGS/SV-ASSIGN/elgs/hsc_elgs.fits', format='fits', overwrite=True)

  
 else:
  snr              = 5.
  
  ##
  fig, axarr = plt.subplots(figsize=(10, 5))
  plt.subplots_adjust(left = 0.05, right = 0.95, hspace=0.6, wspace=0.4, top = 0.925, bottom = 0.05)

  cols             = ['FRANKENZ', 'MORPHTYPE', 'GSNR', 'RSNR', 'ZSNR', 'RZSNR', 'SNR', 'EBV', 'CHI2DIFF', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1',\
                      'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'MIZUKIZ', 'PHOTSYS', 'BFIT', 'BFIT2']
  
  hsc              = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_{}.fits'.format(gcap))[1].data)
  hsc              = hsc[(hsc['FRANKENZ'] >= 0.5) & (hsc['FRANKENZ'] < 1.6)]
  hsc              = hsc[hsc['DEC'] > -30.]
  hsc['MORPHTYPE'] = hsc['TYPE']
  hsc['MIZUKIZ']   = hsc['mizuki_photoz_best']

  hsc.pprint()

  exit(1)
  
  ##
  hsc['PHOTSYS']                      = Column(data=np.array(['N'] * len(hsc)), name='PHOTSYS', dtype='S32')
  hsc['PHOTSYS'][hsc['DEC'] < 32.375] = 'S'
  
  set_photsys(hsc)
    
  ##  S/N should be prior to extinction correction.                                                                                                                                                                         
  hsc['GSNR']  = np.sqrt(hsc['FLUX_IVAR_G']) * hsc['FLUX_G']                                                                                                                               
  hsc['RSNR']  = np.sqrt(hsc['FLUX_IVAR_R']) * hsc['FLUX_R']                                                                                                                                 
  hsc['ZSNR']  = np.sqrt(hsc['FLUX_IVAR_Z']) * hsc['FLUX_Z']                                                                                                                               

  ##                                                                                                                                                                                                                        
  hsc['RZSNR']    = np.sqrt(hsc['RSNR'] ** 2. + hsc['ZSNR'] ** 2.)
  hsc['SNR']      = np.sqrt(hsc['GSNR'] ** 2. + hsc['RSNR'] ** 2. + hsc['ZSNR'] ** 2.)

  hsc['BFIT']     = np.array([np.array(x).max() for x in hsc['DCHISQ']])
  hsc['BFIT2']    = np.array([np.sort(np.array(x))[:-1].max() for x in hsc['DCHISQ']])

  hsc['CHI2DIFF'] = np.clip(np.log10(hsc['BFIT'] - hsc['BFIT2']), a_min=-6., a_max=None)
  
  hsc['BFIT']     = np.clip(np.log10(hsc['BFIT']),  a_min=-6., a_max=None)
  hsc['BFIT2']    = np.clip(np.log10(hsc['BFIT2']), a_min=-6., a_max=None)
  
  ##
  hsc             = hsc[hsc['SNR'] >= snr]
  hsc             = hsc[cols]

  hsc.write(scratch + '/BGS/SV-ASSIGN/elgs/hsc_{}_lite.fits'.format(gcap), format='fits', overwrite=True)
  

  hsc             = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_{}_lite.fits'.format(gcap))[1].data) 
  
  ##  Downsample.
  for i in range(2):
    choice        = np.random.randint(2, size=len(hsc)).astype(bool)
    hsc           = hsc[choice]  
  
  hsc.sort('SNR')
  
  hsc.pprint()

  ##
  mtypes          = ['PSF', 'EXP', 'DEV', 'COMP', 'REX']
  colors          = ['b', 'r', 'g', 'k', 'y']

  photsyss        = ['BMZLS'] ##  np.unique(hsc['PHOTSYS'])
  
  ##
  for photsys in photsyss:
    datbysys        = hsc[hsc['PHOTSYS'] == photsys]

    ##  
    fig, axarr      = plt.subplots(nrows=7, ncols=7, figsize=(20, 20))
    
    for mm, color in zip(mtypes, colors):
      print('Solving for {} ({})'.format(mm, color))

      dat           = datbysys[datbysys['MORPHTYPE'] == mm]
      dat.sort('SNR')
      
      toprint       = dat[:500]
      toprint.pprint(max_lines=-1)

      if len(dat) > 0:       
       ##  np.log10(dat['BFIT']), dat['CHI2DIFF']].   
       data          = np.c_[dat['GSNR'], dat['RSNR'], dat['ZSNR'], dat['SNR'], dat['FRANKENZ'], dat['PSFSIZE_Z'], dat['CHI2DIFF']]

       ##  ['BFIT', r'dChi2'].    
       labels        = [r'gSNR', r'rSNR', r'zSNR', r'SNR', r'FRANKEN-z', 'zPSF', 'CHI2DIFF']
      
       ##  ['BFIT', r'dChi2'].    
       _range        = [(-5., 50.)] * 4 + [(0.4, 1.6)] + [(0.8, 1.5)] + [(-6., 6.)]

       ##  [0.05, 0.1, 0.25, 0.5, 0.75]
       percentiles   = []

       ##  data_kwargs={'alpha': 0.2}.
       fig           = corner.corner(data, plot_contours=False, color=color, labels=labels, range=_range, 
                                     quantiles=percentiles, plot_density=False, fig=fig, hist_kwargs={'log': True},
                                     show_titles=True)
        
    fig.suptitle(r'0.5 $\leq z \leq 1.6$ HSC-PDR1 in {}'.format(photsys), fontsize=20)
    
    ##
    pl.savefig('hsc_{}.png'.format(photsys))
  
  '''
  fig, axarr    = plt.subplots(nrows=1, ncols=3, figsize=(10, 5))

  elgs          = Table(fits.open(scratch + '/BGS/SV-ASSIGN/elgs/hsc_elgs.fits')[1].data)
  elgs['RZSNR'] = np.sqrt(elgs['RSNR'] ** 2. + elgs['ZSNR'] ** 2.)
  
  lowq          = elgs[(elgs['RSNR']  < snr) | (elgs['ZSNR']  < snr)]
  hiq           = elgs[(elgs['RSNR'] >= snr) & (elgs['ZSNR'] >= snr)]

  lowq.pprint(max_width=-1)
  hiq.pprint(max_width=-1)

  ##
  dz         = 0.05
  bins       = np.arange(0.0, 3.5, dz)
  
  axarr[0].hist(lowq['FRANKENZ'], bins=bins, color='darkviolet', label='LOQ')
  axarr[0].hist( hiq['FRANKENZ'], bins=bins, color='royalblue',  label='HIQ')

  axarr[0].set_xlim(0., 3.5)
  
  axarr[0].set_xlabel(r'$z$')
  axarr[0].set_ylabel(r'$\Delta N / \Delta z$')

  axarr[0].legend(frameon=False)
  
  ##
  for ii, [subsample, color, label] in enumerate(zip([lowq, hiq], ['darkviolet', 'royalblue'], ['LOWQ', 'HIQ'])):
   gg = subsample['FLUX_G'] / subsample['MW_TRANSMISSION_G']
   rr =	subsample['FLUX_R'] / subsample['MW_TRANSMISSION_R']
   zz =	subsample['FLUX_Z'] / subsample['MW_TRANSMISSION_Z']
   
   axarr[ii + 1].plot(rr-zz, gg-rr, marker='.', lw=0, c=color, markersize=.5)

   axarr[ii + 1].set_xlim(-2., 0.)
   axarr[ii + 1].set_ylim(-.5, 2.)

   axarr[ii + 1].set_xlabel(r'$(r-z)$') 
   axarr[ii + 1].set_ylabel(r'$(g-r)$')
   
  ##  
  plt.tight_layout()
  
  pl.savefig('hsc_elgs.png')
  '''  

  print('\n\nDone.\n\n')
