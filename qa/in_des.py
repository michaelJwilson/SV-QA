import  os
import  pymangle
import  numpy     as  np
import  healpy    as  hp


def prep_indes(nside=256, fname=None):
    if fname  == None:
      ##  '/global/cscratch1/sd/raichoor/desits/des.ply' 
      fname    = os.getenv('CSCRATCH') + '/BGS/SV-ASSIGN/des/des.ply'

    ##                                                                                                                                                                                                                         
    nside      = np.int(nside)
    npix       = hp.nside2npix(nside)

    # checking hp pixels                                                                                                                                                                                                        
    mng        = pymangle.Mangle(fname)

    theta, phi = hp.pix2ang(nside, np.arange(0, npix, 1), nest=False)
    hpra,hpdec = 180. / np.pi * phi, 90. -180. / np.pi * theta

    allhp      = mng.polyid(hpra, hpdec)

    hpindes    = (allhp != -1).astype(int)

    # pixels with all neighbours in des.                                                                                                                                                                                       
    hpindes_secure  = np.array([i for i in range(npix) if hpindes[i] + hpindes[hp.get_all_neighbours(nside, i)].sum() == 9])

    # pixels with all neighbours outside des.                                                                                                                                                                                  
    hpoutdes_secure = np.array([i for i in range(npix) if hpindes[i] + hpindes[hp.get_all_neighbours(nside, i)].sum() == 0])

    # hpind to be checked                                                                                                                                                                                                      
    tmp                  = np.ones(npix, dtype=bool)
    tmp[hpindes_secure]  = False
    tmp[hpoutdes_secure] = False

    hp_tbc    = np.array(range(npix))[tmp]

    np.savetxt(os.getenv('CSCRATCH') + '/BGS/SV-ASSIGN/des/hpindes_secure.txt', hpindes_secure, fmt='%d')   
    np.savetxt(os.getenv('CSCRATCH') + '/BGS/SV-ASSIGN/des/hp_tbc.txt', hp_tbc, fmt='%d')
    

def in_des(ra, dec, nside=256, fname=None, verbose=True):
  if fname  == None:
    ##  '/global/cscratch1/sd/raichoor/desits/des.ply'  
    fname    = os.getenv('CSCRATCH') + '/BGS/SV-ASSIGN/des/des.ply'

  try:
    if verbose:
      print('Loading pre-calculated pixels.')
        
    hpindes_secure = np.loadtxt(os.getenv('CSCRATCH') + '/BGS/SV-ASSIGN/des/hpindes_secure.txt')
    hp_tbc         = np.loadtxt(os.getenv('CSCRATCH') + '/BGS/SV-ASSIGN/des/hp_tbc.txt')

    mng            = pymangle.Mangle(fname)
        
  except:
    if verbose:
      print('Loading failed;  Re-calculate pixels.')

    ##
    nside      = np.int(nside)
    npix       = hp.nside2npix(nside)
    
    # checking hp pixels
    mng        = pymangle.Mangle(fname)

    theta,phi  = hp.pix2ang(nside, np.arange(0, npix, 1), nest=False)
    hpra,hpdec = 180. / np.pi * phi, 90. -180. / np.pi * theta

    allhp      = mng.polyid(hpra, hpdec)
    
    hpindes    = (allhp != -1).astype(int)

    # print(allhp, np.unique(hpindes))
    
    # pixels with all neighbours in des.
    hpindes_secure  = np.array([i for i in range(npix) if hpindes[i] + hpindes[hp.get_all_neighbours(nside, i)].sum() == 9])

    # pixels with all neighbours outside des.
    hpoutdes_secure = np.array([i for i in range(npix) if hpindes[i] + hpindes[hp.get_all_neighbours(nside, i)].sum() == 0])

    # hpind to be checked
    tmp                  = np.ones(npix, dtype=bool)
    tmp[hpindes_secure]  = False
    tmp[hpoutdes_secure] = False
    
    hp_tbc    = np.array(range(npix))[tmp]

  # now checking indiv. obj. in the tbc pixels
  hppix     = hp.ang2pix(nside, (90. - dec) * np.pi / 180., ra * np.pi / 180., nest=False)
  hpind     = np.unique(hppix)
    
  isdes     = np.zeros(len(ra), dtype=bool)
  isdes[np.in1d(hppix, hpindes_secure)] = True

  tbc       = np.where(np.in1d(hppix, hp_tbc))[0]

  if len(tbc > 0):
    tbcisdes  = (mng.polyid(ra[tbc], dec[tbc]) != -1)
    isdes[tbc][tbcisdes] = True
    
  return  isdes


if __name__ == '__main__':
    ra    = np.array([260., 20., 350., 50.])
    dec   = np.array([20, -10, -40., -40.])

    ##  prep_indes()
    
    isdes = in_des(ra, dec)

    print(isdes)
                     
    print('\n\nDone.\n\n')
