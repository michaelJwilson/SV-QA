import numpy         as     np

from   astropy.table import Column, Table
from   in_des        import in_des


def set_photsys(_dat, verbose=True):
    ##  Rewrite to allow for longer strings.
    cols                  =  _dat.columns
    dat                   =  Table(_dat, names=cols, copy=True)

    dat['PHOTSYS']        =  Column(data=np.array(dat['PHOTSYS']), name='PHOTSYS', dtype='S32')

    isin                  =  dat['PHOTSYS'] == 'N'
    dat['PHOTSYS'][isin]  = 'BMZLS'

    isin                  =  (dat['PHOTSYS'] == 'S') & (dat['RA'].quantity.value < 300.) & (dat['RA'].quantity.value > 100.)
    dat['PHOTSYS'][isin]  = 'DECALS-NGC'

    isin                  =  ~isin & (dat['PHOTSYS'] == 'S')
    dat['PHOTSYS'][isin]  = 'DECALS-SGC'

    indes                 = in_des(dat['RA'].quantity.value, dat['DEC'].quantity.value, verbose=verbose)
    dat['PHOTSYS'][indes] = 'DES'

    return  dat
