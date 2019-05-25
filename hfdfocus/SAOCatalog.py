import sys
import logging
from struct import pack, unpack, calcsize

import numpy as np
#from astropy import units as u
#from astropy.coordinates import SkyCoord
#from astropy.coordinates import FK5
#from astropy.coordinates import Angle

#def write_SAOCatalog_binary(obj, fname):
#    f = open(fname, 'wb')
#    pickle.dump(obj, f)
#    f.close()
#
#def load_SAOCatalog_binary(fname):
#    f = open(fname, 'rb')
#    obj = pickle.load(f)
#    f.close()
#    return obj

#
# Binary format (preliminary!)
#
# nrec      (uint32)
# minmag    (float32)
# maxmag    (float32)
# mindec    (float32)
# epoch     (float32)
# SAO ids   (nrec*int32)
# RA (deg)  (nrec*float32)
# DEC (deg) (nrec*float32)
# Vmag      (nrec*float32)
#
def write_SAOCatalog_binary(saocat, fname):
    f = open(fname, 'wb')

    # make sure all lists the same length!
    nelem = []
    nelem.append(len(saocat.id))
    nelem.append(len(saocat.ra))
    nelem.append(len(saocat.dec))
    nelem.append(len(saocat.vmag))

    # just count first elem of nelem and if it is equal to len
    # then all elem are the same
    if nelem.count(nelem[0]) != len(nelem):
        logging.error('write_SAOCatalog_binary(): Error - not all lists the same length!')
        return None

    nrec = nelem[0]
    f.write(pack('I', nrec))
    f.write(pack('f', saocat.minmag))
    f.write(pack('f', saocat.maxmag))
    f.write(pack('f', saocat.mindec))
    f.write(pack('f', saocat.epoch))
    f.write(pack(f'{nrec}I', *saocat.id))
    f.write(pack(f'{nrec}f', *saocat.ra))
    f.write(pack(f'{nrec}f', *saocat.dec))
    f.write(pack(f'{nrec}f', *saocat.vmag))
    f.close()


def load_SAOCatalog_binary(fname):

    def unpack_next(f, fmt):
        #print(fmt)
        nr = calcsize(fmt)
        b = f.read(nr)
        return unpack(fmt, b)

    saocat = SAOCatalog()

    f = open(fname, 'rb')
    nrec = unpack_next(f, 'I')[0]
    #print(nrec)
    saocat.minmag = unpack_next(f, 'f')[0]
    #print(saocat.minmag)
    saocat.maxmag = unpack_next(f, 'f')[0]
    #print(saocat.maxmag)
    saocat.mindec = unpack_next(f, 'f')[0]
    #print(saocat.mindec)
    saocat.epoch = unpack_next(f, 'f')[0]
    #print(saocat.epoch)
    saocat.id = list(unpack_next(f, f'{nrec}I'))
    saocat.ra = list(unpack_next(f, f'{nrec}f'))
    saocat.dec = list(unpack_next(f, f'{nrec}f'))
    saocat.vmag = list(unpack_next(f, f'{nrec}f'))
    f.close()

    # make sure all lists the same length!
    nelem = []
    nelem.append(len(saocat.id))
    nelem.append(len(saocat.ra))
    nelem.append(len(saocat.dec))
    nelem.append(len(saocat.vmag))
    #print(nelem)

    # just count first elem of nelem and if it is equal to len
    # then all elem are the same
    if nelem.count(nelem[0]) != len(nelem):
        logging.error('load_SAOCatalog_binary(): Error - not all lists the same length!')
        return None
    else:
        return saocat

class SAOCatalog:
    def __init__(self):
        self.id = []
        self.epoch = 2000.0
        self.ra = []
        self.dec = []
        self.vmag = []
        self.maxmag = None
        self.minmag = None
        self.mindec = None

    def find_stars_near_target(self, target, dist, minmag, maxmag, exclude=None):

        t_ra_rad = target.ra.radian
        t_dec_rad = target.dec.radian

        #print(t_ra_rad, t_dec_rad)

        # convert sao coords to radians
        c_ra_rad = np.deg2rad(self.ra)
        c_dec_rad = np.deg2rad(self.dec)

        # equation for distance between two points on sphere
        # try with haversine
        ddec = t_dec_rad - c_dec_rad
        dra = t_ra_rad - c_ra_rad
        a = (np.sin(ddec/2)**2 + np.cos(t_dec_rad)*np.cos(c_dec_rad)*np.sin(dra/2)**2)
        c = 2 * np.arcsin(np.sqrt(a))

        #print(np.max(c), np.min(c), np.median(c))

        #sortargs = np.argsort(c)
        close_idx = np.where(c < np.deg2rad(dist))[0]

        if exclude is not None:
            #logging.info(f'exclude={exclude}')
            #logging.info(f'pre-exclude: {close_idx}')
            close_idx = np.setxor1d(close_idx, exclude)
            #logging.debug(f'post-exclude: {close_idx}')

        logging.debug(f'found {len(close_idx)} within {dist} degrees')

        # filter by mag
        mags = np.array(self.vmag)[close_idx]
        mags_idx = np.where((mags < minmag) & (mags >= maxmag))[0]

        logging.debug(f'found {len(mags_idx)} within {dist} degrees and {maxmag} < mag < {minmag}')
        #logging.debug(f'mags_idx = {mags_idx}')
        ret_idx = close_idx[mags_idx]
        #logging.debug(f'ret_idx = {ret_idx}')
        return ret_idx, c[ret_idx]

