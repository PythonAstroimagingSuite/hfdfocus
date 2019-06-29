import logging
import pickle

import numpy as np
#from astropy import units as u
#from astropy.coordinates import SkyCoord
#from astropy.coordinates import FK5
#from astropy.coordinates import Angle


def write_Tycho2_binary(obj, fname):
    f = open(fname, 'wb')
    pickle.dump(obj, f)
    f.close()

def load_Tycho2_binary(fname):
    f = open(fname, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj

class SAOCatalog:
    def __init__(self):
        self.id = []
        self.radec = []
        self.epoch = 2000.0
        self.ra = []
        self.dec = []
        self.vmag = []
        self.maxmag = None
        self.minmag = None
        self.mindec = None

    def read_catalog(self, maxmag, minmag, mindec):
        f = open('Tycho2_Vizier.tsv')
        first = True
        nread = 0
        self.maxmag = maxmag
        self.minmag = minmag
        self.mindec = mindec
        for l in f.readlines():
            if first:
                first = False
                continue
            fields = l.strip().split(',')
#            print(l.rstrip())
#            print(fields)

            if fields[5] == 'D':
                print('skipping ', fields[0])
                continue

            hid = fields[0]
            ra_deg = fields[1]
            dec_deg = fields[2]
            vmag = fields[3]
            hd = fields[4].strip('"').strip()

            try:
                vmag = float(vmag)
                if vmag >= minmag or vmag < maxmag:
                    continue
                if float(dec_deg) < mindec:
                    continue
            except Exception as err:
                logging.error(f'Error processing #1 |{l.rstrip()}| |{vmag}| {err}')
                continue

            try:
                self.id.append(int(hid))
                self.ra.append(float(ra_deg))
                self.dec.append(float(dec_deg))
                self.vmag.append(float(vmag))
            except Exception as err:
                logging.error(f'Error processing #2 |{l.rstrip()}| |{vmag}| {err}')
                continue

            nread += 1
#            if nread > 100:
#                break

        f.close()
        logging.info(f'# stars loaded = {nread}')

        i = 1
        for id, radec, vmag in zip(self.id, self.radec, self.vmag):
            logging.info(f"{i:5d} {id:10s} {radec.to_string('hmsdms', sep=':'):30s}  {vmag:4.2f}")
            i += 1

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

