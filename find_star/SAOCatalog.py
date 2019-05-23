import logging
import pickle

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.coordinates import Angle


def write_SAOCatalog_binary(obj, fname):
    f = open(fname, 'wb')
    pickle.dump(obj, f)
    f.close()

def load_SAOCatalog_binary(fname):
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
        f = open('SAO_Catalog.dat')
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

    def find_nearby(self, tradec):
        for pradec in self.radec:
            sep = pradec.separation(tradec).degree
            print(tradec, pradec, sep)

