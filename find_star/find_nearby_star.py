# find closest star within mag range
# FIXME DOES NOT CONSIDER MERIDAN FLIP IMPLICATIONS!

import os
import sys
import time
import numpy as np
import argparse
import logging

import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.coordinates import Angle


import matplotlib.pyplot as plt

from SAOCatalog import SAOCatalog, load_SAOCatalog_binary


if __name__ == '__main__':
    logging.basicConfig(filename='find_star.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    parser = argparse.ArgumentParser()
    parser.add_argument('cat', type=str, help='Catalog to search')
    parser.add_argument('ra2000', type=str, help='RA J2000')
    parser.add_argument('dec2000', type=str, help='DEC J2000')
    parser.add_argument('dist', type=float, help='Max distance in degrees')
    parser.add_argument('--minmag', type=float, default=8)
    parser.add_argument('--maxmag', type=float, default=7)

    args = parser.parse_args()

    logging.info(f'args = {args}')

    saocat = load_SAOCatalog_binary(args.cat)

    logging.info(f'Loaded {len(saocat.id)} stars')

    # 'fast' compute distance between points
    # first convert input ra/dec to angles
    target_str = args.ra2000 + " "
    target_str += args.dec2000
    logging.info(f"target_str = {target_str}")

    try:
        target = SkyCoord(target_str, unit=(u.hourangle, u.deg), frame='fk5', equinox='J2000')
    except ValueError:
        logging.error("Invalid RA/DEC POSITION!")
        sys.exit(1)

    # SAO Catalog I have is in J2000 so we can compare directly

    t_ra_rad = target.ra.radian
    t_dec_rad = target.dec.radian

    print(t_ra_rad, t_dec_rad)

    # convert sao coords to radians
    c_ra_rad = np.deg2rad(saocat.ra)
    c_dec_rad = np.deg2rad(saocat.dec)

    # equation for distance between two points on sphere
    # try with haversine
    ddec = t_dec_rad - c_dec_rad
    dra = t_ra_rad - c_ra_rad
    a = (np.sin(ddec/2)**2 + np.cos(t_dec_rad)*np.cos(c_dec_rad)*np.sin(dra/2)**2)
    c = 2 * np.arcsin(np.sqrt(a))

    print(np.max(c), np.min(c), np.median(c))

    #sortargs = np.argsort(c)
    close_idx = np.where(c < np.deg2rad(args.dist))[0]
    print(f'found {len(close_idx)} within {args.dist} degrees')

    # filter by mag
    mags = np.array(saocat.vmag)[close_idx]
    mags_idx = np.where((mags < args.minmag) & (mags >= args.maxmag))

    print(f'found {len(mags_idx)} within {args.dist} degrees and {args.maxmag} < mag < {args.minmag}')

#    print(close_idx)
#    print(mags_idx)

#    ids = []
#    for i in close_idx[:20]:
#        ids.append(int(saocat.id[i]))
#
##    for i in sorted(ids):
##        print(i)

    for i in close_idx[mags_idx]:
        radec = SkyCoord(saocat.ra[i], saocat.dec[i], unit=u.deg, frame='fk5', equinox='J2000')
        logging.info(f"{i:5d} {np.rad2deg(c[i]):5.2f} {saocat.id[i]:10s} {radec.to_string('hmsdms', sep=':'):30s}  {saocat.vmag[i]:4.2f}")
        logging.info(f"{i:5d} {np.rad2deg(c[i]):5.2f} {saocat.id[i]:10s} {c_ra_rad[i]:5.2f} {c_dec_rad[i]:5.2f}  {saocat.vmag[i]:4.2f}")



