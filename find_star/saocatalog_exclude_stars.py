#!/usr/bin/env python3
import os
import sys
import time
import numpy as np
import argparse
import logging

from astropy import units as u
from astropy.coordinates import SkyCoord

from SAOCatalog import load_SAOCatalog_binary

if __name__ == '__main__':
    logging.basicConfig(filename='saocatalog_exclude_stars.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    parser = argparse.ArgumentParser()
    parser.add_argument('cat', type=str, help='Catalog to process')
    parser.add_argument('--excluderad', type=float, default=0.15,
                        help='Exclusion radius in degrees')
    args = parser.parse_args()

    logging.info(f'args = {args}')

    # see if exclusion list already exists so we don't overwrite
    if os.path.isfile('sao_exclude_lst.dat'):
        logging.error(f'The file sao_exclude_lst.dat already exists!')
        logging.error('Move it out of the way before running this command or')
        sys.exit(1)

    saocat = load_SAOCatalog_binary(args.cat)
    logging.info(f'Loaded {len(saocat.id)} stars')

    near_search_rad = 1
    exclusion_rad = 0.15
    exclude_list = []
    ts = time.time()
    ndone = 0
    for cat_idx in range(0, len(saocat.id)):
        logging.debug(f'Evaluating cat index={cat_idx} ' \
                      f'SAO={saocat.id[cat_idx]:>8d}')
        #logging.info("CatIdx    SAO           RA.DEC (J2000)           VMag")

        radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx],
                         unit=u.deg, frame='fk5', equinox='J2000')
        logging.info(f"{cat_idx}   {saocat.id[cat_idx]:>8d} " \
                     f"{radec.to_string('hmsdms', sep=':', precision=3):30s} " \
                     f"{saocat.vmag[cat_idx]:4.2f}")

        # look within 1 degree for other stars - set maxmag so all brighter star considered and minmag to 2 mags fainter
        cand_idx_2, cand_dist_2 = saocat.find_stars_near_target(radec,
                                                                near_search_rad,
                                                                999, -999,
                                                                exclude=[cat_idx])
        if len(cand_idx_2) > 0:
            # figure out closest
            closest_idx = np.argmin(cand_dist_2)
            closest_deg = np.rad2deg(np.min(cand_dist_2))
            logging.info(f"      nstars={len(cand_idx_2)} closest={closest_deg} deg")
            if closest_deg <= exclusion_rad:
                exclude_list.append(cat_idx)

        ndone += 1
        if (ndone % 100) == 0:
            logging.info(f'{ndone} of {len(saocat.id)} complete')

    te = time.time()
    logging.info(f'Excluded {len(exclude_list)} stars in {te-ts:6.2f} seconds')
    f=open('sao_exclude_lst.dat', 'w')
    for exclude_idx in exclude_list:
        f.write(f'{exclude_idx}, {saocat.id[exclude_idx]}\n')
    f.close()

    sys.exit(0)

