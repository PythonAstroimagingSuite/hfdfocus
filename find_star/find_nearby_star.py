# find closest star within mag range
# FIXME DOES NOT CONSIDER MERIDAN FLIP IMPLICATIONS!

import sys
import time
import numpy as np
import argparse
import logging

from astropy import units as u
from astropy.coordinates import SkyCoord

from SAOCatalog import load_SAOCatalog_binary

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
    ch.setLevel(logging.INFO)
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
        target = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                          frame='fk5', equinox='J2000')
    except ValueError:
        logging.error("Invalid RA/DEC POSITION!")
        sys.exit(1)

    # SAO Catalog I have is in J2000 so we can compare directly
    cand_idx, cand_dist = saocat.find_stars_near_target(target, args.dist,
                                                        args.minmag, args.maxmag)

    logging.debug(f'cand_idx = {cand_idx} cand_dist = {cand_dist}')

    if len(cand_idx) == 0:
        logging.error('No stars found!')
        sys.exit(1)

    logging.info('Candidates:')
    logging.info(" CatIdx    Deg   SAO              RA.DEC (J2000)           VMag")
    for i in range(0, len(cand_idx)):
        cat_idx = cand_idx[i]
        radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx], unit=u.deg,
                         frame='fk5', equinox='J2000')
        logging.info(f" {cat_idx}   {np.rad2deg(cand_dist[i]):5.2f} " \
                     f"{saocat.id[cat_idx]:10s} " \
                     f"{radec.to_string('hmsdms', sep=':', precision=3):30s} " \
                     f"{saocat.vmag[cat_idx]:4.2f}")
    logging.info("")

    nnear = []
    nnear_idx = []
    nnear_dist = []
    near_search_rad = 1
    exclusion_rad = 0.15
    for i in range(0, len(cand_idx)):
        cat_idx = cand_idx[i]
        logging.debug(f'Evaluating cat index={cat_idx} SAO={saocat.id[cat_idx]:10s}')
        logging.info("CatIdx    Deg   SAO           RA.DEC (J2000)           VMag")

        radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx], unit=u.deg,
                         frame='fk5', equinox='J2000')
        logging.info(f"{cat_idx}  {np.rad2deg(cand_dist[i]):5.2f} " \
                     f"{saocat.id[cat_idx]:10s} "\
                     f"{radec.to_string('hmsdms', sep=':', precision=3):30s} " \
                     f"{saocat.vmag[cat_idx]:4.2f}")

        # look within 1 degree for other stars - set maxmag so all brighter star considered and minmag to 2 mags fainter
        cand_idx_2, cand_dist_2 = saocat.find_stars_near_target(radec,
                                                                near_search_rad,
                                                                999, -999,
                                                                exclude=[cat_idx])
        logging.debug(f'cand_idx_2 = {cand_idx_2} cand_dist_2 = {cand_dist_2}')

        logging.info('Nearby Candidates:')
        logging.info(" CatIdx    Deg   SAO              RA.DEC (J2000)           VMag")
        for i2 in range(0, len(cand_idx_2)):
            cat_idx_2 = cand_idx_2[i2]
            print(i2, cat_idx_2)
            radec = SkyCoord(saocat.ra[cat_idx_2], saocat.dec[cat_idx_2], unit=u.deg, frame='fk5', equinox='J2000')
            logging.info(f" {cat_idx_2}   {np.rad2deg(cand_dist_2[i2]):5.2f} " \
                         f"{saocat.id[cat_idx_2]:10s} " \
                         f"{radec.to_string('hmsdms', sep=':', precision=3):30s} " \
                         f"{saocat.vmag[cat_idx_2]:4.2f}")

        # figure out closest
        closest_idx = np.argmin(cand_dist_2)
        print(i, closest_idx, cand_idx_2[closest_idx], cand_dist_2[closest_idx])
        closest_deg = np.rad2deg(np.min(cand_dist_2))
        logging.info(f"      nstars={len(cand_idx_2)} closest={closest_deg} deg")

        if closest_deg > exclusion_rad:
            nnear.append(len(cand_idx_2))
            nnear_idx.append(cat_idx)
            nnear_dist.append(cand_dist[i])
        else:
            logging.info(f'Excluding star due to neighbor within {exclusion_rad} degrees.')

        logging.info("")

    logging.debug(f'nnear={nnear} nnear_idx={nnear_idx} nnear_dist={nnear_dist}')


    if len(nnear) < 1:
        logging.error('No star found!')
        sys.exit(1)

    # find cand with least nearby
    least_idx = np.argmin(nnear)
    best_idx = nnear_idx[least_idx]
    radec = SkyCoord(saocat.ra[best_idx], saocat.dec[best_idx],
                     unit=u.deg, frame='fk5', equinox='J2000')
    logging.info(f'BEST STAR = {best_idx} SAO{saocat.id[best_idx]:10s} ' \
                 f'VMag = {saocat.vmag[best_idx]:4.2f} # nearby = {nnear[least_idx]}')
    logging.info("CatIdx    Deg   SAO           RA.DEC (J2000)           VMag")
    logging.info(f"{best_idx}  {np.rad2deg(nnear_dist[least_idx]):5.2f} " \
                 f"{saocat.id[best_idx]:10s} " \
                 f"{radec.to_string('hmsdms', sep=':', precision=3):30s} " \
                 f"{saocat.vmag[best_idx]:4.2f}")


    #logging.info(f"{np.rad2deg(nnear_dist[least_idx]):5.2f} {saocat.id[best_idx]:10s} {radec.to_string('hmsdms', sep=':'):30s}  {saocat.vmag[best_idx]:4.2f}")

    cand_idx_2, cand_dist_2 = saocat.find_stars_near_target(radec, 1, 999, -999, exclude=[best_idx])
    logging.debug(f'cand_idx_2 = {cand_idx_2} cand_dist_2 = {cand_dist_2}')
    logging.info(f'List of stars closest to BEST STAR')
    logging.info(" CatIdx    Deg   SAO              RA.DEC (J2000)           VMag")
    for i2 in range(0, len(cand_idx_2)):
        cat_idx_2 = cand_idx_2[i2]
        radec = SkyCoord(saocat.ra[cat_idx_2], saocat.dec[cat_idx_2],
                         unit=u.deg, frame='fk5', equinox='J2000')
        logging.info(f" {cat_idx_2}   {np.rad2deg(cand_dist_2[i2]):5.2f} " \
                     f"{saocat.id[cat_idx_2]:10s} " \
                     f"{radec.to_string('hmsdms', sep=':', precision=3):30s} " \
                     f"{saocat.vmag[cat_idx_2]:4.2f}")

    distsort_idx = np.argsort(nnear)
    logging.info("")
    logging.info(f'Candidates within {args.dist} degrees of J2000 ' \
                 f'{target.to_string("hmsdms", sep=":"):30s}')
    logging.info('CatIdx    Deg     SAO           RA.DEC (J2000)            ' \
                 'VMag     # near stars')
    for i in range(0, len(distsort_idx)):
        near_idx = distsort_idx[i]
        cat_idx = nnear_idx[near_idx]
        radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx],
                         unit=u.deg, frame='fk5', equinox='J2000')
        logging.info(f"{cat_idx}  {np.rad2deg(nnear_dist[near_idx]):5.2f} " \
                     f"{saocat.id[cat_idx]:>8s}   " \
                     f"{radec.to_string('hmsdms', sep=':', precision=3):30s} "\
                     f"{saocat.vmag[cat_idx]:4.2f}        {nnear[near_idx]}")



