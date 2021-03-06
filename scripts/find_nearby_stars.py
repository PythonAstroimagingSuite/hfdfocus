#!/usr/bin/env python3
#
# find closest star within mag range
#
# Copyright 2020 Michael Fulbright
#
#
#    hfdfocus is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
import os
import sys
import numpy as np
import argparse
import logging
from pathlib import Path

print("find_nearby_stars.py: DISABLED IERS AGE CHECK!")
from astropy.utils.iers import conf
#conf.iers_auto_url = "https://datacenter.iers.org/data/9/finals2000A.all"
conf.iers_auto_url = "ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all"
conf.auto_max_age = None


from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

import hfdfocus
from hfdfocus.SAOCatalog import load_SAOCatalog_binary

# default catalog
CATALOG_NAME = "data/SAO_Catalog_m5_p11_filtered.bin"

if __name__ == '__main__':
    LONG_FORMAT = '%(asctime)s.%(msecs)03d [%(filename)20s:%(lineno)3s - ' \
                  + '%(funcName)20s() ] %(levelname)-8s %(message)s'
    SHORT_FORMAT = '%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s'

    logging.basicConfig(filename='find_star.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format=LONG_FORMAT, #format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    parser = argparse.ArgumentParser()
    parser.add_argument('ra2000', type=str, help='RA J2000 HH:MM:SS')
    parser.add_argument('dec2000', type=str, help='DEC J2000 DD:MM:SS')
    parser.add_argument('dist', type=float, help='Max distance in degrees')
    parser.add_argument('--cat', type=str, help='Override catalog to search')
    parser.add_argument('--minmag', type=float, default=8, help='Faintest mag allowed')
    parser.add_argument('--maxmag', type=float, default=7, help='Brightest mag allowed')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--outfile', type=str, help='Output file with candidates')
    parser.add_argument('--force', action='store_true', help='Overwrite output file')
    parser.add_argument('--lst', type=str, help='Local sidereal time')
    parser.add_argument('--onlyside', type=str, help='EAST or WEST side only')
    parser.add_argument('--currentside', action='store_true',
                        help='Only look on current side of pier')
    parser.add_argument('--meridianthres', type=str, default='00:30:00',
                        help='How close to meridian is allowed (hh:mm:ss)')
    parser.add_argument('--lon', type=float, help='Location longitude')
#    parser.add_argument('--lat', type=float, help='Location latitude')
#    parser.add_argument('--tz', type=str, help='Location tz')
#    parser.add_argument('--allowflip', action='store_true', help='Allow meridian flip to find star')

    args = parser.parse_args()

    logging.debug(f'args = {args}')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter(LONG_FORMAT)
    ch = logging.StreamHandler()
    if args.debug:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO
    ch.setLevel(loglevel)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    # complain if output file exists
    if args.outfile is not None:
        if os.path.isfile(args.outfile):
            if not args.force:
                logging.error(f'Output file {args.outfile} already exists - '
                              'please remove before running')
                sys.exit(1)
            else:
                logging.debug(f'Removing existing output file {args.outfile}')
                os.unlink(args.outfile)

    # first convert input ra/dec to angles
    target_str = args.ra2000 + " "
    target_str += args.dec2000
    logging.debug(f"target_str = {target_str}")

    try:
        target = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                          frame='fk5', equinox='J2000')
    except ValueError:
        logging.error("Invalid RA/DEC POSITION!")
        sys.exit(1)


    # arguments for controlling filtering by which side of meridian are a
    # bit convoluted and conflict
    #
    # to sort it out for now only allow
    #
    # --lst --onlyside
    #
    # OR
    #
    # --lon --onlyside
    #
    meridian_thres = None
    local_sidtime = None
    force_side = None
    if args.lst:
        lst_ok = False
        if args.lst is not None:
            try:
                local_sidtime = Angle(args.lst, unit=u.hour)
            except ValueError:
                logging.error("Invalid LST!")
                sys.exit(1)
            lst_ok = True

        if not lst_ok:
            logging.error('Cannot apply pier side filtering with the supplied LST')
            sys.exit(1)

        try:
            meridian_thres = Angle(args.meridianthres, unit=u.hour)
        except ValueError:
            logging.error("Invalid meridan threshold")
            sys.exit(1)
    else:
        # should have lat/lon/tz
        obs_lon = None

        if args.lon is not None:
            obs_lon = args.lon

        if obs_lon is not None:
            # compute local sidereal time and current pier side
            now = Time.now()
            local_sidtime = now.sidereal_time('apparent', obs_lon * u.degree)
            time_fmt = '%H:%M:%S'
            logging.debug(f'Local sidereal time is {local_sidtime.hms}')

    # if a sid time is given then get rest of args
    if local_sidtime is not None:
        try:
            meridian_thres = Angle(args.meridianthres, unit=u.hour)
        except ValueError:
            logging.error("Invalid meridan threshold")
            sys.exit(1)

        # figure out side requested position is on
        logging.debug(f'{args.onlyside} {args.currentside}')
        if args.onlyside is not None and args.currentside:
            logging.error('Can only specify one of --onlyside and '
                          '--currentside!')
            sys.exit(1)
        if args.onlyside is not None:
            if args.onlyside not in ['EAST', 'WEST']:
                logging.error('--onlyside must specify EAST or WEST only')
                sys.exit(1)
            force_side = args.onlyside
        elif args.currentside is not None:
            # compute current side using HA
            hour_angle = (local_sidtime - target.ra).wrap_at('180d')
            logging.debug(f'target hour angle={hour_angle}')
            if abs(hour_angle.hour) > 6:
                logging.error('Target is too far from meridian! '
                              f'Hour angle = {hour_angle}')
                sys.exit(1)
            elif hour_angle < 0:
                logging.debug('Target is in EAST')
                force_side = 'EAST'
            else:
                logging.debug('Target is in the WEST')
                force_side = 'WEST'

    logging.debug(f'Using pier side constraint of {force_side} ')
    if local_sidtime is not None:
        logging.debug(f'Using LST = {local_sidtime.hms}')

    if meridian_thres is not None:
        logging.debug(f'Using meridian threshold of '
                      '{meridian_thres.to_string(unit=u.hour)}')

    # see if the catalog was overridden
    if args.cat is not None:
        catalog_path = args.cat
    else:
        # sniff location from package location
        catalog_path = Path(hfdfocus.__file__).parent.joinpath(CATALOG_NAME)

    saocat = load_SAOCatalog_binary(catalog_path)
    logging.debug(f'Loaded {len(saocat.id)} stars')

    # 'fast' compute distance between points

    # SAO Catalog I have is in J2000 so we can compare directly
    cand_idx, cand_dist = saocat.find_stars_near_target(target, args.dist,
                                                        args.minmag, args.maxmag)

    logging.debug(f'cand_idx = {cand_idx} cand_dist = {cand_dist}')

    if len(cand_idx) == 0:
        logging.error('No stars found!')
        sys.exit(1)

    if args.verbose:
        logging.debug('Candidates:')
        logging.debug(' CatIdx    Deg   SAO              RA.DEC (J2000)           VMag')
        for i in range(0, len(cand_idx)):
            cat_idx = cand_idx[i]
            radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx], unit=u.deg,
                             frame='fk5', equinox='J2000')
            logging.debug(f' {cat_idx}   {np.rad2deg(cand_dist[i]):5.2f} '
                          f'{saocat.id[cat_idx]:>8d} '
                          f'{radec.to_string("hmsdms", sep=":", precision=3):30s} '
                          f'{saocat.vmag[cat_idx]:4.2f}')
        logging.debug('')

    nnear = []
    nnear_idx = []
    nnear_dist = []
    near_search_rad = 1
    exclusion_rad = 0.15
    nprox_exclude = 0
    npier_exclude = 0
    for i in range(0, len(cand_idx)):
        cat_idx = cand_idx[i]
        if args.verbose:
            logging.debug(f'Evaluating cat index={cat_idx} SAO={saocat.id[cat_idx]:>8d}')
            logging.debug('CatIdx    Deg   SAO           RA.DEC (J2000)           VMag')

        radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx], unit=u.deg,
                         frame='fk5', equinox='J2000')
        if args.verbose:
            logging.debug(f'{cat_idx}  {np.rad2deg(cand_dist[i]):5.2f} '
                          f'{saocat.id[cat_idx]:>8d} '
                          f'{radec.to_string("hmsdms", sep=":", precision=3):30s} '
                          f'{saocat.vmag[cat_idx]:4.2f}')

        # look within 1 degree for other stars - set maxmag so all brighter star considered and minmag to 2 mags fainter
        cand_idx_2, cand_dist_2 = saocat.find_stars_near_target(radec,
                                                                near_search_rad,
                                                                999, -999,
                                                                exclude=[cat_idx])
        #logging.debug(f'cand_idx_2 = {cand_idx_2} cand_dist_2 = {cand_dist_2}')

        if args.verbose:
            logging.debug('Nearby Candidates:')
            logging.debug(' CatIdx    Deg   SAO              RA.DEC (J2000)           VMag')
            for i2 in range(0, len(cand_idx_2)):
                cat_idx_2 = cand_idx_2[i2]
                radec = SkyCoord(saocat.ra[cat_idx_2], saocat.dec[cat_idx_2], unit=u.deg, frame='fk5', equinox='J2000')
                logging.debug(f' {cat_idx_2}   {np.rad2deg(cand_dist_2[i2]):5.2f} '
                              f'{saocat.id[cat_idx_2]:>8d} '
                              f'{radec.to_string("hmsdms", sep=":", precision=3):30s} '
                              f'{saocat.vmag[cat_idx_2]:4.2f}')

        # figure out closest
        closest_idx = np.argmin(cand_dist_2)
        #print(i, closest_idx, cand_idx_2[closest_idx], cand_dist_2[closest_idx])
        closest_deg = np.rad2deg(np.min(cand_dist_2))
        #logging.info(f"      nstars={len(cand_idx_2)} closest={closest_deg} deg")

        star_ok = True
        if closest_deg < exclusion_rad:
            # filter if neighbors too close
            star_ok = False
            nprox_exclude += 1
            if args.verbose:
                logging.debug(f'Excluding star SAO{saocat.id[cat_idx]} due " \
                             f"to neighbor within {exclusion_rad} degrees.')

        elif force_side is not None:
            # filter out if on wrong side of pier from where scope is pointing
            hour_angle = (local_sidtime - radec.ra).wrap_at('180d')
            logging.debug(f'SAO{saocat.id[cat_idx]} {radec.ra} {local_sidtime} '
                          f'{hour_angle.degree} {meridian_thres.degree}')
            # negative hour angle means the object is EAST of the meridian

            if (hour_angle.degree < meridian_thres.degree and force_side == 'WEST'):
                logging.debug('Too far EAST for WEST contraint')
                star_ok = False
            if (hour_angle.degree > -meridian_thres.degree and force_side == 'EAST'):
                logging.debug('Too far WEST for EAST contraint')
                star_ok = False
            if (abs(hour_angle.degree) > 6 * 15):
                logging.debug('More than 6 hours from meridian!')
                star_ok = False

            if not star_ok:
                npier_exclude += 1
                if True or args.verbose:
                    logging.debug(f'Excluding star SAO{saocat.id[cat_idx]} '
                                  'due to pier side or more than 6 hours from '
                                  'meridian.')

        if star_ok:
            nnear.append(len(cand_idx_2))
            nnear_idx.append(cat_idx)
            nnear_dist.append(cand_dist[i])
            if True or args.verbose:
                logging.debug(f'Including star SAO{saocat.id[cat_idx]} ')

        # if args.verbose:
        #     logging.info("")

    logging.debug(f'{nprox_exclude} stars excluded due to having close neighbors')
    logging.debug(f'{npier_exclude} stars excluded due to wrong pier side')

    if True or args.verbose:
        logging.debug(f'nnear={nnear} nnear_idx={nnear_idx} nnear_dist={nnear_dist}')

    if len(nnear) < 1:
        logging.error('No star found!')

        if args.outfile is not None:
            logging.info(f'Writing output file {args.outfile}')
            with open(args.outfile, 'w') as f:
                f.write('nstars = 0\n')
        sys.exit(1)

    #distsort_idx = np.argsort(nnear)
    # sort by # near then by distance
    distsort_idx = np.lexsort((nnear_dist, nnear))

    least_idx = distsort_idx[0]
    best_idx = nnear_idx[least_idx]
    radec = SkyCoord(saocat.ra[best_idx], saocat.dec[best_idx],
                     unit=u.deg, frame='fk5', equinox='J2000')
    logging.info(f'BEST STAR = {best_idx} SAO{saocat.id[best_idx]:>8d} '
                 f'VMag = {saocat.vmag[best_idx]:4.2f} # nearby = {nnear[least_idx]}')
    logging.info('CatIdx    Deg   SAO           RA.DEC (J2000)           VMag')
    logging.info(f'{best_idx}  {np.rad2deg(nnear_dist[least_idx]):5.2f} '
                 f'{saocat.id[best_idx]:>8d} '
                 f'{radec.to_string("hmsdms", sep=":", precision=3):30s} '
                 f'{saocat.vmag[best_idx]:4.2f}')

    #logging.info(f"{np.rad2deg(nnear_dist[least_idx]):5.2f} {saocat.id[best_idx]:10s} {radec.to_string('hmsdms', sep=':'):30s}  {saocat.vmag[best_idx]:4.2f}")

    cand_idx_2, cand_dist_2 = saocat.find_stars_near_target(radec, 1, 999, -999, exclude=[best_idx])
    #logging.debug(f'cand_idx_2 = {cand_idx_2} cand_dist_2 = {cand_dist_2}')

    if args.verbose:
        logging.debug('List of stars closest to BEST STAR')
        logging.debug(' CatIdx    Deg   SAO              RA.DEC (J2000)           VMag')
        for i2 in range(0, len(cand_idx_2)):
            cat_idx_2 = cand_idx_2[i2]
            radec = SkyCoord(saocat.ra[cat_idx_2], saocat.dec[cat_idx_2],
                             unit=u.deg, frame='fk5', equinox='J2000')
            logging.debug(f' {cat_idx_2}   {np.rad2deg(cand_dist_2[i2]):5.2f} '
                          f'{saocat.id[cat_idx_2]:>8d} '
                          f'{radec.to_string("hmsdms", sep=":", precision=3):30s} '
                          f'{saocat.vmag[cat_idx_2]:4.2f}')

    logging.debug(''""'')
    logging.debug(f'{len(distsort_idx)} candidates within {args.dist} degrees of J2000 '
                  f'{target.to_string("hmsdms", sep=":"):30s}')
    logging.debug('CatIdx    Deg     SAO           RA.DEC (J2000)            '
                  'VMag     # near stars')
    for i in range(0, len(distsort_idx)):
        near_idx = distsort_idx[i]
        cat_idx = nnear_idx[near_idx]
        radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx],
                         unit=u.deg, frame='fk5', equinox='J2000')
        logging.debug(f'{cat_idx:>6d}  {np.rad2deg(nnear_dist[near_idx]):5.2f} '
                      f'{saocat.id[cat_idx]:>8d}   '
                      f'{radec.to_string("hmsdms", sep=":", precision=3):30s} '
                      f'{saocat.vmag[cat_idx]:4.2f}        {nnear[near_idx]}')

    # if output file requested create it
    if args.outfile is not None:
        logging.info(f'Writing output file {args.outfile}')
        with open(args.outfile, 'w') as f:
            f.write(f'nstars = {len(distsort_idx)}\n')
            f.write('#CatIdx,Distance (deg), SAO #, RA (J2000)), DEC (J2000),'
                    ' VMag, # near stars\n')
            for i in range(0, len(distsort_idx)):
                near_idx = distsort_idx[i]
                cat_idx = nnear_idx[near_idx]
                radec = SkyCoord(saocat.ra[cat_idx], saocat.dec[cat_idx],
                                 unit=u.deg, frame='fk5', equinox='J2000')

                f.write(f'{cat_idx}, {np.rad2deg(nnear_dist[near_idx])},'
                        f' {saocat.id[cat_idx]},'
                        f' {radec.ra.to_string(unit=u.hour, sep=":", precision=3)},'
                        f' {radec.dec.to_string(unit=u.degree, sep=":", precision=3)},'
                        f' {saocat.vmag[cat_idx]}, {nnear[near_idx]}\n')
