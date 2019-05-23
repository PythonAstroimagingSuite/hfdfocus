# find closest star within mag range
# FIXME DOES NOT CONSIDER MERIDAN FLIP IMPLICATIONS!

import os
import sys
import time
import argparse
import logging

from SAOCatalog import SAOCatalog, write_SAOCatalog_binary

if __name__ == '__main__':
    logging.basicConfig(filename='saocatalog_convert.log',
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
#    parser.add_argument('focus_low', type=int, help='Low end of focus run')
#    parser.add_argument('focus_high', type=int, help='High end of focus run')
    parser.add_argument('outfile', type=str, help='Output filename')
    parser.add_argument('--maxmag', type=float, default=-5)
    parser.add_argument('--minmag', type=float, default=9)
    parser.add_argument('--mindec', type=float, default=-90)

    args = parser.parse_args()

    logging.info(f'args = {args}')

    saocat = SAOCatalog()
    saocat.read_catalog(maxmag=args.maxmag, minmag=args.minmag, mindec=args.mindec)
    #saocat.find_nearby(SkyCoord("15:00:00 +30:00:00", unit=(u.hourangle, u.deg), frame='fk5', equinox='J2000'))

    write_SAOCatalog_binary(saocat, args.outfile)
#    saocat_bin = load_SAOCatalog_binary('SAO_Catalog.bin')
#
#    print(set(saocat.id) == set(saocat_bin.id))
#    print(set(saocat.radec) == set(saocat_bin.radec))
#    print(set(saocat.vmag) == set(saocat_bin.vmag))
#
#    for i in range(0, len(saocat.radec)):
#        logging.info(f"{saocat.id[i]:10s} {saocat.radec[i].to_string('hmsdms', sep=':'):30s}")
#        logging.info(f"{saocat_bin.id[i]:10s} {saocat_bin.radec[i].to_string('hmsdms', sep=':'):30s}")
#        logging.info("")