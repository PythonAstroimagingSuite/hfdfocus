# find closest star within mag range
# FIXME DOES NOT CONSIDER MERIDAN FLIP IMPLICATIONS!

import os
import sys
import time
import argparse
import logging

from SAOCatalog import SAOCatalog, write_SAOCatalog_binary

def numstr(x):
    if x < 0:
        return f'm{abs(x)}'
    else:
        return f'p{x}'

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
    parser.add_argument('outfile', type=str, help='Output filename base')
    parser.add_argument('--mindec', type=float, default=-90)

    args = parser.parse_args()

    logging.info(f'args = {args}')

    # FIXME wasteful to keep reading over and over
    for maxmag, minmag in [(-5, 6), (6,7), (7,8), (8,9), (10,11), (11,12)]:
        logging.info('Running catalog for maxmag = {maxmag} minmag = {minmag}')
        saocat = SAOCatalog()
        saocat.read_catalog(maxmag=maxmag, minmag=minmag, mindec=args.mindec)
        write_SAOCatalog_binary(saocat, f'{args.outfile}_{numstr(maxmag)}_{numstr(minmag)}_{numstr(int(args.mindec))}.bin')
