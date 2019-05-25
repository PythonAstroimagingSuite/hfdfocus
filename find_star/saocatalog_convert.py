# find closest star within mag range
# FIXME DOES NOT CONSIDER MERIDAN FLIP IMPLICATIONS!

import os
import sys
import time
import argparse
import logging

from hfdfocus.SAOCatalog import SAOCatalog, write_SAOCatalog_binary

def read_original_catalog(saocat, maxmag, minmag, mindec):
    f = open('SAO_Catalog.dat')
    first = True
    nread = 0
    saocat.maxmag = maxmag
    saocat.minmag = minmag
    saocat.mindec = mindec
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
            saocat.id.append(int(hid))
            saocat.ra.append(float(ra_deg))
            saocat.dec.append(float(dec_deg))
            saocat.vmag.append(float(vmag))
        except Exception as err:
            logging.error(f'Error processing #2 |{l.rstrip()}| |{vmag}| {err}')
            continue

        nread += 1
        if nread % 1000 == 0:
            logging.info(f'{nread} read')
#            if nread > 100:
#                break

    f.close()
    logging.info(f'# stars loaded = {nread}')

    i = 1
    for id, radec, vmag in zip(saocat.id, saocat.radec, saocat.vmag):
        logging.info(f"{i:5d} {id:10s} {radec.to_string('hmsdms', sep=':'):30s}  {vmag:4.2f}")
        i += 1

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

    read_original_catalog(saocat, maxmag=args.maxmag, minmag=args.minmag, mindec=args.mindec)

    write_SAOCatalog_binary(saocat, args.outfile)
