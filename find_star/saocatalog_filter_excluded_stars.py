#!/usr/bin/env python3
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
import argparse
import logging


from SAOCatalog import load_SAOCatalog_binary, write_SAOCatalog_binary

def delete_entry(cat, idx):
    del cat.id[idx]
    del cat.ra[idx]
    del cat.dec[idx]
    del cat.vmag[idx]


if __name__ == '__main__':
    logging.basicConfig(filename='saocatalog_filter_excluded_stars.log',
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
    parser.add_argument('input_cat', type=str, help='Catalog to process')
    parser.add_argument('output_cat', type=str, help='Output catalog')
    parser.add_argument('exclude_file', type=str, help='Exclusion list')

    args = parser.parse_args()

    logging.info(f'args = {args}')

    # see if exclusion list already exists so we don't overwrite
    if os.path.isfile(args.output_cat):
        logging.error(f'The file {args.output_cat} already exists!')
        logging.error('Move it out of the way before running this command or')
        sys.exit(1)

    saocat = load_SAOCatalog_binary(args.input_cat)
    logging.info(f'Loaded {len(saocat.id)} stars')

    excl_file = open(args.exclude_file, 'r')
    exclude_list = []
    for line in excl_file.readlines():
        catidx, saoid = line.strip().split(',')
        exclude_list.append(int(catidx))

    logging.info(f'Read in {len(exclude_list)} exclusion records')

    # need to delete from end of list backwards so exlcusion idx are still
    # relevant!
    exclude_list = sorted(exclude_list, reverse=True)

    logging.info('Filtering catalog...')
    ndone = 0
    for excl_idx in exclude_list:
        delete_entry(saocat, excl_idx)
        ndone += 1
        if ndone % 10000 == 0:
            logging.info(f'{ndone} of {len(exclude_list)} complete')

    logging.info(f'Writing output catalog {args.output_cat}')
    write_SAOCatalog_binary(saocat, args.output_cat)

    logging.info('Finished!')
