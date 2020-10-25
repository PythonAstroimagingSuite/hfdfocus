# find closest star within mag range
# FIXME DOES NOT CONSIDER MERIDAN FLIP IMPLICATIONS!
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
    write_SAOCatalog_binary(saocat, args.outfile)
