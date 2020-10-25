#
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
import time
import argparse
import logging

import astropy.io.fits as pyfits

import numpy as np
from scipy.stats import siegelslopes

import matplotlib.pyplot as plt

from StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window

if __name__ == '__main__':
    logging.basicConfig(filename='measure_hfd_one_file.log',
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
    parser.add_argument('image', type=str, help='Image name')
    parser.add_argument('--noplot', action='store_true',  help='Do not show plots and do not pause')
    parser.add_argument('--pause', type=float,  help='Pause instead of blocking when done')
    parser.add_argument('--bgmodel', action='store_true',  help='Compute full bg model')
    parser.add_argument('--starmodel', action='store_true',  help='Compute full star model')
    parser.add_argument('--framesize', default=0, type=int,  help='Size of capture frame, 0=full')
    parser.add_argument('--saturation', default=55000, type=int,  help='Saturation level for sensor')
    parser.add_argument('--bgthres', default=50, type=int,  help='Threshold multiplier for star detection')
    args = parser.parse_args()

    logging.info(f'args = {args}')

    if not args.noplot:
        fig = plt.figure(figsize=(4, 3))
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)

    # analyze frame
    bg = 800
    #thres = 10000

    hdu = pyfits.open(args.image)
    starimage_data = hdu[0].data.astype(float)
    hdu.close()

    logging.info(f'Image size is {starimage_data.shape[1]} x {starimage_data.shape[0]}')

    # use subframe
    if args.framesize != 0:
        h, w = starimage_data.shape
        xl = int(w / 2 - args.framesize / 2)
        xw = args.framesize
        yt = int(h / 2 - args.framesize / 2)
        yh = args.framesize
        logging.info(f'Shrinking to framesize = {args.framesize} ({xl}, {yt}) {xw} x {yh}')
        starimage_data = starimage_data[yt:yt + yh, xl:xl + xw]


    xcen, ycen, bg, mad, starmask = find_star(starimage_data, bgmodel=args.bgmodel,
                                    starmodel=args.starmodel, bgfact=args.bgthres,
                                    debugfits=True)

    thres = bg + args.bgthres * mad
    logging.info(f'Using thres = {thres}')

    if np.max(starimage_data[starmask] > args.saturation):
        logging.warning('SATURATED PIXELS DETECTED!')

    win = 300
    xlow = max(0, int(xcen - win / 2))
    xhi = min(starimage_data.shape[0] - 1, int(xcen + win / 2))
    ylow = max(0, int(ycen - win / 2))
    yhi = min(starimage_data.shape[1] - 1, int(ycen + win / 2))
    logging.debug(f'cropping to window={win} x={xlow}:{xhi} y={ylow}:{yhi}')
    crop_data = starimage_data[ylow:yhi, xlow:xhi]

    if not args.noplot:
        fig.clear()
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)
        im = ax_2d.imshow((crop_data - bg).astype(float))
        fig.colorbar(im, ax=ax_2d)
        #plt.show()

    profile = horiz_bin_window(crop_data, bg=bg)

    if not args.noplot:
        ax_1d.plot(profile)

    rc = find_hfd_from_1D(profile, thres=thres)
    if rc is not None:
#        logging.info(rc)
        pass
    else:
        logging.warning('No star found in image')
        plt.show()
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc

    if not args.noplot:
        ax_1d.axvline(scen, color='red')
        if sl is not None and sr is not None:
            ax_1d.axvline(sl, color='green')
            ax_1d.axvline(sr, color='green')
            ax_1d.axvline(hfl, color='blue')
            ax_1d.axvline(hfr, color='blue')
            delta = sr - sl
            ax_1d.set_xlim(sl - delta / 4, sr + delta / 4)
            ax_1d.set_title(f'{hfr - hfl:5.3f}')
        print('drawing plot')
        fig.show()
        if args.pause is not None:
            logging.info(f'Pausing {args.pause} seconds')
            plt.pause(args.pause)
        else:
            plt.show()

    sys.exit(0)