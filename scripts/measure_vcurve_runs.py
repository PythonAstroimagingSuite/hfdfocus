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
import re
import time
import glob
import json
import logging
import argparse
import numpy as np
from scipy.stats import siegelslopes

import astropy.io.fits as pyfits

# for testing
import matplotlib.pyplot as plt

from hfdfocus.StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window


if __name__ == '__main__':
    logging.basicConfig(filename='measure_vcurve.log',
                        filemode='w',
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    #test_1d_with_gaussian()
    parser = argparse.ArgumentParser()
    parser.add_argument('imgdir', type=str, help='Directory containing V Curve Images')
    parser.add_argument('--hfdcutoff', default=10, type=float,  help='Ignore points with HFD less than this value')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')

    args = parser.parse_args()

#    logging.info(f'command args = {args}')

    if args.debugplots:
        fig = plt.figure(figsize=(4.5, 2))
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)
        fig2 = plt.figure(figsize=(4, 3))
        ax_hfd = fig2.add_subplot(111)
        hfd_plot, = ax_hfd.plot([], [], marker='o', fillstyle='none', ls='')

    run = 1
    fit_arr = []
    while True:
        imgfiles = glob.glob(os.path.join(args.imgdir, f'vcurve_focuspos_run{run:03d}*.fit'))

        logging.info(f'imgfiles = {imgfiles}')

        if len(imgfiles) < 1:
            break

        focus_re = re.compile('^.*focuspos_run(?:\d{2,})_(?P<pos>\d{2,}).fit$')

        fpos_arr = []
        hfd_arr = []

        for infile in sorted(imgfiles):

            focus_pos = int(focus_re.match(infile).group(1))
            logging.info(f'infile = {infile} focus_pos = {focus_pos}')

            hdu = pyfits.open(infile)
            image_data = hdu[0].data.astype(float)
            hdu.close()

            bg = 800
            thres = 10000

            xcen, ycen, bg, mad, starmask, alon = find_star(image_data, debugfits=False)

            win = 100
            xlow = int(xcen - win / 2)
            xhi = int(xcen + win / 2)
            ylow = int(ycen - win / 2)
            yhi = int(ycen + win / 2)
            crop_data = image_data[ylow:yhi, xlow:xhi]

            if args.debugplots:
                print('clear')
                fig.clear()
                #mpl.rcParams.update({'axes.labelsize' : 18})
                ax_1d = fig.add_subplot(121)
                ax_2d = fig.add_subplot(122)
                im = ax_2d.imshow((crop_data - bg).astype(float))
                fig.colorbar(im, ax=ax_2d)

            profile = horiz_bin_window(crop_data, bg=bg)

            #print('thres = ', thres)
            #print('profile = ', profile)

            if args.debugplots:
                ax_1d.plot(profile)

            scen, sl, sr, hfl, hfr, tflux = find_hfd_from_1D(profile, thres=thres)

            if args.debugplots:
                ax_1d.axvline(scen, color='red')
                if sl is not None and sr is not None:
                    ax_1d.axvline(sl, color='green')
                    ax_1d.axvline(sr, color='green')
                    ax_1d.axvline(hfl, color='blue')
                    ax_1d.axvline(hfr, color='blue')
                    delta = sr - sl
                    ax_1d.set_xlim(sl-delta/4, sr+delta/4)

            #print('total counts = ', np.sum(crop_data-bg))

            fpos_arr.append(focus_pos)
            hfd_arr.append(hfr - hfl)
            logging.info(f'{fpos_arr} {hfd_arr}')

            if args.debugplots:
                print('fig show')
                fig.show()
                plt.pause(0.05)

        # fit
        fpos_arr = np.array(fpos_arr)
        hfd_arr = np.array(hfd_arr)

        # sort so fpos is increasing
        fargs = fpos_arr.argsort()
        fpos_arr = fpos_arr[fargs]
        hfd_arr = hfd_arr[fargs]

        # find mininum value
        midx = np.argmin(hfd_arr)

        # set new focus center to min
        focus_center = fpos_arr[midx]
        logging.info(f'Set new focus center to {focus_center}')

        fpos_arr_l = np.array(fpos_arr[:midx - 2])
        fpos_arr_r = np.array(fpos_arr[midx + 3:])
        hfd_arr_l = np.array(hfd_arr[:midx - 2])
        hfd_arr_r = np.array(hfd_arr[midx + 3:])

        # apply threshold
        l_hfd_filter = np.where(hfd_arr_l > args.hfdcutoff)
        r_hfd_filter = np.where(hfd_arr_r > args.hfdcutoff)
        fpos_arr_l = fpos_arr_l[l_hfd_filter]
        fpos_arr_r = fpos_arr_r[r_hfd_filter]
        hfd_arr_l = hfd_arr_l[l_hfd_filter]
        hfd_arr_r = hfd_arr_r[r_hfd_filter]

        print('fpos_l', fpos_arr_l)
        print('hfd_l', hfd_arr_l)
        print('fpos_r', fpos_arr_r)
        print('hfd_r', hfd_arr_r)

        siegel_left_fit = siegelslopes(hfd_arr_l, fpos_arr_l)
        siegel_right_fit = siegelslopes(hfd_arr_r, fpos_arr_r)
        siegel_left_zero = -siegel_left_fit[1] / siegel_left_fit[0]
        siegel_right_zero = -siegel_right_fit[1] / siegel_right_fit[0]
        siegel_best_pos = (siegel_left_fit[1] - siegel_right_fit[1]) / (siegel_right_fit[0] - siegel_left_fit[0])
        logging.info(f'siegel left  fit = {siegel_left_fit}')
        logging.info(f'siegel right fit = {siegel_right_fit}')
        logging.info(f'siegel best pos  = {siegel_best_pos}')

        if args.debugplots:
            ax_hfd.plot(fpos_arr_l, hfd_arr_l, marker='+', ls='', color='red')
            ax_hfd.plot(fpos_arr_r, hfd_arr_r, marker='+', ls='', color='red')
            ax_hfd.plot(fpos_arr[midx - 5:], siegel_right_fit[0] * fpos_arr[midx  -5:] + siegel_right_fit[1], color='green')
            ax_hfd.plot(fpos_arr[:midx + 5], siegel_left_fit[0] * fpos_arr[:midx + 5] + siegel_left_fit[1], color='blue')
            ax_hfd.axvline(siegel_best_pos, color='red')
            ax_hfd.set_title(f'Left {siegel_left_fit[0]:7.6f}/{siegel_best_pos - siegel_left_zero:5.3f} Right {siegel_right_fit[0]:7.6f}/{siegel_best_pos - siegel_right_zero:5.3f}')
            ax_hfd.relim()
            ax_hfd.set_ylim(bottom=0)
            ax_hfd.autoscale_view()
            fig2.canvas.draw()
            plt.pause(0.1)

        fit_arr.append((siegel_left_fit[0], siegel_best_pos - siegel_left_zero, siegel_right_fit[0], siegel_best_pos - siegel_right_zero))

        print(fpos_arr[:midx + 5], siegel_left_fit[0] * fpos_arr[:midx + 5] + siegel_left_fit[1])

        logging.info('Left Side:')
        logging.info(f'   slope: {siegel_left_fit[0]}')
        logging.info(f'   inter: {siegel_left_fit[1]}')
        logging.info(f'   yzero: {siegel_left_zero}')
        logging.info(f'   PID  : {siegel_best_pos - siegel_left_zero}')
        logging.info('Right Side:')
        logging.info(f'   slope: {siegel_right_fit[0]}')
        logging.info(f'   inter: {siegel_right_fit[1]}')
        logging.info(f'   yzero: {siegel_right_zero}')
        logging.info(f'   PID  : {siegel_best_pos - siegel_right_zero}')

        f = open('measure_vcurve_runs_fits.json', 'a')
        ls = siegel_left_fit[0]
        lp = siegel_best_pos - siegel_left_zero
        rs = siegel_right_fit[0]
        rp = siegel_best_pos - siegel_right_zero
        tstamp = time.strftime('%Y/%m/%d %H:%M:%S %Z')
        #f.write(f'{tstamp}, {ls}, {lp}, {rs}, {rp}\n')
        j = json.dumps({ 'timestamp' : tstamp, 'rightslope' : rs, 'rightpid' : rp, 'leftslope' : ls, 'leftpid' : lp})
        f.write(j + '\n')
        f.close()

        run = run + 1

#            f = open(os.path.join(args.imgdir, 'hfd.txt'), 'a')
#            f.write(f'{focus_pos}, {hfr-hfl}\n')
#            f.close()
