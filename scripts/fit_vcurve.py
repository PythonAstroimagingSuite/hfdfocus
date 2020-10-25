#
# fit lines to vcurve
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
import logging
import argparse
import numpy as np
from scipy.stats import siegelslopes
import matplotlib.pyplot as plt


if __name__ == '__main__':

    logging.basicConfig(filename='fit_vcurve.log',
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

    parser = argparse.ArgumentParser()
    parser.add_argument('hfd_filename', type=str, nargs='?',
                        default='hfd.txt', help='HFR data file name')
    args = parser.parse_args()

    hfd_file = args.hfd_filename

    f = open(hfd_file, 'r')

    fpos_arr = []
    hfd_arr = []
    for line in f.readlines():
        fpos, hfd = line.strip().split(',')
        fpos = float(fpos)
        hfd = float(hfd)
        print(fpos, hfd)
        fpos_arr.append(fpos)
        hfd_arr.append(hfd)

    fpos_arr = np.array(fpos_arr)
    hfd_arr = np.array(hfd_arr)
    print('fpos', fpos_arr)
    print('hfd', hfd_arr)

    # find mininum value
    midx = np.argmin(hfd_arr)
    fpos_arr_l = np.array(fpos_arr[:midx - 1])
    fpos_arr_r = np.array(fpos_arr[midx + 1:])
    hfd_arr_l = np.array(hfd_arr[:midx - 1])
    hfd_arr_r = np.array(hfd_arr[midx + 1:])

    print('fpos_l', fpos_arr_l)
    print('hfd_l', hfd_arr_l)
    print('fpos_r', fpos_arr_r)
    print('hfd_r', hfd_arr_r)

    fig = plt.figure()
#    ax_1 = fig.add_subplot(131)
#    ax_2 = fig.add_subplot(132)
    ax_3 = fig.add_subplot(111)
#    ax_1.plot(fpos_arr_l, hfd_arr_l)
#    ax_2.plot(fpos_arr_r, hfd_arr_r)
    ax_3.scatter(fpos_arr, hfd_arr, color='green')

    # flip left side around so HFR values INCREASE with index
    print('fit')
#    robust_right = robust_line_fit(fpos_arr_r, hfd_arr_r)
#    robust_left = robust_line_fit(np.flipud(fpos_arr_l), np.flipud(hfd_arr_l))
#    robust_best_pos = (robust_left[1]-robust_right[1])/(robust_right[0]-robust_left[0])
#    print(f"Robust left  -> {robust_left}")
#    print(f"Robust right -> {robust_right}")
#    print(f"Robust intersection/best focus -> {robust_best_pos}")

    siegel_left_fit = siegelslopes(hfd_arr_l, fpos_arr_l)
    siegel_right_fit = siegelslopes(hfd_arr_r, fpos_arr_r)
    siegel_left_zero = -siegel_left_fit[1] / siegel_left_fit[0]
    siegel_right_zero = -siegel_right_fit[1] / siegel_right_fit[0]
    siegel_best_pos = (siegel_left_fit[1] - siegel_right_fit[1]) / (siegel_right_fit[0] - siegel_left_fit[0])
    logging.info(f'siegel left  fit = {siegel_left_fit}')
    logging.info(f'siegel right fit = {siegel_right_fit}')
    logging.info(f'siegel best pos  = {siegel_best_pos}')

    ax_3.plot(fpos_arr[:midx + 5], siegel_left_fit[0] * fpos_arr[:midx + 5]+siegel_left_fit[1])
    ax_3.plot(fpos_arr[midx - 5:], siegel_right_fit[0] * fpos_arr[midx - 5:]+siegel_right_fit[1])
    ax_3.axvline(siegel_best_pos, color='red')

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

    plt.show()
