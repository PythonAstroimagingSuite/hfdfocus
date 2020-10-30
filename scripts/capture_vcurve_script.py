#!/usr/bin/env python3
#
# Capture a Vcurve
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
import os
import sys
import time
import json
import argparse
from datetime import datetime

import astropy.io.fits as pyfits

import numpy as np
from scipy.stats import siegelslopes

import matplotlib.pyplot as plt

from pyastroprofile.AstroProfile import AstroProfile

from hfdfocus.StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window

# for simulator
from hfdfocus.c8_simul_star import C8_F7_Star_Simulator


def parse_commandline():
    parser = argparse.ArgumentParser()
#    parser.add_argument('focus_low', type=int, help='Low end of focus run')
#    parser.add_argument('focus_high', type=int, help='High end of focus run')
    parser.add_argument('focus_center', type=int, help='Center position of focus run')
    parser.add_argument('focus_range', type=int, help='Range of focus run')
    parser.add_argument('focus_nstep', type=int, help='V Curve number of steps')
    parser.add_argument('focus_dir', type=str, help='IN or OUT')
    parser.add_argument('nruns', type=int, help='Number of vcurve runs')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')
    parser.add_argument('--savefits', action='store_true', help='Save all images taken')
    parser.add_argument('--simul', action='store_true', help='Simulate star')
    parser.add_argument('--profile', type=str, help='Name of astro profile')
    parser.add_argument('--backend', type=str,  help='Backend')
    parser.add_argument('--focuser', type=str,  help='Focuser Driver')
    parser.add_argument('--camera', type=str,  help='Camera Driver')
    parser.add_argument('--exposure_start', default=0.25, type=float,  help='Starting exposure value')
    parser.add_argument('--exposure_min', default=0.25, type=float,  help='Minimum exposure value')
    parser.add_argument('--exposure_max', default=8, type=float,  help='Maximum exposure value')
    parser.add_argument('--saturation', default=55000, type=int,  help='Saturation level for sensor')
    parser.add_argument('--starflux_min', default=50000, type=int,  help='Maximum flux in star')
    parser.add_argument('--framesize', default=1000, type=int,  help='Size of capture frame, 0=full')
    parser.add_argument('--runoffset', default=0, type=int,  help='Shift center of run by this amount')
    parser.add_argument('--hfdcutoff', default=10, type=float,  help='Ignore points with HFD less than this value')
    parser.add_argument('--bgthres', default=50, type=int,  help='Threshold multiplier for star detection')
    parser.add_argument('--movedelay', default=0, type=float,  help='Delay in seconds between moves')
    parser.add_argument('--backlash', default=50, type=int,  help='Number of steps of backlash for overshoot method')

    return parser.parse_args()

if __name__ == '__main__':

    log_timestamp = datetime.now()
    logfilename = 'capture_vcurve_script-' + log_timestamp.strftime('%Y%m%d%H%M%S') + '.log'

    LONG_FORMAT = '%(asctime)s.%(msecs)03d [%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'
    SHORT_FORMAT = '%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s'

    logging.basicConfig(filename=logfilename,
                        filemode='w',
                        level=logging.DEBUG,
                        format=LONG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter(SHORT_FORMAT)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    args = parse_commandline()
    logging.info(f'args = {args}')

    # connect focuser and camera
    backend_name = None
    focuser_driver = None
    camera_driver = None

    if not args.simul:
        from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI
        sdi = SDI()

        # use profile if provided
        if args.profile is not None:
            logging.info(f'Setting up device using astro profile {args.profile}')
            ap = AstroProfile()
            ap.read(args.profile)
            prof_backend_name = ap.equipment.backend.name
            logging.info(f'profile backend = {prof_backend_name}')
            if prof_backend_name is not None:
                backend_name = prof_backend_name
            prof_camera_driver = ap.equipment.camera.driver
            logging.info(f'profile camera driver = {prof_camera_driver}')
            if prof_camera_driver is not None:
                camera_driver = prof_camera_driver
            prof_focuser_driver = ap.equipment.focuser.driver
            logging.info(f'profile mount driver = {prof_focuser_driver}')
            if prof_focuser_driver is not None:
                focuser_driver = prof_focuser_driver

        # command line overrides
        if args.backend is not None:
            backend_name = args.backend

        if backend_name is None:
            logging.error('Must specify the backend in profile or on command line!')
            sys.exit(1)

        logging.info(f'Connecting to backend {backend_name}')
        sdi.connect_backend(backend_name)

        if args.focuser is not None:
            focuser_driver = args.focuser

        if focuser_driver is None:
            logging.error('Must specify the focuser driver in profile or on command line!')
            sys.exit(1)

        logging.info(f'Connecting to focuser driver {focuser_driver}')
        focuser = sdi.connect_focuser(focuser_driver)
        logging.info(f'focuser = {focuser}')
        if not focuser:
            logging.error(f'Unabled to connect to focuser driver {focuser_driver}')
            sys.exit(1)

        if args.camera is not None:
            camera_driver = args.camera

        if camera_driver is None:
            logging.error('Must specify the camera driver in profile or on command line!')
            sys.exit(1)

        logging.info(f'Connecting to camera driver {camera_driver}')
        cam = sdi.connect_camera(camera_driver)
        logging.info(f'cam = {cam}')
        if not cam:
            logging.error(f'Unabled to connect to camera driver {camera_driver}')
            sys.exit(1)
    else:
        simul_star = C8_F7_Star_Simulator(companion_offset=(20, 20), companion_flux_ratio=1.0)

    # create output dir
    datestr = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    imagesdir = datestr
    os.mkdir(imagesdir)

    if args.debugplots:
        logging.info('Creating figure')
        fig2 = plt.figure(figsize=(4, 3))
        ax_hfd = fig2.add_subplot(111)
        hfd_plot, = ax_hfd.plot([], [], marker='o', fillstyle='none', ls='')

        fig = plt.figure(figsize=(4, 3))
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)

    fit_arr = []

    focus_center = args.focus_center
    focus_range = args.focus_range
    focus_expos = args.exposure_start
    for iter in range(0, args.nruns):
        if args.debugplots:
            if iter != 0:
                fig2.clear()
                ax_hfd = fig2.add_subplot(111)
                hfd_plot, = ax_hfd.plot([], [], marker='o', ls='')

                fig.clear()
                ax_2d = fig.add_subplot(121)
                ax_1d = fig.add_subplot(122)

        # figure out direction
        backlash = args.backlash
        logging.info(f'Using backlash = {args.backlash}')
        logging.info(f'Shifting focus center by runoffset = {args.runoffset}')
        focus_low = int(focus_center + args.runoffset-focus_range / 2)
        focus_high = focus_low + args.focus_range
        if args.focus_dir == 'OUT':
            # start past desired start and move to it to remove backlash
            focus_init = focus_low - backlash
            focus_start = focus_low
            focus_end = focus_high
            focus_nstep = args.focus_nstep
        elif args.focus_dir == 'IN':
            # start past desired start and move to it to remove backlash
            focus_init = focus_high + backlash
            focus_start = focus_high
            focus_end = focus_low
            focus_nstep = args.focus_nstep
        else:
            logging.error(f'Unknown focus direction {args.focus_dir} - exitting!')
            sys.exit(1)

        # move to init pos
        logging.info(f'Moving to init pos {focus_init}')
        if not args.simul:
            logging.info(f'initial focuser pos = {focuser.get_absolute_position()}')
            if not focuser.move_absolute_position(focus_init):
                logging.error("Focuser error!!")
                sys.exit(1)

            sdi.wait_on_focuser_move(focuser)  # , position=focus_init)

            # while True:
            #     logging.info(f'current focuser pos = {focuser.get_absolute_position()}')
            #     time.sleep(0.25)

            logging.info('Extra wait for focuser 1 sec')
            time.sleep(1)

        focus_step = int((focus_end - focus_start) / (focus_nstep - 1))
        logging.info(f'Focus run from {focus_start} to {focus_end} step {focus_step}')
        fpos_arr = []
        hfd_arr = []
        for focus_pos in range(focus_start, focus_end + focus_step, focus_step):
            logging.info(f'Moving to focus pos {focus_pos}')
            imgname = os.path.join(imagesdir, f'vcurve_focuspos_run{iter+1:03d}_{focus_pos}.fit')
            star_found = False
            while True:
                if not args.simul:
                    if not focuser.move_absolute_position(focus_pos):
                        logging.error("Focuser error!!")
                        sys.exit(1)

                    sdi.wait_on_focuser_move(focuser)  # , position=focus_pos)

                    if args.movedelay > 0:
                        logging.info(f'movedelay {args.movedelay} seconds')
                        time.sleep(args.movedelay)

                    # NOTE roi in INDI not working for some reason (EKOS overriding?)
                    # just take full frame then shrink
                    roi = None

                    if args.framesize != 0:
                        w, h = cam.get_size()
                        xl = int(w / 2 - args.framesize / 2)
                        xw = args.framesize
                        yt = int(h / 2 - args.framesize / 2)
                        yh = args.framesize
                        roi = (xl, yt, xw, yh)
                    else:
                        roi = None

                    logging.info(f'Taking exposure exposure = {focus_expos} seconds roi = {roi}')

                    rc = sdi.take_exposure(cam, focus_expos, imgname, roi=roi, overwrite=True)

                    logging.info(f'exposure result code = {rc}')

                else:
                    tmp_starimage_data = simul_star.get_simul_star_image(focus_pos)
                    pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)

                hdu = pyfits.open(imgname)
                starimage_data = hdu[0].data.astype(float)
                hdu.close()

                if not args.savefits:
                    os.unlink(imgname)

                # use subframe
#                if not args.simul and args.framesize != 0:
#                    w, h = cam.get_size()
#                    xl = int(w/2-args.framesize/2)
#                    xw = args.framesize
#                    yt = int(h/2-args.framesize/2)
#                    yh = args.framesize
#                    logging.info(f'Shrinking to framesize = {args.framesize}')
#                    starimage_data = starimage_data[yt:yt+yh, xl:xl+xw]
#                    print(starimage_data.shape)
#                    pyfits.writeto('starimage_data.fits', starimage_data.astype(float), overwrite=True)

                # analyze frame
                #bg = 800
                #thres = 10000

                result = find_star(starimage_data, bgfact=args.bgthres, debugfits=True)
                if result is None:
                    logging.warning('No star found - moving to next step!')
                    break
                else:
                    xcen, ycen, bg, mad, starmask, alone = result

                #thres = bg + args.bgthres*mad
                thres = 10000

                logging.info(f'Using thres = {thres}')

                if np.max(starimage_data[starmask] > args.saturation):
                    logging.warning(f'SATURATED PIXELS DETECTED!')

                win = 200
                xlow = max(0, int(xcen - win / 2))
                xhi = min(starimage_data.shape[0] - 1, int(xcen + win / 2))
                ylow = max(0, int(ycen - win / 2))
                yhi = min(starimage_data.shape[1] - 1, int(ycen + win / 2))
                #logging.debug(f'cropping to window={win} x={xlow}:{xhi} y={ylow}:{yhi}')
                crop_data = starimage_data[ylow:yhi, xlow:xhi]

                if args.debugplots:
                    fig.clear()
                    ax_2d = fig.add_subplot(121)
                    ax_1d = fig.add_subplot(122)
                    im = ax_2d.imshow((crop_data - bg).astype(float))
                    fig.colorbar(im, ax=ax_2d)

                profile = horiz_bin_window(crop_data, bg=bg)

                if args.debugplots:
                    ax_1d.plot(profile)

                rc = find_hfd_from_1D(profile, thres=thres)
                if rc is not None:
                    scen, sl, sr, hfl, hfr, totflux = rc
                    if totflux >= args.starflux_min:
                        star_found = True
                        break
                    else:
                        logging.warning('Star flux in image {imgname} ({totflux} is less than {args.starflux_min}')
                else:
                    logging.warning(f'No star found in image {imgname}')

                focus_expos += 1
                if args.simul:
                    logging.info(f'Simulated image did not have star - check settings.')
                    sys.exit(1)

                if focus_expos <= args.exposure_max:
                    logging.info(f'Bumping exposure to {focus_expos} and retrying')
                    continue
                else:
                    logging.info(f'Exposure is already at max of {args.exposure_max} seconds - aborting')
                    sys.exit(1)

            if star_found:
                fpos_arr.append(focus_pos)
                hfd_arr.append(hfr - hfl)

                if args.debugplots:
                    hfd_plot.set_data(fpos_arr, hfd_arr)
                    ax_hfd.relim()
                    ax_hfd.autoscale_view()
                    ax_hfd.set_title(f'iter {iter+1} of {args.nruns} pt {len(hfd_arr)+1} of {args.focus_nstep}')

                    fig2.canvas.draw()

                    ax_1d.axvline(scen, color='red')
                    if sl is not None and sr is not None:
                        ax_1d.axvline(sl, color='green')
                        ax_1d.axvline(sr, color='green')
                        ax_1d.axvline(hfl, color='blue')
                        ax_1d.axvline(hfr, color='blue')
                        delta = sr - sl
                        ax_1d.set_xlim(sl - delta / 4, sr + delta / 4)
                        ax_1d.set_title(f'{hfr - hfl:5.3f} {alone}')
                    #print('drawing plot')
                    fig.show()
                    plt.pause(0.05)
                    #plt.show()

        # fit
        if len(fpos_arr) == 0:
            logging.error('No data points so V Curve cannot be fit!')
            sys.exit(1)

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
            ax_hfd.plot(fpos_arr[midx - 5:], siegel_right_fit[0] * fpos_arr[midx - 5:] + siegel_right_fit[1], color='green')
            ax_hfd.plot(fpos_arr[:midx + 5], siegel_left_fit[0] * fpos_arr[:midx + 5] + siegel_left_fit[1], color='blue')
            ax_hfd.axvline(siegel_best_pos, color='red')
            ax_hfd.set_title(f'Left {siegel_left_fit[0]:7.6f}/{siegel_best_pos - siegel_left_zero:5.3f} Right {siegel_right_fit[0]:7.6f}/{siegel_best_pos - siegel_right_zero:5.3f}')
            ax_hfd.relim()
            ax_hfd.set_ylim(bottom=0)
            ax_hfd.autoscale_view()
            fig2.canvas.draw()
            plt.pause(0.1)

        fit_arr.append((siegel_left_fit[0], siegel_best_pos - siegel_left_zero,
                        siegel_right_fit[0], siegel_best_pos - siegel_right_zero))

        print(fpos_arr[:midx+5], siegel_left_fit[0]*fpos_arr[:midx+5]+siegel_left_fit[1])

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

        f = open(os.path.join(imagesdir, 'vcurve_fits.json'), 'a')
        ls = siegel_left_fit[0]
        lp = siegel_best_pos - siegel_left_zero
        rs = siegel_right_fit[0]
        rp = siegel_best_pos - siegel_right_zero
        tstamp = time.strftime('%Y/%m/%d %H:%M:%S %Z')
        #f.write(f'{tstamp}, {ls}, {lp}, {rs}, {rp}\n')
        j = json.dumps({ 'timestamp' : tstamp, 'rightslope' : rs, 'rightpid' : rp, 'leftslope' : ls, 'leftpid' : lp})
        f.write(j + '\n')
        f.close()

    # write out
    f = open(os.path.join(imagesdir, f'hfd_run_{iter+1:03d}.txt'), 'w')
    for (ls, lp, rs, rp) in fit_arr:
        f.write(f'{ls}, {lp}, {rs}, {rp}\n')
    f.close()

    if args.debugplots:
        plt.show()
