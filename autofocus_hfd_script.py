import os
import sys
import time
import argparse
import logging

import astropy.io.fits as pyfits

import numpy as np
from scipy.stats import siegelslopes

import matplotlib.pyplot as plt

ASCOM_FOCUS_DRIVER = 'ASCOM.Simulator.Focuser'

#from pyastrobackend.SimpleDeviceInterface import connect_backend
#from pyastrobackend.SimpleDeviceInterface import connect_focuser
#from pyastrobackend.SimpleDeviceInterface import connect_camera
#from pyastrobackend.SimpleDeviceInterface import take_exposure
#from pyastrobackend.SimpleDeviceInterface import wait_on_focuser_move

from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI


from StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window

# for simulator
from c8_simul_star import C8_F7_Star_Simulator

def take_exposure_and_measure_star(cam, imgname, focus_expos):
    global args, fig, fig2, ax_1d, ax_2d, ax_hfd, ax_hfd

    if not args.simul:
        logging.info(f'Taking exposure exposure = {focus_expos} seconds')
        rc = sdi.take_exposure(cam, focus_expos, imgname)
        logging.info(f'exposure result code = {rc}')

    else:
        tmp_starimage_data = simul_star.get_simul_star_image(focus_pos)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)

    hdu = pyfits.open(imgname)
    starimage_data = hdu[0].data.astype(float)
    hdu.close()

    # analyze frame
    bg = 800
    thres = 10000

    xcen, ycen, bg, mad = find_star(starimage_data, debugfits=True)

    win = 100
    xlow = int(xcen-win/2)
    xhi = int(xcen+win/2)
    ylow = int(ycen-win/2)
    yhi = int(ycen+win/2)
    crop_data = starimage_data[ylow:yhi, xlow:xhi]

    if args.debugplots:
        fig.clear()
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)
        im = ax_2d.imshow((crop_data-bg).astype(float))
        fig.colorbar(im, ax=ax_2d)

    profile = horiz_bin_window(crop_data, bg=bg)

    if args.debugplots:
        ax_1d.plot(profile)

    return find_hfd_from_1D(profile, thres=thres)

def move_focuser(focuser, pos):
    if not focuser.move_absolute_position(pos):
        logging.error("Focuser error!!")
        sys.exit(1)

    sdi.wait_on_focuser_move(focuser)

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('focus_min', type=int, help='Lowest focus travel allowed')
    parser.add_argument('focus_max', type=int, help='Highest focus travel allowed')
    parser.add_argument('focus_dir', type=str, help='IN or OUT')
    parser.add_argument('nruns', type=int, help='Number of vcurve runs')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')
    parser.add_argument('--simul', action='store_true', help='Simulate star')
    parser.add_argument('--focuser_driver', type=str,  help='Focuser Driver')
    parser.add_argument('--camera_driver', type=str,  help='Camera Driver')
    parser.add_argument('--exposure_start', default=1, type=int,  help='Starting exposure value')
    parser.add_argument('--exposure_min', default=1, type=int,  help='Minimum exposure value')
    parser.add_argument('--exposure_max', default=8, type=int,  help='Maximum exposure value')
    parser.add_argument('--starflux_min', default=50000, type=int,  help='Maximum flux in star')

    return parser.parse_args()

if __name__ == '__main__':
    logging.basicConfig(filename='sample_vcurve.log',
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

    args = parse_commandline()
    logging.info(f'args = {args}')

    # connect focuser and camera
    if not args.simul:
        sdi = SDI()

        sdi.connect_backend()

        #focuser = connect_focuser(ASCOM_FOCUS_DRIVER)
        logging.info(f'Connecting to focuser driver {args.focuser_driver}')
        focuser = sdi.connect_focuser(args.focuser_driver)
        logging.info(f'focuser = {focuser}')
        if not focuser:
            logging.error(f'Unabled to connect to focuser driver {args.focuser_driver}')

        logging.info(f'Connecting to camera driver {args.camera_driver}')
        cam = sdi.connect_camera(args.camera_driver)
        logging.info(f'cam = {cam}')
        if not cam:
            logging.error(f'Unabled to connect to camera driver {args.camera_driver}')
    else:
        simul_star = C8_F7_Star_Simulator()

    # create output dir
    datestr = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    imagesdir = datestr
    os.mkdir(imagesdir)

    # create plots if needed
    if args.debugplots:
        logging.info('Creating figure')
        fig2 = plt.figure()
        ax_hfd = fig2.add_subplot(111)
        hfd_plot, = ax_hfd.plot([],[], marker='o', ls='')

        fig = plt.figure()
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)

    focus_expos = args.exposure_start




    focus_step = int((focus_end - focus_start)/(focus_nstep - 1))
    logging.info(f'Focus ex {focus_start} to {focus_end} step {focus_step}')
    fpos_arr = []
    hfd_arr = []
    for focus_pos in range(focus_start, focus_end+focus_step, focus_step):
        logging.info(f'Moving to focus pos {focus_pos}')
        imgname = os.path.join(imagesdir, f'vcurve_focuspos_run{iter+1:03d}_{focus_pos}.fit')

            if rc is not None:
                scen, sl, sr, hfl, hfr, totflux = rc
                if totflux >= args.starflux_min:
                    break
                else:
                    logging.warning('Star flux in image {imgname} ({totflux} is less than {args.starflux_min}')
            else:
                logging.warning(f'No star found in image {imgname}')

            focus_expos += 1
            if not args.simul:
                logging.info(f'Simulated image did not have star - check settings.')
                sys.exit(1)

            if focus_expos <= args.exposure_max:
                logging.info(f'Bumping exposure to {focus_expos} and retrying')
                continue
            else:
                logging.info(f'Exposure is already at max of {args.exposure_max} seconds - aborting')
                sys.exit(1)

        fpos_arr.append(focus_pos)
        hfd_arr.append(hfr-hfl)

        if args.debugplots:
            hfd_plot.set_data(fpos_arr, hfd_arr)
            ax_hfd.relim()
            ax_hfd.autoscale_view()
            fig2.canvas.draw()

            ax_1d.axvline(scen, color='red')
            if sl is not None and sr is not None:
                ax_1d.axvline(sl, color='green')
                ax_1d.axvline(sr, color='green')
                ax_1d.axvline(hfl, color='blue')
                ax_1d.axvline(hfr, color='blue')
                delta = sr-sl
                ax_1d.set_xlim(sl-delta/4, sr+delta/4)
                ax_1d.set_title(f'{hfr-hfl:5.3f}')
            print('drawing plot')
            fig.show()
            plt.pause(0.01)

        # write out
        f = open(os.path.join(imagesdir, f'hfd_run_{iter+1:03d}.txt'), 'w')
        for (ls, lp, rs, rp)  in fit_arr:
            f.write(f'{ls}, {lp}, {rs}, {rp}\n')
        f.close()

        # fit
        fpos_arr = np.array(fpos_arr)
        hfd_arr = np.array(hfd_arr)

        # sort so fpos is increasing
        fargs = fpos_arr.argsort()
        fpos_arr = fpos_arr[fargs]
        hfd_arr = hfd_arr[fargs]

        # find mininum value
        midx = np.argmin(hfd_arr)
        fpos_arr_l = np.array(fpos_arr[:midx-3])
        fpos_arr_r = np.array(fpos_arr[midx+4:])
        hfd_arr_l = np.array(hfd_arr[:midx-3])
        hfd_arr_r = np.array(hfd_arr[midx+4:])

        print('fpos_l', fpos_arr_l)
        print('hfd_l', hfd_arr_l)
        print('fpos_r', fpos_arr_r)
        print('hfd_r', hfd_arr_r)

        siegel_left_fit = siegelslopes(hfd_arr_l, fpos_arr_l)
        siegel_right_fit = siegelslopes(hfd_arr_r, fpos_arr_r)
        siegel_left_zero = -siegel_left_fit[1]/siegel_left_fit[0]
        siegel_right_zero = -siegel_right_fit[1]/siegel_right_fit[0]
        siegel_best_pos = (siegel_left_fit[1]-siegel_right_fit[1])/(siegel_right_fit[0]-siegel_left_fit[0])
        logging.info(f'siegel left  fit = {siegel_left_fit}')
        logging.info(f'siegel right fit = {siegel_right_fit}')
        logging.info(f'siegel best pos  = {siegel_best_pos}')

        if args.debugplots:
            ax_hfd.plot(fpos_arr_l, hfd_arr_l, marker='+', ls='', color='red')
            ax_hfd.plot(fpos_arr_r, hfd_arr_r, marker='+', ls='', color='red')
            ax_hfd.plot(fpos_arr[midx-5:], siegel_right_fit[0]*fpos_arr[midx-5:]+siegel_right_fit[1], color='green')
            ax_hfd.plot(fpos_arr[:midx+5], siegel_left_fit[0]*fpos_arr[:midx+5]+siegel_left_fit[1], color='blue')
            ax_hfd.axvline(siegel_best_pos, color='red')
            ax_hfd.set_title(f'Left {siegel_left_fit[0]:7.6f}/{siegel_best_pos - siegel_left_zero:5.3f} Right {siegel_right_fit[0]:7.6f}/{siegel_best_pos - siegel_right_zero:5.3f}')
            ax_hfd.relim()
            ax_hfd.set_ylim(bottom=0)
            ax_hfd.autoscale_view()
            fig2.canvas.draw()
            plt.pause(0.1)

        fit_arr.append((siegel_left_fit[0], siegel_best_pos - siegel_left_zero, siegel_right_fit[0], siegel_best_pos - siegel_right_zero))

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

    f = open(os.path.join(imagesdir, 'vcurve_fits.txt'), 'w')
    for (ls, lp, rs, rp)  in fit_arr:
        f.write(f'{ls}, {lp}, {rs}, {rp}\n')
    f.close()


    if args.debugplots:
        plt.show()