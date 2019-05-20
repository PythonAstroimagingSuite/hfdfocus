import os
import sys
import time
import argparse
import logging

import astropy.io.fits as pyfits

import numpy as np
from scipy.stats import siegelslopes

import matplotlib as mpl
mpl.rc('font', size=8)
import matplotlib.pyplot as plt

ASCOM_FOCUS_DRIVER = 'ASCOM.Simulator.Focuser'

#from pyastrobackend.SimpleDeviceInterface import connect_backend
#from pyastrobackend.SimpleDeviceInterface import connect_focuser
#from pyastrobackend.SimpleDeviceInterface import connect_camera
#from pyastrobackend.SimpleDeviceInterface import take_exposure
#from pyastrobackend.SimpleDeviceInterface import wait_on_focuser_move

#from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI


from StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window

# for simulator
from c8_simul_star import C8_F7_Star_Simulator

def measure_frame(starimage_data):
    global args, fig, fig2, ax_1d, ax_2d, ax_hfd, ax_hfd

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
        #mpl.rcParams.update({'axes.labelsize' : 18})
        ax_1d = fig.add_subplot(121)
        ax_2d = fig.add_subplot(122)
        im = ax_2d.imshow((crop_data-bg).astype(float))
        fig.colorbar(im, ax=ax_2d)

    profile = horiz_bin_window(crop_data, bg=bg)

    if args.debugplots:
        ax_1d.plot(profile)

    return find_hfd_from_1D(profile, thres=thres)

def take_exposure_and_measure_star(cam, imgname, focus_expos):
    global args, fig, fig2, ax_1d, ax_2d, ax_hfd, ax_hfd

    if not args.simul:
        logging.info(f'Taking exposure exposure = {focus_expos} seconds')
        rc = sdi.take_exposure(cam, focus_expos, imgname)
        logging.info(f'exposure result code = {rc}')

    hdu = pyfits.open(imgname)
    starimage_data = hdu[0].data.astype(float)
    hdu.close()

    return measure_frame(starimage_data)

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
    #parser.add_argument('nruns', type=int, help='Number of vcurve runs')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')
    parser.add_argument('--debugplotsdelay', type=float, default=1, help='Delay (seconds) showing each plot')
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
        from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI
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
#        fig2 = plt.figure()
#        ax_hfd = fig2.add_subplot(111)
#        hfd_plot, = ax_hfd.plot([],[], marker='o', ls='')


        fig = plt.figure(figsize=(4.5,2))
        ax_1d = fig.add_subplot(121)
        ax_2d = fig.add_subplot(122)

    focus_expos = args.exposure_start

    # vcurve parameters
#    vcurve_rs = 0.04462570065893874
#    vcurve_rp = 4.959828107421345
#    vcurve_ls =-0.04565294355776672
#    vcurve_lp =-4.848226361604247

    vcurve_rs = 0.049684986658347155
    vcurve_rp = -6.715524350563101
    vcurve_ls = -0.05181792218702318
    vcurve_lp = 6.4390991317122825

    if not args.simul:
        fpos_1 = focuser.get_absolute_position()
    else:
        fpos_1 = 8000

    imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos_1}.fit')
    if not args.simul:
        rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
    else:
        fpos_1 = 8000
        tmp_starimage_data = simul_star.get_simul_star_image(fpos_1)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    if rc is None:
        logging.error('No star found!')
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc
    hfd_1 = hfr-hfl

    if args.debugplots:
        fig.suptitle(f'First focus {fpos_1} HFD {hfd_1:5.2f}')
        fig.show()
        plt.pause(args.debugplotsdelay)

    # kind of redundant to take it again?
    # take two measurements and see if one the right side of curve
    imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos_1}.fit')
    if not args.simul:
        rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
    else:
        fpos_1 = 8000
        tmp_starimage_data = simul_star.get_simul_star_image(fpos_1)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    if rc is None:
        logging.error('No star found!')
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc
    hfd_1 = hfr-hfl
    logging.info(f'INITIAL FOCUS = {fpos_1}  HFD = {hfd_1}')

    # move out 10 HFD
    nsteps = int(10/vcurve_rs)
    fpos_2 = fpos_1 + nsteps
    logging.info(f'Moving out to {fpos_2}')

    imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos_2}.fit')
    if not args.simul:
        move_focuser(focuser, fpos_2)

        rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
    else:
        tmp_starimage_data = simul_star.get_simul_star_image(fpos_2)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    if rc is None:
        logging.error('No star found!')
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc
    hfd_2 = hfr-hfl
    logging.info(f'SECOND FOCUS = {fpos_2}  HFD = {hfd_2}')

    if args.debugplots:
        fig.suptitle(f'Second focus {fpos_2} HFD {hfd_2:5.2f}')
        fig.show()
        plt.pause(args.debugplotsdelay)

    # make sure hfd got larger
    if hfd_2 < hfd_1:
        logging.error('On wrong side of focus!')
        sys.exit(1)

    # compute location for desired initial HFD
    initial_hfd = 24
    fpos_initial = fpos_2 + int((initial_hfd-hfd_2)/vcurve_rs)
    logging.info(f'Initial HFD = {initial_hfd} pred focus = {fpos_initial}')

    # figure out direction
    backlash = 200
    if args.focus_dir == 'OUT':
        # start past desired start and move to it to remove backlash
        fpos_pre_initial = fpos_initial - backlash
    elif args.focus_dir == 'IN':
        # start past desired start and move to it to remove backlash
        fpos_pre_initial = fpos_initial + backlash
    else:
        logging.error(f'Unknown focus directin {args.focus_dir} - exitting!')
        sys.exit(1)

    imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos_initial}.fit')
    if not args.simul:
        logging.info(f'Moving to pre spot {fpos_pre_initial}')
        move_focuser(focuser, fpos_pre_initial)

        # move out to initial HFD
        logging.info(f'Moving to {fpos_initial}')
        time.sleep(0.5)
        move_focuser(focuser, fpos_initial)

        rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
    else:
        tmp_starimage_data = simul_star.get_simul_star_image(fpos_initial)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    if rc is None:
        logging.error('No star found!')
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc
    hfd_initial = hfr-hfl

    logging.info(f'INITIAL POSITION FOCUS = {fpos_initial}  HFD = {hfd_initial}')

    if args.debugplots:
        fig.suptitle(f'Initial position focus {fpos_initial} HFD {hfd_initial:5.2f}')
        fig.show()
        plt.pause(args.debugplotsdelay)

    # now take several measurements and get average HFD
    avg_initial_hfd = 0
    ninitial = 5
    for i in range(0, ninitial):
        imgname = os.path.join(imagesdir, f'vcurve_focuspos_initial{i:02d}_{fpos_initial}.fit')
        if not args.simul:
            rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
        else:
            tmp_starimage_data = simul_star.get_simul_star_image(fpos_initial)
            pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
            rc = measure_frame(tmp_starimage_data)

        if rc is None:
            logging.error('No star found!')
            sys.exit(1)

        scen, sl, sr, hfl, hfr, totflux = rc
        avg_initial_hfd += hfr-hfl

        if args.debugplots:
            fig.suptitle(f'Initial pos iter #{i+1} focus {fpos_initial} ' \
                         f'HFD {hfr-hfl:5.2f} AVG:{avg_initial_hfd/(i+1):5.2f}')
            fig.show()
            plt.pause(args.debugplotsdelay)

    avg_initial_hfd /= ninitial

    logging.info(f'INITIAL POSITION FOCUS AVERAGE HFD = {avg_initial_hfd}')

    # now based on average at initial compute inner HFD location
    # compute location for desired initial HFD
    inner_hfd = 12
    fpos_inner = fpos_2 + int((inner_hfd-hfd_2)/vcurve_rs)
    logging.info(f'Inner HFD = {inner_hfd} pred focus = {fpos_inner}')

    imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos_inner}.fit')
    if not args.simul:
        # move out to inner HFD
        logging.info(f'Moving to {fpos_inner}')
        time.sleep(0.5)
        move_focuser(focuser, fpos_inner)

        rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
    else:
        tmp_starimage_data = simul_star.get_simul_star_image(fpos_inner)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    if rc is None:
        logging.error('No star found!')
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc
    hfd_inner = hfr-hfl

    logging.info(f'INNER POSITION FOCUS = {fpos_inner}  HFD = {hfd_inner}')

    if args.debugplots:
        fig.suptitle(f'Inner position focus {fpos_inner} HFD {hfd_inner:5.2f}')
        fig.show()
        plt.pause(args.debugplotsdelay)

    # now take several measurements and get average HFD
    avg_inner_hfd = 0
    ninner = 5
    for i in range(0, ninitial):
        imgname = os.path.join(imagesdir, f'vcurve_focuspos_inner{i:02d}_{fpos_initial}.fit')
        if not args.simul:
            rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
        else:
            tmp_starimage_data = simul_star.get_simul_star_image(fpos_inner)
            pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
            rc = measure_frame(tmp_starimage_data)

        if rc is None:
            logging.error('No star found!')
            sys.exit(1)

        scen, sl, sr, hfl, hfr, totflux = rc
        avg_inner_hfd += hfr-hfl

        if args.debugplots:
            fig.suptitle(f'Inner pos iter #{i+1} focus {fpos_inner} ' \
                         f'HFD {hfr-hfl:5.2f} AVG:{avg_inner_hfd/(i+1):5.2f}')
            fig.show()
            plt.pause(args.debugplotsdelay)

    avg_inner_hfd /= ninitial

    logging.info(f'INNER POSITION FOCUS AVERAGE HFD = {avg_inner_hfd}')

    # now compute best focus
    fpos_best = int(fpos_inner - int(avg_inner_hfd/vcurve_rs) + vcurve_rp)

    logging.info(f'BEST FOCUS POSITION = {fpos_best} {fpos_inner} {int(avg_inner_hfd/vcurve_rs)} {vcurve_rp}')

    if fpos_best > args.focus_max or fpos_best > args.focus_min:
        logging.error(f'Best focus position {fpos_best} is outside allowed range {args.focus_min} to {args.focus_max}')

    if not args.simul:
        logging.info(f'Moving to {fpos_best}')
        move_focuser(focuser, fpos_best)
        imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos_best}.fit')

        rc = take_exposure_and_measure_star(cam, imgname, focus_expos)
    else:
        tmp_starimage_data = simul_star.get_simul_star_image(fpos_best)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    if rc is None:
        logging.error('No star found!')
        if args.debugplots:
            plt.show()
        sys.exit(1)

    scen, sl, sr, hfl, hfr, totflux = rc
    best_hfd = hfr-hfl

    logging.info(f'BEST FOCUS POSITION = {fpos_best} HFD = {best_hfd}')

    if args.debugplots:
        fig.suptitle(f'Best pos focus {fpos_best} HFD {best_hfd:5.2f}')
        fig.show()
        plt.show()

#        if args.debugplots:
#            hfd_plot.set_data(fpos_arr, hfd_arr)
#            ax_hfd.relim()
#            ax_hfd.autoscale_view()
#            fig2.canvas.draw()
#
#            ax_1d.axvline(scen, color='red')
#            if sl is not None and sr is not None:
#                ax_1d.axvline(sl, color='green')
#                ax_1d.axvline(sr, color='green')
#                ax_1d.axvline(hfl, color='blue')
#                ax_1d.axvline(hfr, color='blue')
#                delta = sr-sl
#                ax_1d.set_xlim(sl-delta/4, sr+delta/4)
#                ax_1d.set_title(f'{hfr-hfl:5.3f}')
#            print('drawing plot')
#            fig.show()
#            plt.pause(0.01)
#
#        # write out
#        f = open(os.path.join(imagesdir, f'hfd_run_{iter+1:03d}.txt'), 'w')
#        for (ls, lp, rs, rp)  in fit_arr:
#            f.write(f'{ls}, {lp}, {rs}, {rp}\n')
#        f.close()
#
#        # fit
#        fpos_arr = np.array(fpos_arr)
#        hfd_arr = np.array(hfd_arr)
#
#        # sort so fpos is increasing
#        fargs = fpos_arr.argsort()
#        fpos_arr = fpos_arr[fargs]
#        hfd_arr = hfd_arr[fargs]
#
#        # find mininum value
#        midx = np.argmin(hfd_arr)
#        fpos_arr_l = np.array(fpos_arr[:midx-3])
#        fpos_arr_r = np.array(fpos_arr[midx+4:])
#        hfd_arr_l = np.array(hfd_arr[:midx-3])
#        hfd_arr_r = np.array(hfd_arr[midx+4:])
#
#        print('fpos_l', fpos_arr_l)
#        print('hfd_l', hfd_arr_l)
#        print('fpos_r', fpos_arr_r)
#        print('hfd_r', hfd_arr_r)
#
#        siegel_left_fit = siegelslopes(hfd_arr_l, fpos_arr_l)
#        siegel_right_fit = siegelslopes(hfd_arr_r, fpos_arr_r)
#        siegel_left_zero = -siegel_left_fit[1]/siegel_left_fit[0]
#        siegel_right_zero = -siegel_right_fit[1]/siegel_right_fit[0]
#        siegel_best_pos = (siegel_left_fit[1]-siegel_right_fit[1])/(siegel_right_fit[0]-siegel_left_fit[0])
#        logging.info(f'siegel left  fit = {siegel_left_fit}')
#        logging.info(f'siegel right fit = {siegel_right_fit}')
#        logging.info(f'siegel best pos  = {siegel_best_pos}')
#
#        if args.debugplots:
#            ax_hfd.plot(fpos_arr_l, hfd_arr_l, marker='+', ls='', color='red')
#            ax_hfd.plot(fpos_arr_r, hfd_arr_r, marker='+', ls='', color='red')
#            ax_hfd.plot(fpos_arr[midx-5:], siegel_right_fit[0]*fpos_arr[midx-5:]+siegel_right_fit[1], color='green')
#            ax_hfd.plot(fpos_arr[:midx+5], siegel_left_fit[0]*fpos_arr[:midx+5]+siegel_left_fit[1], color='blue')
#            ax_hfd.axvline(siegel_best_pos, color='red')
#            ax_hfd.set_title(f'Left {siegel_left_fit[0]:7.6f}/{siegel_best_pos - siegel_left_zero:5.3f} Right {siegel_right_fit[0]:7.6f}/{siegel_best_pos - siegel_right_zero:5.3f}')
#            ax_hfd.relim()
#            ax_hfd.set_ylim(bottom=0)
#            ax_hfd.autoscale_view()
#            fig2.canvas.draw()
#            plt.pause(0.1)
#
#        fit_arr.append((siegel_left_fit[0], siegel_best_pos - siegel_left_zero, siegel_right_fit[0], siegel_best_pos - siegel_right_zero))
#
#        print(fpos_arr[:midx+5], siegel_left_fit[0]*fpos_arr[:midx+5]+siegel_left_fit[1])
#
#        logging.info('Left Side:')
#        logging.info(f'   slope: {siegel_left_fit[0]}')
#        logging.info(f'   inter: {siegel_left_fit[1]}')
#        logging.info(f'   yzero: {siegel_left_zero}')
#        logging.info(f'   PID  : {siegel_best_pos - siegel_left_zero}')
#        logging.info('Right Side:')
#        logging.info(f'   slope: {siegel_right_fit[0]}')
#        logging.info(f'   inter: {siegel_right_fit[1]}')
#        logging.info(f'   yzero: {siegel_right_zero}')
#        logging.info(f'   PID  : {siegel_best_pos - siegel_right_zero}')
#
#    f = open(os.path.join(imagesdir, 'vcurve_fits.txt'), 'w')
#    for (ls, lp, rs, rp)  in fit_arr:
#        f.write(f'{ls}, {lp}, {rs}, {rp}\n')
#    f.close()
#
#
#    if args.debugplots:
#        plt.show()