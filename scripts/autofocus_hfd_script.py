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

from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI
from pyastroprofile.EquipmentProfile import EquipmentProfile

#ASCOM_FOCUS_DRIVER = 'ASCOM.Simulator.Focuser'

#from pyastrobackend.SimpleDeviceInterface import connect_backend
#from pyastrobackend.SimpleDeviceInterface import connect_focuser
#from pyastrobackend.SimpleDeviceInterface import connect_camera
#from pyastrobackend.SimpleDeviceInterface import take_exposure
#from pyastrobackend.SimpleDeviceInterface import wait_on_focuser_move

#from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI


from hfdfocus.StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window

# for simulator
from hfdfocus.c8_simul_star import C8_F7_Star_Simulator

def measure_frame(starimage_data):
    # analyze frame
    bg = 800
    thres = 10000

    xcen, ycen, bg, mad, starmask, alone = find_star(starimage_data, debugfits=True)

    if np.max(starimage_data[starmask] > args.saturation):
        logging.warning(f'SATURATED PIXELS DETECTED!')

    win = args.winsize
    xlow = max(0, int(xcen-win/2))
    xhi = min(starimage_data.shape[0]-1, int(xcen+win/2))
    ylow = max(0, int(ycen-win/2))
    yhi = min(starimage_data.shape[1]-1, int(ycen+win/2))
    logging.debug(f'cropping to window={win} x={xlow}:{xhi} y={ylow}:{yhi}')
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

    rc = find_hfd_from_1D(profile, thres=thres)

    if rc is not None:
        scen, sl, sr, hfl, hfr, totflux = rc
        return hfr-hfl
    else:
        return None

def take_exposure_and_measure_star(imgname, focus_expos):
    if not args.simul:
        if args.framesize != 0:
            w, h = cam.get_size()
            xl = int(w/2-args.framesize/2)
            xw = args.framesize
            yt = int(h/2-args.framesize/2)
            yh = args.framesize
            roi = (xl, yt, xw, yh)
        else:
            roi = None

        logging.info(f'Taking exposure exposure = {focus_expos} seconds roi = {roi}')

        rc = sdi.take_exposure(cam, focus_expos, imgname, roi=roi)

#        logging.info(f'Taking exposure exposure = {focus_expos} seconds')
#        rc = sdi.take_exposure(cam, focus_expos, imgname)

        logging.info(f'exposure result code = {rc}')

    hdu = pyfits.open(imgname)
    starimage_data = hdu[0].data.astype(float)
    hdu.close()

    return measure_frame(starimage_data)

def move_focuser(pos):
    if pos > args.focus_max or pos < args.focus_min:
        logging.error(f'Requested focus position {pos} is outside allowed range {args.focus_min} to {args.focus_max}!!!')
        sys.exit(1)

    if not focuser.move_absolute_position(pos):
        logging.error("Focuser error!!")
        sys.exit(1)

    sdi.wait_on_focuser_move(focuser)

def measure_at_focus_pos(fpos, focus_expos):
    if not args.simul:
        if focuser.get_absolute_position() != fpos:
            logging.info(f'Moving focuser to {fpos}')
            move_focuser(fpos)
            time.sleep(args.focusdelay)

    imgname = os.path.join(imagesdir, f'vcurve_focuspos_{fpos}.fit')
    if not args.simul:
        rc = take_exposure_and_measure_star(imgname, focus_expos)
    else:
        tmp_starimage_data = simul_star.get_simul_star_image(fpos)
        pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
        rc = measure_frame(tmp_starimage_data)

    return rc

def average_measure_at_focus_pos(fpos, focus_expos, niter, tag=''):
    if not args.simul:
        if focuser.get_absolute_position() != fpos:
            logging.info(f'Moving focuser to {fpos}')
            move_focuser(fpos)
            time.sleep(args.focusdelay)

    avg_hfd = 0
    ncap = 0
    for i in range(0, niter):
        imgname = os.path.join(imagesdir, f'vcurve_focuspos_{tag}{i:02d}_{fpos}.fit')
        if not args.simul:
            hfd = take_exposure_and_measure_star(imgname, focus_expos)
        else:
            tmp_starimage_data = simul_star.get_simul_star_image(fpos)
            pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
            hfd= measure_frame(tmp_starimage_data)

        if hfd is None:
            logging.error('No star found - continuing!')
            continue

        ncap += 1
        avg_hfd += hfd

        if args.debugplots:
            fig.suptitle(f'{tag} pos iter #{i+1} focus {fpos} ' \
                         f'HFD {hfd:5.2f} AVG:{avg_hfd/(ncap):5.2f}')
            fig.show()
            plt.pause(args.debugplotsdelay)

    if ncap > 0:
        return avg_hfd/ncap
    else:
        return None

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('focus_min', type=int, help='Lowest focus travel allowed')
    parser.add_argument('focus_max', type=int, help='Highest focus travel allowed')
    parser.add_argument('focus_dir', type=str, help='IN or OUT')
    parser.add_argument('--focus_start', type=int, help='Starting focus pos')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')
    parser.add_argument('--debugplotsdelay', type=float, default=1, help='Delay (seconds) showing each plot')
    parser.add_argument('--simul', action='store_true', help='Simulate star')
    parser.add_argument('--stayopen', action='store_true', help='stay open when done')
    parser.add_argument('--profile', type=str, help='Name of equipment profile')
    parser.add_argument('--focuser', type=str,  help='Focuser Driver')
    parser.add_argument('--camera', type=str,  help='Camera Driver')
    parser.add_argument('--exposure_start', default=1, type=int,  help='Starting exposure value')
    parser.add_argument('--exposure_min', default=1, type=int,  help='Minimum exposure value')
    parser.add_argument('--exposure_max', default=8, type=int,  help='Maximum exposure value')
    parser.add_argument('--starflux_min', default=50000, type=int,  help='Maximum flux in star')
    parser.add_argument('--saturation', default=55000, type=int,  help='Saturation level for sensor')
    parser.add_argument('--framesize', default=800, type=int,  help='Size of capture frame, 0=full')
    parser.add_argument('--winsize', default=250, type=int,  help='Size of window used to analyze star')
    parser.add_argument('--focusdelay', default=0.5, type=float,  help='Delay (seconds) after focus moves')
    parser.add_argument('--numaverage', default=5, type=int,  help='Number of images to average')

    args = parser.parse_args()
    return args


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

        # load profile
        if args.profile is not None:
            logging.info(f'Using equipment profile {args.profile}')
            equip_profile = EquipmentProfile('astroprofiles/equipment', args.profile)
            equip_profile.read()
            camera_driver = equip_profile.camera_driver
            focuser_driver = equip_profile.focuser_driver
        else:
            focuser_driver = args.focuser
            camera_driver = args.camera

        sdi = SDI()

        sdi.connect_backend()

        #focuser = connect_focuser(ASCOM_FOCUS_DRIVER)
        logging.info(f'Connecting to focuser driver {focuser_driver}')
        focuser = sdi.connect_focuser(focuser_driver)
        logging.info(f'focuser = {focuser}')
        if not focuser:
            logging.error(f'Unabled to connect to focuser driver {focuser_driver}')
            sys.exit(-1)
        logging.info(f'Connecting to camera driver {camera_driver}')
        cam = sdi.connect_camera(camera_driver)
        logging.info(f'cam = {cam}')
        if not cam:
            logging.error(f'Unabled to connect to camera driver {camera_driver}')
            sys.exit(-1)
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
    logging.info(f'Starting exposure is {focus_expos} seconds')

    # vcurve parameters
#    vcurve_rs = 0.04462570065893874
#    vcurve_rp = 4.959828107421345
#    vcurve_ls =-0.04565294355776672
#    vcurve_lp =-4.848226361604247

    # from 2019/05/20 C8 @ f/7 - seemed to work well with L filter (only tested)
    vcurve_rs = 0.04472660490995652
    vcurve_rp = 5.333072819832182
    vcurve_ls = -0.045806498798410436
    vcurve_lp = -5.22113594480723

    if args.focus_start is None:
        if not args.simul:
            fpos_1 = focuser.get_absolute_position()
        else:
            fpos_1 = 8000
    else:
        fpos_1 = args.focus_start

    if args.focus_dir == 'OUT':
        fdir = 1
        vslope = vcurve_ls
        vpid = vcurve_lp
    elif args.focus_dir == 'IN':
        fdir = -1
        vslope = vcurve_rs
        vpid = vcurve_rp
    else:
        logging.error(f'Unknown focus directin {args.focus_dir} - exitting!')
        sys.exit(1)

    ntries = 0
    right_side = False
    while ntries < 3:
        hfd_1 = measure_at_focus_pos(fpos_1, focus_expos)

        if hfd_1 is None:
            logging.error('No star found!')
            if args.debugplots:
                plt.show()
            sys.exit(1)

        if args.debugplots:
            fig.suptitle(f'First focus {fpos_1} HFD {hfd_1:5.2f}')
            fig.show()
            plt.pause(args.debugplotsdelay)

        logging.info(f'INITIAL FOCUS = {fpos_1}  HFD = {hfd_1}')

        # move out 10 HFD
        nsteps = int(abs(10/vslope))
        fpos_2 = fpos_1 - fdir*nsteps

        hfd_2 = measure_at_focus_pos(fpos_2, focus_expos)

        if hfd_2 is None:
            logging.error('No star found!')
            if args.debugplots:
                plt.show()
            sys.exit(1)

        logging.info(f'SECOND FOCUS = {fpos_2}  HFD = {hfd_2}')

        if args.debugplots:
            fig.suptitle(f'Second focus {fpos_2} HFD {hfd_2:5.2f}')
            fig.show()
            plt.pause(args.debugplotsdelay)

        # make sure hfd got larger
        if hfd_2 < hfd_1:
            logging.error('On wrong side of focus!')

            #  bump to near zero point
            fpos_1 = fpos_2 - fdir*int(abs(hfd_2/vslope))
            ntries += 1
            continue

        right_side = True
        break

    if not right_side:
        logging.error('Could not find right side for focusing')
        sys.exit(1)

    # compute location for desired initial HFD
    initial_hfd = 24
    fpos_initial = fpos_2 - fdir*int(abs((initial_hfd-hfd_2)/vslope))
    logging.info(f'Initial HFD = {initial_hfd} pred focus = {fpos_initial}')

    # figure out direction
    backlash = 200
    fpos_pre_initial = fpos_initial - fdir*backlash

    logging.info(f'Moving to pre-initial pos {fpos_pre_initial}')
    if not args.simul:
        move_focuser(fpos_pre_initial)
        time.sleep(0.5)

#    hfd_initial = measure_at_focus_pos(fpos_pre_initial, focus_expos)
#
#    if hfd_initial is None:
#        logging.error('No star found!')
#        if args.debugplots:
#            plt.show()
#        sys.exit(1)
#
#    logging.info(f'INITIAL POSITION FOCUS = {fpos_initial}  HFD = {hfd_initial}')
#
#    if args.debugplots:
#        fig.suptitle(f'Initial position focus {fpos_initial} HFD {hfd_initial:5.2f}')
#        fig.show()
#        plt.pause(args.debugplotsdelay)

    # now take several measurements and get average HFD
    avg_initial_hfd = average_measure_at_focus_pos(fpos_initial, focus_expos, args.numaverage, tag='initial')

    logging.info(f'INITIAL POSITION FOCUS AVERAGE HFD = {avg_initial_hfd}')

    # now based on average at initial compute inner HFD location
    # compute location for desired initial HFD
    inner_hfd = 12
    fpos_inner = fpos_initial + fdir*int(abs((inner_hfd-avg_initial_hfd)/vslope))
    logging.info(f'Inner HFD = {inner_hfd} pred focus = {fpos_inner}')

    hfd_inner = measure_at_focus_pos(fpos_inner, focus_expos)

    if hfd_inner is None:
        logging.error('No star found!')
        if args.debugplots:
            plt.show()
        sys.exit(1)

    logging.info(f'INNER POSITION FOCUS = {fpos_inner}  HFD = {hfd_inner}')

    if args.debugplots:
        fig.suptitle(f'Inner position focus {fpos_inner} HFD {hfd_inner:5.2f}')
        fig.show()
        plt.pause(args.debugplotsdelay)

    # now take several measurements and get average HFD
    avg_inner_hfd = average_measure_at_focus_pos(fpos_inner, focus_expos, args.numaverage, tag='inner')

    logging.info(f'INNER POSITION FOCUS AVERAGE HFD = {avg_inner_hfd}')

    # now compute best focus
    fpos_best = int(fpos_inner + fdir*int(abs(avg_inner_hfd/vslope)) + vpid)

    logging.info(f'BEST FOCUS POSITION = {fpos_best} {fpos_inner} {int(avg_inner_hfd/vslope)} {vpid}')

    if fpos_best > args.focus_max or fpos_best < args.focus_min:
        logging.error(f'Best focus position {fpos_best} is outside allowed range {args.focus_min} to {args.focus_max}')
    else:
        best_hfd = measure_at_focus_pos(fpos_best, focus_expos)

        if best_hfd is None:
            logging.error('No star found!')
            if args.debugplots:
                plt.show()
            sys.exit(1)

        logging.info(f'BEST FOCUS POSITION = {fpos_best} HFD = {best_hfd}')

        if args.debugplots:
            fig.suptitle(f'Best pos focus {fpos_best} HFD {best_hfd:5.2f}')

    # keep plots up until keypress
    if args.debugplots:
        if args.stayopen:
            plt.show()
        else:
            plt.pause(5)
