import os
import sys
import math
import time
import argparse
import logging
import tempfile
from datetime import datetime

import astropy.io.fits as pyfits

import numpy as np
#from scipy.stats import siegelslopes

import matplotlib as mpl
#mpl.use("Qt5agg")
mpl.use("TkAgg")
mpl.rcParams['toolbar'] = 'None'
mpl.rc('font', size=8)
import matplotlib.pyplot as plt

from pyastrobackend.SimpleDeviceInterface import SimpleDeviceInterface as SDI
from pyastroprofile.AstroProfile import AstroProfile

#from pyastroprofile.EquipmentProfile import EquipmentProfile

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

FOCUSER_MIN_POS = None
FOCUSER_MAX_POS = None
FOCUSER_DIR = None


def show_fig_and_wait(fig, wait):
    fig.show()
    plt.pause(wait)
#    fig.canvas.flush_events()
#    time.sleep(wait)

def measure_frame(starimage_data):
    # analyze frame
    bg = 800
    thres = 10000

    rc_find = find_star(starimage_data, debugfits=False)

    if rc_find is None:
        logging.error('measure_frame: find_star failed!')
        return None

    xcen, ycen, bg, mad, starmask, alon = rc_find

    satur = False
    logging.info(f'Star Max Pixel Value = {np.max(starimage_data[starmask])}')
    if np.max(starimage_data[starmask] > args.saturation):
        logging.warning(f'SATURATED PIXELS DETECTED! Value = {np.max(starimage_data[starmask])}')
        satur = True

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

    rc = find_hfd_from_1D(profile, thres=thres, debugplots=False)

    if rc is not None:
        scen, sl, sr, hfl, hfr, totflux = rc
        if args.debugplots:
            #hfw = 2*(hfr-hfl)
            #hfc = (hfr+hfl)/2
            #ax_1d.set_xlim(hfc-hfw, hfc+hfw)

            ax_1d.axvline(scen, color='red')
            if sl is not None and sr is not None:
                ax_1d.axvline(sl, color='green')
                ax_1d.axvline(sr, color='green')
                ax_1d.axvline(hfl, color='blue')
                ax_1d.axvline(hfr, color='blue')
                delta = sr-sl
                ax_1d.set_xlim(sl-delta/4, sr+delta/4)

        return hfr-hfl, satur, rc
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

        # FIXME remove file if exists - had problem with taking multiple
        #       exposures at same focus pos because overwrite didnt work
        try:
            os.remove(imgname)
        except:
            pass

        rc = sdi.take_exposure(cam, focus_expos, imgname, roi=roi)

#        logging.info(f'Taking exposure exposure = {focus_expos} seconds')
#        rc = sdi.take_exposure(cam, focus_expos, imgname)

        logging.info(f'exposure result code = {rc}')

    hdu = pyfits.open(imgname)
    starimage_data = hdu[0].data.astype(float)
    hdu.close()

    return measure_frame(starimage_data)

def move_focuser(pos):
    if pos > FOCUSER_MAX_POS or pos < FOCUSER_MIN_POS:
        logging.error(f'Requested focus position {pos} is outside allowed range {FOCUSER_MIN_POS} to {FOCUSER_MAX_POS}!!!')
        sys.exit(1)

    if not focuser.move_absolute_position(pos):
        logging.error("Focuser error!!")
        sys.exit(1)

    sdi.wait_on_focuser_move(focuser)

def measure_at_focus_pos(fpos, focus_expos):
    if not args.simul or args.forcehw:
        if focuser.get_absolute_position() != fpos:
            logging.info(f'Moving focuser to {fpos}')
            move_focuser(fpos)
            time.sleep(args.focusdelay)

    use_expos = focus_expos
    for i in range(0, 4):

        imgname = os.path.join(IMAGESDIR, f'vcurve_focuspos_{fpos}.fit')
        if not args.simul:
            rc = take_exposure_and_measure_star(imgname, use_expos)
        else:
            tmp_starimage_data = simul_star.get_simul_star_image(fpos)
            pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
            rc = measure_frame(tmp_starimage_data)

        if rc is None:
            use_expos *= 2
            logging.info(f'No star - bumping exposusure to {use_expos}')
        else:
            satur = rc[1]  # hfr-hfl, satur, rc
            if satur:
                use_expos /= 2
                logging.info(f'Saturated star - dropping exposusure to {use_expos}')
            else:
                break


    return rc

def average_measure_at_focus_pos(fpos, focus_expos, niter, tag=''):
    if not args.simul or args.forcehw:
        if focuser.get_absolute_position() != fpos:
            logging.info(f'Moving focuser to {fpos}')
            move_focuser(fpos)
            time.sleep(args.focusdelay)

    avg_hfd = 0
    ncap = 0
    for i in range(0, niter):
        imgname = os.path.join(IMAGESDIR, f'vcurve_focuspos_{tag}{i:02d}_{fpos}.fit')
        if not args.simul:
            rc_focus = take_exposure_and_measure_star(imgname, focus_expos)
            if rc_focus is None:
                return None
            hfd, satur, focus_result = rc_focus
        else:
            tmp_starimage_data = simul_star.get_simul_star_image(fpos)
            pyfits.writeto(imgname, tmp_starimage_data.astype(float), overwrite=True)
            hfd, satur, focus_result = measure_frame(tmp_starimage_data)

        if hfd is None:
            logging.error('No star found - continuing!')
            continue

        ncap += 1
        avg_hfd += hfd

        msg = f'{tag} pos iter #{i+1} focus {fpos} ' \
              f'HFD {hfd:5.2f} AVG:{avg_hfd/(ncap):5.2f}'
        logging.info(msg)
        if args.debugplots:
            fig.suptitle(msg)
#            fig.show()
#            plt.pause(args.debugplotsdelay)
            show_fig_and_wait(fig, args.debugplotsdelay)

    if ncap > 0:
        return avg_hfd/ncap
    else:
        return None

def determine_final_hfd(fpos_best, best_expos):
    best_tries = 0
    s_max = 1.0
    s_min = 0.01
    first = True
    final_hfd = None
    while True:
        s_fact = (s_max + s_min)/2
        logging.info(f'scale min/max/cur: {s_min} {s_max} {s_fact}')
        rc = measure_at_focus_pos(fpos_best, s_fact*best_expos)
        if rc is None:
            if first:
                return None
#                logging.error('No star found!')
#                if not args.simul or args.forcehw:
#                    restore_focus_pos(starting_focus_pos)
#                cleanup_files()
#                logging.info('Exitting!')
#                sys.exit(1)

            logging.error(f'Measurement failed!')
            logging.warning(f'Setting s_min to {s_fact}')
            s_min = s_fact
        else:
            # show result
            if args.debugplots:
                plt.pause(0.05)

            best_hfd, satur, fresult = rc
            scen, sl, sr, hfl, hfr, totflux = fresult
            logging.info(f'Flux  : {totflux}')
            logging.info(f'HFD   : {best_hfd}')
            logging.info(f'Satur : {satur}')

            final_hfd = best_hfd

            if args.debugplots:
                if satur:
                    s_str = 'SATUR'
                else:
                    s_str = 'UNSAT'
                fig.suptitle(f'Test best focus {fpos_best} HFD {best_hfd:5.2f} {s_str}')
                plt.pause(0.05)

            if not satur:
                final_hfd = best_hfd
                break
            else:
                logging.warning(f'Saturated setting s_max to {s_fact}')
                s_max = s_fact

            first = False

        best_tries += 1
        if best_tries > 5:
            logging.error(f'Unble to take best focus exposure after {best_tries-1} tries!.')
            break

    return final_hfd, satur


def measure_file_only(image_file, debug_scale_factor=1.0):
    m_hdu = pyfits.open(image_file)
    m_starimage_data = m_hdu[0].data.astype(float)
    m_hdu.close()
    m_med = np.median(m_starimage_data)

    rc = measure_frame(m_starimage_data)
    logging.info(f'Measure frame for {args.measurefile} returned {rc}')
    if rc is not None:
        hfd, satur, fresult = rc
        scen, sl, sr, hfl, hfr, totflux = fresult
        logging.info(f'Measurement report:')
        logging.info(f'Median: {m_med}')
        logging.info(f'Flux  : {totflux}')
        logging.info(f'HFD   : {hfd}')
    else:
        logging.error(f'Measurement failed!')
        return False
    return True

def test_exposure_scaling(image_file):
    m_hdu = pyfits.open(image_file)
    m_starimage_data = m_hdu[0].data.astype(float)
    m_hdu.close()
    m_med = np.median(m_starimage_data)

    #
    # testing - simulate more/less exposure
    #
    # guess based on ONTC focused frames
    s_min = 0.01
    s_max = 16.0
    e_min = 0.1
    while True:
        scale_factor = (s_min+s_max)/2

        s_starimage_data = (m_starimage_data-m_med)*scale_factor + m_med

        rc = measure_frame(s_starimage_data)
        logging.info(f'Measure frame for {args.measurefile} returned {rc}')
        logging.info(f'Measurement report:')
        logging.info(f'Scale Factor: {s_min} {s_max} {scale_factor}')
        logging.info(f'Median: {m_med}')

        if rc is not None:
            hfd, satur, fresult = rc
            scen, sl, sr, hfl, hfr, totflux = fresult

            logging.info(f'Flux  : {totflux}')
            logging.info(f'HFD   : {hfd}')
            logging.info(f'Flux/PI/HFD^2/4: {totflux/hfd/hfd/math.pi/4}')

            if not satur:
                return True
            else:
                logging.warning(f'Saturated setting s_max to {scale_factor}')
                s_max = scale_factor
        else:
            logging.error(f'Measurement failed!')
            logging.warning(f'Setting s_min to {scale_factor}')
            s_min = scale_factor


def cleanup_files():
    logging.info(f'Cleaning up temp dir {IMAGESDIROBJ}')
    IMAGESDIROBJ.cleanup()
#    if args.keepfiles:
#        logging.debug('Keeping focus image files.')
#        return
#
#    if IMAGESDIR is not None and IMAGESDIR != '':
#        if os.path.isdir(IMAGESDIR):
#            logging.debug('Removing focus files')
#            fit_files = glob.glob(os.path.join(IMAGESDIR, '*.fit'))
#            fits_files = glob.glob(os.path.join(IMAGESDIR, '*.fits'))
#            log_files = glob.glob(os.path.join(IMAGESDIR, '*.log'))
#
#            files = list(set(fit_files + fits_files + log_files))
#            for f in files:
#                logging.info(f'Would have removed {f}')
#
#            logging.info(f'would have rmdir {IMAGESDIR}')

def restore_focus_pos(pos):
    logging.info(f'Returning focuser to initial position {pos}')
    return move_focuser(pos)

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('--focus_min', type=int, help='Lowest focus travel allowed')
    parser.add_argument('--focus_max', type=int, help='Highest focus travel allowed')
    parser.add_argument('--focus_dir', type=str, help='IN or OUT')
    parser.add_argument('--focus_start', type=int, help='Starting focus pos')
    parser.add_argument('--debug', action='store_true', help='Show debug output on console')
    parser.add_argument('--debugplots', action='store_true', help='Show debug plots')
    parser.add_argument('--debugplotsdelay', type=float, default=0.25, help='Delay (seconds) showing each plot')
    parser.add_argument('--simul', action='store_true', help='Simulate star')
    parser.add_argument('--stayopen', action='store_true', help='stay open when done')
    parser.add_argument('--profile', type=str, help='Name of equipment profile')
    parser.add_argument('--focuser', type=str,  help='Focuser Driver')
    parser.add_argument('--camera', type=str,  help='Camera Driver')
    parser.add_argument('--exposure_start', default=3, type=float,  help='Starting exposure value')
    parser.add_argument('--exposure_min', default=1, type=float,  help='Minimum exposure value')
    parser.add_argument('--exposure_max', default=8, type=float,  help='Maximum exposure value')
    parser.add_argument('--starflux_min', default=50000, type=int,  help='Maximum flux in star')
    parser.add_argument('--saturation', default=55000, type=int,  help='Saturation level for sensor')
    parser.add_argument('--framesize', default=800, type=int,  help='Size of capture frame, 0=full')
    parser.add_argument('--winsize', default=250, type=int,  help='Size of window used to analyze star')
    parser.add_argument('--focusdelay', default=0.5, type=float,  help='Delay (seconds) after focus moves')
    parser.add_argument('--numaverage', default=5, type=int,  help='Number of images to average')
    parser.add_argument('--simuldatadir', type=str,  help='Location of simulated star data')
    parser.add_argument('--vcurve_rs', type=float, help='VCurve Right Slope')
    parser.add_argument('--vcurve_ls', type=float, help='VCurve Left Slope')
    parser.add_argument('--vcurve_rp', type=float, help='VCurve Right PID')
    parser.add_argument('--vcurve_lp', type=float, help='VCurve Left PID')
    parser.add_argument('--start_hfd', type=float, help='Starting (outer) HFD value)')
    parser.add_argument('--near_hfd', type=float, help='Nearer (inner) HFD value)')
    parser.add_argument('--backlash', type=float, default=0, help='Steps of backlash')
    parser.add_argument('--forcehw', action='store_true', help='Force connecting to hw in simul mode')
    parser.add_argument('--keepfiles', action='store_true', help='Keep images taken during focusing')
    parser.add_argument('--measurefile', type=str,  help='Measure file and exit')
    parser.add_argument('--testfinal', action='store_true', help='Test final HFD determination')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # FIXME assumes tz is set properly in system?
    now = datetime.now()
    logfilename = 'autofocus_hfd_script-' + now.strftime('%Y%m%d%H%M%S') + '.log'

#    FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
#    FORMAT = '[%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'
    #FORMAT = '%(asctime)s [%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'

    LONG_FORMAT = '%(asctime)s.%(msecs)03d [%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'
    SHORT_FORMAT = '%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s'

    logging.basicConfig(filename=logfilename,
                        filemode='a',
                        level=logging.DEBUG,
                        format=LONG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    CONSOLE_FORMAT = '[%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'

    log = logging.getLogger()
    formatter = logging.Formatter(CONSOLE_FORMAT)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    args = parse_commandline()
    logging.debug(f'args = {args}')

    if args.debug:
        ch.setLevel(logging.DEBUG)

    # create plots if needed
    if args.debugplots:
        logging.debug('Creating figure')
#        fig2 = plt.figure()
#        ax_hfd = fig2.add_subplot(111)
#        hfd_plot, = ax_hfd.plot([],[], marker='o', ls='')
        fig = plt.figure(figsize=(4.5,2))
        ax_1d = fig.add_subplot(121)
        ax_2d = fig.add_subplot(122)
        plt.pause(0.01)

    if args.measurefile is not None:
        rc = measure_file_only(args.measurefile)
        #rc = test_exposure_scaling(args.measurefile)
        # keep plots up until keypress
        if args.debugplots:
            if args.stayopen:
                plt.show()
            else:
                plt.pause(args.debugplotsdelay)
        if rc:
            sys.exit(0)
        else:
            sys.exit(1)

    # initialize some variables
    backlash = 0

    start_time = time.time()

    # connect focuser and camera
    if not args.simul or args.forcehw:
        # load profile
        if args.profile is not None:
            logging.info(f'Using astro profile {args.profile}')
            ap = AstroProfile()
            ap.read(args.profile)
            logging.debug(f'profile = {ap.equipment}')
            backend_name = ap.equipment.backend.get('name')
            camera_driver = ap.equipment.camera.driver
            #print(dir(ap.equipment.focuser))
            logging.debug(f'ap.equipment.backend = {ap.equipment.backend}')
            logging.debug(f'ap.equipment.camera = {ap.equipment.camera}')
            logging.debug(f'ap.equipment.focuser = {ap.equipment.focuser}')
            logging.debug(f'ap.settings.autofocus = {ap.settings.autofocus}')
            focuser_driver = ap.equipment.focuser.get('driver')
            FOCUSER_MIN_POS = ap.equipment.focuser.get('minpos', None)
            FOCUSER_MAX_POS = ap.equipment.focuser.get('maxpos', None)
            FOCUSER_DIR = ap.settings.autofocus.get('focus_dir', None)
            # vcurve_rs = 0.04472660490995652
            # vcurve_rp = 5.333072819832182
            # vcurve_ls = -0.045806498798410436
            # vcurve_lp = -5.22113594480723
            vcurve_rs = ap.settings.autofocus.get('vcurve_rs', None)
            vcurve_rp = ap.settings.autofocus.get('vcurve_rp', None)
            vcurve_ls = ap.settings.autofocus.get('vcurve_ls', None)
            vcurve_lp = ap.settings.autofocus.get('vcurve_lp', None)
            logging.debug(f'vcurve_rs,rp,ls,lp = {vcurve_rs} {vcurve_rp} {vcurve_ls} {vcurve_lp}')
            start_hfd = ap.settings.autofocus.get('start_hfd', None)
            near_hfd = ap.settings.autofocus.get('near_hfd', None)
            max_hfd = ap.settings.autofocus.get('max_hfd', None)
            backlash = ap.settings.autofocus.get('backlash', 0)
        else:
            focuser_driver = args.focuser
            camera_driver = args.camera
            FOCUSER_MIN_POS = args.focus_min
            FOCUSER_MAX_POS = args.focus_max
            FOCUSER_DIR = args.focus_dir
            vcurve_rs = args.vcurve_rs
            vcurve_rp = args.vcurve_rp
            vcurve_ls = args.vcurve_ls
            vcurve_lp = args.vcurve_lp
            start_hfd = args.start_hfd
            near_hfd = args.near_hfd
            backlash = args.backlash
            backend_name = None

        sdi = SDI()

        sdi.connect_backend(backend_name)

        #focuser = connect_focuser(ASCOM_FOCUS_DRIVER)
        logging.debug(f'Connecting to focuser driver {focuser_driver}')
        focuser = sdi.connect_focuser(focuser_driver)
        logging.debug(f'focuser = {focuser}')
        if not focuser:
            logging.error(f'Unabled to connect to focuser driver {focuser_driver}')
            sys.exit(-1)
        logging.debug(f'Connecting to camera driver {camera_driver}')
        cam = sdi.connect_camera(camera_driver)
        logging.debug(f'cam = {cam}')
        if not cam:
            logging.error(f'Unabled to connect to camera driver {camera_driver}')
            sys.exit(-1)

    if args.simul:
        if args.simuldatadir is not None:
            starfile = os.path.join(args.simuldatadir, 'C8_Simul_Defocus_Star.fit')
            bgfile = os.path.join(args.simuldatadir, 'C8_Simul_BG.fit')
            simul_star = C8_F7_Star_Simulator(starimage_name=starfile,
                                              bgimage_name=bgfile)
        else:
            simul_star = C8_F7_Star_Simulator()

        # from 2019/05/20 C8 @ f/7 - seemed to work well with L filter (only tested)
        vcurve_rs = 0.04472660490995652
        vcurve_rp = 5.333072819832182
        vcurve_ls = -0.045806498798410436
        vcurve_lp = -5.22113594480723

        if args.focus_dir is not None:
            FOCUSER_DIR = args.focus_dir
        else:
            logging.info('Since --simul mode active focus dir of IN since none specified.')
            FOCUSER_DIR = 'IN'

        if args.focus_min is not None:
            FOCUSER_MIN_POS = args.focus_min
        else:
            logging.info('Since --simul mode active FOCUSER_MIN_POS set to 5000 since none specified.')
            FOCUSER_MIN_POS = 5000

        if args.focus_max is not None:
            FOCUSER_MAX_POS = args.focus_max
        else:
            logging.info('Since --simul mode active FOCUSER_MAX_POS set to 12000 since none specified.')
            FOCUSER_MAX_POS = 12000

        if args.start_hfd is not None:
            start_hfd = args.start_hfd
        else:
            logging.info('Since --simul mode active start_hfd set to 25 since none specified.')
            start_hfd = 25

        if args.near_hfd is not None:
            near_hfd = args.near_hfd
        else:
            logging.info('Since --simul mode active near_hfd set to 12 since none specified.')
            near_hfd = 12

    # some sanity checks
    if None in [near_hfd, start_hfd]:
        logging.error('Must specify both start and near HFD values.')
        sys.exit(1)

    if None in [vcurve_rs, vcurve_rp, vcurve_ls, vcurve_lp]:
        logging.error('Must specify VCurve Parameters!')
        sys.exit(1)

    # since we're on real hardware need a min/max pos
    if FOCUSER_MIN_POS is None or FOCUSER_MAX_POS is None:
        logging.error('Must specify a min and max allowed pos for focuser!')
        sys.exit(1)

    if FOCUSER_DIR is None:
        logging.error('Must specify the preferred focus direction!')
        sys.exit(1)


    logging.debug(f'Using focus min/max/dir of {FOCUSER_MIN_POS} {FOCUSER_MAX_POS} {FOCUSER_DIR}')

    if FOCUSER_MAX_POS <= FOCUSER_MIN_POS:
        logging.error('Focus max position must be less than focus min position!')
        sys.exit(1)

    # create output dir
    datestr = time.strftime("%Y%m%d_%H%M%S", time.localtime())

    IMAGESDIROBJ = tempfile.TemporaryDirectory(prefix='autofocus_hfd_script_'+datestr+'_')
    IMAGESDIR = IMAGESDIROBJ.name

    logging.info(f'Using temporary directory {IMAGESDIR}')

    #IMAGESDIR = datestr
    #os.mkdir(IMAGESDIR)

    # save initial focus position
    focus_expos = args.exposure_start
    logging.info(f'Starting exposure is {focus_expos} seconds')

    # now get from profile or command line only!
    # here for reference only
    # from 2019/05/20 C8 @ f/7 - seemed to work well with L filter (only tested)
    # vcurve_rs = 0.04472660490995652
    # vcurve_rp = 5.333072819832182
    # vcurve_ls = -0.045806498798410436
    # vcurve_lp = -5.22113594480723

    if args.focus_start is None:
        if not args.simul:
            fpos_1 = focuser.get_absolute_position()
        else:
            fpos_1 = 8000
    else:
        fpos_1 = args.focus_start

    # test code just measure HFD and adjust exposure
    # so it is not saturated
    if args.testfinal:
        final_hfd = determine_final_hfd(fpos_1, focus_expos)
        logging.info(f'Final HFD = {final_hfd}')
        sys.exit(0)

    # save starting position for if we need to restore
    starting_focus_pos = fpos_1

    if FOCUSER_DIR == 'OUT':
        fdir = 1
        vslope = vcurve_ls
        vpid = vcurve_lp
    elif FOCUSER_DIR == 'IN':
        fdir = -1
        vslope = vcurve_rs
        vpid = vcurve_rp
    else:
        logging.error(f'Unknown focus direction {FOCUSER_DIR} - exitting!')
        cleanup_files()
        sys.exit(1)

    ntries = 0
    right_side = False
    while ntries < 3:
        rc_fpos1 = measure_at_focus_pos(fpos_1, focus_expos)
        if rc_fpos1 is None:
            hfd_1 = None
        else:
            hfd_1, satur, freturn = rc_fpos1

        if hfd_1 is None:
            logging.error('No star found!')
            #if args.debugplots:
            #    plt.show()
            cleanup_files()
            logging.info('Exitting!')
            sys.exit(1)

        if args.debugplots:
            fig.suptitle(f'First focus {fpos_1} HFD {hfd_1:5.2f}')
#            fig.show()
#            plt.pause(args.debugplotsdelay)
            show_fig_and_wait(fig, args.debugplotsdelay)

        logging.info(f'INITIAL FOCUS = {fpos_1}  HFD = {hfd_1}')

        # move out 10 HFD
        nsteps = int(abs(10/vslope))
        fpos_2 = fpos_1 - fdir*nsteps

        rc_fpos2 = measure_at_focus_pos(fpos_2, focus_expos)
        if rc_fpos2 is None:
            hfd_2 = None
        else:
            hfd_2, satur, freturn = rc_fpos2

        if hfd_2 is None:
            logging.error('No star found!')
            #if args.debugplots:
            #    plt.show()
            if not args.simul or args.forcehw:
                restore_focus_pos(starting_focus_pos)
            cleanup_files()
            logging.info('Exitting!')
            sys.exit(1)

        logging.info(f'SECOND FOCUS = {fpos_2}  HFD = {hfd_2}')

        if args.debugplots:
            fig.suptitle(f'Second focus {fpos_2} HFD {hfd_2:5.2f}')
#            fig.show()
#            plt.pause(args.debugplotsdelay)
            show_fig_and_wait(fig, args.debugplotsdelay)

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
        cleanup_files()
        if not args.simul or args.forcehw:
            restore_focus_pos(starting_focus_pos)
        sys.exit(1)

    # compute location for desired Start HFD
    #initial_hfd = 24
    fpos_start = fpos_2 - fdir*int(abs((start_hfd-hfd_2)/vslope))
    logging.info(f'Start HFD = {start_hfd} pred focus = {fpos_start}')

    # figure out direction
    #backlash = 200
    fpos_pre_start = fpos_start - fdir*backlash

    logging.info(f'Moving to pre-start pos {fpos_pre_start}')
    if not args.simul or args.forcehw:
        move_focuser(fpos_pre_start)
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
    avg_start_hfd = average_measure_at_focus_pos(fpos_start, focus_expos, args.numaverage, tag='start')

    logging.info(f'START POSITION FOCUS AVERAGE HFD = {avg_start_hfd}')

    # now based on average at start compute near HFD location
    # compute location for desired start HFD
    #inner_hfd = 12

    if avg_start_hfd is not None:
        fpos_near = fpos_start + fdir*int(abs((near_hfd-avg_start_hfd)/vslope))
        logging.info(f'Near HFD = {near_hfd} pred focus = {fpos_near}')
    else:
        fpos_near = None

    if fpos_near is None:
        logging.error('Near HFD - No star found!')
        #if args.debugplots:
        #    plt.show()
        if not args.simul or args.forcehw:
            restore_focus_pos(starting_focus_pos)
        cleanup_files()
        logging.info('Exitting!')
        sys.exit(1)

    rc_near = measure_at_focus_pos(fpos_near, focus_expos)
    if rc_near is None:
        hfd_near = None
    else:
        hfd_near, satur, freturn = rc_near

    if hfd_near is None:
        logging.error('No star found!')
        #if args.debugplots:
        #    plt.show()
        if not args.simul or args.forcehw:
            restore_focus_pos(starting_focus_pos)
        cleanup_files()
        logging.info('Exitting!')
        sys.exit(1)

    logging.info(f'NEAR POSITION FOCUS = {fpos_near}  HFD = {hfd_near}')

    if args.debugplots:
        fig.suptitle(f'Near position focus {fpos_near} HFD {hfd_near:5.2f}')
        #fig.show()
        #plt.pause(args.debugplotsdelay)
        show_fig_and_wait(fig, args.debugplotsdelay)

    # now take several measurements and get average HFD
    avg_near_hfd = average_measure_at_focus_pos(fpos_near, focus_expos, args.numaverage, tag='near')

    logging.info(f'NEAR POSITION FOCUS AVERAGE HFD = {avg_near_hfd}')

    # now compute best focus
    if avg_near_hfd is not None:
        fpos_best = int(fpos_near + fdir*int(abs(avg_near_hfd/vslope)) + vpid)
    else:
        fpos_best = None

    logging.info(f'BEST FOCUS POSITION = {fpos_best} {fpos_near} {int(avg_near_hfd/vslope)} {vpid}')

    # when using the internal 'focus simulator' we don't impose a min/max position limit
    if (not args.simul or args.forcehw) and (fpos_best > FOCUSER_MAX_POS or fpos_best < FOCUSER_MIN_POS):
        logging.error(f'Best focus position {fpos_best} is outside allowed range {FOCUSER_MIN_POS} to {FOCUSER_MAX_POS}')
        restore_focus_pos(starting_focus_pos)
    else:
        tries = 0
        best_hfd = None
        # old method tried several times - not doing this for now
        # while tries < 4:
            # final_hfd, satur = determine_final_hfd(fpos_best, focus_expos)
            # logging.info(f'Final HFD = {final_hfd}')
            # if final_hfd is not None:
                # best_hfd = final_hfd
                # break
            # tries += 1

        # just try once for now
        # try to get autoexposed non-saturated star profile and measure
        best_hfd, satur = determine_final_hfd(fpos_best, focus_expos)
        logging.info(f'Best HFD = {best_hfd}')

        if best_hfd is None:
            logging.error('Could not determine final HFD!')
            #if not args.simul or args.forcehw:
            #    restore_focus_pos(starting_focus_pos)
            #cleanup_files()
            #logging.info('Exitting!')
            #sys.exit(1)
        else:
            logging.info(f'BEST FOCUS POSITION = {fpos_best} HFD = {best_hfd}')

            if args.debugplots:
                if satur:
                    s_str = 'SATUR'
                else:
                    s_str = 'UNSAT'
                fig.suptitle(f'Best pos focus {fpos_best} HFD {best_hfd:5.2f} {s_str}')

    # keep plots up until keypress
    if args.debugplots:
        if args.stayopen:
            plt.show()
        else:
            plt.pause(5)

    cleanup_files()

    logging.info(f'Focus run took {time.time() - start_time} seconds.')

    if args.debugplots:
        plt.close('all')

    logging.info('Returning with rc of 0')
    #sys.exit(0)
    os._exit(0)

