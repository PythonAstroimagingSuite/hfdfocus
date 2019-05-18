import os
import sys
import time
import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt

ASCOM_FOCUS_DRIVER = 'ASCOM.Simulator.Focuser'
#ASCOM_FOCUS_DRIVER = 'ASCOM.MoonliteDRO.Focuser'

def get_backend_for_os():
    import os
    # chose an implementation, depending on os
    if os.name == 'nt': #sys.platform == 'win32':
        return 'ASCOM'
    elif os.name == 'posix':
        return 'INDI'
    else:
        raise Exception("Sorry: no implementation for your platform ('%s') available" % os.name)

BACKEND = get_backend_for_os()

# debugging override with simulator
BACKEND = 'SIMULATOR'

import astropy.io.fits as pyfits
from c8_simul_star import C8_F7_Star_Simulator

class FocusSimlatorDriver:
    def __init__(self):
        self.destpos = 8000

    def move_absolute_position(self, pos):
        self.destpos = pos
        return True

    def get_absolute_position(self):
        return self.destpos

    def is_moving(self):
        return False

    def __getattr__(self, name):
        def method(*args):
            print("tried to handle unknown method " + name)
            if args:
                print("it had arguments: " + str(args))
            return True

        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return method

# camera driver needs to access simuluated focus driver
# in order to create star image scaled
#
# NOTE NOTE NOTE
# focus center of 7983 is based on the
class CameraSimlatorDriver:
    def __init__(self, focus_driver):
        self.focus_driver = focus_driver
        self.simul = C8_F7_Star_Simulator()
        self.focus_center = self.simul.get_best_focus_pos()
        self.image_shape = self.simul.get_simul_star_image(0, 0).shape

    def get_size(self):
        #print(self.image_shape)
        return self.image_shape

    def save_image_data(self, fname):
        data = self.simul.get_simul_star_image(self.focus_driver.get_absolute_position(), self.focus_center)
        pyfits.writeto(fname, data.astype(float), overwrite=False)
        return True

    def __getattr__(self, name):
        def method(*args):
            print("tried to handle unknown method " + name)
            if args:
                print("it had arguments: " + str(args))
            return True

        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return method

# camera driver needs to access simuluated focus driver
# in order to create star image scaled
class BackendSimlatorDriver:
    def __init__(self):
        pass

    def __getattr__(self, name):
        def method(*args):
            print("tried to handle unknown method " + name)
            if args:
                print("it had arguments: " + str(args))
            return True
        return method

if BACKEND == 'ASCOM':
    from pyastrobackend.ASCOMBackend import DeviceBackend as Backend
    from pyastrobackend.ASCOM.Focuser import Focuser
elif BACKEND == 'INDI':
    from pyastrobackend.INDIBackend import DeviceBackend as Backend
    from pyastrobackend.INDIBackend import Focuser
elif BACKEND == 'SIMULATOR':
    Backend = BackendSimlatorDriver()
    SimFocuser = FocusSimlatorDriver()
else:
    raise Exception(f'Unknown backend {BACKEND} - choose ASCOM or INDI in BackendConfig.py')


if BACKEND == 'ASCOM':
    from pyastrobackend.MaximDL.Camera import Camera as MaximDL_Camera
elif BACKEND == 'INDI':
    from pyastrobackend.INDIBackend import Camera as INDI_Camera
elif BACKEND == 'SIMULATOR':
    SimCamera = CameraSimlatorDriver(SimFocuser)
else:
    raise Exception(f'Unknown backend {BACKEND} - choose ASCOM or INDI in BackendConfig.py')



def connect_focuser(backend):
    global SimFocuser

    if BACKEND == 'ASCOM':
        focuser = Focuser()
    elif BACKEND == 'INDI':
        focuser = Focuser(backend)
    elif BACKEND == 'SIMULATOR':
        focuser = SimFocuser

    rc = focuser.connect(ASCOM_FOCUS_DRIVER)

    if rc:
        return focuser
    else:
        return None

def wait_on_focuser_move(focuser, timeout=60):
    ts = time.time()
    while (time.time()-ts) < timeout:
        logging.info(f'waiting on focuser move - curpos = {focuser.get_absolute_position()}')
        if not focuser.is_moving():
            break
        time.sleep(0.5)
    time.sleep(0.5) # just be sure its done


# FIXME INDI stuff is broken!!!!
def connect_camera(backend):
    global SimCamera

    if BACKEND == 'ASCOM':
        driver = 'MaximDL'
        cam = MaximDL_Camera()
    elif BACKEND == 'INDI':
        driver = 'INDICamera'
        cam = INDI_Camera(backend)
    elif BACKEND == 'SIMULATOR':
        driver = 'Simulator'
        cam = SimCamera

    logging.info(f'connect_camera: driver = {driver}')

    rc = None
    if driver == 'INDICamera':
        if ':' in self.settings.camera_driver:
            indi_cam_driver = self.settings.camera_driver.split(':')[1]
            rc = self.cam.connect(indi_cam_driver)
        else:
            logging.error('connect_camera(): Must configure INDI camera driver first!')
            return None
    else:
        rc = cam.connect(driver)

    if not rc:
        logging.error('connect_camera(): Unable to connect to camera!')
        return None
    return cam

# take exposure
def take_exposure(focus_expos, output_filename):

    focus_expos = 1

    # reset frame to full sensor
    cam.set_binning(1, 1)
    width, height = cam.get_size()
    cam.set_frame(0, 0, width, height)
    cam.start_exposure(focus_expos)

    # give things time to happen (?) I get Maxim not ready errors so slowing it down
    time.sleep(0.25)

    elapsed = 0
    while not cam.check_exposure():
        logging.info(f"Taking image with camera {elapsed} of {focus_expos} seconds")
        time.sleep(0.5)
        elapsed += 0.5
        if elapsed > focus_expos:
            elapsed = focus_expos

    logging.info('Exposure complete')
    # give it some time seems like Maxim isnt ready if we hit it too fast
    time.sleep(0.5)

    ff = os.path.join(os.getcwd(), output_filename)

    retries = 0
    while True:
        logging.info(f"Going to save {ff}")

        # FIXME we only call this because the
        # MaximDL backend needs it to save to disk
        # RPC backend already has saved it to disk by now
        if BACKEND == 'INDI':
            # FIXME need better way to handle saving image to file!
            image_data = cam.get_image_data()
            # this is an hdulist
            image_data.writeto(ff, overwrite=True)
            result = True
        else:
            result = cam.save_image_data(ff)

        if result is True:
            logging.info("Save successful")
            break

        retries += 1
        if retries > 3:
            logging.error(f"Failed to save {ff}!! Aborting!")
            return False

        logging.error(f"Failed to save {ff}!! Trying again")
        time.sleep(1)

    return True

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('focus_low', type=int, help='Low end of focus run')
    parser.add_argument('focus_high', type=int, help='High end of focus run')
    parser.add_argument('focus_nstep', type=int, help='V Curve number of steps')
    parser.add_argument('focus_dir', type=str, help='IN or OUT')

    #    parser.add_argument('--debuggraphs', action='store_true', help="Display debug graphs")

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

    # connect focuser
    focuser = connect_focuser(Backend)
    logging.info(f'focuser = {focuser}')

    # connect camera
    cam = connect_camera(Backend)
    logging.info(f'cam = {cam}')
    logging.info(f'cam size = {cam.get_size()}')

    # create output dir
    datestr = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    imagesdir = datestr
    os.mkdir(imagesdir)

    # figure out direction
    backlash = 200
    if args.focus_dir == 'OUT':
        # start past desired start and move to it to remove backlash
        focus_init = args.focus_low - backlash
        focus_start = args.focus_low
        focus_end = args.focus_high
        focus_nstep = args.focus_nstep
    elif args.focus_dir == 'IN':
        # start past desired start and move to it to remove backlash
        focus_init = args.focus_high + backlash
        focus_start = args.focus_high
        focus_end = args.focus_low
        focus_nstep = args.focus_nstep
    else:
        logging.error(f'Unknown focus directin {args.focus_dir} - exitting!')
        sys.exit(1)

    # move to init pos
    logging.info(f'Moving to init pos {focus_init}')
    if not focuser.move_absolute_position(focus_init):
        logging.error("Focuser error!!")
        sys.exit(1)

    wait_on_focuser_move(focuser)

    focus_expos = 1

    focus_step = int((focus_end - focus_start)/(focus_nstep - 1))
    logging.info(f'Focus run from {focus_start} to {focus_end} step {focus_step}')
    for focus_pos in range(focus_start, focus_end+focus_step, focus_step):
        logging.info(f'Moving to focus pos {focus_pos}')
        if not focuser.move_absolute_position(focus_pos):
            logging.error("Focuser error!!")
            sys.exit(1)

        wait_on_focuser_move(focuser)

        logging.info('Taking exposure')
        rc = take_exposure(focus_expos, os.path.join(imagesdir, f'vcurve_focuspos_{focus_pos}.fit'))
        logging.info(f'exposure result code = {rc}')

