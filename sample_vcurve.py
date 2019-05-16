import sys
import time
import logging
import numpy as np
import matplotlib.pyplot as plt

print(sys.path)

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

if BACKEND == 'ASCOM':
    from pyastrobackend.MaximDL.Camera import Camera as MaximDL_Camera
elif BACKEND == 'INDI':
    from pyastrobackend.INDIBackend import Camera as INDI_Camera
else:
    raise Exception(f'Unknown backend {BACKEND} - choose ASCOM or INDI in BackendConfig.py')




def connect_camera():
    if BACKEND == 'ASCOM':
        driver = 'MaximDL'
        cam = MaximDL_Camera()
    elif BACKEND == 'INDI':
        driver = 'INDICamera'
        cam = INDI_Camera(self.backend)

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

    # connect camera
    cam = connect_camera()
    logging.info(f'cam = {cam}')

    # take exposure



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
        logging(f"Taking image with camera {elapsed} of {focus_expos} seconds")
        time.sleep(0.5)
        elapsed += 0.5
        if elapsed > focus_expos:
            elapsed = focus_expos

    logging.info('Exposure complete')
    # give it some time seems like Maxim isnt ready if we hit it too fast
    time.sleep(0.5)
