import sys
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




def connect_camera(self):
    if BACKEND == 'ASCOM':
        driver = 'MaximDL'
        self.cam = MaximDL_Camera()
    elif BACKEND == 'INDI':
        driver = 'INDICamera'
        self.cam = INDI_Camera(self.backend)

    logging.info(f'connect_camera: driver = {driver}')

    rc = False
    if driver == 'INDICamera':
        if ':' in self.settings.camera_driver:
            indi_cam_driver = self.settings.camera_driver.split(':')[1]
            rc = self.cam.connect(indi_cam_driver)
        else:
            logging.error('connect_camera(): Must configure INDI camera driver first!')
            return False
    else:
        rc = self.cam.connect(driver)



    if not rc:
        logging.error('connect_camera(): Unable to connect to camera!')
        return False
    return True


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