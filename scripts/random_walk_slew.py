import sys
import time
from datetime import datetime
import logging

import numpy as np

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.coordinates import Angle

from pyastrobackend.BackendConfig import get_backend_for_os

BACKEND = get_backend_for_os()

from pyastrobackend.BackendConfig import get_backend_for_os

BACKEND = get_backend_for_os()

if BACKEND == 'ASCOM':
    from pyastrobackend.ASCOMBackend import DeviceBackend as Backend
elif BACKEND == 'INDI':
    from pyastrobackend.INDIBackend import DeviceBackend as Backend
else:
    raise Exception(f'Unknown backend {BACKEND} - choose ASCOM or INDI in BackendConfig.py')

if BACKEND == 'ASCOM':
    from pyastrobackend.ASCOM.Mount import Mount as MountClass
elif BACKEND == 'INDI':
    from pyastrobackend.INDIBackend import Mount as MountClass
else:
    raise Exception(f'Unknown backend {BACKEND} - choose ASCOM or INDI in BackendConfig.py')


class Telescope(MountClass):
    def __init__(self, *args):
        super().__init__(*args)

#        self.tel = None
        self.connected = False

        # FIXME Currently INDI uses backend when creating devices and ASCOM
        #       doesnt!
        #       Assume if init has argument it is backend for INDI
        if len(args) > 0:
            self.backend = args[0]

    @staticmethod
    def precess_J2000_to_JNOW(pos_J2000):
        """Precess J2000 coordinates to JNOW

        Parameters
        ----------
        pos_J2000 - SkyCoord
            J2000 sky coordinate to precess

        Returns
        -------
        pos_JNOW : SkyCoord
            JNow coordinate
        """
        time_now = Time(datetime.utcnow(), scale='utc')
        return pos_J2000.transform_to(FK5(equinox=Time(time_now.jd, format='jd', scale='utc')))

    @staticmethod
    def precess_JNOW_to_J2000(pos_JNOW):
        """Precess J2000 coordinates to JNOW

        Parameters
        ----------
        pos_JNOW - SkyCoord
            JNow sky coordinate to precess

        Returns
        -------
        pos_J2000 : SkyCoord
            J2000 coordinate
        """
        return pos_JNOW.transform_to(FK5(equinox='J2000'))

    def connect_to_telescope(self, driver):
        if self.connected:
            logging.warning('connect_to_telescope: already connected!')

        logging.info(f"Connect to telescope driver {driver}")

        if not super().connect(driver):
            return False

        self.connected = True
        return True

    def is_connected(self):
        return self.connected

    def get_position_jnow(self):
        if not self.connected:
            return None
        time_now = Time(datetime.utcnow(), scale='utc')
        ra_now, dec_now = super().get_position_radec()

        return SkyCoord(ra=ra_now*u.hour, dec=dec_now*u.degree, frame='fk5', equinox=Time(time_now.jd, format="jd", scale="utc"))

    def get_position_j2000(self):
        if not self.connected:
            return None
        pos_jnow = self.get_position_jnow()
        return self.precess_JNOW_to_J2000(pos_jnow)

    def sync(self, pos):
        if not self.connected:
            return False

        logging.info(f"Syncing to {pos.ra.hour}  {pos.dec.degree}")
        try:
           super().sync(pos.ra.hour, pos.dec.degree)
        except Exception as e:
            logging.error('sync() Exception ->', exc_info=True)
            return False

        return True

    def goto(self, pos):
        if not self.connected:
            return False
        logging.info(f"Goto to {pos.ra.hour}  {pos.dec.degree}")
        super().slew(pos.ra.hour, pos.dec.degree)
        return True

if __name__ == '__main__':
    logging.basicConfig(filename='random_walk_slew.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    LOG = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    CH = logging.StreamHandler()
    CH.setLevel(logging.DEBUG)
    CH.setFormatter(formatter)
    LOG.addHandler(CH)

    backend = Backend()

    rc = backend.connect()
    if not rc:
        logging.error('Failed to connect to backend!')
        sys.exit(-1)

    tel = Telescope(backend)

    rc = tel.connect_to_telescope("Telescope Simulator")
    if not rc:
        logging.error('Unable to connect to telescope!')

    curpos = tel.get_position_jnow()
    logging.info(f"Current position is {curpos.to_string('hmsdms', sep=':')}")

    # compute random offset
    walk_dist = Angle(5*np.random.rand()*u.degree)
    walk_angle = Angle(2*np.pi*np.random.rand()*u.rad)

    logging.info(f'walk_dist = {walk_dist} degrees at angle {walk_angle.degree}')
    newpos = curpos.directional_offset_by(walk_angle, walk_dist)
    logging.info(f"Current position is {newpos.to_string('hmsdms', sep=':')}")

    rc = tel.goto(newpos)

    if not rc:
        logging.error('slew failed')
        sys.exit(1)

    while True:
        isslew = tel.is_slewing()
        logging.info(f'is_slewing = {isslew}')
        if not isslew:
            break
        time.sleep(3)

    sys.exit(0)

