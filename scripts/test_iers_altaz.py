# 2019/12/06
#
# Due to some naval obs sites being down this has big delay at start as it
# times out trying to update IERS data
import logging

# this uses a different site to download it
from astropy.utils import iers
from astroplan import download_IERS_A
#iers.IERS_A_URL = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'

print('Downloading IERS')
download_IERS_A()

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation
from astroplan import Observer

def get_altaz_from_radec(radec, observer, obstime):
    altaz = observer.altaz(obstime, target=radec)
#    logging.debug(f'get_mount_altaz: computerd alt = {altaz.alt.degree}')
#    logging.debug(f'get_mount_altaz: computerd az  = {altaz.az.degree}')
    return altaz.alt.degree, altaz.az.degree

if __name__ == '__main__':

    LONG_FORMAT = '%(asctime)s.%(msecs)03d [%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'
    SHORT_FORMAT = '%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s'

    logging.basicConfig(filename='test_iers_altaz.log',
                        filemode='a',
                        level=logging.DEBUG,
                        format=LONG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter(LONG_FORMAT)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(ch)



cur_radec = SkyCoord("1:12:43.2 +1:12:43", unit=(u.hourangle, u.deg))
print(cur_radec)

location = EarthLocation.from_geodetic(-272*u.deg, 36*u.deg, 140*u.m)

observer = Observer(location=location, name='TestSite', timezone='US/Eastern')
alt, az = get_altaz_from_radec(cur_radec, observer, Time.now())
print(f'Current mount alt/az is {alt:3.1f} {az:4.1f}')