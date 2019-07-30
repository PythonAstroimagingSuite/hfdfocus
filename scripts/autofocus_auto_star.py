#
#
# currently relies on being run from hfdfocus/scripts directory
# and having pyastrometry_cli checked out in same tree
#

import os
import sys
import json
import time
import shlex
import argparse
import logging
import subprocess
import tempfile
from datetime import datetime

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pyastroprofile.AstroProfile import AstroProfile

# FIXME this should be handled automaticall!

#if os.name == 'nt':
#    PYTHON_EXE_PATH = 'C:\\Users\\msf\\Anaconda3\\envs\\AstronomyUtilitiesPipNumpy\\python'
#elif os.name == 'posix':
#    PYTHON_EXE_PATH = '/home/msf/anaconda3/envs/pyastro37/bin/python3 '
#else:
#    logging.error('PYTHON_EXE_PATH NOT SET')
#    sys.exit(1)

# find interpretter we're running under and use it?
PYTHON_EXE_PATH = sys.executable

def run_program(cmd_line, label='', trimlog=True):
    """ runs exec_path with cmd_line returning the return code rc and stdout
        and stderr output as output """

    # seems like we have to do things differently for Windows and Linux
    if os.name == 'nt':
        cmd_val = cmd_line
    elif os.name == 'posix':
        cmd_val = shlex.split(cmd_line)

    logging.debug(f'run_program() command value = |{cmd_val}|')

    with subprocess.Popen(cmd_val,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True) as ps_proc:

        if label == '':
            label = 'run_program()'
        logging.debug('ps_proc output:')
        output = ''
        for l in ps_proc.stdout:
            #logging.debug(l.strip())
            if trimlog:
                # get rid of first 3 'words' which are the
                # logging info from program
                words = l.strip().split(' ')
                out = ' '.join(words[3:])
            else:
                out = l.strip()
            logging.info(f'{label}: {out}')
            output += l
        logging.debug('end of output')

#    poll_value = None
#    while True:
#        poll_value = ps_proc.poll()
#
#        if poll_value is not None:
#            break

    # check return code
    rc = ps_proc.returncode

    logging.debug(f'run_program: return code was {rc}!')
    return rc, output


def run_platesolve():
    # parse out device info
    # FIXME seems alot of duplication need a better way to represent
    # device info like a config file

    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of astro profile')
    parser.add_argument('--pixelscale', type=float, help='Pixel scale for plate solving')
    dev_args, unknown = parser.parse_known_args(sys.argv)

    logging.debug(f'run_platesolve: dev_args, unknown = {dev_args} {unknown}')

    #result_fname = './autofocus_auto_origpos.json'

    tmp_fd, tmp_result_fname = tempfile.mkstemp(suffix='.json')
    os.close(tmp_fd)
    result_fname = tmp_result_fname
    logging.debug(f'Using solution json tmp file {result_fname}')

    cmd_line = PYTHON_EXE_PATH + ' '
    script = 'pyastrometry_cli_main.py'
    if PYASTROMETRY_SCRIPT_PATH is not None:
        script = os.path.join(PYASTROMETRY_SCRIPT_PATH, script)
    cmd_line += f'{script} solvepos '
    cmd_line += f'--outfile {result_fname} '
    if dev_args.pixelscale is not None:
        cmd_line += f'--pixelscale {dev_args.pixelscale} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    # unlink previous solve if any
    if os.path.isfile(result_fname):
        os.unlink(result_fname)

    rc, output = run_program(cmd_line, label='platesolve')

    # FIXME following might  leave temporary file around if there is
    #       an exception or it returns before reaching end of function
    if rc < 0:
        logging.error(f'run_platesolve: return code was {rc}!')
        return None

    try:
        result_file = open(result_fname, 'r')
    except OSError as err:
        logging.error(f"Error opening results file: {err}")
        return None

    try:
        solve_dict = json.load(result_file)
        target_str = solve_dict['ra2000'] + " "
        target_str += solve_dict['dec2000']
        logging.debug(f"target_str = {target_str}")

        radec = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                          frame='fk5', equinox='J2000')

    except Exception as err:
        logging.error(f"Error converting solve results! {err}")
        radec = None

    result_file.close()

    try:
        os.unlink(result_fname)
    except:
        logging.error(f'Error cleaning up tmp file {result_fname}', exc_info=True)

    if radec is not None:
        logging.debug(f'radec = {radec.to_string()}')

    return radec

def run_getpos():
    # parse out device info
    # FIXME seems alot of duplication need a better way to represent
    # device info like a config file

    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of astro profile')
    dev_args, unknown = parser.parse_known_args(sys.argv)

    logging.debug(f'run_getpos: dev_args, unknown = {dev_args} {unknown}')

    #result_fname = './autofocus_auto_origpos.json'

    tmp_fd, tmp_result_fname = tempfile.mkstemp(suffix='.json')
    os.close(tmp_fd)
    result_fname = tmp_result_fname
    logging.debug(f'Using solution json tmp file {result_fname}')

    cmd_line = PYTHON_EXE_PATH + ' '
    script = 'pyastrometry_cli_main.py'
    if PYASTROMETRY_SCRIPT_PATH is not None:
        script = os.path.join(PYASTROMETRY_SCRIPT_PATH, script)
    cmd_line += f'{script} getpos '
    cmd_line += f'--outfile {result_fname} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    # unlink previous solve if any
    if os.path.isfile(result_fname):
        os.unlink(result_fname)

    rc, output = run_program(cmd_line, label='platesolve')

    if rc < 0:
        logging.error(f'run_platesolve: return code was {rc}!')
        return None

    try:
        result_file = open(result_fname, 'r')
    except OSError as err:
        logging.error(f"Error opening results file: {err}")
        return None

    try:
        solve_dict = json.load(result_file)
        target_str = solve_dict['ra2000'] + " "
        target_str += solve_dict['dec2000']
        logging.debug(f"target_str = {target_str}")

        radec = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                          frame='fk5', equinox='J2000')

    except Exception as err:
        logging.error(f"Error converting solve results! {err}")
        radec = None

    if radec is not None:
        logging.debug(f'radec = {radec.to_string()}')

    result_file.close()

    try:
        os.unlink(result_fname)
    except:
        logging.error(f'Error cleaning up tmp file {result_fname}', exc_info=True)

    return radec


def run_findstars(curpos, args, lon=None):
    logging.info(f'Finding nearby stars within {args.dist} deg around mag {float(args.mag)}')

    #result_fname = './autofocus_auto_starlist.dat'

    tmp_fd, tmp_result_fname = tempfile.mkstemp(prefix='autofocus_auto_starlist_', suffix='.dat')
    os.close(tmp_fd)
    result_fname = tmp_result_fname
    logging.debug(f'Using solution json tmp file {result_fname}')

    import time
    time.sleep(10)

    cmd_line = PYTHON_EXE_PATH + ' '
    if AUTOFOCUS_SCRIPT_PATH is not None:
        cmd_line += f'{AUTOFOCUS_SCRIPT_PATH}/'
    cmd_line += 'find_nearby_stars.py '
    if AUTOFOCUS_DATA_PATH is not None:
        cmd_line += f'{AUTOFOCUS_DATA_PATH}/'
    #cmd_line += '../data/SAO_Catalog_m5_p11_filtered.bin '
    cmd_line += 'SAO_Catalog_m5_p11_filtered.bin '
    rastr = curpos.ra.to_string(u.hour, sep=":", pad=True)
    cmd_line += rastr + ' '
    decstr = curpos.dec.to_string(alwayssign=True, sep=":", pad=True)
    cmd_line += f'" {decstr}" '
    cmd_line += f'{args.dist} '
    cmd_line += f'--outfile {result_fname} '
    mag = float(args.mag)
    cmd_line += f'--minmag {mag+0.5} '
    cmd_line += f'--maxmag {mag-0.5} '
    if args.lst is not None:
        cmd_line += f'--lst {args.lst} '
    if args.lon is not None:
        cmd_line += f'--lon {args.lon} '
    elif lon is not None:
        cmd_line += f'--lon {lon} '
    if args.onlyside is not None:
        cmd_line += f'--onlyside {args.onlyside} '
    if args.meridianthres is not None:
        cmd_line += f' --meridianthres {args.meridianthres}'

    logging.debug(cmd_line)

    # unlink previous solve if any
    if os.path.isfile(result_fname):
        os.unlink(result_fname)

    rc, output = run_program(cmd_line, label='findstars')

    if rc != 0:
        logging.error(f'run_findstars: return code was {rc}!')
        return None

    try:
        result_file = open(result_fname, 'r')
    except OSError as err:
        logging.error(f"Error opening results file: {err}")
        return None

    star_list = []
    try:
        nlines = 0
        for l in result_file.readlines():
            #print(l.strip())
            if nlines == 0:
                xxx, field = l.strip().split('=')
                nstars = int(field)
                logging.info(f'# candidate stars = {nstars}')
                nlines += 1
                continue
            elif nlines == 1:
                # skip over headers
                nlines += 1
                continue

            catidx, distdeg, sao, rastr, decstr, vmag, nneigh = l.strip().split(',')
            target_str = rastr + " " + decstr
            logging.debug(f"star target_str = {target_str}")
            radec = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                              frame='fk5', equinox='J2000')
            star_list.append((sao, radec))

    except Exception as err:
        logging.error(f"Error converting solve results! {err}")
        star_list = None

    #logging.debug(f'star_list = {star_list}')

    result_file.close()

    try:
        os.unlink(result_fname)
    except:
        logging.error(f'Error cleaning up tmp file {result_fname}', exc_info=True)

    return star_list

def run_precise_slew(target, args, extra_args):
    # parse out device info
    # FIXME seems alot of duplication need a better way to represent
    # device info like a config file

    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of astro profile')
    parser.add_argument('--mount', type=str, help='Name of mount driver')
    parser.add_argument('--camera', type=str, help='Name of camera driver')
    parser.add_argument('--exposure', type=float, default=5, help='Exposure time')
    parser.add_argument('--binning', type=int, default=2, help='Camera binning')
    parser.add_argument('--pixelscale', type=float, help='Pixel scale (arcsec/pixel)')
#    parser.add_argument('--framesize', default=0, type=int,  help='Size of capture frame, 0=full')
    dev_args, unknown = parser.parse_known_args(sys.argv)

    logging.debug(f'dev_args, unknown = {dev_args} {unknown}')

    cmd_line = PYTHON_EXE_PATH + ' '
    script = 'pyastrometry_cli_main.py'
    if PYASTROMETRY_SCRIPT_PATH is not None:
        script = os.path.join(PYASTROMETRY_SCRIPT_PATH, script)
    cmd_line += f'{script} slewsolve '
    rastr = target.ra.to_string(u.hour, sep=":", pad=True)
    cmd_line += rastr + ' '
    decstr = target.dec.to_string(alwayssign=True, sep=":", pad=True)
    cmd_line += f'" {decstr}" '
    if dev_args.pixelscale:
        cmd_line += f'--pixelscale {dev_args.pixelscale} '
#    if dev_args.framesize is not None:
#        cmd_line += f'--framesize {dev_args.framesize} '

    # use json to handle double quotes in camera and mount
    # arguments which might be passed as something like "CCD Simulator"
    # and the shlex.split() will split them
    if dev_args.camera is not None:
        cmd_line += f'--camera {json.dumps(dev_args.camera)} '
    if dev_args.mount is not None:
        cmd_line += f'--mount {json.dumps(dev_args.mount)} '
    cmd_line += f'--exposure {dev_args.exposure} '
    cmd_line += f'--binning {dev_args.binning} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    rc, output = run_program(cmd_line, label='precise_slew')

    if rc != 0:
        logging.error(f'run_precise_slew: return code was {rc}!')
        return None

    return True

def run_slew(target, args, extra_args):
    # parse out device info
    # FIXME seems alot of duplication need a better way to represent
    # device info like a config file

    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of astro profile')
    parser.add_argument('--mount', type=str, help='Name of mount driver')

    dev_args, unknown = parser.parse_known_args(sys.argv)

    logging.debug(f'dev_args, unknown = {dev_args} {unknown}')

    cmd_line = PYTHON_EXE_PATH + ' '
    script = 'pyastrometry_cli_main.py'
    if PYASTROMETRY_SCRIPT_PATH is not None:
        script = os.path.join(PYASTROMETRY_SCRIPT_PATH, script)
    cmd_line += f'{script} slew '
    rastr = target.ra.to_string(u.hour, sep=":", pad=True)
    cmd_line += rastr + ' '
    decstr = target.dec.to_string(alwayssign=True, sep=":", pad=True)
    cmd_line += f'" {decstr}" '


    # use json to handle double quotes in camera and mount
    # arguments which might be passed as something like "CCD Simulator"
    # and the shlex.split() will split them
    if dev_args.mount is not None:
        cmd_line += f'--mount {json.dumps(dev_args.mount)} '

    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    rc, output = run_program(cmd_line, label='precise_slew')

    if rc != 0:
        logging.error(f'run_precise_slew: return code was {rc}!')
        return None

    return True

def run_autofocus(args, extra_args):
    # FIXME need to add parameters for autofocus
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of equipment profile')
    parser.add_argument('--focuser', type=str,  help='Focuser Driver')
    parser.add_argument('--camera', type=str, help='Name of camera driver')
    parser.add_argument('--simul', action='store_true', help='Simulate star')
    parser.add_argument('--simuldatadir', type=str,  help='Location of simulated star data')
    parser.add_argument('--debugplots', action='store_true', help='Show plots')
    parser.add_argument('--framesize', default=0, type=int,  help='Size of capture frame, 0=full')
    parser.add_argument('--focusexposure', type=float, help='Focus exposure time')
    dev_args, unknown = parser.parse_known_args(sys.argv)
    logging.debug(f'dev_args, unknown = {dev_args} {unknown}')

    cmd_line = PYTHON_EXE_PATH + ' '
    if AUTOFOCUS_SCRIPT_PATH is not None:
        cmd_line += f'{AUTOFOCUS_SCRIPT_PATH}/'
    cmd_line += 'autofocus_hfd_script.py '
    #FIXME defaults for C8
#    cmd_line += f'{args.focusmin} {args.focusmax} '
#    cmd_line += f'{args.focusdir} '
    if dev_args.debugplots:
        cmd_line += '--debugplots '
    if dev_args.simul:
        cmd_line += '--simul '
        if dev_args.simuldatadir is not None:
            cmd_line += f'--simuldatadir {dev_args.simuldatadir} '
    if dev_args.framesize:
        cmd_line += f'--framesize {dev_args.framesize} '
    # use json to handle double quotes in camera and mount
    # arguments which might be passed as something like "CCD Simulator"
    # and the shlex.split() will split them
    if dev_args.camera is not None:
        cmd_line += f'--camera {json.dumps(dev_args.camera)} '
    if dev_args.focuser is not None:
        cmd_line += f'--focuser {json.dumps(dev_args.focuser)} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile} '
    if dev_args.focusexposure:
        cmd_line += f'--exposure_start {dev_args.focusexposure}'

    logging.debug(f'run_autofocus cmd_line = {cmd_line}')

    rc, output = run_program(cmd_line, label='autofocus')

    logging.debug(f'autofocus rc = {rc}')
    if rc != 0:
        logging.error(f'run_autofocus: return code was {rc}!')
        #sys.exit(1)
        return False

    return True

def get_altaz_from_radec(radec, observer, obstime):
    altaz = observer.altaz(obstime, target=radec)
    logging.debug(f'get_mount_altaz: computerd alt = {altaz.alt.degree}')
    logging.debug(f'get_mount_altaz: computerd az  = {altaz.az.degree}')
    return altaz.alt.degree, altaz.az.degree

if __name__ == '__main__':

    # FIXME assumes tz is set properly in system?
    now = datetime.now()
    logfilename = 'autofocus_auto_star-' + now.strftime('%Y%m%d%H%M%S') + '.log'

#    FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
    FORMAT = '[%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'

    logging.basicConfig(filename=logfilename,
                        filemode='a',
                        level=logging.DEBUG,
                        format=FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    start_time = time.time()

    parser = argparse.ArgumentParser()
#    parser.add_argument('focusmin', type=int, help='Min allowed focuser pos')
#    parser.add_argument('focusmax', type=int, help='Max allowed focuser pos')
#    parser.add_argument('focusdir', type=str, help='Focus IN or OUT')

    # see if focusonly is present and get the minimum args we need
    parser.add_argument('--focusonly', action='store_true', help='Focus at current position - no slewing')
    parser.add_argument('--profile', type=str, help='Name of astro profile')
    parser.add_argument('--usedebugpaths', action='store_true', help='Run auxilary programs from checked out sources')
    parser.add_argument('--maxtries', type=int, default=3, help='Number of tries before giving up')
    args, extra_args = parser.parse_known_args()

    # if not using focusonly then add in the rest of the args and reparse
    if not args.focusonly:
        parser.add_argument('--dist', type=float, help='Max distance in degrees')
        parser.add_argument('--mag', type=float, help='Desired mag focus star')
        parser.add_argument('--lst', type=str, help='Local sidereal time')
        parser.add_argument('--onlyside', type=str, help='EAST or WEST side only')
        parser.add_argument('--lon', type=float, help='Location longitude')
        parser.add_argument('--meridianthres', type=str, default='00:30:00',
                            help='How close to meridian is allowed (hh:mm:ss)')
        #parser.add_argument('--noplatesolve', action='store_true', help='Just slew do not improve accuracy with plate solving')
        parser.add_argument('--preciseslewstar', action='store_true', help='Use precise slew to star')
        parser.add_argument('--preciseslewreturn', action='store_true', help='Use precise slew returning from star')
        parser.add_argument('--nolower', type=float, default=None, help='How many degrees lower star can be')
        parser.add_argument('--nolowerthres', type=float, default=None, help='Minimum altitude to enforce --nolower')
        parser.add_argument('--nousehorizon', action='store_true', help='Ignore horizon in astroprofile when choosing stars')
        parser.add_argument('--errorsimul', type=int, help='Randomly have autofocus fail for testing - give percentage fail rate')

        args, extra_args = parser.parse_known_args()

    logging.debug(f'args = {args}')
    logging.debug(f'extra_args = {extra_args}')

    # FIXME setup DEBUG paths if we're running all from git
    if args.usedebugpaths:
        # find paths to other modules
        from importlib.util import find_spec
        s = find_spec('hfdfocus')
        head, rest = os.path.split(s.origin)
        AUTOFOCUS_SCRIPT_PATH = os.path.normpath(head + '/../scripts')
        AUTOFOCUS_DATA_PATH = os.path.normpath(head + '/../data')
        s = find_spec('pyastrometry')
        head, rest = os.path.split(s.origin)
        PYASTROMETRY_SCRIPT_PATH = os.path.normpath(head + '/../scripts')
    else:
        AUTOFOCUS_SCRIPT_PATH = None
        AUTOFOCUS_DATA_PATH = None
        PYASTROMETRY_SCRIPT_PATH = None

    logging.debug(f'AUTOFOCUS_SCRIPT_PATH = {AUTOFOCUS_SCRIPT_PATH}')
    logging.debug(f'AUTOFOCUS_DATA_PATH = {AUTOFOCUS_DATA_PATH}')
    logging.debug(f'PYASTROMETRY_SCRIPT_PATH = {PYASTROMETRY_SCRIPT_PATH}')

    # get astro profile if specified
    if args.profile is not None:
        logging.info(f'Using astro profile {args.profile}')
        astro_profile = AstroProfile()
        astro_profile.read(args.profile)
        lon = astro_profile.observatory.location.get('longitude', None)
    else:
        astro_profile = None

    # if just running autofocus here then get to it
    if args.focusonly:
        logging.info('Running autofocus without finding a focus star as requested.')
        retries = args.maxtries
        while retries >= 0:
            focus_result = run_autofocus(args, extra_args)
            if focus_result:
                logging.info('Autofocus successful!')
                sys.exit(0)
            else:
                logging.error(f'Autofocus failed - {retries} retries left!')
                retries = retries - 1
            sys.exit(1)

    # from here on down we are going to find a focus star, slew to it,
    # focus and then return to starting point

    # need --dist --mag if not doing focus only
    if args.dist is None:
        logging.error('Must specify --dist!')
        sys.exit(1)
    if args.mag is None:
        logging.error('Must specify --mag!')
        sys.exit(1)

    if args.preciseslewreturn or args.preciseslewstar:
        logging.info('Getting current position by plate solving')
        cur_radec = run_platesolve()
    else:
        logging.info('Getting current position from mount')
        cur_radec = run_getpos()

    if cur_radec is None:
        logging.error('Could not determine current position!')
        sys.exit(1)

    logging.info(f'Original location is {cur_radec.to_string("hmsdms", sep=":")}')

    star_list = run_findstars(cur_radec, args, lon=lon)
    #logging.debug(f'Star list = {star_list}')

    if star_list is None or len(star_list) < 1:
        logging.error('No star candidate available!')
        sys.exit(1)

    # MAX_NOLOWER_ALT sets an upper limit on the check for the focus star
    # being at a higher altitude that current position.  Otherwise as the
    # scope gets closer and closer to the zenith it would not be possible
    # to find a higher star!
    MAX_NOLOWER_ALT = 65
    result = None
    stars_tried = 0
    for sao, radec in star_list:
        logging.info(f'Trying SAO{sao} at {radec.to_string("hmsdms", sep=":")}')

        # check alititude version current altitude if we have a location
        if astro_profile is not None:
            alt, az = get_altaz_from_radec(cur_radec,
                                           astro_profile.observatory.observer,
                                           Time.now())
            logging.info(f'Current mount alt/az is {alt:3.1f} {az:4.1f}')
            star_alt, star_az = get_altaz_from_radec(radec,
                                                     astro_profile.observatory.observer,
                                                     Time.now())
            logging.info(f'Star alt/az is {star_alt:3.1f} {star_az:4.1f}')
            logging.info(f'args.nousehorizon = {args.nousehorizon}')
            if not args.nousehorizon:
                logging.debug(f'Checking star against horizion definition')
                hor_alt = astro_profile.observatory.horizon.get_alt(star_az)
                logging.info(f'Horiizon alt at az={star_az:4.1f} is {hor_alt:3.1f} ')
                if star_alt <= hor_alt:
                    logging.info('Star is below horizon - skipping!')
                    continue
            logging.debug(f'args.nolower = {args.nolower}')
            if args.nolower is not None:
                logging.info(f'Checking alt of star')

                logging.info(f'nolowerthres = {args.nolowerthres}')
                logging.info(f'MAX_NOLOWER_ALT = {MAX_NOLOWER_ALT}')
                if alt < MAX_NOLOWER_ALT and (args.nolowerthres is None or alt < args.nolowerthres):

#                    star_alt, star_az = get_altaz_from_radec(radec,
#                                                             astro_profile.observatory.observer,
#                                                             Time.now())
#                    logging.info(f'Star alt/az is {star_alt:3.1f} {star_az:4.1f}')

                    if star_alt < alt:
                        logging.info(f'Star SAO{sao} is below current position so skipping!')
                        continue

        if args.preciseslewstar:
            slew_result = run_precise_slew(radec, args, extra_args)
        else:
            slew_result = run_slew(radec, args, extra_args)

        logging.debug(f'slew_result = {slew_result}')

        if slew_result:
            retries = args.maxtries
            result = False
            while retries >= 0:
                focus_result = run_autofocus(args, extra_args)
                if args.errorsimul is not None:
                    import random
                    if int(random.SystemRandom().random()*100) < args.errorsimul:
                        logging.info('Simulating failure because of --simulerrors!')
                        focus_result = False
                if not focus_result:
                    logging.error(f'Autofocus failed - {retries} retries left!')
                    retries = retries - 1
                    continue
                else:
                    logging.info(f'Autofocus succeeded - result = {focus_result}')
                    result = focus_result
                    break
            if result:
                # break out of loop over star candidates
                break
            stars_tried = stars_tried + 1
            if stars_tried < 4:
                logging.error(f'Autofocus failed - trying next candidate!')
                continue
            logging.error(f'Tried {stars_tried} stars unsuccessfully - quitting!')
            sys.exit(1)
        else:
            stars_tried = stars_tried + 1
            if stars_tried < 4:
                logging.error(f'Precise slew failed - trying next candidate!')
                continue
            logging.error(f'Tried {stars_tried} stars unsuccessfully - quitting!')
            sys.exit(1)

    if not result:
        logging.error('Could not autofocus sucessfully')
    else:
        logging.info('Autofocus sucessful!')
        logging.info(f'Autofocus auto star task took {time.time() - start_time} seconds')

    # return to original position
    logging.info(f'Returning to original position {cur_radec.to_string("hmsdms", sep=":")}')
    #return_result = run_precise_slew(cur_radec, args, extra_args)
    if args.preciseslewreturn:
        return_result = run_precise_slew(cur_radec, args, extra_args)
    else:
        return_result = run_slew(cur_radec, args, extra_args)
    logging.debug(f'return_result = {return_result}')
    if not return_result:
        logging.error('Failed to return to original position!')

        sys.exit(1)
    else:
        logging.info('Return to original position suceeded!')

    if not result:
        sys.exit(1)
    else:
        sys.exit(0)


