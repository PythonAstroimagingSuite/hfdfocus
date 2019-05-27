#
#
# currently relies on being run from hfdfocus/scripts directory
# and having pyastrometry_cli checked out in same tree
#

import os
import sys
import json
import shlex
import argparse
import logging
import subprocess
import numpy as np

# FIXME this should be handled automaticall!

if os.name == 'nt':
    PYTHON_EXE_PATH = 'C:\\Users\\msf\\Anaconda3\\envs\\AstronomyUtilitiesPipNumpy\\python'
elif os.name == 'posix':
    PYTHON_EXE_PATH = '/home/msf/anaconda3/envs/pyastro37/bin/python3 '
else:
    logging.error('PYTHON_EXE_PATH NOT SET')
    sys.exit(1)

def run_platesolve(pixelscale):
    # parse out device info
    # FIXME seems alot of duplication need a better way to represent
    # device info like a config file

    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of equipment profile')
    dev_args, unknown = parser.parse_known_args(sys.argv)

    logging.debug(f'run_platesolve: dev_args, unknown = {dev_args} {unknown}')

    result_fname = './autofocus_auto_origpos.json'
    cmd_line = PYTHON_EXE_PATH + ' '
    cmd_line += '../../pyastrometry/scripts/pyastrometry_cli_main.py solvepos '
    cmd_line += f'--outfile {result_fname} '
    cmd_line += f'--pixelscale {pixelscale} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    cmd_args = shlex.split(cmd_line)

    logging.debug(f'run_platesolve() command line = |{cmd_line}|')
    logging.debug(f'run_platesolve() command args = |{cmd_args}|')

    # unlink previous solve if any
    if os.path.isfile(result_fname):
        os.unlink(result_fname)

    if os.name == 'nt':
        cmd_val = cmd_line
    elif os.name == 'posix':
        cmd_val = cmd_args
        
    ps_proc = subprocess.Popen(cmd_val,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)

    logging.debug('ps_proc output:')
    for l in ps_proc.stdout.readlines():
        logging.debug(f'ps_proc: {l.strip()}')
    logging.debug('end of output')

    poll_value = None
    while True:
        poll_value = ps_proc.poll()

        if poll_value is not None:
            break

    # check return code
    rc = ps_proc.returncode

    if rc < 0:
        logging.error(f'run_platesolve: return code was {rc}!')
        return None

    try:
        result_file = open(result_fname, 'r')
    except OSError as err:
        print(f"Error opening results file: {err}")
        return None

    try:
        solve_dict = json.load(result_file)
        target_str = solve_dict['ra2000'] + " "
        target_str += solve_dict['dec2000']
        logging.info(f"target_str = {target_str}")

        radec = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                          frame='fk5', equinox='J2000')

    except Exception as err:
        print(f"Error converting solve results! {err}")
        return None

    logging.info(f'radec = {radec.to_string()}')
    return radec

def run_findstars(curpos, args):
    result_fname = './autofocus_auto_starlist.dat'
    cmd_line = PYTHON_EXE_PATH + ' '
    #cmd_line = 'python '
    cmd_line += 'find_nearby_stars.py '
    cmd_line += '../data/SAO_Catalog_m5_p11_filtered.bin '
    rastr = curpos.ra.to_string(u.hour, sep=":", pad=True)
    cmd_line += rastr + ' '
    decstr = curpos.dec.to_string(alwayssign=True, sep=":", pad=True)
    cmd_line += f'" {decstr}" '
    cmd_line += f'{args.dist} '
    cmd_line += f'--outfile {result_fname} '
    mag = float(args.mag)
    cmd_line += f'--minmag {mag+0.5} '
    cmd_line += f'--maxmag {mag-0.5} '
    cmd_line += f'--lst {args.lst} '
    cmd_line += f'--onlypierside {args.onlypierside}'
    if args.meridianthres is not None:
        cmd_line += f' --meridianthres {args.meridianthres}'

    cmd_args = shlex.split(cmd_line)

    logging.debug(f'run_findstars() command line = |{cmd_line}|')
    logging.debug(f'run_findstars() command args = |{cmd_args}|')

    # unlink previous solve if any
    if os.path.isfile(result_fname):
        os.unlink(result_fname)

    if os.name == 'nt':
        cmd_val = cmd_line
    elif os.name == 'posix':
        cmd_val = cmd_args
        
    ps_proc = subprocess.Popen(cmd_val,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)

    logging.debug('ps_proc output:')
    for l in ps_proc.stdout.readlines():
        logging.debug(f'ps_proc: {l.strip()}')
    logging.debug('end of output')

    poll_value = None
    while True:
        poll_value = ps_proc.poll()

        if poll_value is not None:
            break

    # check return code
    rc = ps_proc.returncode

    if rc != 0:
        logging.error(f'run_findstars: return code was {rc}!')
        return None

    try:
        result_file = open(result_fname, 'r')
    except OSError as err:
        print(f"Error opening results file: {err}")
        return None

    star_list = []
    try:
        nlines = 0
        for l in result_file.readlines():
            print(l.strip())
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
            logging.info(f"star target_str = {target_str}")
            radec = SkyCoord(target_str, unit=(u.hourangle, u.deg),
                              frame='fk5', equinox='J2000')
            star_list.append((sao, radec))

    except Exception as err:
        print(f"Error converting solve results! {err}")
        return None

    logging.info(f'star_list = {star_list}')
    return star_list

def run_precise_slew(target, args, extra_args):
    # parse out device info
    # FIXME seems alot of duplication need a better way to represent
    # device info like a config file

    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', type=str, help='Name of equipment profile')
    parser.add_argument('--telescope', type=str, help='Name of telescope driver')
    parser.add_argument('--camera', type=str, help='Name of camera driver')
    parser.add_argument('--exposure', type=float, default=5, help='Exposure time')
    parser.add_argument('--binning', type=int, default=2, help='Camera binning')
    parser.add_argument('--framesize', default=0, type=int,  help='Size of capture frame, 0=full')
    dev_args, unknown = parser.parse_known_args(sys.argv)

    logging.debug(f'dev_args, unknown = {dev_args} {unknown}')

    cmd_line = PYTHON_EXE_PATH + ' '
    cmd_line += '../../pyastrometry/scripts/pyastrometry_cli_main.py slewsolve '
    rastr = target.ra.to_string(u.hour, sep=":", pad=True)
    cmd_line += rastr + ' '
    decstr = target.dec.to_string(alwayssign=True, sep=":", pad=True)
    cmd_line += f'" {decstr}" '
    cmd_line += f'--pixelscale {args.pixelscale} '
    if dev_args.framesize is not None:
        cmd_line += f'--framesize {dev_args.framesize} '
        
    # use json to handle double quotes in camera and telescope
    # arguments which might be passed as something like "CCD Simulator"
    # and the shlex.split() will split them
    if dev_args.camera is not None:
        cmd_line += f'--camera {json.dumps(dev_args.camera)} '
    if dev_args.telescope is not None:
        cmd_line += f'--telescope {json.dumps(dev_args.telescope)} '
    cmd_line += f'--exposure {dev_args.exposure} '
    cmd_line += f'--binning {dev_args.binning} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    cmd_args = shlex.split(cmd_line)

    logging.debug(f'run_precise_slew() command line = |{cmd_line}|')
    logging.debug(f'run_precise_slew() command args = |{cmd_args}|')

    # unlink previous solve if any
#    if os.path.isfile(result_fname):
#        os.unlink(result_fname)

    if os.name == 'nt':
        cmd_val = cmd_line
    elif os.name == 'posix':
        cmd_val = cmd_args
        
    ps_proc = subprocess.Popen(cmd_val,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)

    logging.debug('ps_proc output:')
    for l in ps_proc.stdout.readlines():
        logging.debug(f'ps_proc: {l.strip()}')
    logging.debug('end of output')

    poll_value = None
    while True:
        poll_value = ps_proc.poll()

        if poll_value is not None:
            break

    # check return code
    rc = ps_proc.returncode

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
    parser.add_argument('--debugplots', action='store_true', help='Show plots')
    parser.add_argument('--framesize', default=0, type=int,  help='Size of capture frame, 0=full')
    dev_args, unknown = parser.parse_known_args(sys.argv)
    logging.debug(f'dev_args, unknown = {dev_args} {unknown}')

    cmd_line = PYTHON_EXE_PATH + ' '
    cmd_line += 'autofocus_hfd_script.py '
    #FIXME defaults for C8
    cmd_line += '6000 12000 IN '
    if dev_args.debugplots:
        cmd_line += '--debugplots '
    if dev_args.simul:
        cmd_line += '--simul '
    if dev_args.framesize:
        cmd_line += f'--framesize {dev_args.framesize} '
    # use json to handle double quotes in camera and telescope
    # arguments which might be passed as something like "CCD Simulator"
    # and the shlex.split() will split them
    if dev_args.camera is not None:
        cmd_line += f'--camera {json.dumps(dev_args.camera)} '
    if dev_args.focuser is not None:
        cmd_line += f'--focuser {json.dumps(dev_args.focuser)} '
    if dev_args.profile is not None:
        cmd_line += f'--profile {dev_args.profile}'

    cmd_args = shlex.split(cmd_line)

    logging.debug(f'autofocus command line = |{cmd_line}|')
    logging.debug(f'autofocus command args = |{cmd_args}|')

    # unlink previous solve if any
#    if os.path.isfile(result_fname):
#        os.unlink(result_fname)

    if os.name == 'nt':
        cmd_val = cmd_line
    elif os.name == 'posix':
        cmd_val = cmd_args
        
    ps_proc = subprocess.Popen(cmd_val,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True)

    logging.debug('ps_proc output:')
    for l in ps_proc.stdout.readlines():
        logging.debug(f'ps_proc: {l.strip()}')
    logging.debug('end of output')

    poll_value = None
    while True:
        poll_value = ps_proc.poll()

        if poll_value is not None:
            break

    # check return code
    rc = ps_proc.returncode
    logging.debug(f'autofocus rc = {rc}')
    if rc != 0:
        logging.error(f'run_autofocus: return code was {rc}!')
        return None

    return True


from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

if __name__ == '__main__':
    logging.basicConfig(filename='autofocus_auto_star.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    log.addHandler(ch)

    parser = argparse.ArgumentParser()
    parser.add_argument('dist', type=float, help='Max distance in degrees')
    parser.add_argument('mag', type=float, help='Desired mag focus star')
    parser.add_argument('lst', type=str, help='Local sidereal time')
    parser.add_argument('onlypierside', type=str, help='EAST or WEST side only')
    parser.add_argument('pixelscale', type=float, help='Pixel scale for plate solving')
    parser.add_argument('--meridianthres', type=str, default='00:30:00',
                        help='How close to meridian is allowed (hh:mm:ss)')

    args, extra_args = parser.parse_known_args()

    logging.info(f'args = {args}')
    logging.info(f'extra_args = {extra_args}')

    cur_radec = run_platesolve(args.pixelscale)

    logging.info(f'Original location is {cur_radec.to_string("hmsdms", sep=":")}')

    star_list = run_findstars(cur_radec, args)

    logging.info(f'Star list = {star_list}')

    if star_list is None or len(star_list) < 1:
        logging.error('No star candidate available!')
        sys.exit(1)

    result = None
    for sao, radec in star_list:
        logging.info(f'Trying SAO{sao} at {radec.to_string("hmsdms", sep=":")}')

        slew_result = run_precise_slew(radec, args, extra_args)

        logging.debug(f'slew_result = {slew_result}')

        if slew_result:
            focus_result = run_autofocus(args, extra_args)
            if not focus_result:
                logging.error('Autofocus failed - trying next candidate!')
                continue
            else:
                logging.info(f'Autofocus suceeded - result = {focus_result}')
                result = focus_result
                break
        else:
            logging.error('Precise slew failed!')
            continue

    if not result:
        logging.error('Could not autofocus sucessfully')

    # return to original position
    logging.info(f'Returning to original position {cur_radec.to_string("hmsdms", sep=":")}')
    return_result = run_precise_slew(cur_radec, args, extra_args)
    logging.debug(f'return_result = {return_result}')
    if not return_result:
        logging.error('Failed to return to original position!')

        sys.exit(1)

    logging.info('Autofocus and return to original position suceeded!')
    sys.exit(0)


