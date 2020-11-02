#
# Server that listens for requests through stdin and then runs a star fit
# and returns result over stdout.
#
# Copyright 2020 Michael Fulbright
#
#
#    hfdfocus is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
import sys
import time
import json
import logging
from queue import Empty, Queue
from threading import Thread
from multiprocessing import Process  # , Queue
from marshmallow import fields, Schema, ValidationError
from marshmallow.validate import Range

import astropy.io.fits as pyfits

from hfdfocus.StarFitHFR_RadialProfile import star_fit_hfr_radial_profile
from hfdfocus.MultipleStarFitHFD import StarFitResult

def _setup_logging(debug=False):
    LONG_FORMAT = '%(asctime)s.%(msecs)03d [%(filename)20s:%(lineno)3s ' \
                  + '%(funcName)20s() ] %(levelname)-8s %(message)s'
    SHORT_FORMAT = '%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s'

    logging.basicConfig(filename='StarFitHFR_RadialProfile.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # # add to screen as well
    # log = logging.getLogger()
    # formatter = logging.Formatter(SHORT_FORMAT)
    # ch = logging.StreamHandler()
    # if debug:
    #     ch.setLevel(logging.DEBUG)
    # else:
    #     ch.setLevel(logging.INFO)
    # ch.setFormatter(formatter)
    # log.addHandler(ch)

def field_maxstars(*args, **kwargs):
    """ Fields descriptor for maxstars """
    return fields.Integer(*args,
                          **kwargs,
                          validate=Range(min=0, min_inclusive=True,
                                         max=10000, max_inclusive=True))
def field_window(*args, **kwargs):
    """ Fields descriptor for window """
    return fields.Integer(*args,
                          **kwargs,
                          validate=Range(min=0, min_inclusive=True,
                                         max=30, max_inclusive=True))

def field_bgfact(*args, **kwargs):
    """ Fields descriptor for bgfact """
    return fields.Integer(*args,
                          **kwargs,
                          validate=Range(min=0, min_inclusive=True,
                                         max=1000, max_inclusive=True))

class StarFitResultEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, StarFitResult):
            #print(obj)
            rdict = {}
            for key in ['star_cx', 'star_cy', 'star_r1', 'star_r2',
                        'star_r', 'star_angle', 'star_f']:
                val = getattr(obj, key, [])
                logging.debug(f'StarFitResultEncoder: {key} = {val}')
                rdict[key] = list(val)

            for key in ['nstars', 'bgest', 'noiseest', 'width', 'height']:
                logging.debug(f'StarFitResultEncoder: {key} = {val}')
                rdict[key] = getattr(obj, key, None)

            #print(rdict)
            return rdict

# def return_result(status, msg, result=None):
#     rdict = {'Result': status, 'Message': msg}
#     if result is not None:
#         rdict['Value'] = result

#     sys.stdout.write(json.dumps(rdict, cls=StarFitResultEncoder))
#     sys.stdout.write('\n')
#     sys.stdout.flush()


class MeasureHFRServer:

    def __init__(self, q, debug=False):
        #setup_logging(debug)
        self.queue = q
        self.debug = debug

    def return_result(self, status, msg, result=None):
        rdict = {'Result': status, 'Message': msg}
        if result is not None:
            rdict['Value'] = result

        retstr = json.dumps(rdict, cls=StarFitResultEncoder)
        logging.debug(f'retstr = {retstr}')

        self.queue.put(retstr)

    def run(self, request):
        #print('MeasureHFRServer run starting')

        schema_d = dict(filename=fields.String(required=True),
                        request_id=fields.Integer(required=True),
                        maxstars=field_maxstars(missing=100),
                        bgfact=field_bgfact(missing=50),
                        window=field_window(missing=7))

        schema = Schema.from_dict(schema_d)

        if True:
            # for request in sys.stdin:
            #     request = request.strip()
                logging.debug(f'request = |{request}|')

                try:
                    jdict = json.loads(request)
                except json.JSONDecodeError:
                    logging.error(f'Invalid JSON passed: {request}')
                    self.return_result('Error', 'Invalid JSON')
                    return

                # requests are JSON formatted
                # {
                #    'filename' : <image file to be analyzed>,
                #    'request_id' : <integer unique to this request>,
                #    'maxstars' : <max # stars to analyze, optional>,
                #    'bgfact' : <background factor, optional>,
                #    'window' : <window size in pixels for analyzing stars, optional>,
                # }

                logging.debug(jdict)

                result = None
                msg = ''
                try:
                    try_result = schema().load(jdict)
                except ValidationError:
                    msg = 'Failed to validate'
                    logging.error(msg, exc_info=True)
                except ValueError:
                    msg = 'Invalid values'
                    logging.error(msg, exc_info=True)
                else:
                    result = try_result

                if result is None:
                    self.return_result('Error', msg)
                    return
                    #continue

                logging.debug(f'MeasureHFRServer: result={result}')

                try:
                    logging.debug('MeasureHFRServer: opening fits {result["filename"]}')
                    hdu = pyfits.open(result['filename'])
                    logging.debug('MeasureHFRServer: reading data')
                    image_data = hdu[0].data.astype(float)
                    logging.debug('MeasureHFRServer: closing')
                    hdu.close()
                    logging.debug('MeasureHFRServer: done file i/o')
                except:
                    logging.error('Error opening image!', exc_info=True)
                    self.return_result('Error', 'Unable to open image!')
                    return

                try:
                    logging.debug(f'MeasureHFRServer: called starfit')
                    starfit = star_fit_hfr_radial_profile(image_data,
                                                          max_stars=result['maxstars'],
                                                          bgfact=result['bgfact'],
                                                          window=result['window'],
                                                          debugplots=False,
                                                          debugfits=False)
                    logging.debug(f'MeasureHFRServer: startfit = {starfit}')

                except:
                    logging.error('Error analyzing image!', exc_info=True)
                    self.return_result('Error', 'Unable to analyze image!')
                    return

                logging.debug('MeasureHFRServer: Done!')
                self.return_result('Success', 'Completed', result=starfit)
                return

class MeasureHFRClient(Thread):
    def __init__(self, job_queue, result_queue):
        """
        Create client to call the server to measure hfr on an image.

        Expects a message on job_queue which is a dict with the following
        fields:
            filename: name of image file
            request_id: unique integer value for this job
            maxstars: max stars to prcess
            bgfact: BG factor - higher means higher noise rejection
            window: size in pixels of window used to analyze a star

        :param job_queue: Queue for requests.
        :type queue: Multiprocessing.Queue
        :param result_queue: Queue for results.
        :type queue: Multiprocessing.Queue
        """
        super().__init__()
        self.job_queue = job_queue
        self.result_queue = result_queue

    def run(self):
        logging.debug('MeasureHFRCient: start MeasureHFRClient run()')
        while True:
            logging.debug('waiting on job')
            job_dict = self.job_queue.get()
            logging.debug(f'job data = {job_dict}')

            # rdict = dict(filename='SH2-157-m14.0C-gain200-bin_1-300s-Ha-Light-010.fits',
            #               request_id=10,
            #               maxstars=10)
            rstr = json.dumps(job_dict)
            server = MeasureHFRServer(self.result_queue, False)
            #logging.debug(f'rstr = {rstr}')
            p = Process(target=server.run, args=(rstr,))

            logging.debug('start process')
            p.start()
            logging.debug('join process')
            p.join() # this blocks until the process terminates
            logging.debug('process done')
            logging.debug(f'process exitcode = {p.exitcode}')

class MeasureHFRClientStdin:
    def __init__(self):
        """
        Create client to call the server to measure hfr on an image.

        Note: Only runs the first request and exits.

        Listens to stdin for a JSON formatted request with the following fields:
            filename: name of image file
            request_id: unique integer value for this job
            maxstars: max stars to prcess
            bgfact: BG factor - higher means higher noise rejection
            window: size in pixels of window used to analyze a star

        Will output a JSON response with format:
            'Result':  'Success' or 'Error'
            'Message': If error will have details.
            'Value': JSON encoded results for stars detected.
        """
        super().__init__()

    def run(self):
        #logging.debug('MeasureHFRCientStdin run()')
        #sys.stdout.write('waiting\n')
        #sys.stdout.flush()
        while True:
            #logging.debug('waiting on job')
            sys.stdout.write('waiting\n')
            sys.stdout.flush()
            for request in sys.stdin:
                request = request.strip()
                #logging.debug(f'request = |{request}|')
                #sys.stdout.write(f'request = |{request}|\n')
                #sys.stdout.flush()

                if request == "exit":
                    sys.stdout.write('exiting\n')
                    sys.stdout.flush()
                    sys.exit(0)

                try:
                    jdict = json.loads(request)
                except json.JSONDecodeError:
                    #logging.error(f'Invalid JSON passed: {request}')
                    sys.stdout.write(json.dumps({'Result': 'Error',
                                                 'Message': 'Invalid JSON'}))
                    sys.stdout.write('\n')
                    sys.stdout.flush()
                    sys.exit(1)

                #logging.debug(f'job data = {jdict}')
                #sys.stdout.write(f'job data = {jdict}\n')
                #sys.stdout.flush()

                rstr = json.dumps(jdict)
                #logging.debug(f'rstr = {rstr}')
                result_queue = Queue()
                s = MeasureHFRServer(result_queue, debug=False)
                #logging.debug('start process')
                s.run(rstr)
                stars = result_queue.get()
                sys.stdout.write(f'{stars}\n')
                sys.stdout.flush()

                time.sleep(10)
                sys.stdout.write('done\n')
                sys.stdout.flush()

                sys.exit(0)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str, help='Target')
    parser.add_argument('--bgfact', type=int, default=50, help='BG Factor')
    parser.add_argument('--maxstars', type=int, default=500,
                        help='Max stars detected')
    parser.add_argument('--window', type=int, default=7,
                        help='Star detection window (pixels)')
    parser.add_argument('--debug', action='store_true',
                        help='show debug info')

    args = parser.parse_args()

    print('reading')
    hdu = pyfits.open(args.infile)
    data = hdu[0].data
    hdu.close()
    print('done reading')

   # sys.exit(0)



    _setup_logging((args.debug))

    job_queue = Queue()
    result_queue = Queue()
    logging.info('test')

    #client = MeasureHFRClient(job_queue, result_queue).start()
    #client.run()
    time.sleep(1)

    print('reading')
    hdu = pyfits.open(args.infile)
    data = hdu[0].data
    hdu.close()
    print('done reading')

    #sys.exit(0)

    rdict = dict(filename=args.infile, #filename='SH2-157-m14.0C-gain200-bin_1-300s-Ha-Light-010.fits',
                  request_id=10,
                  maxstars=10)

    s = MeasureHFRServer(result_queue, debug=False)
    s.run(json.dumps(rdict))

    #sys.exit(0)

    print('trying 1st file')
    #job_queue.put(rdict)
    print('waiting on result')
    while True:
        print('checking queue')
        try:
            result = result_queue.get(False)
        except Empty:
            print('empty')
        else:
            print('not empty')
            break
        time.sleep(0.1)
    print('done waiting in main')
    print('reading result queue')
    print('stars = ', result)

    sys.exit(0)

    # rdict = dict(filename='SH2-157-m14.0C-gain200-bin_1-300s-Ha-Light-015.fits',
    #               request_id=10,
    #               maxstars=10)

    # print('trying 2nd file')
    # job_queue.put(rdict)
    # print('waiting on result')
    # while True:
    #     print('checking queue')
    #     try:
    #         result = result_queue.get(False)
    #     except Empty:
    #         print('empty')
    #     else:
    #         break
    #     time.sleep(0.1)
    sys.exit(0)

    # print('starting thread')
    # x = Thread(target=wait_thread, args=(queue, 1))
    # x.start()
    # print('waiting on result')
    # while True:
    #     print('checking queue')
    #     try:
    #         result = queue.get(False)
    #     except Empty:
    #         print('empty')
    #     else:
    #         break
    #     time.sleep(0.1)
    # print('result=', result)

    # sys.exit(0)



    infile = args.infile

    hdu = pyfits.open(infile)
    #print(hdu[0].data)
    image_data = hdu[0].data.astype(float)
    hdu.close()

    starfit = star_fit_hfr_radial_profile(image_data,
                                          max_stars=args.maxstars,
                                          bgfact=args.bgfact,
                                          window=args.window,
                                          debugplots=False,
                                          debugfits=False)

    print(starfit)
