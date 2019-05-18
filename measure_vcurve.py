import os
import re
import glob
import logging
import argparse
import numpy as np

import astropy.io.fits as pyfits

# for testing
import matplotlib.pyplot as plt

from StarFitHFD import find_hfd_from_1D, find_star, horiz_bin_window


if __name__ == '__main__':
    logging.basicConfig(filename='measure_vcurve.log',
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

    #test_1d_with_gaussian()
    parser = argparse.ArgumentParser()
    parser.add_argument('imgdir', type=str, help='Directory containing V Curve Images')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')

    args = parser.parse_args()

#    logging.info(f'command args = {args}')

    imgfiles = glob.glob(os.path.join(args.imgdir, '*.fit'))

    logging.info(f'imgfiles = {imgfiles}')

    focus_re = re.compile('^.*focuspos_(?P<pos>\d{2,}).fit$')

    for infile in sorted(imgfiles):

        focus_pos = focus_re.match(infile).group(1)
        logging.info(f'infile = {infile} focus_pos = {focus_pos}')

        hdu = pyfits.open(infile)
        image_data = hdu[0].data.astype(float)
        hdu.close()

        bg = 800
        thres = 10000

        xcen, ycen, bg, mad = find_star(image_data, debugfits=True)

        win = 100
        xlow = int(xcen-win/2)
        xhi = int(xcen+win/2)
        ylow = int(ycen-win/2)
        yhi = int(ycen+win/2)
        crop_data = image_data[ylow:yhi, xlow:xhi]

        if args.debugplots:
            fig = plt.figure()
            ax_2d = fig.add_subplot(121)
            ax_1d = fig.add_subplot(122)
            im = ax_2d.imshow((crop_data-bg).astype(float))
            fig.colorbar(im, ax=ax_2d)

        profile = horiz_bin_window(crop_data, bg=bg)

        print('thres = ', thres)
        print('profile = ', profile)

        if args.debugplots:
            ax_1d.plot(profile)

        scen, sl, sr, hfl, hfr = find_hfd_from_1D(profile, thres=thres)

        if args.debugplots:
            ax_1d.axvline(scen, color='red')
            if sl is not None and sr is not None:
                ax_1d.axvline(sl, color='green')
                ax_1d.axvline(sr, color='green')
                ax_1d.axvline(hfl, color='blue')
                ax_1d.axvline(hfr, color='blue')
                delta = sr-sl
                ax_1d.set_xlim(sl-delta/4, sr+delta/4)

        print('total counts = ', np.sum(crop_data-bg))

        if args.debugplots:
            plt.show()

        f = open(os.path.join(args.imgdir, 'hfd.txt'), 'a')
        f.write(f'{focus_pos}, {hfr-hfl}\n')
        f.close()
