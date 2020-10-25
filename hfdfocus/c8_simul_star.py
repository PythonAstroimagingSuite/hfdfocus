# simulate C8 @ f/7 star based on v curve data
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

#
#Fit:
#
#2019-05-16 18:51:29,897 INFO     siegel left  fit = (-0.04791842143243222, 382.9895493357121)
#2019-05-16 18:51:29,898 INFO     siegel right fit = (0.04235856430687786, -337.733122948905)
#2019-05-16 18:51:29,898 INFO     siegel best pos  = 7983.459642370256
#2019-05-16 18:51:29,901 INFO     Left Side:
#2019-05-16 18:51:29,901 INFO        slope: -0.04791842143243222
#2019-05-16 18:51:29,901 INFO        inter: 382.9895493357121
#2019-05-16 18:51:29,901 INFO        yzero: 7992.532681314426
#2019-05-16 18:51:29,901 INFO        PID  : -9.073038944169639
#2019-05-16 18:51:29,901 INFO     Right Side:
#2019-05-16 18:51:29,901 INFO        slope: 0.04235856430687786
#2019-05-16 18:51:29,901 INFO        inter: -337.733122948905
#2019-05-16 18:51:29,901 INFO        yzero: 7973.195703756808
#2019-05-16 18:51:29,901 INFO        PID  : 10.263938613447863

import argparse
import logging
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
#import skimage.transform as tf
import scipy.ndimage as ndimage
from hfdfocus.StarFitHFD import find_brightest_star_HFD

# to silence scipy warning about zoom()
import warnings

def shrink_star(starimage_data, bgimage_data, reduction):
    """
    Given a 2D image shrink by the reduction factor and pad the outer areas
    with the median of the image data.

    :param starimage_data: 2D numpy image data to be reduced.
    :param bgimage_data: 2D numpy image data for a background region to be sampled.
    :param reduction: Scaling factor (0-1).
    :return: Scaled down image:
    """

    ht, wd = starimage_data.shape
#    shrunk_star_data = tf.resize(starimage_data,
#                                    [int(ht*reduction), int(wd*reduction)],
#                                    order=0, mode='constant', anti_aliasing=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shrunk_star_data = ndimage.zoom(starimage_data,
                                        [reduction, reduction],
                                        order=0, mode='constant')#, anti_aliasing=True)
    # scale it up so flux is the same
    bgmed = np.median(bgimage_data)
    shrunk_star_data = (shrunk_star_data - bgmed) / reduction / reduction + bgmed
    sh_ht, sh_wd = shrunk_star_data.shape

    # composite shrunk data centered on bg data
    bg_ht, bg_wd = bgimage_data.shape
    lx = int((bg_wd - sh_wd) / 2)
    hx = lx + sh_wd
    ly = int((bg_ht - sh_ht) / 2)
    hy = ly + sh_ht
    #print(bg_ht, bg_wd, sh_ht, sh_wd)
    #print(lx,hx,ly,hy)
    shrunk_image = np.copy(bgimage_data)
    shrunk_image[ly:hy, lx:hx] = np.maximum(shrunk_image[ly:hy, lx:hx], shrunk_star_data)
    #pyfits.writeto(f'shrunk_data.fits', shrunk_data.astype(float), overwrite=True)

    return shrunk_image


class C8_F7_Star_Simulator:
    """
    Simulates the defocused star from a Celestron C8 8 inch f/10 SCT reduced
    to f/7.  Uses actual image data and scales down the image depending on the
    focus position requested.

    :param starimage_name: Name of FITS image data for defocused star.
    :param bgimage_name: Name of FITS image containing background only.
    :param companion_offset: Optional offset (x, y) in pixels for another simulated star.
    :param companion_ratio: Ratio (0-1) of brightness of companion.
    """

    def __init__(self, starimage_name='../data/C8_Simul_Defocus_Star.fit',
                 bgimage_name='../data/C8_Simul_BG.fit',
                 companion_offset=None, companion_flux_ratio=1.0):

        # load background image
        #bgimage_name = 'data/C8_Simul_BG.fit'
        hdu = pyfits.open(bgimage_name)
        self.bgimage_data = hdu[0].data.astype(float)
        hdu.close()

        # load star image
        #starimage_name = 'data/C8_Simul_Defocus_Star.fit'
        hdu = pyfits.open(starimage_name)
        self.starimage_data = hdu[0].data.astype(float)
        hdu.close()

        # will compute when required so we don't
        # create logging and other computations on creation
        self.ref_hfd = None

        # simulate a nearby star
        self.companion_offset = companion_offset
        self.companion_flux_ratio = companion_flux_ratio

    # measured best focus position from sampled V curve
    def get_best_focus_pos(self):
        """
        Returns the focus position for simulated model of smallest star size.

        :returns: Focus position of smallest star size.
        """
        return 7983

    # based on data measured 2019/05/13 on C8 @ f/7
    def simul_hfd_size(self, focus_pos, focus_cen):
        """
        Compute HFD size at focus position based on a model.

        :param focus_pos: Focus position for computation.
        :param focus_cen: Focus position of best focus.
        :return: HFD at requested focus position.
        """
        # load this on demand
        if self.ref_hfd is None:
            # measure star size
            scen, sl, sr, hfl, hfr, totflux, alone = find_brightest_star_HFD(self.starimage_data)
            self.ref_hfd = hfr - hfl
            logging.info(f'Reference star HFD = {self.ref_hfd}')


        # equation is for fit to left side of vcurve
        # we will just mirror it to right side
        # df needs to be negative because of how left side fit was done
        #
        df = -abs(focus_pos - focus_cen)
        hfd = 2.03638808511032E-10*df**4 + 3.38300076148509E-07*df**3 + \
              0.00018332706661*df**2 - 0.011952045599155*df + 2.38830504454329
        return hfd

    # given a desired 'best focus' position focus_cen and
    # a current focuser position return a scaled image of
    # a defocused star based on sampled V curve data
    def get_simul_star_image(self, focus_pos, focus_cen=None):
        """
        Return simulated star image for requested focus position.

        :param focus_pos: Focus position for computation.
        :param focus_cen: Focus position of best focus.
        :return: 2D simulated star image.
        """
        if focus_cen is None:
            focus_cen = self.get_best_focus_pos()
        hfd = self.simul_hfd_size(focus_pos, focus_cen)
        red = min(1.0, hfd / self.ref_hfd)
        shrunk_image = shrink_star(self.starimage_data, self.bgimage_data, red)

        # make another star offset
        if self.companion_offset is not None:
#            tform = tf.SimilarityTransform(scale=1, translation=self.companion_offset)
#            comp_image = tf.warp(shrunk_image, tform, mode='constant',
#                                 cval=np.median(shrunk_image))
            # flip x/y offsets
            ox, oy = self.companion_offset[1], self.companion_offset[0]
            comp_image = ndimage.shift(shrunk_image, (ox, oy), mode='constant',
                                 cval=np.median(shrunk_image))
            shrunk_image = (shrunk_image+self.companion_flux_ratio*comp_image) / (1 + self.companion_flux_ratio)

        return shrunk_image


def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('focus_cen', type=int, help='Focus position for best focus')
    parser.add_argument('focus_pos', type=int, help='Current focus position')
    parser.add_argument('focus_slope', type=float, help='V curve slope')
    parser.add_argument('focus_pid', type=float, help='V curve PID')

    #    parser.add_argument('--debuggraphs', action='store_true', help="Display debug graphs")

    return parser.parse_args()


if __name__ == '__main__':
    logging.basicConfig(filename='c8_simul_star.log',
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

    # load background image
    bgimage_name = 'data/C8_Simul_BG.fit'
    hdu = pyfits.open(bgimage_name)
    bgimage_data = hdu[0].data.astype(float)
    hdu.close()

    # load star image
    starimage_name = 'data/C8_Simul_Defocus_Star.fit'
    hdu = pyfits.open(starimage_name)
    starimage_data = hdu[0].data.astype(float)
    hdu.close()

    # measure star size
    scen, sl, sr, hfl, hfr = find_brightest_star_HFD(starimage_data, debugplots=True)
    ref_hfd = hfr - hfl
    logging.info(f'Reference star HFD = {ref_hfd}')

    shrunk_image = shrink_star(starimage_data, bgimage_data, 0.5)

    scen, sl, sr, hfl, hfr = find_brightest_star_HFD(shrunk_image, debugplots=True, debugfits=True)
    shrunk_hfd = hfr - hfl
    logging.info(f'50% shrunk star HFD = {shrunk_hfd}')

#    focus_cen = 7983
#    focus_ref = 7350
#    focus_slope = 0.04235856430687786
#    focus_pid =  10.263938613447863
#    f = open('c8_simul_hfd.txt', 'w')
#    for fpos in range(7350, 8500, 50):
#        df = fpos - focus_cen
#        hfd = simul_hfd_size(fpos, focus_cen)
#        red = min(1.0, hfd/ref_hfd)
#        shrunk_image = shrink_star(starimage_data, bgimage_data, red)
#        scen, sl, sr, hfl, hfr = find_brightest_star_HFD(shrunk_image, debugplots=True, debugfits=True)
#        shrunk_hfd = hfr-hfl
#        f.write(f'{fpos}, {shrunk_hfd}\n')
#    f.close()

    c8simul = C8_F7_Star_Simulator('data/C8_Simul_Defocus_Star.fit', 'data/C8_Simul_BG.fit')
    focus_cen = 7983
    f = open('c8_simul_hfd_2.txt', 'w')
    for fpos in range(7350, 8500, 50):
        shrunk_image = c8simul.get_simul_star_image(fpos, focus_cen)
        scen, sl, sr, hfl, hfr = find_brightest_star_HFD(shrunk_image, debugplots=True, debugfits=True)
        shrunk_hfd = hfr-hfl
        f.write(f'{fpos}, {shrunk_hfd}\n')
    f.close()




