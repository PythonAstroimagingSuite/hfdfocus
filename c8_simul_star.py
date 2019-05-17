# simulate C8 @ f/7 star based on v curve data
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
from skimage.transform import resize as image_resize
from StarFitHFD import find_brightest_star_HFD

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('focus_cen', type=int, help='Focus position for best focus')
    parser.add_argument('focus_pos', type=int, help='Current focus position')
    parser.add_argument('focus_slope', type=float, help='V curve slope')
    parser.add_argument('focus_pid', type=float, help='V curve PID')

    #    parser.add_argument('--debuggraphs', action='store_true', help="Display debug graphs")

    return parser.parse_args()

def shrink_star(starimage_data, bgimage_data, reduction):
    ht, wd = starimage_data.shape
    shrunk_star_data = image_resize(starimage_data,
                                    [int(ht*reduction), int(wd*reduction)],
                                    order=0, mode='constant', anti_aliasing=True)

    # scale it up so flux is the same
    bgmed = np.median(bgimage_data)
    shrunk_star_data = (shrunk_star_data-bgmed)/reduction/reduction + bgmed
    sh_ht, sh_wd = shrunk_star_data.shape

    # composite shrunk data centered on bg data
    bg_ht, bg_wd = bgimage_data.shape
    lx = int((bg_wd-sh_wd)/2)
    hx = lx + sh_wd
    ly = int((bg_ht-sh_ht)/2)
    hy = ly + sh_ht
    print(bg_ht, bg_wd, sh_ht, sh_wd)
    print(lx,hx,ly,hy)
    shrunk_image = np.copy(bgimage_data)
    shrunk_image[ly:hy, lx:hx] = np.maximum(shrunk_image[ly:hy, lx:hx], shrunk_star_data)
    #pyfits.writeto(f'shrunk_data.fits', shrunk_data.astype(float), overwrite=True)

    return shrunk_image


class C8_F7_Star_Simulator:
    def __init__(self, starimage_name, bgimage_name):

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

        # measure star size
        scen, sl, sr, hfl, hfr = find_brightest_star_HFD(self.starimage_data)
        self.ref_hfd = hfr-hfl
        logging.info(f'Reference star HFD = {self.ref_hfd}')


    # based on data measured 2019/05/13 on C8 @ f/7
    def simul_hfd_size(self, focus_pos, focus_cen):
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
    def get_simul_star_image(self, focus_pos, focus_cen):
        hfd = self.simul_hfd_size(focus_pos, focus_cen)
        red = min(1.0, hfd/self.ref_hfd)
        shrunk_image = shrink_star(self.starimage_data, self.bgimage_data, red)
        return shrunk_image


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
    ref_hfd = hfr-hfl
    logging.info(f'Reference star HFD = {ref_hfd}')

    shrunk_image = shrink_star(starimage_data, bgimage_data, 0.5)

    scen, sl, sr, hfl, hfr = find_brightest_star_HFD(shrunk_image, debugplots=True, debugfits=True)
    shrunk_hfd = hfr-hfl
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




