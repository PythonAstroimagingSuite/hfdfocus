#
# star fitting routines
#
# Copyright 2019 Michael Fulbright
#
#
#    This program is free software: you can redistribute it and/or modify
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

import sys
import math
import time
import logging
import argparse
import numpy as np
from scipy import signal, ndimage
from scipy.interpolate import interp1d
from scipy.integrate import romb

import astropy.io.fits as pyfits

# for testing
import matplotlib.pyplot as plt


#
# using lots of window x window samples compute median
#
def compute_bg_model(image, window):
    ht, wd = image.shape

    bgmodel = np.empty_like(image)
    for y in range(0, ht, window):
        yl = y
        ym = min(ht-1, yl+window)

        for x in range(0, wd, window):
            xl = x
            xm = min(wd-1, xl+window)

#            print(f'sample range y={yl}:{ym}  x={xl}:{xm} median=', np.median(image[yl:ym, xl:xm]))

            bgmodel[yl:ym, xl:xm] = np.median(image[yl:ym, xl:xm])

    return bgmodel

# assumes gray data
def compute_noise_level(data):
    mad = np.median(np.abs(data - np.median(data)))
    return mad

def compute_median(data):
    return np.median(data)

def find_centroid(image_data, thres):

    image_data = image_data - thres
    image_data[image_data < 0] = 0

    total = (image_data).sum()
    iidx, jidx = np.indices(image_data.shape)

    iidx = iidx + 0.5
    jidx = jidx + 0.5

    wsum_i =(iidx*(image_data)).sum()
    wsum_j = (jidx*(image_data)).sum()

    ci = wsum_i/total
    cj = wsum_j/total

    return (ci, cj)

def find_star(image_data, bgfact=50, satur=50000, window=100,
              starmodel=False, bgmodel=False, debugfits=False):
    logging.info(f'find_star: debugfits = {debugfits}')

    logging.info(f'find_stars start: bgfact={bgfact} window={window}')

    if debugfits:
        pyfits.writeto('thres_test.fits', image_data.astype(float), overwrite=True)

    ttot_s = time.time()
    if bgmodel:
        logging.info('compute_bg_model START')
        ts = time.time()
        bgmodel = compute_bg_model(image_data, 100)
        te = time.time()
        logging.info(f'compute_bg_model DONE took {te-ts} seconds')
    else:
        bg = compute_median(image_data)
        logging.info(f'using constant bg model = {bg}')
        bgmodel = np.full(image_data.shape, bg)

    if debugfits:
        pyfits.writeto('bgmodel_test.fits', bgmodel.astype(float), overwrite=True)

    # compute image with bg removed
    bgrem_data = image_data.astype(float) - bgmodel.astype(float)

    if debugfits:
        pyfits.writeto('bgrem_test.fits', bgrem_data.astype(float), overwrite=True)

    # find otsu threshold for each window
    ht, wd = bgrem_data.shape
    bg_window = 100

    logging.info('computing star_model START')
    ts = time.time()
    star_model = np.zeros_like(bgrem_data)
    if starmodel:
        for y in range(0, ht, bg_window):
            yl = y
            ym = min(ht-1, yl+bg_window)

            for x in range(0, wd, bg_window):
                xl = x
                xm = min(wd-1, xl+bg_window)

                data = bgrem_data[yl:ym, xl:xm]
                data_med = np.median(data)
                data_mad = compute_noise_level(data)

                thres = data_med + bgfact*data_mad
                star_model[yl:ym, xl:xm] = bgrem_data[yl:ym, xl:xm] > thres
    else:
        thres = compute_median(bgrem_data) + bgfact*compute_noise_level(bgrem_data)
        logging.info(f'Using constant thres for star model = {thres}')
        star_model = bgrem_data > thres

    te = time.time()
    logging.info(f'computing star_model DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'star_model_bgfact{bgfact}.fits', star_model.astype(float), overwrite=True)

    # use dilation to stregthen any structures
    ndilate = 3
    tmp_image = star_model
    logging.info(f'dilating star_model {ndilate} times')
    ts = time.time()
    for i in range(0, ndilate):
        #tmp_image = ndimage.grey_dilation(tmp_image, size=(3,3))
        tmp_image = ndimage.binary_dilation(tmp_image, np.ones((3,3)))
    te = time.time()
    logging.info(f'computing dilation took {te-ts} seconds')

    dilate_star_model = tmp_image
    if debugfits:
        pyfits.writeto(f'dilate_tar_model_bgfact{bgfact}.fits', dilate_star_model.astype(float), overwrite=True)

    # find structures
    star_label, nstars = ndimage.measurements.label(dilate_star_model)

    #print(star_label, nstars)
    if debugfits:
        pyfits.writeto(f'star_label.fits', star_label.astype(float), overwrite=True)

    logging.info('Finding stars with >= 9 pix START')
    star_boxes = np.zeros_like(star_model)
    max_pix = 0
    max_l = None
    for l in ndimage.find_objects(star_label):
        #print(f'l={l}')
        npix = (l[0].stop-l[0].start)*(l[1].stop-l[1].start)
        #print(bgrem_data[l].shape, npix)

        if npix < 9:
            print('not enuf pix')
            continue

#        if (bgrem_data[l] > satur).sum() > 0:
#            print('bgrem too high')
#            continue

        if npix > max_pix:
            max_pix = npix
            max_l = l

    if max_l is None:
        logging.error('find_star(): no object found')
        return None

    star_boxes[max_l] = 1

    if debugfits:
        pyfits.writeto(f'star_boxes.fits', star_boxes.astype(float), overwrite=True)

    cy, cx = ndimage.measurements.center_of_mass(star_boxes)
    logging.info(f'COM cx, cy = {cx}, {cy}')

    # compute background near star
    #
    yl = max(0, int(cy)-bg_window)
    ym = min(ht-1, yl+bg_window)
    xl = max(0, int(cx)-bg_window)
    xm = min(wd-1, xl+bg_window)
    #print('xl/xm/yl/xm=', xl, xm, yl, ym)
    data = bgmodel[yl:ym, xl:xm]
    bglevel = np.median(data)
    data = bgrem_data[yl:ym, xl:xm]
    bgmad = compute_noise_level(data)
    logging.info(f'bg = {bglevel} mad = {bgmad}')

    ttot_e = time.time()

    logging.info(f'find_star took {ttot_e-ttot_s} seconds')

    return cx, cy, bglevel, bgmad, star_boxes


# horizontal bin data to form a 1D star profile
# compute radial profile - returns (rad[], val[])
def horiz_bin_window(data, bg=0):
    #ts = time.time()

    ni, nj = data.shape

    profile = np.sum(data-bg, axis=0)

#    print(data)
#    print(profile)

    return profile

# find left/right limits of star disk
def find_star_limits_robust(profile, thres=0):

    # search from left side (idx=0) for first pixel over thres
    #idx = np.where(profile > thres)[0]

    left = None
    right = None
    idx = 0
    while idx < len(profile)-2:
        #print('l', idx, profile[idx], profile[idx+1], thres)
        if profile[idx+1] >= thres and profile[idx] < thres:
            left = idx
            break
        idx += 1

    if left is None:
        return None, None

    idx = len(profile)-2
    while idx > left+2:
        #print('r', idx, profile[idx], profile[idx-1], thres)
        if profile[idx-1] >= thres and profile[idx] < thres:
            right = idx-1
            break
        idx -= 1

    return left, right

# more robust find left/right limits of star disk
def find_star_limits(profile, thres=0):

    # start from left and find adjacent pixels where leftmost
    # is less than thres and rightmost is greater than thres

    # search from left side (idx=0) for first pixel over thres
    idx = np.where(profile > thres)[0]
    #print('where -> ', idx)
    if len(idx) < 2:
        return None, None

    left = np.min(idx)
    right = np.max(idx)

    return left, right


# find HFD from 1D profile
def find_hfd_from_1D(profile, thres=0, debugplots=False):

    # compute flux within values lx and rx in profile.
    # proffunc should be an scipy interp1d object representing
    # the 1D profile of the star which has been background
    # subtracted
#    def flux(proffunc, lx, rx):
#        return quad(proffunc, lx, rx, epsabs=1, limit=100)[0]

    # different way
    def flux2(proffunc, rlo, rhi):
        nr = 2**7+1
        rvals, dr = np.linspace(rlo, rhi, num=nr, retstep=True)
        fvals = proffunc(rvals)
        return romb(fvals, dr)

    # compute center of profile with weighted average
    plen = profile.shape[0]

    #print('profile len = ', plen)

    idx = np.arange(0, plen, 1)
    #print(idx)

    #lidx, ridx = find_star_limits(profile, thres)
    lidx, ridx = find_star_limits_robust(profile, thres)
    logging.info(f'left, right = {lidx}, {ridx}')
    if lidx is None or ridx is None:
        return None

    cidx = np.sum(idx[lidx:ridx]*profile[lidx:ridx])/np.sum(profile[lidx:ridx])
    logging.info(f'center = {cidx}')

    # find left, right limits
    # invert profile as y as a function of x on left and right
    # sides of the profile
    lidx_low = int(lidx)
    lidx_hi = int(lidx+2)
    ridx_low = int(ridx)
    ridx_hi = int(ridx+2)
#    print(lidx, lidx_low, lidx_hi, ridx, ridx_low, ridx_hi)
#    print(profile[lidx_low:lidx_hi], idx[lidx_low:lidx_hi])
#    print(profile[ridx_low:ridx_hi], idx[ridx_low:ridx_hi])
    lfunc = interp1d(profile[lidx_low:lidx_hi], idx[lidx_low:lidx_hi])
    rfunc = interp1d(profile[ridx_low:ridx_hi], idx[ridx_low:ridx_hi])

    lx = lfunc(thres)
    rx = rfunc(thres)
    proffunc = interp1d(idx, profile)

    if debugplots:
        fig = plt.figure()
        ax_1 = fig.add_subplot(131)
        ax_2 = fig.add_subplot(132)
        ax_3 = fig.add_subplot(133)
        ax_1.plot(profile[lidx_low:lidx_hi], lfunc(profile[lidx_low:lidx_hi]))
        ax_1.axhline(lx, color='green')
        ax_2.plot(profile[ridx_low:ridx_hi], rfunc(profile[ridx_low:ridx_hi]))
        ax_2.axhline(rx, color='green')

#    print('lx=', lx, thres, proffunc(lx))
#    print('rx=', rx, thres, proffunc(rx))

    # make sure inversion worked
    INVERSE_FRAC = 0.15
    if abs(thres-proffunc(lx))/thres > INVERSE_FRAC or abs(thres-proffunc(rx))/thres > INVERSE_FRAC:
        logging.error('inversion failed!')
        logging.error(f'lx={lx} thres={thres} fit={proffunc(lx)}')
        logging.error(f'rx={rx} thres={thres} fit={proffunc(rx)}')
        if debugplots:
            fig2 = plt.figure()
            ax_4 = fig2.add_subplot(111)
            ax_4.plot(np.arange(lidx_low,ridx_hi), profile[lidx_low:ridx_hi])
            plt.show()
        return None

    # now find flux inside left/right
    simple_totflux = np.sum([profile[lidx:ridx]])
    logging.info(f'flux in star (simple) = {simple_totflux}')

    totflux = flux2(proffunc, lx, rx)
    #print('integrated flux in star = ', flux(proffunc, lidx, ridx))
    logging.info(f'integrated flux in star = {totflux}')

    # compute flux as function of distance from center
    r_max = min(cidx-lx, rx-cidx)

    r_arr = []
    #flux_vs_r = []
    flux2_vs_r = []
    for r in range(0, int(r_max)):
        r_arr.append(r)
        rlo = cidx - r
        rhi = cidx + r
        #flux_vs_r.append(flux(proffunc, rlo, rhi))

        # different way
#        nr = 2**7+1
#        rvals, dr = np.linspace(rlo, rhi, num=nr, retstep=True)
#        fvals = proffunc(rvals)
#
#        flux2_vs_r.append(romb(fvals, dr))

        flux2_vs_r.append(flux2(proffunc, rlo, rhi))

    #print(f'flux2_vs_r, r = {flux2_vs_r} {r_arr}')

    # make inverse function of r vs flux to find half flux
    flux_inv = interp1d(flux2_vs_r, r_arr)

    # uncomment to see raw profile
    if debugplots:
        ax_3.plot(r_arr, flux2_vs_r)
    #plt.show()

    half_flux_r = flux_inv(totflux/2)

    logging.info(f'half flux rad = {half_flux_r}')
    #print(r_arr)
    #print(flux_vs_r)
    #print(flux2_vs_r)
    #ax_3.plot(r_arr, flux_vs_r)
    if debugplots:
        ax_3.axvline(half_flux_r, color='green')

    return (cidx, lx, rx, cidx-half_flux_r, cidx+half_flux_r, totflux)

def find_brightest_star_HFD(image_data, thres=10000, win=100, debugplots=False, debugfits=False):
    logging.info(f'find_brightest_star_HFD: debugfits = {debugfits}')
    xcen, ycen, bg, noise, starmask = find_star(image_data, debugfits=debugfits)

    img_ht, img_wd = image_data.shape

    xlow = max(0, int(xcen-win/2))
    xhi = min(img_wd-1, int(xcen+win/2))
    ylow = max(0, int(ycen-win/2))
    yhi = min(img_ht-1, int(ycen+win/2))
    crop_data = image_data[ylow:yhi, xlow:xhi]

    if debugplots:
        fig = plt.figure()
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)
        im = ax_2d.imshow((crop_data-bg).astype(float))
        fig.colorbar(im, ax=ax_2d)

    profile = horiz_bin_window(crop_data, bg=bg)

    logging.info(f'thres = {thres}')
    logging.info(f'profile = {profile}')

    if debugplots:
        ax_1d.plot(profile)
        plt.draw()

    return find_hfd_from_1D(profile, thres=thres)


def test_1d_with_gaussian():

    # generate 2D gaussian sample star
    sigma_x = 10
    sigma_y = sigma_x
    size = 45
    peak = 200000/sigma_x/sigma_y/2/np.pi
    bg = 800
    bgnoise = 50
    thres = math.sqrt(size)*5*bgnoise

    x = np.linspace(-size/2, size/2, size)
    y = np.linspace(-size/2, size/2, size)

    x, y = np.meshgrid(x, y)
    #z = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-(x**2/(2*sigma_x**2)
    #     + y**2/(2*sigma_y**2))))
    z = (peak * np.exp(-(x**2/(2*sigma_x**2)
         + y**2/(2*sigma_y**2))))
    z[np.where((x**2+y**2) < 3)] = 0
    z += np.random.normal(loc=bg, scale=bgnoise, size=z.shape)

    #plt.contourf(x, y, z, cmap='Blues')
    fig = plt.figure()
    ax_2d = fig.add_subplot(121)
    ax_1d = fig.add_subplot(122)
    im = ax_2d.imshow(z-peak)
    fig.colorbar(im, ax=ax_2d)

    profile = horiz_bin_window(z, bg=bg)

    print('thres = ', thres)
    print('profile = ', profile)

    scen, sl, sr = find_hfd_from_1D(profile, thres=thres)

    ax_1d.plot(profile)

    ax_1d.axvline(scen, color='red')
    if sl is not None and sr is not None:
        ax_1d.axvline(sl, color='green')
        ax_1d.axvline(sr, color='green')

    print('total counts = ', np.sum(z-bg))

    plt.show()

    sys.exit(1)





if __name__ == '__main__':
    logging.basicConfig(filename='StarFitHFD.log',
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
    parser.add_argument('infile', type=str, help='Target')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')

    args = parser.parse_args()

#    logging.info(f'command args = {args}')

    infile = args.infile

    hdu = pyfits.open(infile)
    print(hdu[0].data)
    image_data = hdu[0].data.astype(float)
    hdu.close()

    bg = 800
    thres = 10000

    xcen, ycen = find_star(image_data, debugfits=True)

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

    scen, sl, sr, hfl, hfr = find_hfr_from_1D(profile, thres=thres)

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

    f = open('hfd.txt', 'a')
    f.write(f'{infile[26:30]}, {hfr-hfl}\n')
    f.close()



