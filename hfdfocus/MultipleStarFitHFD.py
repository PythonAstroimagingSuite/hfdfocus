#
# Star fitting routines for mutiple stars - based on StarFitHFD.py but
# generalized to multiple stars instead of a single star.
#
# NOTE: Code is duplicated and refactored from StarFitHFD.py - decided it was
#       better to star over and design it better and leave the existing
#       single star code alone for now as it is used only by autofocus and
#       is working really well!
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
from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.integrate import romb

import astropy.io.fits as pyfits

# for testing
import matplotlib as mpl
#mpl.use("Qt5agg")
print(mpl.get_backend())
mpl.use("TkAgg")
print(mpl.get_backend())
#mpl.rcParams['toolbar'] = 'None'
#mpl.rc('font', size=8)
import matplotlib.pyplot as plt


@dataclass
class StarInfo:
    """
    Container for information about a single detected star.

    Attributes:
        cx:                Centroid X
        cy:                Centroid Y
        bglevel:           Background level at star
        bgmad:             Background MAD at star
        estsize:           Estimated size in pixels of detected blob - NOT HFD measurement!
        estflux:           Estimated flux from detected blob - NOT precise!
        hfd:               HFD
        flux:              HFD derived flux measurement
        alone:             Whether another star was nearby
    """
    cx: float
    cy: float
    bglevel: float
    bgmad: float
    estsize: float
    estflux: float
    hfd: float
    flux: float
    alone: bool

@dataclass
class DetectedStars:
    """
    Container for information generated when star detection is run on an image.

    Attributes:
        bg_window:          Background around star computed in bg_window*2 x bg_window*2 box
        bgmodel:           Modeled background level at pixel
        bgrem_data:        Image with background subtracted
        star_model:        True if pixel was above detection threshold
        dilate_star_model: Dilated version of star model
        star_label:        Labels for objects in star_model
        nlabels:           Number of labels created
        stars:             List of StarInfo objects for detected stars in order of contained flux
    """
    bg_window: int
    bgmodel: np.ndarray
    bgrem_data: np.ndarray
    star_model: np.ndarray
    dilate_star_model: np.ndarray
    star_label: np.ndarray
    nlabels: int
    stars: list

#
# using lots of window x window samples compute median
#
def compute_bg_model(image, window):
    """
    Compute a model of the background of an image by using a grid of sample
    boxes each (window x window) pixels and replacing pixels within with the
    median value of the box.

    :param image: Numpy 2D array image data.
    :param window: Size of sampling window in pixels.

    :return: Image containing background model.
    """
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
    """
    Compute the median absolute deviation (MAD) of image data.

    :param data: Numpy 2D array image data.
    :return: MAD of image data.
    :rtype: float
    """
    mad = np.median(np.abs(data - np.median(data)))
    return mad

def compute_median(data):
    """
    Compute the median of image data.

    :param data: Numpy 2D array image data.
    :return: Median of image data.
    :rtype: float
    """
    return np.median(data)

# from https://alyssaq.github.io/2015/computing-the-axes-or-orientation-of-a-blob/
def raw_moment(data, i_order, j_order):
    """
    Compute the raw image and central moments.
    Based on from https://alyssaq.github.io/2015/computing-the-axes-or-orientation-of-a-blob/.

    :param data: Numpy 2D array image data.
    :param i_order, j_order: Indices of moment matrix to compute
    :return: Request moment.
    :rtype: float
    """
    nrows, ncols = data.shape
    y_indices, x_indicies = np.mgrid[:nrows, :ncols]
    return (data * x_indicies**i_order * y_indices**j_order).sum()

def moments_cov(data):
    """
    Compute the second order central moments and covariance matrix.
    Based on from https://alyssaq.github.io/2015/computing-the-axes-or-orientation-of-a-blob/.

    :param data: Numpy 2D array image data.
    :return: Covariance matrix.
    """
    data_sum = data.sum()
    m10 = raw_moment(data, 1, 0)
    m01 = raw_moment(data, 0, 1)
    x_centroid = m10 / data_sum
    y_centroid = m01 / data_sum
    u11 = (raw_moment(data, 1, 1) - x_centroid * m01) / data_sum
    u20 = (raw_moment(data, 2, 0) - x_centroid * m10) / data_sum
    u02 = (raw_moment(data, 0, 2) - y_centroid * m01) / data_sum
    cov = np.array([[u20, u11], [u11, u02]])
    return cov

def get_major_minor_axes(data):
    """
    Compute the vector for major and minor axes.
    Based on from https://alyssaq.github.io/2015/computing-the-axes-or-orientation-of-a-blob/.

    :param data: Numpy 2D array image data.
    :returns: Tuple of two tuples, each containing the X, Y components of vector and its norm.
    :rtype: ((float, float, float), (float, float, float))
    """
    cov = moments_cov(data)
    evals, evecs = np.linalg.eig(cov)

    sort_indices = np.argsort(evals)[::-1]
    x_v1, y_v1 = evecs[:, sort_indices[0]]  # Eigenvector with largest eigenvalue
    x_v2, y_v2 = evecs[:, sort_indices[1]]

    n_v1 = evals[sort_indices][0]
    n_v2 = evals[sort_indices][1]
#    print(evals)
#    print(evals[sort_indices])
    logging.debug(f'Major axis = {x_v1}, {y_v1}, {n_v1}')
    logging.debug(f'Minor axis = {x_v2}, {y_v2}, {n_v2}')

    return ((x_v1, y_v1, n_v1), (x_v2, y_v2, n_v2))


# end of moment code

def find_centroid(image_data, thres):
    """
    Compute centroid of image data.

    :param image_data: Numpy 2D array image data.
    :param thres: Ignore pixels less than this value, also subtracted before calculation.
    :returns: Tuple of centroid X, Y values.
    :rtype: (float, float)
    """
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

def detect_stars(image_data, bgfact=50, satur=50000, window=100,
              starmodel=False, bgmodel=False, debugfits=False):
    """
    Find the stars in given image data.  All pixels less than bgfact*MAD over
    the background are ignored.  A region of (window x window) pixels is used to
    search of stars.  If this value is too small then large defocused stars will be
    missed.

    :param image_data: Numpy 2D image data.
    :param bgfact:  Used to set threshold for rejecting background pixels.
    :param satur: Any pixel over this value is considered saturated.
    :param window: Size of square window used for calculating potential star parameters.
    :param starmodel: Whether a model of star pixels is used - usually disabled.
    :param bgmodel: Whether a model of background is computed or just median used.
    :param debugfits: If True then FITS files of various stages of computation are output.
    :returns:
        DetectedStars object
    """
    logging.debug(f'detect_stars start: bgfact={bgfact} window={window}')

    if debugfits:
        pyfits.writeto('thres_test.fits', image_data.astype(float), overwrite=True)

    ttot_s = time.time()
    if bgmodel:
        logging.debug('compute_bg_model START')
        ts = time.time()
        bgmodel = compute_bg_model(image_data, 100)
        te = time.time()
        logging.debug(f'compute_bg_model DONE took {te-ts} seconds')
    else:
        bg = compute_median(image_data)
        logging.debug(f'using constant bg model = {bg}')
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

    logging.debug('computing star_model START')
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
        logging.debug(f'Using constant thres for star model = {thres}')
        star_model = bgrem_data > thres

    te = time.time()
    logging.debug(f'computing star_model DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'star_model_bgfact{bgfact}.fits', star_model.astype(float), overwrite=True)

    # use dilation to stregthen any structures
    ndilate = 3
    tmp_image = star_model
    logging.debug(f'dilating star_model {ndilate} times')
    ts = time.time()
    for i in range(0, ndilate):
        #tmp_image = ndimage.grey_dilation(tmp_image, size=(3,3))
        tmp_image = ndimage.binary_dilation(tmp_image, np.ones((3,3)))
    te = time.time()
    logging.debug(f'computing dilation took {te-ts} seconds')

    dilate_star_model = tmp_image
    if debugfits:
        pyfits.writeto(f'dilate_tar_model_bgfact{bgfact}.fits', dilate_star_model.astype(float), overwrite=True)

    # find structures
    star_label, nlabels = ndimage.measurements.label(dilate_star_model)

    #print(star_label, nlabels)
    if debugfits:
        pyfits.writeto(f'star_label.fits', star_label.astype(float), overwrite=True)

    logging.debug('Finding stars with >= 9 pix START')
#    max_pix = 0
#    max_l = None
#    for l in ndimage.find_objects(star_label):
#        #print(f'l={l}')
#        npix = (l[0].stop-l[0].start)*(l[1].stop-l[1].start)
#        #print(bgrem_data[l].shape, npix)
#
#        if npix < 9:
#            logging.warning('not enuf pix')
#            continue
#
##        if (bgrem_data[l] > satur).sum() > 0:
##            print('bgrem too high')
##            continue
#
#        if npix > max_pix:
#            max_pix = npix
#            max_l = l
#
#    if max_l is None:
#        logging.error('find_star(): no object found')
#        return None

    # create return object now we know there is at least one star
    detected = DetectedStars( bg_window, bgmodel, bgrem_data,
                             star_model, dilate_star_model,
                             star_label, nlabels, [])

    # now loop over labels and measure stars
    debug_star_boxes = np.zeros_like(star_model)
    star_results = []
    for cur_l in ndimage.find_objects(star_label):
        star_boxes = np.zeros_like(star_model)
        npix_cur_l = (cur_l[0].stop-cur_l[0].start)*(cur_l[1].stop-cur_l[1].start)

        if npix_cur_l < 9:
            #logging.warning('not enuf pix')
            continue

        # see if another object is within window centered on main object
        alone = True
        logging.debug('looking for nearby stars')
        for l in ndimage.find_objects(star_label):
            npix = (l[0].stop-l[0].start)*(l[1].stop-l[1].start)
            logging.debug(f'{cur_l}, {l}, {npix_cur_l}, {npix}')
            if l[0].start == cur_l[0].start and l[1].start == cur_l[1].start and \
               l[0].stop == cur_l[0].stop and l[1].stop == cur_l[1].stop:
                   logging.debug('skipping is original label')
                   continue

            if npix > npix_cur_l/2:
                alone = False
                break

        if not alone:
            logging.warning('ERROR: find_star() Another star is nearby!')

        star_boxes[cur_l] = 1
        debug_star_boxes[cur_l] = 1

        axes = get_major_minor_axes(dilate_star_model[cur_l])
        majax = axes[0][2]
        minax = axes[1][2]
        ecc = math.sqrt(1.0-(minax/majax)**2)
        logging.debug(f'Major/Minor/Ecc = {majax} {minax} {ecc}')

        if ecc > 0.75:
            logging.warning(f'fERROR find_star(): Eccentricity {ecc} is outside acceptable range - probably not alone')
            alone = False

        cy, cx = ndimage.measurements.center_of_mass(star_boxes)
        logging.debug(f'COM cx, cy = {cx}, {cy}')
    #    logging.debug(f'box={max_l}')
    #    logging.debug(f'boxperim = {2*pw+2*ph} predperim={np.pi*(pw+ph)/2} actperim = {perim}')

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

        estsize = math.sqrt((cur_l[0].stop-cur_l[0].start)*(cur_l[0].stop-cur_l[0].start)+
                            (cur_l[1].stop-cur_l[1].start)*(cur_l[1].stop-cur_l[1].start))

        estflux = np.sum(bgrem_data[cur_l])
        logging.debug(f'cx = {cx} cy = {cy} bg = {bglevel} mad = {bgmad} '
                      f'estsize = {estsize} estflux = {estflux}')

        res = StarInfo(cx, cy, bglevel, bgmad, estsize, estflux, None, None, alone)

        star_results.append(res)

    if debugfits:
        pyfits.writeto(f'multiple_star_boxes.fits', debug_star_boxes.astype(float), overwrite=True)

    detected.stars = star_results

    ttot_e = time.time()

    logging.info(f'detect_stars took {ttot_e-ttot_s} seconds')

    return detected




# horizontal bin data to form a 1D star profile
# compute radial profile - returns (rad[], val[])
def horiz_bin_window(data, bg=0):
    """
    Bin image data in 1D.

    :param data: Numpy 2D image data.
    :param bg: Background value to subtract from data before binning.
    :return: 1D numpy profile.
    """
    #ts = time.time()

    ni, nj = data.shape

    profile = np.sum(data-bg, axis=0)

#    print(data)
#    print(profile)

    return profile

# find left/right limits of star disk
def find_star_limits_robust_from_edges(profile, thres=0):
    """
    Given a 1D star profile search to right and left of peak for edges of profile.

    :param profile: 1D star profile (numpy array)
    :param thres: Ignore all pixels below this value.
    :returns:
        Tuple containing index for left and right extremes of star profile.
    :rtype: (int, int)
    """
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

# find left/right limits of star disk
# start from center of profile
def find_star_limits_robust_from_center(profile, thres=0):
    """
    Given a 1D star profile search to right and left of peak for edges of profile.

    :param profile: 1D star profile (numpy array)
    :param thres: Ignore all pixels below this value.
    :returns:
        Tuple containing index for left and right extremes of star profile.
        Left limit is FIRST pixel to left of 'star' BELOW thres
        Right limit is LAST pixel to right of 'star' ABOVE thres
    :rtype: (int, int)
    """
    # search from left side (idx=0) for first pixel over thres
    #idx = np.where(profile > thres)[0]

    #assume 'peak' is in center of profile

    left = None
    right = None

    # check to left first
    #logging.debug('Finding left')
    idx = int(len(profile)/2)
    while idx > 1:
        #logging.debug(f'{idx} {thres} {profile[idx]} {profile[idx-1]}')
        if profile[idx] >= thres and profile[idx-1] < thres:
            left = idx-1
            break
        idx -= 1

    if left is None:
        return None, None

    idx = int(len(profile)/2)
    #logging.debug('Finding right')
    while idx < len(profile)-2:
        #logging.debug(f'{idx} {thres} {profile[idx]} {profile[idx+1]}')
        if profile[idx] >= thres and profile[idx+1] < thres:
            right = idx
            break
        idx += 1

    return left, right

# find left/right limits of star disk
# start from center of profile
def find_star_limits_robust_from_max_region(profile, thres=0):
    """
    Given a 1D star profile search to right and left of peak for edges of profile.

    :param profile: 1D star profile (numpy array)
    :param thres: Ignore all pixels below this value.
    :returns:
        Tuple containing index for left and right extremes of star profile.
        Left limit is FIRST pixel to left of 'star' BELOW thres
        Right limit is LAST pixel to right of 'star' ABOVE thres
    :rtype: (int, int)
    """
    # search from left side (idx=0) for first pixel over thres
    #idx = np.where(profile > thres)[0]

    #assume 'peak' is in center of profile

    left = None
    right = None

    from scipy.ndimage import find_objects, label
    profile_arr = np.array(profile)
    profile_thres = (profile_arr > thres).astype(int)
    logging.debug(f'profile_thres = {profile_thres}')
    profile_label, label_count = label(profile_thres)
    reg = find_objects(profile_label)
    logging.debug(f'region = {reg}')

    if reg is None:
        return None, None

    maxreg = None
    maxsum = 0
    for a in reg:
        rsum = np.sum(profile[a])
        logging.debug(f'{a[0].start} {a[0].stop} {rsum}')
        if rsum > maxsum:
            maxreg = a[0]
            maxsum = rsum

    logging.debug(f'max region {maxreg.start} {maxreg.stop} {maxsum}')

    left = maxreg.start-1
    right = maxreg.stop-1

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
    """
    Given a 1D star profile compute Half Flux Diameter (HFD).

    :param profile: 1D star profile (numpy array)
    :param thres: Ignore all pixels below this value.
    :param debugplots: If True then matplotlib plots of progress will be generated.
    :return: HFD of star profile in pixels.
    :rtype: (float)
    """
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
    lidx, ridx = find_star_limits_robust_from_edges(profile, thres)
    logging.debug(f'from edge left, right = {lidx}, {ridx}')
    lidx, ridx = find_star_limits_robust_from_center(profile, thres)
    logging.debug(f'from center left, right = {lidx}, {ridx}')
    lidx, ridx = find_star_limits_robust_from_max_region(profile, thres)
    logging.debug(f'from max region left, right = {lidx}, {ridx}')
    if lidx is None or ridx is None:
        return None

    cidx = np.sum(idx[lidx:ridx]*profile[lidx:ridx])/np.sum(profile[lidx:ridx])
    logging.debug(f'center = {cidx}')

    # find left, right limits
    # invert profile as y as a function of x on left and right
    # sides of the profile
    lidx_low = int(lidx)
    lidx_hi = int(lidx+2)
    ridx_low = int(ridx)
    ridx_hi = int(ridx+2)
    print(lidx, lidx_low, lidx_hi, ridx, ridx_low, ridx_hi)
    print(profile[lidx_low:lidx_hi], idx[lidx_low:lidx_hi])
    print(profile[ridx_low:ridx_hi], idx[ridx_low:ridx_hi])
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
    logging.debug(f'flux in star (simple) = {simple_totflux}')

    totflux = flux2(proffunc, lx, rx)
    #print('integrated flux in star = ', flux(proffunc, lidx, ridx))
    logging.debug(f'integrated flux in star = {totflux}')

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
    try:
        flux_inv = interp1d(flux2_vs_r, r_arr)

        # uncomment to see raw profile
        if debugplots:
            ax_3.plot(r_arr, flux2_vs_r)
        #plt.show()

        half_flux_r = flux_inv(totflux/2)
    except:
        logging.error('Error finding half flux rad!', exc_info=True)
        return None

    logging.debug(f'half flux rad = {half_flux_r}')
    #print(r_arr)
    #print(flux_vs_r)
    #print(flux2_vs_r)
    #ax_3.plot(r_arr, flux_vs_r)
    if debugplots:
        ax_3.axvline(half_flux_r, color='green')
        #plt.pause(15)

    return (cidx, lx, rx, cidx-half_flux_r, cidx+half_flux_r, totflux)


def find_star_HFD(star, image_data, thres=10000, win=100, debugplots=False, debugfits=False):
    # yuck for testing!
    global fig, ax_2d, ax_1d

    #xcen, ycen, bg, noise, starmask, alone = find_star(image_data, debugfits=debugfits)

    xcen = star.cx
    ycen = star.cy
    bg = star.bglevel
    mad = star.bgmad
    estsize = star.estsize/1.5
    logging.info(f'cx: {xcen} cy: {ycen} - using estsize of {estsize} for profile window')

    img_ht, img_wd = image_data.shape

    xlow = max(0, int(xcen-estsize))
    xhi = min(img_wd-1, int(xcen+estsize))
    ylow = max(0, int(ycen-estsize))
    yhi = min(img_ht-1, int(ycen+estsize))
    crop_data = image_data[ylow:yhi, xlow:xhi]

    if debugplots:
        fig = plt.figure()
        ax_2d = fig.add_subplot(121)
        ax_1d = fig.add_subplot(122)
        im = ax_2d.imshow((crop_data-bg).astype(float))
        fig.colorbar(im, ax=ax_2d)

    profile = horiz_bin_window(crop_data, bg=bg)

    logging.debug(f'thres = {thres}')
    logging.debug(f'profile = {profile}')

    if debugplots:
        ax_1d.plot(profile)
        ax_1d.axhline(bg, color='green')
        ax_1d.axhline(bg+25*mad, color='green')

        plt.draw()

    return find_hfd_from_1D(profile, thres=bg+10*mad)

def measure_stars_in_image(image_data):

    bg = 800
    thres = 10000

    #xcen, ycen = find_star(image_data, debugfits=False)
    detected = detect_stars(image_data, debugfits=True)
    print(detected)

    if detected is not None and len(detected.stars) > 0:
        for star in detected.stars:

            rc = find_star_HFD(star, image_data, debugplots=True, debugfits=False)

            if rc is None:
                continue

            cidx, lx, rx, lr, rr, totflux = rc
            hfd = rr - lr

            if args.debugplots:
                #hfw = 2*(hfr-hfl)
                #hfc = (hfr+hfl)/2
                #ax_1d.set_xlim(hfc-hfw, hfc+hfw)
                fig.suptitle(f'HFD {hfd:5.2f}')

                ax_1d.axvline(cidx, color='red')
                if lx is not None and rx is not None:
                    ax_1d.axvline(lx, color='green')
                    ax_1d.axvline(rx, color='green')
                    ax_1d.axvline(lr, color='blue')
                    ax_1d.axvline(rr, color='blue')
                    delta = rx-lx
                    ax_1d.set_xlim(lx-delta/4, rx+delta/4)

            fig.show()
            plt.pause(1)

        if args.debugplots:
            plt.show(block=True)

if __name__ == '__main__':

    LONG_FORMAT = '%(asctime)s.%(msecs)03d [%(filename)20s:%(lineno)3s - %(funcName)20s() ] %(levelname)-8s %(message)s'
    SHORT_FORMAT = '%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s'

    logging.basicConfig(filename='MultipleStarFitHFD.log',
                        filemode='w',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # add to screen as well
    log = logging.getLogger()
    formatter = logging.Formatter(LONG_FORMAT) #'%(asctime)s %(levelname)-8s %(message)s')
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


    measure_stars_in_image(image_data)

    sys.exit(0)

    bg = 800
    thres = 10000

    #xcen, ycen = find_star(image_data, debugfits=False)
    detected = detect_stars(image_data, debugfits=True)
    print(detected)

    if detected is not None and len(detected.stars) > 0:
        for star in detected.stars:

            rc = find_star_HFD(star, image_data, debugplots=True, debugfits=False)

            if rc is None:
                continue

            cidx, lx, rx, lr, rr, totflux = rc
            hfd = rr - lr

            if args.debugplots:
                #hfw = 2*(hfr-hfl)
                #hfc = (hfr+hfl)/2
                #ax_1d.set_xlim(hfc-hfw, hfc+hfw)
                fig.suptitle(f'HFD {hfd:5.2f}')

                ax_1d.axvline(cidx, color='red')
                if lx is not None and rx is not None:
                    ax_1d.axvline(lx, color='green')
                    ax_1d.axvline(rx, color='green')
                    ax_1d.axvline(lr, color='blue')
                    ax_1d.axvline(rr, color='blue')
                    delta = rx-lx
                    ax_1d.set_xlim(lx-delta/4, rx+delta/4)

            fig.show()
            plt.pause(1)

        if args.debugplots:
            plt.show(block=True)




#            win = 100
#            xlow = int(s.cx-win/2)
#            xhi = int(s.cx+win/2)
#            ylow = int(s.cy-win/2)
#            yhi = int(s.cy+win/2)
#            crop_data = image_data[ylow:yhi, xlow:xhi]
#
#            print(xlow, xhi, ylow, yhi)
#
#            if args.debugplots:
#                fig = plt.figure()
#                ax_2d = fig.add_subplot(121)
#                ax_1d = fig.add_subplot(122)
#                im = ax_2d.imshow((crop_data-bg).astype(float))
#                fig.colorbar(im, ax=ax_2d)
#
#            plt.show()
#
#            profile = horiz_bin_window(crop_data, bg=bg)
#
#            print('thres = ', thres)
#            print('profile = ', profile)
#
#            if args.debugplots:
#                ax_1d.plot(profile)
#
#            scen, sl, sr, hfl, hfr = find_hfd_from_1D(profile, thres=thres)
#
#            if args.debugplots:
#                ax_1d.axvline(scen, color='red')
#                if sl is not None and sr is not None:
#                    ax_1d.axvline(sl, color='green')
#                    ax_1d.axvline(sr, color='green')
#                    ax_1d.axvline(hfl, color='blue')
#                    ax_1d.axvline(hfr, color='blue')
#                    delta = sr-sl
#                    ax_1d.set_xlim(sl-delta/4, sr+delta/4)
#
#            print('total counts = ', np.sum(crop_data-bg))
#
#            if args.debugplots:
#                plt.show()
#
#            f = open('hfd.txt', 'a')
#            f.write(f'{infile[26:30]}, {hfr-hfl}\n')
#            f.close()



