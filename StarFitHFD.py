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
from scipy.integrate import quad, romb

import astropy.io.fits as pyfits

# for testing
import matplotlib.pyplot as plt

class StarFitResult:
    def __init__(self, star_cx, star_cy, star_r, star_f,
                 nstars, bgest, noiseest, width, height):
        self.star_cx = star_cx
        self.star_cy = star_cy
        self.star_r = star_r
        self.star_f = star_f
        self.nstars = nstars
        self.bgest = bgest
        self.noiseest = noiseest
        self.width = width
        self.height = height

    # compute split of stars in different annuli in percentage of image diagonal

    def compute_hfr_in_out(self, r_ex=0, r_in=0.35, r_gap=0.1, r_out=0.9):

        hfr_result_in = self.filter_range(r_low=r_ex, r_high=r_in)
        hfr_result_out = self.filter_range(r_low=r_in+r_gap, r_high=r_out)

        if hfr_result_in.nstars > 0:
            hfr_in = np.median(hfr_result_in.star_r)
        else:
            hfr_in = np.nan

        if hfr_result_out.nstars > 0:
            hfr_out = np.median(hfr_result_out.star_r)
        else:
            hfr_out = np.nan

        self.hfr_in = hfr_in
        self.hfr_out = hfr_out

        return

        diag2 = (self.width**2+self.height**2)/4

        halfwid = self.width/2
        halfht  = self.height/2

        thres_ex2 = diag2*r_ex**2
        thres_in2 = diag2*r_in**2
        thres_gap2 = diag2*(r_in+r_gap)**2
        thres_out2 = diag2*r_out**2

        hfr_in = []
        hfr_out = []

        for x, y, r in zip(self.star_cx, self.star_cy, self.star_r):
            rad2 = (x-halfwid)**2+(y-halfht)**2
            if rad2 < thres_ex2:
                continue
            if rad2 < thres_in2:
                hfr_in.append(r)
                continue
            if rad2 > thres_gap2 and rad2 < thres_out2:
                hfr_out.append(r)

        self.hfr_in = np.median(hfr_in)
        self.hfr_out = np.median(hfr_out)

    # return StarFitResult object with only stars in specified range
    # radius is specified as fraction (0 to 1) of 1/2 diagonal distance
    def filter_range(self, r_low=0, r_high=1):
        diag2 = (self.width**2+self.height**2)/4

        halfwid = self.width/2
        halfht  = self.height/2

        thres_low = diag2*r_low**2
        thres_high = diag2*r_high**2

        rad2 = (self.star_cx-halfwid)**2 + (self.star_cy-halfht)**2

        valid_r_idx = np.where(np.logical_and(rad2 >= thres_low, rad2 <= thres_high))

        logging.info(f'valid_r_idx = {valid_r_idx[0].shape[0]}')

        return StarFitResult(self.star_cx[valid_r_idx],
                             self.star_cy[valid_r_idx],
                             self.star_r[valid_r_idx],
                             self.star_f[valid_r_idx],
                             valid_r_idx[0].shape[0],
                             self.bgest,
                             self.noiseest,
                             self.width,
                             self.height)

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


# compute radial profile - returns (rad[], val[])
def compute_radial_profile(data, ci, cj, bg):
    #ts = time.time()

    ni, nj = data.shape

    rad = np.zeros(ni*nj+1)
    val = np.zeros(ni*nj+1)
    rad[0] = 0.0
    val[0] = 0.0

    # build array containing X indices and Y indices
    ij = np.mgrid[0:ni, 0:nj]
    i=ij[0]
    j=ij[1]

    di = i - ci
    dj = j - cj
    rad  = np.sqrt(di*di + dj*dj).flatten()
    val = (data - bg).flatten()

    # first sort by rad
    sortidx = rad.argsort()
    rad = np.concatenate(([0], rad[sortidx]))
    val = np.concatenate(([0], val[sortidx]))

    #print "compute_radial_profile took", time.time() - ts

    return (rad, val)

def find_hfr_from_radial_profile(rad, val, window, extras=False):
    # compute integrated flux - limit to just the circular area of diameter window
    fluxrad = np.zeros(len(rad))
    fluxtot = 0.0
    nneg = 0
    for i in range(0, len(rad)):

        if rad[i] > window:
            break

        # OLD CODE TO STOP INTEGRATION3
        if nneg < 100:
            fluxtot += val[i]
            fluxrad[i] = fluxtot
        else:
            fluxrad[i] = fluxrad[i-1]

        if i > 0:
            if fluxrad[i] < fluxrad[i-1]:
                nneg += 1
        # END OLD CODE

        # NEW CODE USING MAXIMUM
        fluxtot += val[i]
        fluxrad[i] = fluxtot
        # END NEW CODE

# DEBUG
        if False: #or DEBUG_FITSTAR_LOGGING:
            if i == 0:
                print(rad[i], val[i], fluxrad[i])
            else:
                print(rad[i], val[i], fluxrad[i], (fluxrad[i]-fluxrad[i-1])/(fluxrad[i]))

# END DEBUG

    # capture length of filled entries
    radlen = i

    # NEW CODE FOR MAXIMUM CONT'D
    fluxmaxidx = np.argmax(fluxrad)
    fluxmax = fluxrad[fluxmaxidx]

    #Dprint "fluxmax at idx ",fluxmaxidx, " value=", fluxmax

    # set rest of array to max value
    fluxrad[fluxmaxidx:] = fluxmax

    #DEBUG RADPROF
    # fidx = 0
    # while True and fidx < 20:
        # try:
            # f=open("radprof-%03d.csv" % (fidx,), "w")
            # for i in range(0, len(rad)):
                # f.write("%f,%f,%f\n" % (rad[i], val[i], fluxrad[i]))
            # f.close()
            # break
        # except:
            # fidx += 1
            # print "fidx = ",fidx

    #END DEBUG RADPROF

    # did we detect a star?
    if fluxtot <= 0:
        return None

    # find half flux
    rfit = interp1d(rad, fluxrad)

    rmin = 0.0
    rmax = np.amax(rad)
    ftarg = fluxmax/2.0
    while True:
        rmid = (rmin+rmax)/2.0
        f = rfit(rmid)
        #print(ftarg, f, rmin, rmid, rmax)
        if f < ftarg:
            rmin = rmid
        elif f > ftarg:
            rmax = rmid
        else:
            break

# DEGUG
#        if DEBUG_FITSTAR_LOGGING:
#            print(rmid, rmin, rmax, ftarg, f)
# END DEBUG

        if (rmax-rmin) < 0.001:
            break

    # compute contrast (ratio of center value divided by hfr value
    #
    proffit = interp1d(rad, val)
    #contrast = val[1]/proffit(rmid)
    # FIXME contrast doesnt work for defocused stars and we're getting really small denominators (or 0!)
    contrast = 1

    # return position relative to corner of crop window
    if not extras:
        return (rmid, 0, fluxtot, contrast)
    else:
        return ((rmid, 0, fluxtot, contrast), (rad, val, fluxrad, fluxtot))

def star_fit_hfd(image_data, max_stars=100, bgfact=50, satur=50000, window=7,
                 debugfits=False):

    logging.info(f'star_fit_hfr start: max_stars={max_stars} bgfact={bgfact} window={window}')

    if debugfits:
        pyfits.writeto('thres_test.fits', image_data.astype(float), overwrite=True)

    logging.info('compute_bg_model START')
    ttot_s = time.time()
    ts = time.time()
    bgmodel = compute_bg_model(image_data, 100)
    te = time.time()
    logging.info(f'compute_bg_model DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto('bgmodel_test.fits', bgmodel.astype(float), overwrite=True)

    # smooth
#    smkern = np.array([[1,1,1], [1,0,1], [1, 1, 1]])
#    print(smkern.shape, smkern)
#    bgmodel_sm = signal.convolve2d(bgmodel, smkern)
#    pyfits.writeto('bgmodel_test.fits', bgmodel.astype(float), overwrite=True)

    # compute image with bg removed
    bgrem_data = image_data.astype(float) - bgmodel.astype(float)

    if debugfits:
        pyfits.writeto('bgrem_test.fits', bgrem_data.astype(float), overwrite=True)

    # find otsu threshold for each window
    ht, wd = bgrem_data.shape
    bg_window = 100

#    bgfact = 10.0
#    satur = 50000

    logging.info('computing star_model START')
    ts = time.time()
    star_model = np.zeros_like(bgrem_data)
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
    te = time.time()
    logging.info(f'computing star_model DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'star_model_bgfact{bgfact}.fits', star_model.astype(float), overwrite=True)

    # compute number if isolated pixels
    logging.info('neighbor computation START')
    ts = time.time()
    neigh_kern = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    neigh_map = signal.convolve2d(star_model, neigh_kern, mode='same')
    niso = (neigh_map == 1).sum()
    te = time.time()
    logging.info(f'neighbor computation DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'neigh_map.fits', neigh_map.astype(float), overwrite=True)

    # remove pixels with 1 neighbor or less
    logging.info('Removing stragglers START')
    ts = time.time()
    iso_mask = (neigh_map <= 1)
    star_model[iso_mask] = 0
    te = time.time()
    logging.info(f'Removing stragglers DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'star_model_bgfact{bgfact}_cleaned.fits', star_model.astype(float), overwrite=True)

    #print(bgfact, niso)

    # find structures
    logging.info('Finding structures START')
    ts = time.time()
    star_label, nstars = ndimage.measurements.label(star_model)
    te = time.time()
    logging.info(f'Finding structures DONE took {te-ts} seconds')

    logging.info(f'possible nstars = {nstars}')

    if debugfits:
        pyfits.writeto(f'star_label.fits', star_label.astype(float), overwrite=True)

    #print(star_label, star_label.dtype)

    # filter for stars with at least 9 pixels
    logging.info('Finding stars with >= 9 pix START')
    ts = time.time()
    star_obj = []
    star_boxes = np.zeros_like(star_model)
    for l in ndimage.find_objects(star_label):
        npix = (l[0].stop-l[0].start)*(l[1].stop-l[1].start)
        #print(bgrem_data[l].shape, npix)
        if npix < 9:
            continue

        if (bgrem_data[l] > satur).sum() > 0:
            continue

        star_obj.append(l)
        star_boxes[l] = 1
#        else:
#            star_label[l] = 0
    te = time.time()
    logging.info(f'Finding stars with >= 9 pix DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'star_boxes.fits', star_boxes.astype(float), overwrite=True)
        pyfits.writeto(f'star_label_cleaned.fits', star_label.astype(float), overwrite=True)

    logging.info(f'final nstars = {len(star_obj)}')

    # sort by brightness
    logging.info('Sort by brighness START')
    ts = time.time()
    lablist = []
    for l in star_obj:
        lablist.append((l, bgrem_data[l].sum()))
    lablist_sorted = sorted(lablist, key=lambda x: x[1], reverse=True)
    te = time.time()
    logging.info(f'Sort by brighness DONE took {te-ts} seconds')

    # keep only top 100
    #lablist_sorted = lablist_sorted[:max_stars]

    # now find centroid of each box
    logging.info('Find centroids START')
    ts = time.time()
    centroids = np.zeros_like(star_model)
    ni, nj = bgrem_data.shape
#    hfr_result = []
    #for l in ndimage.find_objects(star_label):
    star_cx = []
    star_cy = []
    star_r = []
    star_f = []
    nstars = 0
#    rexclude_pix = math.sqrt((ht**2+wd**2)/4)*(float(rexclude)/100)
#    logging.info(f'rexclude = {rexclude} pix = {rexclude_pix}')
#    rexclude_pix2 = rexclude_pix**2
    for l, s in lablist_sorted:
        data = bgrem_data[l]
        ci, cj = find_centroid(data, 0)

        # if exclusion zone exists filter
#        if rexclude > 0:
#            cr2 = (l[0].start+ci-ht/2)**2+(l[1].start+cj-wd/2)**2
#            print(ci, cj, ht/2, wd/2, math.sqrt(cr2), rexclude_pix, cr2, rexclude_pix2)
#
#            if cr2 < rexclude_pix2:
#                print('too close')
#                continue

        istart = max(0,  int(l[0].start+ci+0.5)-window)
        iend   = min(ni, int(l[0].start+ci+0.5)+window)
        jstart = max(0,  int(l[1].start+cj+0.5)-window)
        jend   = min(nj, int(l[1].start+cj+0.5)+window)

        crop_data = bgrem_data[istart:iend, jstart:jend]
        nci, ncj = find_centroid(crop_data, 0)

        #print(window, l, ci, cj, nci, ncj, istart, iend, jstart, jend)

        r, v = compute_radial_profile(crop_data, nci, ncj, np.median(crop_data))

        #print(r,v)

        ci = nci + istart
        cj = ncj + jstart

        hfr = find_hfr_from_radial_profile(r, v, window)

        if hfr is not None:
            rmid, iii, fluxtot, contrast = hfr

        if hfr is not None:
            #print(l, ci, cj, rmid, fluxtot, contrast)
            #hfr_result.append((ci, cj, rmid, fluxtot))
            star_cx.append(cj)
            star_cy.append(ci)
            star_r.append(rmid)
            star_f.append(fluxtot)
            nstars += 1

        centroids[int(ci), int(cj)] = 1
        if nstars == max_stars:
            break

    star_cx = np.array(star_cx)
    star_cy = np.array(star_cy)
    star_r = np.array(star_r)
    star_f = np.array(star_f)
    te = time.time()
    logging.info(f'Find centroids DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'centroids.fits', centroids.astype(float), overwrite=True)

    # need to figure out a bgest and noiseest to return, just 0 for now
    result = StarFitResult(star_cx, star_cy, star_r, star_f, len(star_cx),
                           0, 0, wd, ht)

    ttot_e = time.time()

    logging.info(f'Total star fitting took {ttot_e-ttot_s} seocnds')

    return result


def find_star(image_data, bgfact=50, satur=50000, window=100, debugfits=False):

    logging.info(f'find_stars start: bgfact={bgfact} window={window}')

    if debugfits:
        pyfits.writeto('thres_test.fits', image_data.astype(float), overwrite=True)

    logging.info('compute_bg_model START')
    ttot_s = time.time()
    ts = time.time()
    bgmodel = compute_bg_model(image_data, 100)
    te = time.time()
    logging.info(f'compute_bg_model DONE took {te-ts} seconds')

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
    te = time.time()
    logging.info(f'computing star_model DONE took {te-ts} seconds')

    if debugfits:
        pyfits.writeto(f'star_model_bgfact{bgfact}.fits', star_model.astype(float), overwrite=True)

    # use dilation to stregthen any structures
    ndilate = 3
    tmp_image = star_model
    logging.info(f'dilating star_model {ndilate} times')
    for i in range(0, ndilate):
        tmp_image = ndimage.grey_dilation(tmp_image, size=(3,3))

    dilate_star_model = tmp_image
    if debugfits:
        pyfits.writeto(f'dilate_tar_model_bgfact{bgfact}.fits', dilate_star_model.astype(float), overwrite=True)

    # find structures
    star_label, nstars = ndimage.measurements.label(dilate_star_model)

    print(star_label, nstars)
    if debugfits:
        pyfits.writeto(f'star_label.fits', star_label.astype(float), overwrite=True)

    logging.info('Finding stars with >= 9 pix START')
    star_boxes = np.zeros_like(star_model)
    max_pix = 0
    max_l = None
    for l in ndimage.find_objects(star_label):
        npix = (l[0].stop-l[0].start)*(l[1].stop-l[1].start)
        #print(bgrem_data[l].shape, npix)
        if npix < 9:
            continue

        if (bgrem_data[l] > satur).sum() > 0:
            continue

        if npix > max_pix:
            max_pix = npix
            max_l = l

    star_boxes[max_l] = 1

    if debugfits:
        pyfits.writeto(f'star_boxes.fits', star_boxes.astype(float), overwrite=True)

    cy, cx = ndimage.measurements.center_of_mass(star_boxes)
    logging.info(f'COM cx, cy = {cx}, {cy}')

    ttot_e = time.time()

    logging.info(f'find_star took {ttot_s-ttot_e} seconds')

    return cx, cy


# horizontal bin data to form a 1D star profile
# compute radial profile - returns (rad[], val[])
def horiz_bin_window(data, bg=0):
    #ts = time.time()

    ni, nj = data.shape

    profile = np.sum(data-bg, axis=0)

    print(data)
    print(profile)

    return profile

# find left/right limits of star disk
def find_star_limits(profile, thres=0):

    # search from left side (idx=0) for first pixel over thres
    idx = np.where(profile > thres)[0]
    print('where -> ', idx)
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
    def flux(proffunc, lx, rx):
        return quad(proffunc, lx, rx, epsabs=1, limit=100)[0]

    # different way
    def flux2(proffunc, rlo, rhi):
        nr = 2**7+1
        rvals, dr = np.linspace(rlo, rhi, num=nr, retstep=True)
        fvals = proffunc(rvals)
        return romb(fvals, dr)

    # compute center of profile with weighted average
    plen = profile.shape[0]

    print('profile len = ', plen)

    idx = np.arange(0, plen, 1)
    print(idx)

    cidx = np.sum(idx*profile)/np.sum(profile)
    print('center = ', cidx)

    lidx, ridx = find_star_limits(profile, thres)
    print('left, right = ', lidx, ridx)

    # find left, right limits
    # invert profile as y as a function of x on left and right
    # sides of the profile
    lidx_low = int(lidx-1)
    lidx_hi = int(lidx+2)
    ridx_low = int(ridx-1)
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

    print('lx=', lx, thres, proffunc(lx))
    print('rx=', rx, thres, proffunc(rx))

    # make sure inversion worked
    INVERSE_FRAC = 0.05
    if abs(thres-proffunc(lx))/thres > INVERSE_FRAC or abs(thres-proffunc(rx))/thres > INVERSE_FRAC:
        print('inversion failed!')
        return None, None, None

    # now find flux inside left/right
    simple_totflux = np.sum([profile[lidx:ridx]])
    print('flux in star (simple) = ', simple_totflux)

    totflux = flux2(proffunc, lx, rx)
    #print('integrated flux in star = ', flux(proffunc, lidx, ridx))
    print('integrated flux in star = ', flux)

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

    # make inverse function of r vs flux to find half flux
    flux_inv = interp1d(flux2_vs_r, r_arr)

    half_flux_r = flux_inv(totflux/2)

    logging.info(f'half flux rad = {half_flux_r}')
    print(r_arr)
    #print(flux_vs_r)
    print(flux2_vs_r)
    #ax_3.plot(r_arr, flux_vs_r)
    if debugplots:
        ax_3.plot(r_arr, flux2_vs_r)
        ax_3.axvline(half_flux_r, color='green')

    return cidx, lx, rx, cidx-half_flux_r, cidx+half_flux_r

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

    scen, sl, sr = find_hfr_from_1D(profile, thres=thres)

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



