#
# Determine HFR of stars using radial profie
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

import math
import time
import logging
import argparse
import numpy as np
from scipy import signal, ndimage, interpolate

import astropy.io.fits as pyfits

# for testing
import matplotlib as mpl
mpl.use("TkAgg")
mpl.rcParams['toolbar'] = 'None'
#mpl.rc('font', size=8)
import matplotlib.pyplot as plt

from hfdfocus.MultipleStarFitHFD import StarFitResult, find_centroid, measure_stars_in_image

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
    rfit = interpolate.interp1d(rad, fluxrad)

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
    proffit = interpolate.interp1d(rad, val)
    #contrast = val[1]/proffit(rmid)
    # FIXME contrast doesnt work for defocused stars and we're getting really small denominators (or 0!)
    contrast = 1

    # return position relative to corner of crop window
    if not extras:
        return (rmid, 0, fluxtot, contrast)
    else:
        return ((rmid, 0, fluxtot, contrast), (rad, val, fluxrad, fluxtot))


def find_star_HFR_Radial_Profle(star, image_data, thres=10000, win=100,
                                debugplots=False, debugfits=False):
    # yuck for testing!
   # global fig, ax_2d, ax_1d

    #xcen, ycen, bg, noise, starmask, alone = find_star(image_data, debugfits=debugfits)

    xcen = star.cx
    ycen = star.cy
    bg = star.bglevel
    #mad = star.bgmad
    estsize = star.estsize
    logging.debug(f'cx: {xcen} cy: {ycen} - using estsize of {estsize} for profile window')

    img_ht, img_wd = image_data.shape

    boxsize = 15

    xlow = max(0, int(xcen-boxsize))
    xhi = min(img_wd-1, int(xcen+boxsize))
    ylow = max(0, int(ycen-boxsize))
    yhi = min(img_ht-1, int(ycen+boxsize))
    logging.debug(f'xlow:xhi={xlow}:{xhi}  ylow:yhi={ylow}:{yhi}')
    crop_data = image_data[ylow:yhi, xlow:xhi]

    if debugplots:
        plt.xticks(fontsize=7)
        fig = plt.figure(figsize=(4,3), dpi=200)
        #fig.subplots_adjust(wspace=1.5)
        ax_2d = fig.add_subplot(131)
        ax_2d.tick_params(axis='both', which='major', labelsize=7)
        ax_2d.tick_params(axis='both', which='minor', labelsize=7)
        ax_1da = fig.add_subplot(132)
        ax_1da.tick_params(axis='both', which='major', labelsize=7)
        ax_1da.tick_params(axis='both', which='minor', labelsize=7)
        ax_1db = fig.add_subplot(133)
        ax_1db.tick_params(axis='both', which='major', labelsize=7)
        ax_1db.tick_params(axis='both', which='minor', labelsize=7)
        im = ax_2d.imshow((crop_data-bg).astype(float))
        cb = fig.colorbar(im, ax=ax_2d)
        cb.ax.tick_params(labelsize=7)

    nci, ncj = find_centroid(crop_data, 0)

#        print(window, l, ci, cj, nci, ncj, istart, iend, jstart, jend)

    r, v = compute_radial_profile(crop_data, nci, ncj, np.median(crop_data))

        #print(r,v)
#        for rp,vp in zip(r,v):
#            print(rp, vp)

    if debugplots:
        ax_1da.plot(r, v)

    # guess window from estsize
    window = 2*estsize

    hfr = find_hfr_from_radial_profile(r, v, window, extras=True)

    if hfr is not None:
        #rmid, iii, fluxtot, contrast = hfr

        rmid, iii, fluxtot, contrast = hfr[0]
        rad, val, fluxrad, fluxtot = hfr[1]

        if debugplots:
            ax_1da.axvline(rmid, color='green')
            ax_1db.plot(rad, fluxrad)
            ax_1db.axvline(rmid, color='green')
            plt.pause(0.1)

            fig.suptitle(f'idx={star.idx} cx={star.cx:5.3f} cy={star.cy:5.3f} '
                         f'R1 {rmid:5.2f} R2 {rmid:5.2f} A 0')

        return rmid, rmid, 0, fluxtot
    else:
        return None

def star_fit_hfr_radial_profile(image_data, max_stars=100, bgfact=50, satur=50000, window=7,
                                debugplots=False, debugfits=False):

    detected = measure_stars_in_image(image_data, find_star_HFR_Radial_Profle,
                                      max_stars=max_stars, window=window,
                                      bgfact=bgfact,
                                      debugplots=debugplots,
                                      debugfits=debugfits)

    if detected is not None:
        star_cx = np.array([x.cx for x in detected.stars])
        star_cy = np.array([x.cy for x in detected.stars])
        star_r1 = np.array([x.r1 for x in detected.stars])
        star_r2 = np.array([x.r2 for x in detected.stars])
        star_angle = np.array([x.angle for x in detected.stars])
        star_f = np.array([x.flux for x in detected.stars])
        print(star_cx)
        ht, wd = image_data.shape
        result = StarFitResult(star_cx, star_cy, star_r1, star_r2, star_angle,
                               star_f, len(star_cx), 0, 0, wd, ht)
    else:
        result = None

    return result