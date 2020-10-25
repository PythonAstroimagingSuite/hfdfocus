# average vcurve runs
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
import argparse
import numpy as np


# for testing
import matplotlib as mpl
mpl.rc('font', size=8)
import matplotlib.pyplot as plt

class VCurve:
    def __init__(self):
        pass

    def __repr__(self):
        return f'Time={self.t} rs={self.rs} rp={self.rp} ls={self.ls} lp={self.lp}'

if __name__ == '__main__':
    logging.basicConfig(filename='average_vcurve_runs.log',
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
    parser.add_argument('vcurve_file', type=str, help='V Curve filename')
    parser.add_argument('--debugplots', action='store_true', help='show debug plots')

    args = parser.parse_args()

    f = open(args.vcurve_file, 'r')

    vcurves = []
    tstamps = []
    rs = []
    rp = []
    ls = []
    lp = []
    for l in f.readlines():
#        v = VCurve()
#        v.t, v.ls, v.lp, v.rs, v.rp = l.strip().split(',')
#        v.t = time.strptime(v.t, '%Y/%m/%d %H:%M:%S %Z')
#        for k in ['rs', 'rp', 'ls', 'lp']:
#            v.__dict__[k] = float(v.__dict__[k])
#        vcurves.append(v)
        v = json.loads(l)
        vcurves.append(v)
    f.close()

    for v in vcurves:
        #print('v:', v)
        tstamps.append(time.strptime(v['timestamp'], '%Y/%m/%d %H:%M:%S %Z'))
        rs.append(float(v['rightslope']))
        rp.append(float(v['rightpid']))
        ls.append(float(v['leftslope']))
        lp.append(float(v['leftpid']))

    ls = np.array(ls)
    lp = np.array(lp)
    rs = np.array(rs)
    rp = np.array(rp)

    lfit = np.polyfit(ls, lp, 1)
    rfit = np.polyfit(rs, rp, 1)
    ls_avg = np.average(ls)
    rs_avg = np.average(rs)
    lp_avg = np.poly1d(lfit)(ls_avg)
    rp_avg = np.poly1d(rfit)(rs_avg)

    j = json.dumps({'leftslope' : ls_avg, 'leftpid' : lp_avg,
                    'rightslope' : rs_avg, 'rightpid' : rp_avg})
    #sys.stdout.write('\n\nAverage Fit Params\n')
    sys.stdout.write(j + '\n')

    # format for YAML
    sys.stdout.write('\nYAML\n')
    sys.stdout.write(f'  vcurve_ls: {ls_avg}\n')
    sys.stdout.write(f'  vcurve_lp: {lp_avg}\n')
    sys.stdout.write(f'  vcurve_rs: {rs_avg}\n')
    sys.stdout.write(f'  vcurve_rp: {rp_avg}\n')

    if args.debugplots:
        fig = plt.figure()
        ax_l = fig.add_subplot(221)
        ax_r = fig.add_subplot(222)
        ax_s = fig.add_subplot(223)
        ax_p = fig.add_subplot(224)

        ax_l.set_xlabel('ls')
        ax_l.set_ylabel('lp')
        ax_r.set_xlabel('rs')
        ax_r.set_ylabel('rp')
        ax_s.set_xlabel('ls')
        ax_s.set_ylabel('rs')
        ax_p.set_xlabel('lp')
        ax_p.set_ylabel('rpq')

        ax_l.plot(ls, lp, marker='+', ls='')
        ax_r.plot(rs, rp, marker='+', ls='')
        ax_s.plot(ls, rs, marker='.', ls='')
        ax_p.plot(lp, rp, marker='.', ls='')

        ax_l.plot(np.unique(ls), np.poly1d(lfit)(np.unique(ls)))
        ax_r.plot(np.unique(rs), np.poly1d(rfit)(np.unique(rs)))

        ax_l.plot([ls_avg], [lp_avg], marker='*', color='red')
        ax_r.plot([rs_avg], [rp_avg], marker='*', color='red')
        fig.suptitle(f'ls={ls_avg} lp={lp_avg:5.2f} rs={rs_avg} rs={rp_avg:5.2f}')

        fig.show()
        plt.show()

