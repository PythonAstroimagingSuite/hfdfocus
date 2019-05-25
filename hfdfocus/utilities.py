import numpy as np

def find_stars_near_target(target, dist, minmag, maxmag, exclude=None):

    t_ra_rad = target.ra.radian
    t_dec_rad = target.dec.radian

    #print(t_ra_rad, t_dec_rad)

    # convert sao coords to radians
    c_ra_rad = np.deg2rad(saocat.ra)
    c_dec_rad = np.deg2rad(saocat.dec)

    # equation for distance between two points on sphere
    # try with haversine
    ddec = t_dec_rad - c_dec_rad
    dra = t_ra_rad - c_ra_rad
    a = (np.sin(ddec/2)**2 + np.cos(t_dec_rad)*np.cos(c_dec_rad)*np.sin(dra/2)**2)
    c = 2 * np.arcsin(np.sqrt(a))

    #print(np.max(c), np.min(c), np.median(c))

    #sortargs = np.argsort(c)
    close_idx = np.where(c < np.deg2rad(dist))[0]

    if exclude is not None:
        #logging.info(f'exclude={exclude}')
        #logging.info(f'pre-exclude: {close_idx}')
        close_idx = np.setxor1d(close_idx, exclude)
        logging.debug(f'post-exclude: {close_idx}')

    logging.debug(f'found {len(close_idx)} within {args.dist} degrees')

    # filter by mag
    mags = np.array(saocat.vmag)[close_idx]
    mags_idx = np.where((mags < minmag) & (mags >= maxmag))[0]

    logging.debug(f'found {len(mags_idx)} within {dist} degrees and {maxmag} < mag < {minmag}')
    #logging.debug(f'mags_idx = {mags_idx}')
    ret_idx = close_idx[mags_idx]
    #logging.debug(f'ret_idx = {ret_idx}')
    return ret_idx, c[ret_idx]
