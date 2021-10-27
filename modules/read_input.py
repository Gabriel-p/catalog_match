
import os.path
from os import listdir
from os.path import isfile, join
import numpy as np
import configparser
from astropy import units as u
from astropy.io import ascii
import logging


def readIniFile():
    """
    Read .ini config file
    """
    in_params = configparser.ConfigParser()

    # Local file takes precedence
    if os.path.isfile('params_local.ini'):
        in_params.read('params_local.ini')
    else:
        in_params.read('params.ini')

    data_in = in_params['Data input']
    mag, ra, dec = data_in.get('m_in'), data_in.get('ra'), data_in.get('dec')
    data_cols = [mag, ra, dec]
    ra_hour = data_in.getboolean('ra_hour')

    q_mode = in_params['Query mode']
    cat_mode, m_qry, ra_qry, dec_qry, box_s = q_mode.get('mode'),\
        q_mode.get('m_q'), q_mode.get('ra_q'), q_mode.get('dec_q'),\
        q_mode.get('box_s')

    m_filts = in_params['Match filters']
    max_arcsec, max_mag_delta = m_filts.getfloat('max_arcsec'),\
        m_filts.getfloat('max_mag_delta')

    return data_cols, ra_hour, cat_mode, m_qry, ra_qry, dec_qry, box_s,\
        max_arcsec, max_mag_delta


def get_files():
    """
    Store the paths and names of all the input clusters stored in the
    input folder.
    """
    cl_files = [join('input/', f) for f in listdir('input/') if
                isfile(join('input/', f))]

    # Remove readme file it is still there.
    try:
        cl_files.remove('input/README.md')
    except ValueError:
        pass
    # Remove any '_query.dat' file.
    rm_idx = []
    for i, _ in enumerate(cl_files):
        if _.endswith('_query.dat'):
            rm_idx.append(i)
    cl_files = np.delete(cl_files, rm_idx).tolist()

    return cl_files


def in_data(clust_file, data_cols, ra_hour):
    """
    Read the file that stores the photometric data for all stars.
    """
    # Loads the data in 'clust_file' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    # try:
    inp_data = ascii.read(
        clust_file,  # format='commented_header',
        fill_values=[('', '0'), ('NA', '0'), ('INDEF', '0')])
    # except ascii.core.InconsistentTableError:

    m_obs_nam, ra_nam, dec_nam = data_cols
    m_obs, ra_obs, dec_obs = inp_data[m_obs_nam], inp_data[ra_nam],\
        inp_data[dec_nam]

    # Convert RA coordinates from hours to degrees
    if ra_hour is True:
        ra_obs = ra_obs * u.hourangle.to(u.deg)

    N_obs = len(ra_obs)
    logging.info("N (all stars) = {}".format(N_obs))

    ra_mid, dec_mid = .5 * (min(ra_obs) + max(ra_obs)),\
        .5 * (max(dec_obs) + min(dec_obs))
    ra_rang = (max(ra_obs) - min(ra_obs)) * np.cos(np.deg2rad(dec_mid))
    dec_rang = max(dec_obs) - min(dec_obs)

    logging.info('ra (min, max): {:.4f}, {:.4f}'.format(
        min(ra_obs), max(ra_obs)))
    logging.info('dec (min, max): {:.4f}, {:.4f}'.format(
        min(dec_obs), max(dec_obs)))
    logging.info('Range (ra, dec): {:.4f}, {:.4f}'.format(ra_rang, dec_rang))
    logging.info('Centre (ra, dec): {:.4f}, {:.4f}'.format(ra_mid, dec_mid))

    return inp_data, m_obs, ra_obs, dec_obs, N_obs, ra_mid, dec_mid, ra_rang,\
        dec_rang, m_obs_nam
