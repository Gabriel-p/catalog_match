
import numpy as np
from astroquery.vizier import Vizier
# from astroquery.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
import logging


def in_data(clust_file, data_mode, data_cols, ra_hour):
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

    if data_mode == 'num':
        m_idx, ra_idx, dec_idx = map(int, data_cols)
        m_obs, ra_obs, dec_obs = inp_data.columns[m_idx],\
            inp_data.columns[ra_idx], inp_data.columns[dec_idx]
        m_obs_nam = None
    elif data_mode == 'nam':
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


def cat_query(
    clust_name, N_obs, ra_mid, dec_mid, ra_rang, dec_rang, box_s, cat_mode,
        catalog):
    """
    Query selected catalog or read from file a previous stored version.
    """
    if cat_mode == 'query':
        logging.info("\nFetching data from {} catalog.".format(catalog))
        txt = 'queried'
        cent = SkyCoord(
            ra=ra_mid * u.degree, dec=dec_mid * u.degree, frame='icrs')

        if str(box_s) == 'auto':
            width = max(ra_rang, dec_rang) * u.deg
            logging.info("Using auto width={:.3f}".format(width))
        else:
            width = box_s * u.deg
            logging.info("Using manual width={:.3f}".format(width))

        # Vizier query
        # Unlimited rows, all columns
        v = Vizier(row_limit=-1, columns=['all'])
        query = v.query_region(SkyCoord(
            cent, frame='icrs'), width=width, catalog=[catalog])

        try:
            query = query[0]
        except IndexError:
            raise ValueError("Queried catalog came back empty")

        # If the last row is empty (not sure why it happens), remove it
        if 'END' in query[-1]:
            print("Removing last empty row from queried data")
            query = query[:-1]

        # # Irsa query
        # # Set maximum limit for retrieved stars.
        # Irsa.ROW_LIMIT = int(4 * N_obs)
        # query = Irsa.query_region(
        #     cent, catalog=catalog, spatial='Box', width=ra_rang * u.deg)

        out_file = clust_name + '_query.dat'
        logging.info("Writing queried data to '{}' file.".format(out_file))
        # TODO waiting for fix to
        # https://github.com/astropy/astropy/issues/7744
        ascii.write(
            query, 'output/' + out_file, overwrite=True)  # format='csv',

        catalog = "Queried catalog {}".format(catalog)

    elif cat_mode == 'read':
        logging.info("\nReading input catalog")
        txt = 'read'
        q_file = 'input/' + clust_name + '_query.dat'
        query = ascii.read(q_file)
        catalog = "Read catalog {}".format(catalog)

    logging.info("N ({} catalog): {}".format(txt, len(query)))

    return query, catalog
