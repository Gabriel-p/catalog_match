
from os.path import exists
from os import makedirs
from os import listdir
from os.path import isfile, join
from pathlib import Path

import numpy as np
from astropy.coordinates import Angle
from astropy import units as u
from astropy.io import ascii

# from astroquery.xmatch import XMatch
import logging
import time

from modules import read_input
from modules import match_cats
from modules import make_plots
from modules import out_data


def main():
    """
    Observed catalog matcher.
    """

    data_mode, data_cols, ra_hour, cat_mode, catalog_n, m_qry, ra_qry,\
        dec_qry, box_s, max_arcsec, max_mag_delta, out_fig = params_input()

    # Generate output dir if it doesn't exist.
    if not exists('output'):
        makedirs('output')

    # Process all files inside 'input/' folder.
    cl_files = get_files()
    if not cl_files:
        print("No input cluster files found")

    mypath = Path().absolute()
    for clust_file in cl_files:
        clust_name = clust_file.split('/')[1].split('.')[0]

        # Set up logging module
        level = logging.INFO
        frmt = ' %(message)s'
        handlers = [
            logging.FileHandler(
                join(mypath, 'output', clust_name + '.log'), mode='w'),
            logging.StreamHandler()]
        logging.basicConfig(level=level, format=frmt, handlers=handlers)

        logging.info(time.strftime("\n%Y-%m-%d, %H:%M"))
        logging.info("\nProcessing: {}\n".format(clust_name))
        # Get input data from file.
        keep_going = True
        try:
            inp_data, m_obs, ra_obs, dec_obs, N_obs, ra_mid, dec_mid, ra_rang,\
                dec_rang, m_obs_nam = read_input.in_data(
                    clust_file, data_mode, data_cols, ra_hour)
        except ascii.core.InconsistentTableError as err:
            logging.info("{}\n\nERROR: could not read data file {}".format(
                err, clust_file))
            keep_going = False

        if keep_going is True:
            # Query catalog.
            query, catalog = read_input.cat_query(
                clust_name, N_obs, ra_mid, dec_mid, ra_rang, dec_rang, box_s,
                cat_mode, catalog_n)

            # Match catalogs.
            logging.info("\nMatching catalogs...")
            match_c1_ids_all, no_match_c1_all, match_d2d_all,\
                no_match_d2d_all, match_c2_ids_all, q_rjct_mks,\
                mag_filter_data = match_cats.main(
                    ra_qry, dec_qry, max_arcsec, ra_obs, dec_obs, query, m_obs,
                    m_qry, max_mag_delta)

            # Store match and no match separations as 'Angle' objects.
            match_d2d_all, no_match_d2d_all = Angle(match_d2d_all * u.degree),\
                Angle(no_match_d2d_all * u.degree)

            if out_fig == 'y':

                # Unique stars from observed catalog.
                m_unq, ra_unq, dec_unq =\
                    m_obs[match_c1_ids_all], ra_obs[match_c1_ids_all],\
                    dec_obs[match_c1_ids_all]

                if m_qry is not None:
                    query_mag = query[m_qry]
                    # Unique/rejected magnitudes from queried catalog.
                    m_unq_q = query_mag[match_c2_ids_all]
                    # Rejected magnitudes from queried catalog.
                    m_rjct_q = query_mag[q_rjct_mks]
                else:
                    query_mag = [1. for _ in range(len(query))]
                    m_unq_q, m_rjct_q = [], [_ for _ in range(sum(q_rjct_mks))]

                # Differences in matched coordinates.
                ra_unq_delta, dec_unq_delta =\
                    Angle(ra_unq - query[ra_qry][match_c2_ids_all],
                          unit='deg').arcsec,\
                    Angle(dec_unq - query[dec_qry][match_c2_ids_all],
                          unit='deg').arcsec
                # Observed stars with no match.
                m_rjct, ra_rjct, dec_rjct = m_obs[no_match_c1_all],\
                    ra_obs[no_match_c1_all], dec_obs[no_match_c1_all]

                if match_c1_ids_all.any():
                    logging.info('\nCreating output plots.')
                    make_plots.main(
                        clust_name, m_qry, catalog, max_arcsec, m_obs,
                        m_obs_nam, ra_obs, dec_obs, query_mag,
                        query[ra_qry], query[dec_qry], m_unq, ra_unq, dec_unq,
                        m_unq_q, m_rjct_q, match_d2d_all, no_match_d2d_all,
                        ra_unq_delta, dec_unq_delta, m_rjct, ra_rjct, dec_rjct,
                        max_mag_delta, mag_filter_data)
                    logging.info("Output figure created.")
                else:
                    logging.info("  ERROR: no matches to plot.")
            else:
                logging.info("No output figure created.")

            logging.info('\nCreating output files.')
            out_data.main(
                clust_name, inp_data, cat_mode, ra_qry, dec_qry, query,
                match_c1_ids_all, match_c2_ids_all, match_d2d_all,
                no_match_c1_all, no_match_d2d_all)

            logging.info("\nProcessed {}".format(clust_name))

    logging.info("\nEnd.")


def params_input():
    """
    Read input parameters from 'params_input.dat' file.
    """
    nn_lst = ('n', 'N', 'none', 'None', 'NONE')
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                if reader[0] == 'DC':
                    data_mode = reader[1]
                    data_cols = [_.replace('\_', ' ') for _ in reader[2:-1]]
                    ra_hour = reader[-1]
                    ra_hour = True if ra_hour in ('y', 'Y', 'True') else False
                if reader[0] == 'CA':
                    cat_mode = reader[1]
                    catalog_n = reader[2]
                    m_qry = None if reader[3] in nn_lst else reader[3]
                    ra_qry = reader[4]
                    de_qry = reader[5]
                    box_s = reader[6] if reader[6] == 'auto' else\
                        float(reader[6])
                if reader[0] == 'MA':
                    max_arcsec = float(reader[1])
                    max_mag_delta = float(reader[2])
                if reader[0] == 'FI':
                    out_fig = reader[1]

    return data_mode, data_cols, ra_hour, cat_mode, catalog_n, m_qry, ra_qry,\
        de_qry, box_s, max_arcsec, max_mag_delta, out_fig


def get_files():
    '''
    Store the paths and names of all the input clusters stored in the
    input folder.
    '''
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


if __name__ == '__main__':
    # To see available catalogs:
    if False:
        from astroquery.vizier import Vizier
        # catalog_list = Vizier.find_catalogs('Pan-STARRS')
        catalog_list = Vizier.find_catalogs('ALLWISE')
        catalogs = Vizier.get_catalogs(catalog_list.keys())
        print(catalogs)
    main()
