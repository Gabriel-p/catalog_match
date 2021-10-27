
from os.path import exists
from os import makedirs
from os.path import join
from pathlib import Path

from astropy.coordinates import Angle
from astropy import units as u
from astropy.io import ascii

# from astroquery.xmatch import XMatch
import logging
import time

from modules import read_input
from modules import cat_query
from modules import match_cats
from modules import make_plots
from modules import out_data


def main():
    """
    Catalog cross-match.
    """

    data_cols, ra_hour, cat_mode, m_qry, ra_qry, dec_qry, box_s, max_arcsec,\
        max_mag_delta = read_input.readIniFile()

    # Process all files inside 'input/' folder.
    cl_files = read_input.get_files()
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
                    clust_file, data_cols, ra_hour)
        except ascii.core.InconsistentTableError as err:
            logging.info("{}\n\nERROR: could not read data file {}".format(
                err, clust_file))
            keep_going = False

        if keep_going is True:
            # Query catalog.
            query, catalog = cat_query.get(
                clust_name, N_obs, ra_mid, dec_mid, ra_rang, dec_rang, box_s,
                cat_mode)

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

            logging.info('\nCreating output files.')
            out_data.main(
                clust_name, inp_data, cat_mode, ra_qry, dec_qry, query,
                match_c1_ids_all, match_c2_ids_all, match_d2d_all,
                no_match_c1_all, no_match_d2d_all)

            logging.info("\nProcessed {}".format(clust_name))

    logging.info("\nEnd.")


if __name__ == '__main__':
    # # To see available catalogs:
    # from astroquery.vizier import Vizier
    # # catalog_list = Vizier.find_catalogs('Pan-STARRS')
    # catalog_list = Vizier.find_catalogs('ALLWISE')
    # catalogs = Vizier.get_catalogs(catalog_list.keys())
    # print(catalogs)
    # Generate output dir if it doesn't exist.
    if not exists('output'):
        makedirs('output')
    main()
