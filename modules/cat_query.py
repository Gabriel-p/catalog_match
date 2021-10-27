
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
import logging


def get(
    clust_name, N_obs, ra_mid, dec_mid, ra_rang, dec_rang, box_s,
        cat_mode):
    """
    Query selected catalog or read from file a previous stored version.
    """
    if cat_mode != 'read':
        logging.info("\nFetching data from {} catalog.".format(cat_mode))
        txt = 'queried'
        cent = SkyCoord(
            ra=ra_mid * u.degree, dec=dec_mid * u.degree, frame='icrs')

        if str(box_s) == 'auto':
            width = max(ra_rang, dec_rang)
            width = width + width * .05
            width = width * u.deg
            logging.info("Using auto width={:.3f}".format(width))
        else:
            width = box_s * u.deg
            logging.info("Using manual width={:.3f}".format(width))

        # Vizier query
        # Unlimited rows, all columns
        v = Vizier(row_limit=-1, columns=['all'])
        query = v.query_region(SkyCoord(
            cent, frame='icrs'), width=width, catalog=[cat_mode])

        try:
            query = query[0]
        except IndexError:
            raise ValueError("Queried catalog came back empty")

        # If the last row is empty (not sure why this happens), remove it
        if 'END' in query[-1] or "EN" in query[-1]:
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

        catalog = "Queried catalog {}".format(cat_mode)

    elif cat_mode == 'read':
        logging.info("\nReading input catalog")
        txt = 'read'
        q_file = 'input/' + clust_name + '_query.dat'
        query = ascii.read(q_file)
        catalog = "Read catalog"

    logging.info("N ({} catalog): {}".format(txt, len(query)))

    return query, catalog
