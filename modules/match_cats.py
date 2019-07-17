
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u
import warnings
import logging
from . import LinRegressionCI


def cat_match(ra_obs, dec_obs, ra_qry, dec_qry):
    """
    Match catalogs.
    idx are indices into c2 that are the closest objects to each of the
    coordinates in c1
    d2d are the on-sky distances between them,
    d3d are the 3-dimensional distances
    """
    # Catch 'RuntimeWarning' issued because we are inserting 'nan' values
    # into the queried coordinates.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Define observed and queried catalogs.
        c1 = SkyCoord(ra_obs, dec_obs, unit=(u.degree, u.degree))
        c2 = SkyCoord(ra_qry, dec_qry, unit=(u.degree, u.degree))
        idx, d2d, d3d = c1.match_to_catalog_sky(c2)

    return idx, d2d.deg


def match_filter(c1_ids, c2_ids, d2d, max_deg):
    """
    Reject duplicated matches and matched stars with distances beyond
    the maximum separation defined.
    """
    match_c1_ids, dupl_c1_ids, no_match_c1 = [], [], []
    match_c2_ids, match_d2d, no_match_d2d = [], [], []
    # Indexes of matched star in both catalogs and distance between them.
    for c1_i, c2_i, d in zip(*[c1_ids, c2_ids, d2d]):
        # Filter by maximum allowed match distance.
        if d <= max_deg:
            # Check if this queried star in c2 was already stored.
            if c2_i in match_c2_ids:
                # Index of this stored c2 star.
                j = match_c2_ids.index(c2_i)
                # Only replace with this match if the previous match had a
                # larger distance.
                if match_d2d[j] > d:
                    # Store rejected replaced star.
                    dupl_c1_ids.append(match_c1_ids[j])
                    # Replace star
                    match_c1_ids[j] = c1_i
                    match_d2d[j] = d
                else:
                    dupl_c1_ids.append(c1_i)
            else:
                match_c1_ids.append(c1_i)
                match_c2_ids.append(c2_i)
                match_d2d.append(d)
        else:
            # If this observed star has no queried star closer than the max
            # distance allowed, discard.
            no_match_c1.append(c1_i)
            no_match_d2d.append(d)

    return match_c1_ids, match_c2_ids, no_match_c1, np.array(dupl_c1_ids),\
        match_d2d, no_match_d2d


def main(
    ra_qry, dec_qry, max_arcsec, ra_obs, dec_obs, query, m_obs, m_qry,
        max_mag_delta):
    """
    Catalog matching module.
    """
    logging.info("Maximum match distance: {:.1f} [arcsec]".format(max_arcsec))
    # Maximum separation in arcsec, convert to decimal degrees.
    max_deg = Angle(max_arcsec * u.arcsec).deg

    # Initial full list of observed and queried catalogs.
    c1_ids = np.arange(len(ra_obs))
    c2_ids = np.arange(len(query))
    ra_q, dec_q = query[ra_qry][c2_ids], query[dec_qry][c2_ids]

    # Values used to replace already matched stars. Use values far away from
    # the frame's limits.
    ra_max = max(ra_obs.max(), ra_q.max())
    dec_max = max(dec_obs.max(), dec_q.max())
    ra_offset = ra_max + .5 * ra_max
    dec_offset = dec_max + .5 * dec_max

    # Store all unique matches, and observed stars with no match.
    match_c1_ids_all, match_c2_ids_all, no_match_c1_all, match_d2d_all,\
        no_match_d2d_all = [], [], [], [], []
    # Continue until no more duplicate matches exist.
    j = 0
    while c1_ids.any():
        j += 1

        # Match observed and queried catalogs.
        c2_ids_dup, c1c2_d2d = cat_match(
            ra_obs[c1_ids], dec_obs[c1_ids], ra_q, dec_q)

        logging.info('')
        logging.info("{}. Matched catalogs".format(j))

        # Return unique ids for matched stars between catalogs, ids of
        # observed stars with no match found, and ids of observed stars
        # with a duplicated match that will be re-processed ('c1_ids').
        match_c1_ids, match_c2_ids, no_match_c1, c1_ids, match_d2d,\
            no_match_d2d = match_filter(c1_ids, c2_ids_dup, c1c2_d2d, max_deg)
        logging.info("Unique stars filtered")

        # Store all unique solutions and no match solutions.
        match_c1_ids_all += match_c1_ids
        match_c2_ids_all += match_c2_ids
        no_match_c1_all += no_match_c1
        match_d2d_all += match_d2d
        no_match_d2d_all += no_match_d2d

        logging.info("Unique matches: {}".format(len(match_c1_ids)))
        if match_d2d_all:
            logging.info("Average match separation: {:.2f} arcsec".format(
                np.mean(Angle(match_d2d_all * u.deg).arcsec)))
        logging.info("No match: {}".format(len(no_match_c1)))
        logging.info("Re match: {}".format(len(c1_ids)))

        if c1_ids.any():
            # Queried stars that were not matched to an observed star.
            c2_ids_r = np.setdiff1d(c2_ids, match_c2_ids_all)
            logging.info("Queried stars for re-match: {}".format(
                len(c2_ids_r)))
            # To avoid messing with the indexes, change the coordinates
            # of already matched queried stars so that they can not
            # possibly be matched again.
            ra_q[match_c2_ids_all] = ra_offset
            dec_q[match_c2_ids_all] = dec_offset

    logging.info('')
    logging.info('Observed stars matched: {}'.format(len(match_c1_ids_all)))
    logging.info('Observed stars not matched: {}'.format(len(no_match_c1_all)))

    # Magnitude filter to reject 'magnitude outliers'.
    logging.info("Magnitude delta filter (<{})".format(max_mag_delta))
    mag_filter_data = LinRegressionCI.main(
        m_obs[match_c1_ids_all], query[m_qry][match_c2_ids_all], max_mag_delta)

    mag_msk = np.array(mag_filter_data[0])
    # Update not-matched arrays
    no_match_c1_all += list(np.array(match_c1_ids_all)[~mag_msk])
    no_match_d2d_all += list(np.array(match_d2d_all)[~mag_msk])
    # Update matched arrays
    match_c1_ids_all = np.array(match_c1_ids_all)[mag_msk]
    match_d2d_all = np.array(match_d2d_all)[mag_msk]
    match_c2_ids_all = np.array(match_c2_ids_all)[mag_msk]

    logging.info(
        'Observed stars matched after mag filter: {}'.format(
            len(match_c1_ids_all)))
    logging.info(
        'Observed stars not matched after mag filter: {}'.format(
            len(no_match_c1_all)))

    # Rejected magnitudes from queried catalog.
    q_rjct_mks = np.ones(len(query), dtype=bool)
    q_rjct_mks[match_c2_ids_all] = False
    logging.info('Queried stars not matched: {}'.format(sum(q_rjct_mks)))

    return match_c1_ids_all, no_match_c1_all, match_d2d_all, no_match_d2d_all,\
        match_c2_ids_all, q_rjct_mks, mag_filter_data
