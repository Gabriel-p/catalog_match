
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u


def cat_match(ra_obs, dec_obs, ra_qry, dec_qry):
    """
    Match catalogs.
    idx are indices into c2 that are the closest objects to each of the
    coordinates in c1
    d2d are the on-sky distances between them,
    d3d are the 3-dimensional distances
    """
    # Define observed and queried catalogs.
    c1 = SkyCoord(ra_obs, dec_obs, unit=(u.degree, u.degree))
    c2 = SkyCoord(ra_qry, dec_qry, unit=(u.degree, u.degree))

    idx, d2d, d3d = c1.match_to_catalog_sky(c2)
    # print('Unique matches:', len(set(idx)))
    # print('Filtered matches:', np.sum([d2d.arcsec < max_arcsec]))

    return idx, d2d


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
        if d.deg <= max_deg:
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


def main(ra_qry, de_qry, max_arcsec, ra_obs, dec_obs, query):
    """
    Catalog matching module.
    """

    print("\nMaximum match distance: {:.1f} [arcsec]".format(max_arcsec))
    # Maximum separation in arcsec, convert to decimal degrees.
    max_deg = Angle(max_arcsec * u.arcsec).deg

    # Initial full list of observed and queried catalogs.
    c1_ids = np.arange(len(ra_obs))
    c2_ids = np.arange(len(query))
    ra_q, dec_q = query[ra_qry][c2_ids], query[de_qry][c2_ids]
    # Store all unique matches, and observed stars with no match.
    match_c1_ids_all, match_c2_ids_all, no_match_c1_all, match_d2d_all,\
        no_match_d2d_all = [], [], [], [], []
    # Continue until no more duplicate matches exist.
    print("\nMatching catalogs.")
    j = 0
    while c1_ids.any():

        print("\n  {}.".format(j + 1))
        j += 1

        # Match observed and queried catalogs.
        c2_ids_dup, c1c2_d2d = cat_match(
            ra_obs[c1_ids], dec_obs[c1_ids], ra_q, dec_q)

        # Return unique ids for matched stars between catalogs, ids of
        # observed stars with no match found, and ids of observed stars
        # with a duplicated match that will be re-processed ('c1_ids').
        match_c1_ids, match_c2_ids, no_match_c1, c1_ids, match_d2d,\
            no_match_d2d = match_filter(c1_ids, c2_ids_dup, c1c2_d2d, max_deg)

        # Store all unique solutions and no match solutions.
        match_c1_ids_all += match_c1_ids
        match_c2_ids_all += match_c2_ids
        no_match_c1_all += no_match_c1
        match_d2d_all += match_d2d
        no_match_d2d_all += no_match_d2d

        print("  Unique matches:", len(match_c1_ids))
        if match_d2d_all:
            print("  Average match separation: {:.2f} arcsec".format(
                np.mean(Angle(match_d2d_all).arcsec)))
        print("  No match:", len(no_match_c1))
        print("  Re match:", len(c1_ids))

        if c1_ids.any():
            # Queried stars that were not matched to an observed star.
            c2_ids_r = np.setdiff1d(c2_ids, match_c2_ids_all)
            print("  Queried stars for re-match:", len(c2_ids_r))
            # To avoid messing with the indexes, change the coordinates
            # of already matched queried stars so that they can not
            # possibly be matched again.
            # TODO not sure '0.' is a value to be used.
            ra_q[match_c2_ids_all] = 0.
            dec_q[match_c2_ids_all] = 0.

    print('\nObserved stars matched:', len(match_c1_ids_all))
    print('Observed stars not matched:', len(no_match_c1_all))

    # Rejected magnitudes from queried catalog.
    q_rjct_mks = np.ones(len(query), dtype=bool)
    q_rjct_mks[match_c2_ids_all] = False
    print('Queried stars not matched:', sum(q_rjct_mks))

    return match_c1_ids_all, no_match_c1_all, match_d2d_all, no_match_d2d_all,\
        match_c2_ids_all, q_rjct_mks
