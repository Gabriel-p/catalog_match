
from os.path import exists
from os import makedirs
from os import listdir
from os.path import isfile, join
import numpy as np

from astroquery.irsa import Irsa
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord

# from astroquery.xmatch import XMatch

from modules import read_input
from modules import make_plots
from modules import out_data


def params_input():
    """
    Read data from file.
    """
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                if reader[0] == 'DC':
                    data_cols = map(int, reader[1:])
                if reader[0] == 'CA':
                    catalog = reader[1]
                    m_cat = reader[2]
                if reader[0] == 'MA':
                    max_arcsec = float(reader[1])

    return data_cols, catalog, m_cat, max_arcsec


def get_files():
    '''
    Store the paths and names of all the input clusters stored in the
    input folder.
    '''
    cl_files = [join('input/', f) for f in listdir('input/') if
                isfile(join('input/', f))]

    # Remove readme file it is still there.
    cl_files.remove('input/README.md')

    return cl_files


def cat_query(ids, ra_obs, dec_obs, catalog, m_cat):
    """
    """
    ra_rang, dec_rang = max(ra_obs) - min(ra_obs),\
        max(dec_obs) - min(dec_obs)
    ra_mid, dec_mid = .5 * (min(ra_obs) + max(ra_obs)),\
        .5 * (max(dec_obs) + min(dec_obs))

    print('ra (min, max):', min(ra_obs), max(ra_obs))
    print('dec (min, max):', min(dec_obs), max(dec_obs))
    print('Range (ra, dec):', ra_rang, dec_rang)
    print('Centre (ra, dec):', ra_mid, dec_mid)

    cent = SkyCoord(ra=ra_mid * u.degree, dec=dec_mid * u.degree, frame='icrs')

    # Set maximum limit for retrieved stars.
    Irsa.ROW_LIMIT = int(4 * len(ids))
    query = Irsa.query_region(
        cent, catalog=catalog, spatial='Box', width=dec_rang * u.deg)
    print("N (queried catalog):", len(query))

    m_qry, ra_qry, dec_qry = query[m_cat], query['ra'], query['dec']

    return query, m_qry, ra_qry, dec_qry


def cat_match(ra_obs, dec_obs, ra_qry, dec_qry):
    """
    Match catalogs.
    idx are indices into c2 that are the closest objects to each of the
    coordinates in c1
    d2d are the on-sky distances between them,
    d3d are the 3-dimensional distances
    """
    # Define observed and queried catalogs.
    c1 = SkyCoord(ra=ra_obs * u.degree, dec=dec_obs * u.degree)
    c2 = SkyCoord(ra_qry, dec_qry, unit=(u.degree, u.degree))

    idx, d2d, d3d = c1.match_to_catalog_sky(c2)
    # print('Unique matches:', len(set(idx)))
    # print('Filtered matches:', np.sum([d2d.arcsec < max_arcsec]))

    return idx, d2d


def match_filter(c1_ids, c2_ids, d2d, max_arcsec):
    """
    Reject duplicated matches and matched stars with distances beyond
    the maximum separation defined.
    """
    # Maximum separation in arcsec, convert to decimal degrees.
    max_deg = Angle(max_arcsec * u.arcsec).deg

    unique_c1_ids, dupl_c1_ids, no_match_c1 = [], [], []
    unique_c2_ids, unique_d2d = [], []
    # Indexes of matched star in both catalogs and distance between them.
    for c1_i, c2_i, d in zip(*[c1_ids, c2_ids, d2d]):
        # Filter by maximum allowed match distance.
        if d.deg <= max_deg:
            # Check if this queried star in c2 was already stored.
            if c2_i in unique_c2_ids:
                # Index of this stored c2 star.
                j = unique_c2_ids.index(c2_i)
                # Only replace with this match if the previous match had a
                # larger distance.
                if unique_d2d[j] > d:
                    # Store rejected replaced star.
                    dupl_c1_ids.append(unique_c1_ids[j])
                    # Replace star
                    unique_c1_ids[j] = c1_i
                    unique_d2d[j] = d
                else:
                    dupl_c1_ids.append(c1_i)
            else:
                unique_c1_ids.append(c1_i)
                unique_c2_ids.append(c2_i)
                unique_d2d.append(d)
        else:
            # If this observed star has no queried star closer than the max
            # distance allowed, discard.
            no_match_c1.append(c1_i)

    return unique_c1_ids, unique_c2_ids, no_match_c1, dupl_c1_ids


def main():
    """
    Observed catalog matcher.
    """
    # Irsa.print_catalogs()

    data_cols, catalog, m_cat, max_arcsec = params_input()

    # Process all files inside 'input/' folder.
    cl_files = get_files()
    for clust_file in cl_files:

        clust_name = clust_file.split('/')[1].split('.')[0]
        print("\nProcessing: {}\n".format(clust_name))

        # Get input data from file.
        ids, m_obs, ra_obs, dec_obs = read_input.main(clust_file, data_cols)

        # Query catalog.
        query, m_qry, ra_qry, dec_qry = cat_query(
            ids, ra_obs, dec_obs, catalog, m_cat)

        # Initial full list of observed and queried catalogs.
        c1_ids = [_ for _ in range(len(ids))]
        c2_ids = [_ for _ in range(len(query))]
        ra_q, dec_q = query['ra'][c2_ids], query['dec'][c2_ids]
        # Store all unique matches, and observed stars with no match.
        match_c1_ids_all, match_c2_ids_all, no_match_c1_all = [], [], []
        # Continue until no more duplicate matches exist.
        first_pass = True
        while c1_ids:

            # Match observed and queried catalogs.
            c2_ids_dup, c1c2_d2d = cat_match(
                ra_obs[c1_ids], dec_obs[c1_ids], ra_q, dec_q)
            # Store for plotting
            if first_pass:
                d2d = c1c2_d2d
                first_pass = False

            # Return unique ids for matched stars between catalogs, ids of
            # observed stars with no match found, and ids of observed stars
            # with a duplicated match that will be re-processed ('c1_ids').
            match_c1_ids, match_c2_ids, no_match_c1, c1_ids =\
                match_filter(c1_ids, c2_ids_dup, c1c2_d2d, max_arcsec)
            # Store all unique solutions and no match solutions.
            match_c1_ids_all += match_c1_ids
            match_c2_ids_all += match_c2_ids
            no_match_c1_all += no_match_c1

            print("  Match:", len(match_c1_ids))
            print("  Average match separation: {:.2f} arcsec".format(
                np.mean(c1c2_d2d).arcsec))
            print("  No match:", len(no_match_c1))
            print("  Re match:", len(c1_ids))

            if c1_ids:
                # Queried stars that were not matched to an observed star.
                c2_ids_r = [_ for _ in c2_ids if _ not in match_c2_ids_all]
                print("  Queried stars for re-match:", len(c2_ids_r))
                # print(len(set(c2_ids_r)))
                ra_q, dec_q = [], []
                for c2_i in c2_ids:
                    if c2_i in match_c2_ids_all:
                        # TODO
                        ra_q.append(0.)
                        dec_q.append(0.)
                    else:
                        ra_q.append(query['ra'][c2_i])
                        dec_q.append(query['dec'][c2_i])

        # Unique stars from observed catalog c1.
        id_unq, m_unq, ra_unq, dec_unq =\
            list(map(int, ids[match_c1_ids_all])),\
            m_obs[match_c1_ids_all], ra_obs[match_c1_ids_all],\
            dec_obs[match_c1_ids_all]

        print('Total stars matched:', len(m_unq))
        # Unique magnitudes from queried catalog c2.
        m_unq_q = m_qry[match_c2_ids_all]
        # Differences in matched coordinates.
        ra_unq_delta, dec_unq_delta =\
            Angle(ra_unq - ra_qry[match_c2_ids_all], unit='deg').arcsec,\
            Angle(dec_unq - dec_qry[match_c2_ids_all], unit='deg').arcsec
        # Observed stars with no match.
        id_rjct, m_rjct, ra_rjct, dec_rjct = ids[no_match_c1_all],\
            m_obs[no_match_c1_all], ra_obs[no_match_c1_all],\
            dec_obs[no_match_c1_all]
        print('Total stars not matched:', len(m_rjct))

        # Generate output dir if it doesn't exist.
        if not exists('output'):
            makedirs('output')

        make_plots.main(
            clust_name, m_cat, catalog, max_arcsec, m_obs, ra_obs, dec_obs,
            m_qry, ra_qry, dec_qry, m_unq, ra_unq, dec_unq, m_unq_q,
            d2d, ra_unq_delta, dec_unq_delta, m_rjct, ra_rjct, dec_rjct)

        out_data.main(
            clust_name, query, id_unq, ra_unq, dec_unq, match_c2_ids_all,
            id_rjct, m_rjct, ra_rjct, dec_rjct)

    print("End.")


if __name__ == '__main__':
    main()
