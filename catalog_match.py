
from os.path import exists
from os import makedirs
from os import listdir
from os.path import isfile, join
import numpy as np

from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord

# from astroquery.xmatch import XMatch

from modules import read_input
from modules import make_plots
from modules import out_data


def params_input():
    """
    Read input parameters from 'params_input.dat' file.
    """
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):
            if not line.startswith("#") and line.strip() != '':
                reader = line.split()
                if reader[0] == 'DC':
                    data_mode = reader[1]
                    data_cols = reader[2:]
                if reader[0] == 'CA':
                    cat_mode = reader[1]
                    catalog = reader[2]
                    m_cat = reader[3]
                if reader[0] == 'MA':
                    max_arcsec = float(reader[1])
                if reader[0] == 'FI':
                    out_fig = reader[1]
                if reader[0] == 'OF':
                    out_format = reader[1]
                    out_cols = reader[2:]

    return data_mode, data_cols, cat_mode, catalog, m_cat, max_arcsec,\
        out_fig, out_format, out_cols


def get_files():
    '''
    Store the paths and names of all the input clusters stored in the
    input folder.
    '''
    cl_files = [join('input/', f) for f in listdir('input/') if
                isfile(join('input/', f))]

    # Remove readme file it is still there.
    cl_files.remove('input/README.md')
    # Remove any '_query.dat' file.
    for _ in cl_files:
        if _.endswith('_query.dat'):
            cl_files.remove(_)

    return cl_files


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


def match_filter(c1_ids, c2_ids, d2d, max_arcsec):
    """
    Reject duplicated matches and matched stars with distances beyond
    the maximum separation defined.
    """
    # Maximum separation in arcsec, convert to decimal degrees.
    max_deg = Angle(max_arcsec * u.arcsec).deg

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

    return match_c1_ids, match_c2_ids, no_match_c1, dupl_c1_ids, match_d2d,\
        no_match_d2d


def main():
    """
    Observed catalog matcher. To obtain list of supported catalogs:

    Irsa.print_catalogs()
    """

    data_mode, data_cols, cat_mode, catalog, m_cat, max_arcsec, out_fig,\
        out_format, out_cols = params_input()

    # Process all files inside 'input/' folder.
    cl_files = get_files()
    for clust_file in cl_files:

        clust_name = clust_file.split('/')[1].split('.')[0]
        print("\nProcessing: {}\n".format(clust_name))

        # Get input data from file.
        m_obs, ra_obs, dec_obs, N_obs, ra_mid, dec_mid, ra_rang,\
            dec_rang = read_input.in_data(clust_file, data_mode, data_cols)

        # Query catalog.
        query = read_input.cat_query(
            clust_name, N_obs, ra_mid, dec_mid, ra_rang, dec_rang, cat_mode,
            catalog, m_cat)

        # Initial full list of observed and queried catalogs.
        c1_ids = [_ for _ in range(len(ra_obs))]
        c2_ids = [_ for _ in range(len(query))]
        ra_q, dec_q = query['ra'][c2_ids], query['dec'][c2_ids]
        # Store all unique matches, and observed stars with no match.
        match_c1_ids_all, match_c2_ids_all, no_match_c1_all, match_d2d_all,\
            no_match_d2d_all = [], [], [], [], []
        # Continue until no more duplicate matches exist.
        while c1_ids:

            # Match observed and queried catalogs.
            c2_ids_dup, c1c2_d2d = cat_match(
                ra_obs[c1_ids], dec_obs[c1_ids], ra_q, dec_q)

            # Return unique ids for matched stars between catalogs, ids of
            # observed stars with no match found, and ids of observed stars
            # with a duplicated match that will be re-processed ('c1_ids').
            match_c1_ids, match_c2_ids, no_match_c1, c1_ids, match_d2d,\
                no_match_d2d = match_filter(
                    c1_ids, c2_ids_dup, c1c2_d2d, max_arcsec)
            # Store all unique solutions and no match solutions.
            match_c1_ids_all += match_c1_ids
            match_c2_ids_all += match_c2_ids
            no_match_c1_all += no_match_c1
            match_d2d_all += match_d2d
            no_match_d2d_all += no_match_d2d

            print("  Match:", len(match_c1_ids))
            print("  Average match separation: {:.2f} arcsec".format(
                np.mean(c1c2_d2d).arcsec))
            print("  No match:", len(no_match_c1))
            print("  Re match:", len(c1_ids))

            if c1_ids:
                # Queried stars that were not matched to an observed star.
                c2_ids_r = [_ for _ in c2_ids if _ not in match_c2_ids_all]
                print("  Queried stars for re-match:", len(c2_ids_r))
                # To avoid messing with the indexes, change the coordinates
                # of already matched queried stars so that they can not
                # possibly be matched again.
                ra_q, dec_q = [], []
                for c2_i in c2_ids:
                    if c2_i in match_c2_ids_all:
                        # TODO
                        ra_q.append(0.)
                        dec_q.append(0.)
                    else:
                        ra_q.append(query['ra'][c2_i])
                        dec_q.append(query['dec'][c2_i])

        print('Total stars matched:', len(match_c1_ids_all))
        print('Total stars not matched:', len(no_match_c1_all))

        # Unique stars from observed catalog.
        m_unq, ra_unq, dec_unq =\
            m_obs[match_c1_ids_all], ra_obs[match_c1_ids_all],\
            dec_obs[match_c1_ids_all]

        # Generate output dir if it doesn't exist.
        if not exists('output'):
            makedirs('output')

        # Store match and no match separations as 'Angle' objects.
        match_d2d_all, no_match_d2d_all = Angle(match_d2d_all * u.degree),\
            Angle(no_match_d2d_all * u.degree)

        if out_fig == 'y':
            m_qry, ra_qry, dec_qry = query[m_cat], query['ra'], query['dec']
            # Unique magnitudes from queried catalog.
            m_unq_q = m_qry[match_c2_ids_all]
            # Differences in matched coordinates.
            ra_unq_delta, dec_unq_delta =\
                Angle(ra_unq - ra_qry[match_c2_ids_all], unit='deg').arcsec,\
                Angle(dec_unq - dec_qry[match_c2_ids_all], unit='deg').arcsec
            # Observed stars with no match.
            m_rjct, ra_rjct, dec_rjct = m_obs[no_match_c1_all],\
                ra_obs[no_match_c1_all], dec_obs[no_match_c1_all]
            make_plots.main(
                clust_name, m_cat, catalog, max_arcsec, m_obs, ra_obs, dec_obs,
                m_qry, ra_qry, dec_qry, m_unq, ra_unq, dec_unq, m_unq_q,
                match_d2d_all, no_match_d2d_all, ra_unq_delta,
                dec_unq_delta, m_rjct, ra_rjct, dec_rjct)
        else:
            print("No output figure created.")

        out_data.main(
            clust_name, clust_file, out_format, out_cols, query,
            match_c1_ids_all, match_c2_ids_all, match_d2d_all, no_match_c1_all,
            no_match_d2d_all)

    print("End.")


if __name__ == '__main__':
    main()
