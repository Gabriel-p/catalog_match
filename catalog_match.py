
from os.path import exists
from os import makedirs
from os import listdir
from os.path import isfile, join

from astroquery.irsa import Irsa
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord

# from astroquery.xmatch import XMatch
# from astropy.coordinates import ICRS
# import astropy.coordinates as coord

from modules import read_input
from modules import make_plots


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

    return cl_files


def match_filter(idx, d2d, max_arcsec):
    """
    Reject duplicated matches and matched stars with distances beyond
    the maximum separation defined.
    """
    # Maximum separation in arcsec, convert to decimal degrees.
    max_deg = Angle(max_arcsec * u.arcsec).deg

    unique_c1_ids, unique_c1_d2d, unique_c2_ids, reject_c1_ids =\
        [], [], [], []
    # Indexes of matched stars in both catalogs.
    for c1_i, c2_i in enumerate(idx):
        # Filter by maximum match distance.
        if d2d[c1_i].deg <= max_deg:
            # Check if this star in c2 was already processed.
            if c2_i in unique_c2_ids:
                # Index of this processed c2 star.
                j = unique_c2_ids.index(c2_i)
                # Only replace with this match if the previous match had a
                # larger distance.
                if unique_c1_d2d[j] > d2d[c1_i]:
                    # Store rejected replaced star
                    reject_c1_ids.append(unique_c1_ids[j])
                    # Replace star
                    unique_c1_ids[j] = c1_i
                    unique_c1_d2d[j] = d2d[c1_i]
                    unique_c2_ids[j] = c2_i
                else:
                    reject_c1_ids.append(c1_i)
            else:
                unique_c1_ids.append(c1_i)
                unique_c1_d2d.append(d2d[c1_i])
                unique_c2_ids.append(c2_i)
        else:
            reject_c1_ids.append(c1_i)

    return unique_c1_ids, unique_c2_ids, reject_c1_ids


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

        ra_rang, dec_rang = max(ra_obs) - min(ra_obs),\
            max(dec_obs) - min(dec_obs)
        ra_mid, dec_mid = .5 * (min(ra_obs) + max(ra_obs)),\
            .5 * (max(dec_obs) + min(dec_obs))

        print('ra (min, max):', min(ra_obs), max(ra_obs))
        print('dec (min, max):', min(dec_obs), max(dec_obs))
        print('Range (ra, dec):', ra_rang, dec_rang)
        print('Centre (ra, dec):', ra_mid, dec_mid)

        cent = SkyCoord(
            ra=ra_mid * u.degree, dec=dec_mid * u.degree, frame='icrs')

        # Set maximum limit for retrieved stars.
        Irsa.ROW_LIMIT = int(4 * len(ids))
        query = Irsa.query_region(
            cent, catalog=catalog, spatial='Box', width=dec_rang * u.deg)
        # import pickle
        # with open('query.pkl', 'wb') as f:
        #     pickle.dump(query, f)
        # with open('query.pkl', 'rb') as f:
        #     query = pickle.load(f)
        # print('ra:', query['ra'].min(), query['ra'].max())
        # print('dec:', query['dec'].min(), query['dec'].max())
        # print(query)

        m_qry, ra_qry, dec_qry = query[m_cat], query['ra'], query['dec']

        # Define observed and queried catalogs.
        c1 = SkyCoord(ra=ra_obs * u.degree, dec=dec_obs * u.degree)
        c2 = SkyCoord(query['ra'], query['dec'], unit=(u.degree, u.degree))

        # Match catalogs.
        # idx are indices into c2 that are the closest objects to each of the
        # coordinates in c1
        # d2d are the on-sky distances between them,
        # d3d are the 3-dimensional distances
        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
        print("Queried catalog:", len(query))
        N_unq_no_filter = len(set(idx))
        print('Unique matches:', N_unq_no_filter)
        # print('Filtered matches:', np.sum([d2d.arcsec < max_arcsec]))
        # Store all matched magnitudes for plotting.
        m_match_all = query[idx][m_cat]

        unique_c1_ids, unique_c2_ids, reject_c1_ids = match_filter()

        # Unique stars from observed c1.
        m_unq, ra_unq, dec_unq = m_obs[unique_c1_ids], ra_obs[unique_c1_ids],\
            dec_obs[unique_c1_ids]
        print('Unique filtered matches:', len(m_unq))
        # Unique magnitudes from queried c2.
        m_unq_q = m_qry[unique_c2_ids]
        # Differences in matched coordinates.
        ra_unq_delta, dec_unq_delta =\
            Angle(ra_unq - ra_qry[unique_c2_ids], unit='deg').arcsec,\
            Angle(dec_unq - dec_qry[unique_c2_ids], unit='deg').arcsec
        # Rejected observed stars.
        m_rjct, ra_rjct, dec_rjct = m_obs[reject_c1_ids],\
            ra_obs[reject_c1_ids], dec_obs[reject_c1_ids]
        print('Observed stars with no match:', len(m_rjct))

        # Generate output dir if it doesn't exist.
        if not exists('output'):
            makedirs('output')

        make_plots.main(
            clust_name, m_cat, catalog, max_arcsec, m_obs, ra_obs, dec_obs,
            m_qry, ra_qry, dec_qry, m_unq, ra_unq, dec_unq, m_unq_q,
            N_unq_no_filter, m_match_all, d2d, ra_unq_delta, dec_unq_delta,
            m_rjct, ra_rjct, dec_rjct)

    print("End.")


if __name__ == '__main__':
    main()
