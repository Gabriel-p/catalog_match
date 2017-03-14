
from astroquery.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii


def in_data(clust_file, data_mode, data_cols):
    '''
    Read the file that stores the photometric data for all stars.
    '''
    # Loads the data in 'clust_file' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    data = ascii.read(clust_file, fill_values=[(ascii.masked, '99.999')])

    if data_mode == 'num':
        m_idx, ra_idx, dec_idx = map(int, data_cols)
        m_obs, ra_obs, dec_obs = data.columns[m_idx], data.columns[ra_idx],\
            data.columns[dec_idx]
    elif data_mode == 'nam':
        m_nam, ra_nam, dec_nam = data_cols
        m_obs, ra_obs, dec_obs = data[m_nam], data[ra_nam], data[dec_nam]

    N_obs = len(ra_obs)
    print("N (all stars) = {}".format(N_obs))

    ra_rang, dec_rang = max(ra_obs) - min(ra_obs),\
        max(dec_obs) - min(dec_obs)
    ra_mid, dec_mid = .5 * (min(ra_obs) + max(ra_obs)),\
        .5 * (max(dec_obs) + min(dec_obs))

    print('ra (min, max):', min(ra_obs), max(ra_obs))
    print('dec (min, max):', min(dec_obs), max(dec_obs))
    print('Range (ra, dec):', ra_rang, dec_rang)
    print('Centre (ra, dec):', ra_mid, dec_mid)

    return m_obs, ra_obs, dec_obs, N_obs, ra_mid, dec_mid, ra_rang, dec_rang


def cat_query(clust_name, N_obs, ra_mid, dec_mid, ra_rang, dec_rang, cat_mode,
              catalog, m_cat):
    """
    Query selected catalog or read from file a previous stored version.
    """

    if cat_mode == 'query':
        txt = 'queried'
        cent = SkyCoord(
            ra=ra_mid * u.degree, dec=dec_mid * u.degree, frame='icrs')
        # Set maximum limit for retrieved stars.
        Irsa.ROW_LIMIT = int(4 * N_obs)
        query = Irsa.query_region(
            cent, catalog=catalog, spatial='Box', width=dec_rang * u.deg)
    elif cat_mode == 'read':
        txt = 'read'
        q_file = 'input/' + clust_name + '_query.dat'
        query = ascii.read(q_file)
        # TODO these two lines below are almost certainly not needed.
        query['ra'] = query['ra'] * u.deg
        query['dec'] = query['dec'] * u.deg

    print("N ({} catalog): {}".format(txt, len(query)))

    return query
