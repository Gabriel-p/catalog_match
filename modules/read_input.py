
from astroquery.vizier import Vizier
# from astroquery.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii


def in_data(clust_file, data_mode, data_cols):
    '''
    Read the file that stores the photometric data for all stars.
    '''
    # Loads the data in 'clust_file' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    # try:
    inp_data = ascii.read(
        clust_file, #format='commented_header',
        fill_values=[(ascii.masked, '99.999')])
    # except ascii.core.InconsistentTableError:

    if data_mode == 'num':
        m_idx, ra_idx, dec_idx = map(int, data_cols)
        m_obs, ra_obs, dec_obs = inp_data.columns[m_idx],\
            inp_data.columns[ra_idx], inp_data.columns[dec_idx]
    elif data_mode == 'nam':
        m_nam, ra_nam, dec_nam = data_cols
        m_obs, ra_obs, dec_obs = inp_data[m_nam], inp_data[ra_nam],\
            inp_data[dec_nam]

    N_obs = len(ra_obs)
    print("N (all stars) = {}".format(N_obs))

    ra_rang, dec_rang = max(ra_obs) - min(ra_obs),\
        max(dec_obs) - min(dec_obs)
    ra_mid, dec_mid = .5 * (min(ra_obs) + max(ra_obs)),\
        .5 * (max(dec_obs) + min(dec_obs))

    print('ra (min, max): {:.4f}, {:.4f}'.format(min(ra_obs), max(ra_obs)))
    print('dec (min, max): {:.4f}, {:.4f}'.format(min(dec_obs), max(dec_obs)))
    print('Range (ra, dec): {:.4f}, {:.4f}'.format(ra_rang, dec_rang))
    print('Centre (ra, dec): {:.4f}, {:.4f}'.format(ra_mid, dec_mid))

    return inp_data, m_obs, ra_obs, dec_obs, N_obs, ra_mid, dec_mid, ra_rang,\
        dec_rang


def cat_query(
    clust_name, N_obs, ra_mid, dec_mid, dec_rang, box_s, cat_mode,
        catalog):
    """
    Query selected catalog or read from file a previous stored version.
    """

    if cat_mode == 'query':
        print("\nFetching data from {} catalog.".format(catalog))
        txt = 'queried'
        cent = SkyCoord(
            ra=ra_mid * u.degree, dec=dec_mid * u.degree, frame='icrs')

        if str(box_s) == 'auto':
            width = dec_rang * u.deg
            print("Using auto width={:.3f}".format(width))
        else:
            width = box_s * u.deg
            print("Using manual width={:.3f}".format(width))

        # Vizier query
        # Unlimited rows, all columns
        v = Vizier(row_limit=-1, columns=['all'])
        query = v.query_region(SkyCoord(
            cent, frame='icrs'), width=width, catalog=[catalog])

        try:
            query = query[0]
        except IndexError:
            raise ValueError("Queried catalog came back empty")

        # # Irsa query
        # # Set maximum limit for retrieved stars.
        # Irsa.ROW_LIMIT = int(4 * N_obs)
        # query = Irsa.query_region(
        #     cent, catalog=catalog, spatial='Box', width=ra_rang * u.deg)

        out_file = clust_name + '_query.dat'
        print("Writing queried data to '{}' file.".format(out_file))
        # TODO waiting for fix to https://github.com/astropy/astropy/issues/7744
        ascii.write(
            query, 'output/' + out_file, #format='csv',
            overwrite=True)

    elif cat_mode == 'read':
        print("\nReading input catalog")
        txt = 'read'
        q_file = 'input/' + clust_name + '_query.dat'
        query = ascii.read(q_file)

    print("N ({} catalog): {}".format(txt, len(query)))

    return query
