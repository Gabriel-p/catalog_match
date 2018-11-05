
from os.path import exists
from os import makedirs
from os import listdir
from os.path import isfile, join

from astropy.coordinates import Angle
from astropy import units as u
from astropy.io import ascii

# from astroquery.xmatch import XMatch

from modules import read_input
from modules import match_cats
from modules import make_plots
from modules import out_data


def main():
    """
    Observed catalog matcher.
    """

    data_mode, data_cols, cat_mode, catalog, m_qry, ra_qry, de_qry, box_s,\
        max_arcsec, out_fig, out_format, out_cols = params_input()

    # Generate output dir if it doesn't exist.
    if not exists('output'):
        makedirs('output')

    # Process all files inside 'input/' folder.
    cl_files = get_files()
    for clust_file in cl_files:

        clust_name = clust_file.split('/')[1].split('.')[0]
        print("\nProcessing: {}\n".format(clust_name))

        # Get input data from file.
        keep_going = True
        try:
            inp_data, m_obs, ra_obs, dec_obs, N_obs, ra_mid, dec_mid, ra_rang,\
                dec_rang = read_input.in_data(clust_file, data_mode, data_cols)
        except ascii.core.InconsistentTableError as err:
            print("{}\n\nERROR: could not read data file {}".format(
                err, clust_file))
            keep_going = False

        if keep_going is True:
            # Query catalog.
            query = read_input.cat_query(
                clust_name, N_obs, ra_mid, dec_mid, dec_rang, box_s,
                cat_mode, catalog)

            # Match catalogs.
            match_c1_ids_all, no_match_c1_all, match_d2d_all,\
                no_match_d2d_all, match_c2_ids_all = match_cats.main(
                    ra_qry, de_qry, max_arcsec, ra_obs, dec_obs, query)

            # Store match and no match separations as 'Angle' objects.
            match_d2d_all, no_match_d2d_all = Angle(match_d2d_all * u.degree),\
                Angle(no_match_d2d_all * u.degree)

            if out_fig == 'y':

                # Unique stars from observed catalog.
                m_unq, ra_unq, dec_unq =\
                    m_obs[match_c1_ids_all], ra_obs[match_c1_ids_all],\
                    dec_obs[match_c1_ids_all]

                # Unique magnitudes from queried catalog.
                m_unq_q = query[m_qry][match_c2_ids_all]
                # Differences in matched coordinates.
                ra_unq_delta, dec_unq_delta =\
                    Angle(ra_unq - query[ra_qry][match_c2_ids_all],
                          unit='deg').arcsec,\
                    Angle(dec_unq - query[de_qry][match_c2_ids_all],
                          unit='deg').arcsec
                # Observed stars with no match.
                m_rjct, ra_rjct, dec_rjct = m_obs[no_match_c1_all],\
                    ra_obs[no_match_c1_all], dec_obs[no_match_c1_all]

                if match_c1_ids_all:
                    make_plots.main(
                        clust_name, m_qry, catalog, max_arcsec, m_obs, ra_obs,
                        dec_obs, query[m_qry], query[ra_qry], query[de_qry],
                        m_unq, ra_unq, dec_unq, m_unq_q, match_d2d_all,
                        no_match_d2d_all, ra_unq_delta, dec_unq_delta, m_rjct,
                        ra_rjct, dec_rjct)
                    print("Output figure created.")
                else:
                    print("  ERROR: no matches to plot.")
            else:
                print("No output figure created.")

            out_data.main(
                clust_name, inp_data, cat_mode, ra_qry, de_qry, query,
                out_format, out_cols, match_c1_ids_all, match_c2_ids_all,
                match_d2d_all, no_match_c1_all, no_match_d2d_all)

    print("\nEnd.")


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
                    m_qry = reader[3]
                    ra_qry = reader[4]
                    de_qry = reader[5]
                    box_s = reader[6] if reader[6] == 'auto' else\
                        float(reader[6])
                if reader[0] == 'MA':
                    max_arcsec = float(reader[1])
                if reader[0] == 'FI':
                    out_fig = reader[1]
                if reader[0] == 'OF':
                    out_format = reader[1]
                    out_cols = reader[2:]

    return data_mode, data_cols, cat_mode, catalog, m_qry, ra_qry, de_qry,\
        box_s, max_arcsec, out_fig, out_format, out_cols


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
    for _ in cl_files:
        if _.endswith('_query.dat'):
            cl_files.remove(_)

    return cl_files


if __name__ == '__main__':
    main()
