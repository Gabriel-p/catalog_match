
from astropy.io import ascii
from astropy.table import Column
from astropy.table import hstack
import logging


def main(
        clust_name, inp_data, cat_mode, ra_qry, de_qry, query, out_format,
        out_cols, match_c1_ids_all, match_c2_ids_all, match_d2d_all,
        no_match_c1_all, no_match_d2d_all):
    """
    Write output data to files.
    """

    if cat_mode == 'query':
        # Store full queried catalog.
        f_out = 'output/' + clust_name + '_query.dat'
        ascii.write(query, output=f_out, overwrite=True)  # format='csv'

    # Filter matched observed stars only.
    in_data_match = inp_data[match_c1_ids_all]
    # Generate column with separation data and add it to the table.
    d2d_match = Column(match_d2d_all.arcsec, name='d_arcsec')
    in_data_match.add_column(d2d_match)
    # Define table using only matched stars from queried catalog.
    t_match_c2 = query[match_c2_ids_all]

    # Filter *not* matched observed stars only.
    in_data_no_match = inp_data[no_match_c1_all]
    # Generate column with separation data and add it to the table.
    d2d_no_match = Column(no_match_d2d_all.arcsec, name='d_arcsec')
    in_data_no_match.add_column(d2d_no_match)

    # Output data file name for matched and not matched stars.
    f_match = 'output/' + clust_name + '_match.dat'
    f_no_match = 'output/' + clust_name + '_no_match.dat'

    if out_format == 'all':
        # Write matched stars to file.
        # Combine input data with queried data for matched stars.
        comb_dat = hstack([in_data_match, t_match_c2])
        ascii.write(
            comb_dat, output=f_match, overwrite=True,  # format='csv',
            formats={'d_arcsec': '%.4f'})
        logging.info("Data for all matched stars written to file.")

    elif out_format == 'man':
        in_data_match.add_column(t_match_c2[ra_qry])
        in_data_match.add_column(t_match_c2[de_qry])
        # Add selected columns.
        for col in out_cols:
            in_data_match.add_column(t_match_c2[col])

        # Write matched stars to file.
        ascii.write(
            in_data_match, output=f_match, overwrite=True,  # format='csv',
            formats={'d_arcsec': '%.4f'})

        logging.info("Data for all matched stars written to file.")

    # Write *not* matched stars to file.
    ascii.write(
        in_data_no_match, output=f_no_match, overwrite=True,  # format='csv',
        formats={'d_arcsec': '%.4f'})
    logging.info("Data for stars with no match written to file.")
