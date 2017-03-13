
from astropy.io import ascii
from astropy.table import Column
from astropy.table import hstack


def main(
        clust_name, clust_file, out_format, out_cols, query, match_c1_ids_all,
        match_c2_ids_all, match_d2d_all, no_match_c1_all, no_match_d2d_all):
    """
    Write output data to files.
    """
    # Store full queried catalog.
    f_out = 'output/' + clust_name + '_query.dat'
    ascii.write(
        query, output=f_out, overwrite=True, format='fixed_width',
        delimiter=' ', fill_values=[(ascii.masked, '--')],
        formats={'ra': '%4.9f', 'dec': '%4.9f'})

    # Read input data.
    in_data = ascii.read(
        clust_file, fill_values=[(ascii.masked, '99.999')])

    # Filter matched observed stars only.
    in_data_match = in_data[match_c1_ids_all]
    # Generate column with separation data and add it to the table.
    d2d_match = Column(match_d2d_all.arcsec, name='d_arcsec')
    in_data_match.add_column(d2d_match)
    # Define table using only matched stars from queried catalog.
    t_match_c2 = query[match_c2_ids_all]

    # Filter *not* matched observed stars only.
    in_data_no_match = in_data[no_match_c1_all]
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
            comb_dat, output=f_match, overwrite=True, format='fixed_width',
            delimiter=' ', fill_values=[(ascii.masked, '--')],
            formats={'d_arcsec': '%.5f'})
        print("Data for all matched stars written to file.")

    elif out_format == 'man':
        in_data_match.add_column(t_match_c2['ra'])
        in_data_match.add_column(t_match_c2['dec'])
        # Add selected columns.
        for col in out_cols:
            in_data_match.add_column(t_match_c2[col])

        # Write matched stars to file.
        ascii.write(
            in_data_match, output=f_match, overwrite=True,
            format='fixed_width', delimiter=' ',
            fill_values=[(ascii.masked, '--')], formats={'d_arcsec': '%.5f'})

        print("Data for all matched stars written to file.")

    # Write *not* matched stars to file.
    ascii.write(
        in_data_no_match, output=f_no_match, overwrite=True,
        format='fixed_width', delimiter=' ',
        fill_values=[(ascii.masked, '--')], formats={'d_arcsec': '%.5f'})
    print("Data for stars with no match written to file.")
