
from astropy.io import ascii
from astropy.table import Column


def main(
        clust_name, query, id_unq, ra_unq, dec_unq, match_c2_ids_all,
        id_rjct, m_rjct, ra_rjct, dec_rjct):
    """
    """
    # Define table using only matched stars from queried catalog.
    t_match_c2 = query[match_c2_ids_all]
    # Add columns with observed data.
    obs_ids = Column(id_unq, name='IDs_Obs')
    ra_obs = Column(ra_unq, name='ra_Obs')
    dec_obs = Column(dec_unq, name='dec_Obs')
    # Insert before the first table column
    t_match_c2.add_column(obs_ids, index=0)
    t_match_c2.add_column(ra_obs, index=1)
    t_match_c2.add_column(dec_obs, index=2)

    # Write matched stars to file.
    f_out = 'output/' + clust_name + '_match.dat'
    ascii.write(
        t_match_c2, output=f_out, overwrite=True, format='fixed_width',
        delimiter=' ', fill_values=[(ascii.masked, '--')])

    print("Data for all matched stars written to file.")

    # Write not matched stars to file.
    f_name = 'output/' + clust_name + '_no_match.dat'
    with open(f_name, 'w') as f_out:
        f_out.write("#ID           mag       ra_obs      dec_obs\n")
    with open(f_name, "a") as f_out:
        for line_f in zip(*[map(int, id_rjct), m_rjct, ra_rjct, dec_rjct]):
            f_out.write('''{:<10} {:>6} {:>12} {:>12}'''.format(*line_f))
            f_out.write('\n')

    print("Data for stars with no match written to file.")
