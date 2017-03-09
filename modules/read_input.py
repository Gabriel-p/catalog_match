
import numpy as np


def main(data_file, data_cols):
    '''
    Read the file that stores the photometric data for all stars.
    '''
    # Loads the data in 'data_file' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    data = np.genfromtxt(data_file, dtype=float, filling_values=99.999,
                         unpack=True)
    id_idx, m_idx, ra_idx, dec_idx = data_cols
    ids, m_obs, ra_obs, dec_obs = data[id_idx], data[m_idx], data[ra_idx],\
        data[dec_idx]
    print("N (all stars) = {}".format(len(ids)))

    return ids, m_obs, ra_obs, dec_obs
