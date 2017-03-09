
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.offsetbox as offsetbox
from astropy.visualization import ZScaleInterval
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.interpolate


def star_size(mag, zmin, zmax):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    mag = np.array(mag).clip(zmin, zmax)
    sizes = .1 + 100. * (10 ** ((np.array(mag) - zmin) / -2.5))
    return sizes


def reject_outliers(data, m=2.):
    """
    Reject outliers, http://stackoverflow.com/a/16562028/1391441
    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else 0.
    return data[s < m]


def main(
    clust_name, m_cat, catalog, max_arcsec, m_obs, ra_obs, dec_obs, m_qry,
        ra_qry, dec_qry, m_unq, ra_unq, dec_unq, m_unq_q, N_unq_no_filter,
        m_match_all, d2d, ra_unq_delta, dec_unq_delta, m_rjct, ra_rjct,
        dec_rjct):
    """
    Generate final plots.
    """
    print('Creating output plots.')
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(30, 30))
    gs = gridspec.GridSpec(12, 12)

    ax = plt.subplot(gs[0:4, 0:4])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Observed chart ({})".format(len(m_obs)), fontsize=16)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    # Plot all observed stars.
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(np.array(m_obs))
    st_sizes_arr = star_size(m_obs, zmin, zmax)
    plt.scatter(ra_obs, dec_obs, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    ax.set_aspect('equal')

    ax = plt.subplot(gs[0:4, 4:8])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title(
        "Queried catalog {} ({})".format(catalog, len(m_qry)), fontsize=16)
    plt.xlabel(r'$\alpha_{queried}$', fontsize=18)
    plt.ylabel(r'$\delta_{queried}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    zmin, zmax = interval.get_limits(m_obs)
    st_sizes_arr = star_size(m_qry, zmin, zmax)
    plt.scatter(ra_qry, dec_qry, marker='o', c='k', s=st_sizes_arr,
                zorder=4)
    ax.invert_xaxis()
    ax.set_aspect('equal')

    ax = plt.subplot(gs[0:4, 8:12])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Matched chart ({})".format(len(m_unq)), fontsize=16)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    zmin, zmax = interval.get_limits(np.array(m_obs))
    st_sizes_arr = star_size(m_unq, zmin, zmax)
    plt.scatter(ra_unq, dec_unq, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    ax.set_aspect('equal')

    ax = plt.subplot(gs[4:8, 0:4])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Stars with no match ({})".format(len(m_rjct)), fontsize=16)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    zmin, zmax = interval.get_limits(np.array(m_obs))
    st_sizes_arr = star_size(m_rjct, zmin, zmax)
    plt.scatter(ra_rjct, dec_rjct, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    ax.set_aspect('equal')

    # Plot density map.
    ax = plt.subplot(gs[4:8, 4:8])
    ax.set_title("Difference between observed and queried " + r"$\alpha$",
                 fontsize=16)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    xi, yi = np.linspace(ra_unq.min(), ra_unq.max(), 50),\
        np.linspace(dec_unq.min(), dec_unq.max(), 50)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(
        ra_unq, dec_unq, ra_unq_delta, function='linear')
    zi = rbf(xi, yi)
    zero_pt = 1. - zi.max() / (zi.max() - zi.min())
    cmap = LinearSegmentedColormap.from_list(
        'mycmap', [(0, 'blue'), (zero_pt, 'white'), (1, 'red')])
    im = plt.imshow(
        zi, vmin=zi.min(), vmax=zi.max(), origin='lower',
        extent=[ra_unq.min(), ra_unq.max(), dec_unq.min(), dec_unq.max()],
        cmap=cmap)
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(r'$\Delta\,\alpha\,(arcsec)$', fontsize=12)
    ax.invert_xaxis()

    # Plot density map.
    ax = plt.subplot(gs[4:8, 8:12])
    ax.set_title("Difference between observed and queried " + r"$\delta$",
                 fontsize=16)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    xi, yi = np.linspace(ra_unq.min(), ra_unq.max(), 50),\
        np.linspace(dec_unq.min(), dec_unq.max(), 50)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    rbf = scipy.interpolate.Rbf(
        ra_unq, dec_unq, dec_unq_delta, function='linear')
    zi = rbf(xi, yi)
    zero_pt = 1. - zi.max() / (zi.max() - zi.min())
    cmap = LinearSegmentedColormap.from_list(
        'mycmap', [(0, 'blue'), (zero_pt, 'white'), (1, 'red')])
    im = plt.imshow(
        zi, vmin=zi.min(), vmax=zi.max(), origin='lower',
        extent=[ra_unq.min(), ra_unq.max(), dec_unq.min(), dec_unq.max()],
        cmap=cmap)
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(r'$\Delta\,\delta\,(arcsec)$', fontsize=12)
    ax.invert_xaxis()

    ax = plt.subplot(gs[8:10, 0:4])
    ax.set_title("Differences between selected magnitudes", fontsize=16)
    plt.xlabel(r'$m_{obs}$', fontsize=18)
    plt.ylabel(r'$\Delta\, (m_{{obs}} - $' + '{})'.format(m_cat), fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='k', linestyle='--', lw=.5,
            zorder=1)
    plt.scatter(m_obs, m_obs - m_match_all, marker='o', c='grey', s=30,
                zorder=2)
    plt.scatter(m_unq, m_unq - m_unq_q, marker='o', c='b', s=20, lw=.5,
                edgecolors='k', zorder=4)
    ax.invert_xaxis()
    ax.invert_yaxis()

    ax = plt.subplot(gs[8:10, 4:8])
    ax.set_title("Distances between matched stars", fontsize=16)
    plt.xlabel(r'$d\,(arcsec)$', fontsize=18)
    plt.hist(d2d.arcsec, 50)
    ax.axvline(max_arcsec, color='k', linestyle='--')
    text = 'N = {}'.format(len(d2d)) + '\n' +\
        'N(<{}) = {}'.format(max_arcsec, len(m_unq))
    ob = offsetbox.AnchoredText(text, loc=1, prop=dict(size=14))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)

    ax = plt.subplot(gs[8:10, 8:12])
    ax.set_title("Number of stars in different instances", fontsize=16)
    ax.set_xlim(-0.2, 4.2)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='k', linestyle=':', lw=.5,
            zorder=1)
    x = np.arange(5)
    y = np.array(
        [len(m_obs), len(m_qry), N_unq_no_filter, len(m_unq), len(m_rjct)])
    ax.bar(x, y, align='center', width=0.2, color='g', zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(
        ['Observed', 'Queried', 'Unique match', 'Unique + filter', 'No match'])
    ax.tick_params(axis='x', which='major', labelsize=12)

    #
    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '.png', dpi=150, bbox_inches='tight')
