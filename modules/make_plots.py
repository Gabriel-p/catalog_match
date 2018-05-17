
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
        ra_qry, dec_qry, m_unq, ra_unq, dec_unq, m_unq_q,
        match_d2d_all, no_match_d2d_all, ra_unq_delta, dec_unq_delta, m_rjct,
        ra_rjct, dec_rjct):
    """
    Generate final plots.
    """
    print('Creating output plots.')
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(20, 40))
    gs = gridspec.GridSpec(24, 12)

    interval = ZScaleInterval()
    try:
        zmin_obs, zmax_obs = interval.get_limits(m_obs[~m_obs.mask])
        m_obs.fill_value = np.nanmax(m_obs)
        m_obs = m_obs.filled()
    except AttributeError:
        zmin_obs, zmax_obs = interval.get_limits(m_obs)
    try:
        zmin_qry, zmax_qry = interval.get_limits(m_qry[~m_qry.mask])
        m_qry.fill_value = np.nanmax(m_qry)
        m_qry = m_qry.filled()
    except AttributeError:
        zmin_qry, zmax_qry = interval.get_limits(m_qry)

    ax = plt.subplot(gs[0:6, 0:6])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Observed stars ({})".format(len(m_obs)), fontsize=14)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    # Plot all observed stars.
    st_sizes_arr = star_size(m_obs, zmin_obs, zmax_obs)
    plt.scatter(ra_obs, dec_obs, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    # ax.set_aspect('equal')

    ax = plt.subplot(gs[0:6, 6:12])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title(
        "Queried catalog {} ({})".format(catalog, len(m_qry)), fontsize=14)
    plt.xlabel(r'$\alpha_{queried}$', fontsize=18)
    plt.ylabel(r'$\delta_{queried}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    st_sizes_arr = star_size(m_qry, zmin_qry, zmax_qry)
    plt.scatter(ra_qry, dec_qry, marker='o', c='k', s=st_sizes_arr,
                zorder=4)
    ax.invert_xaxis()
    # ax.set_aspect('equal')

    ax = plt.subplot(gs[6:12, 0:6])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Observed stars with no match ({})".format(len(m_rjct)),
                 fontsize=14)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    st_sizes_arr = star_size(m_rjct, zmin_obs, zmax_obs)
    plt.scatter(ra_rjct, dec_rjct, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    # ax.set_aspect('equal')

    ax = plt.subplot(gs[6:12, 6:12])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Matched stars ({})".format(len(m_unq)), fontsize=14)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    st_sizes_arr = star_size(m_unq, zmin_obs, zmax_obs)
    plt.scatter(ra_unq, dec_unq, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    # ax.set_aspect('equal')

    # Plot density map.
    ax = plt.subplot(gs[12:18, 0:6])
    ax.set_title("Difference between observed and queried " + r"$\alpha$",
                 fontsize=14)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    xi, yi = np.linspace(ra_unq.min(), ra_unq.max(), 200),\
        np.linspace(dec_unq.min(), dec_unq.max(), 200)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    vals = np.array([ra_unq, dec_unq]).T
    zi = scipy.interpolate.griddata(
        vals, ra_unq_delta, (xi, yi), method='linear')

    zero_pt = 1. - zi.max() / (zi.max() - zi.min())
    cmap = LinearSegmentedColormap.from_list(
        'mycmap', [(0, 'blue'), (zero_pt, 'white'), (1, 'red')])
    im = plt.imshow(
        zi, vmin=zi.min(), vmax=zi.max(), origin='lower',
        extent=[ra_unq.min(), ra_unq.max(), dec_unq.min(), dec_unq.max()],
        cmap=cmap, aspect='auto')
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(r'$\Delta\,\alpha\,(arcsec)$', fontsize=12)
    ax.invert_xaxis()

    # Plot density map.
    ax = plt.subplot(gs[12:18, 6:12])
    ax.set_title("Difference between observed and queried " + r"$\delta$",
                 fontsize=14)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    xi, yi = np.linspace(ra_unq.min(), ra_unq.max(), 200),\
        np.linspace(dec_unq.min(), dec_unq.max(), 200)
    xi, yi = np.meshgrid(xi, yi)
    # Interpolate
    vals = np.array([ra_unq, dec_unq]).T
    zi = scipy.interpolate.griddata(
        vals, dec_unq_delta, (xi, yi), method='linear')

    zero_pt = 1. - zi.max() / (zi.max() - zi.min())
    cmap = LinearSegmentedColormap.from_list(
        'mycmap', [(0, 'blue'), (zero_pt, 'white'), (1, 'red')])
    im = plt.imshow(
        zi, vmin=zi.min(), vmax=zi.max(), origin='lower',
        extent=[ra_unq.min(), ra_unq.max(), dec_unq.min(), dec_unq.max()],
        cmap=cmap, aspect='auto')
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(r'$\Delta\,\delta\,(arcsec)$', fontsize=12)
    ax.invert_xaxis()

    ax = plt.subplot(gs[18:21, 0:6])
    ax.set_title("Separation between stars", fontsize=14)
    plt.xlabel(r'$d\,(arcsec)$', fontsize=18)
    ax.hist(match_d2d_all.arcsec, color='green', alpha=0.5,
            label='Match ({})'.format(len(match_d2d_all)))
    ax.hist(no_match_d2d_all.arcsec, color='red', alpha=0.5,
            label='No match ({})'.format(len(no_match_d2d_all)))
    ax.axvline(max_arcsec, color='k', linestyle='--')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, loc='upper right')
    leg.get_frame().set_alpha(0.5)

    ax = plt.subplot(gs[18:21, 6:12])
    ax.set_title("Number of stars in different instances", fontsize=14)
    ax.set_xlim(-0.2, 3.2)
    ax.grid(b=True, which='major', color='k', linestyle=':', lw=.5,
            zorder=1)
    x = np.arange(4)
    y = np.array([len(m_obs), len(m_qry), len(m_unq), len(m_rjct)])
    up = max(y) * .03
    ax.set_ylim(0, max(y) + 3 * up)
    ax.bar(x, y, align='center', width=0.2, color='g', zorder=4)
    for xi, yi, l in zip(*[x, y, list(map(str, y))]):
        ax.text(xi - len(l) * .02, yi + up, l,
                bbox=dict(facecolor='w', edgecolor='w', alpha=.5))
    ax.set_xticks(x)
    ax.set_xticklabels(
        ['Observed', 'Queried', 'Match', 'No match'])
    ax.tick_params(axis='x', which='major', labelsize=12)

    ax = plt.subplot(gs[21:24, 0:6])
    ax.set_title("Differences between matched selected magnitudes",
                 fontsize=14)
    plt.xlabel(r'$m_{obs}$', fontsize=18)
    plt.ylabel(r'$\Delta\, (m_{{obs}} - $' + '{})'.format(m_cat), fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='k', linestyle='--', lw=.5,
            zorder=1)
    plt.scatter(m_unq, m_unq - m_unq_q, marker='o', c='b', s=20, lw=.5,
                edgecolors='k', zorder=4)
    # ax.invert_xaxis()
    # ax.invert_yaxis()

    #
    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '.png', dpi=150, bbox_inches='tight')
