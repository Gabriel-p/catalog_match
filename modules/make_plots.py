
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
    clust_name, m_cat, catalog, max_arcsec, m_obs, m_obs_nam, ra_obs, dec_obs,
        m_qry, ra_qry, dec_qry, m_unq, ra_unq, dec_unq, m_unq_q, m_rjct_q,
        match_d2d_all, no_match_d2d_all, ra_unq_delta, dec_unq_delta, m_rjct,
        ra_rjct, dec_rjct, max_mag_delta, mag_filter_data):
    """
    Generate final plots.
    """
    print('\nCreating output plots.')
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(40, 25))
    gs = gridspec.GridSpec(15, 24)

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
    plt.xlim(min(ra_qry), max(ra_qry))
    plt.ylim(min(dec_qry), max(dec_qry))
    ax.set_title("{} ({})".format(catalog, len(m_qry)), fontsize=14)
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

    ax = plt.subplot(gs[0:6, 12:18])
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

    ax = plt.subplot(gs[0:6, 18:24])
    plt.xlim(min(ra_obs), max(ra_obs))
    plt.ylim(min(dec_obs), max(dec_obs))
    ax.set_title("Observed stars with no match ({}; ~{:.1f}%)".format(
        len(m_rjct), (100. * len(m_rjct)) / float(len(m_obs))), fontsize=14)
    plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
    plt.ylabel(r'$\delta_{obs}$', fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='-', lw=.5,
            zorder=1)
    st_sizes_arr = star_size(m_rjct, zmin_obs, zmax_obs)
    plt.scatter(ra_rjct, dec_rjct, marker='o', c='k', s=st_sizes_arr, zorder=4)
    ax.invert_xaxis()
    # ax.set_aspect('equal')

    def densMap(gs_x, strID, radec_unq_delta):
        # Plot density map.
        ax = plt.subplot(gs[6:12, gs_x[0]:gs_x[1]])
        ax.set_title(r'$\Delta\,$' + strID + ' (observed - queried) [arcsec]',
                     fontsize=14)
        plt.xlabel(r'$\alpha_{obs}$', fontsize=18)
        if gs_x[0] == 0:
            plt.ylabel(r'$\delta_{obs}$', fontsize=18)
        else:
            # plt.yticks([])
            ax.yaxis.set_ticklabels([])
        xi, yi = np.linspace(ra_unq.min(), ra_unq.max(), 50),\
            np.linspace(dec_unq.min(), dec_unq.max(), 50)
        xi, yi = np.meshgrid(xi, yi)
        # Interpolate
        vals = np.array([ra_unq, dec_unq]).T
        zi = scipy.interpolate.griddata(
            vals, radec_unq_delta, (xi, yi), method='linear')

        zero_pt = 1. - np.nanmax(zi) / (np.nanmax(zi) - np.nanmin(zi))
        cmap = LinearSegmentedColormap.from_list(
            'mycmap', [(0, 'blue'), (zero_pt, 'white'), (1, 'red')])
        im = plt.imshow(
            zi, vmin=np.nanmin(zi), vmax=np.nanmax(zi), origin='lower',
            extent=[ra_unq.min(), ra_unq.max(), dec_unq.min(), dec_unq.max()],
            cmap=cmap)  # , aspect='square')
        # Colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        # cbar.set_label(r'$\Delta\,$' + strID + ' (arcsec)', fontsize=12)
        ax.invert_xaxis()

    # RA, DEC difference dens maps
    densMap([0, 6], r"$\alpha$", ra_unq_delta)
    densMap([6, 12], r"$\delta$", dec_unq_delta)

    m_obs_str = m_obs_nam if m_obs_nam is not None else 'm_{obs}'

    ax = plt.subplot(gs[6:12, 12:18])
    ax.set_title(r"Matched magnitudes ($\delta_{{mag}}<{}$)".format(
        max_mag_delta), fontsize=14)
    plt.xlabel(r'${}$'.format(m_obs_str), fontsize=18)
    plt.ylabel(r'${}$'.format(m_cat), fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='k', linestyle='--', lw=.5,
            zorder=0)
    plt.scatter(m_unq, m_unq_q, marker='o', c='g', s=20, lw=.5,
                edgecolors='k', zorder=1, label="Matched stars")
    msk_in, fit, fit_l, fit_u, x_out, y_out = mag_filter_data
    plt.scatter(
        x_out, y_out, c='r', edgecolors='k', s=10, lw=.2, zorder=1,
        label="Rejected matches")
    oldmin, oldmax = 1000., 0.
    for dd in (m_unq, m_unq_q, x_out, y_out):
        if dd.any():
            xymin = min(oldmin, np.min(dd))
            xymax = max(oldmax, np.max(dd))
    if xymin == xymax:
        xymin, xymax = ax.get_xlim()
    p_x = np.linspace(xymin, xymax, 10)
    p_y = fit(p_x)
    lower = fit_l(p_x)
    upper = fit_u(p_x)
    plt.plot(p_x, p_y, 'b--', zorder=4, label="Matches fit")
    plt.plot(p_x, lower, 'r--', zorder=4, label="Magnitude limits")
    plt.plot(p_x, upper, 'r--', zorder=4)
    plt.plot(
        [xymin, xymax], [xymin, xymax], c='k', ls='--', zorder=4,
        label="1:1 line")
    plt.legend(fontsize=12)
    plt.xlim(xymin, xymax)
    plt.ylim(xymin, xymax)

    ax = plt.subplot(gs[6:12, 18:24])
    ax.set_title("Magnitudes for not matched stars", fontsize=14)
    plt.xlabel(
        r"${}\;$".format(m_obs_str) + r"$|\;{}$".format(m_cat), fontsize=18)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='k', linestyle='--', lw=.5,
            zorder=1)
    if len(m_rjct) > 20 and len(m_rjct_q) > 20:
        dens = True
        plt.ylabel(r'$N\;(norm)$', fontsize=18)
    else:
        dens = False
        plt.ylabel(r'$N$', fontsize=18)
    plt.hist(m_rjct, bins=20, alpha=.5, density=dens, label="Observed")
    plt.hist(m_rjct_q, bins=20, alpha=.5, density=dens, label="Queried")
    plt.legend()

    ax = plt.subplot(gs[12:15, 0:6])
    ax.set_title("Separation between stars (match<{:.1f} [arcsec])".format(
        max_arcsec), fontsize=14)
    plt.xlabel(r'${}$'.format(m_obs_str), fontsize=18)
    plt.ylabel(r'$d\,[arcsec]$', fontsize=18)
    plt.scatter(
        m_unq, match_d2d_all.arcsec, c='g', edgecolors='w', lw=.5,
        label='Match ({})'.format(len(match_d2d_all)))
    plt.scatter(
        m_rjct, no_match_d2d_all.arcsec, c='r', edgecolors='w', lw=.5,
        label='No match ({})'.format(len(no_match_d2d_all)))
    ax.axhline(max_arcsec, color='k', linestyle='--')
    dmean = np.mean(match_d2d_all.arcsec)
    ax.axhline(
        dmean, color='orange', linestyle='--',
        label="Mean={:.3f} arcsec".format(dmean))
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, loc='upper left')
    leg.get_frame().set_alpha(0.5)
    plt.ylim(0., max_arcsec * 2)

    ax = plt.subplot(gs[12:15, 6:12])
    ax.set_title("Distance between matches", fontsize=14)
    plt.hist(match_d2d_all.arcsec, bins=20, density=True)
    # plt.hist(no_match_d2d_all.arcsec, bins=20, density=True)
    plt.xlabel(r'$d\,[arcsec]$')
    plt.ylabel('N (norm)')
    ax.axvline(dmean, color='orange', linestyle='--')
    ax.axvline(max_arcsec, color='k', linestyle='--')

    ax = plt.subplot(gs[12:15, 12:18])
    ax.set_title("Number of stars in different instances", fontsize=14)
    ax.set_xlim(-0.2, 4.2)
    ax.grid(b=True, which='major', color='k', linestyle=':', lw=.5,
            zorder=1)
    x = np.arange(5)
    y = np.array(
        [len(m_obs), len(m_qry), len(m_unq), len(m_rjct), len(m_rjct_q)])
    up = max(y) * .03
    ax.set_ylim(0, max(y) + 4 * up)
    barlst = ax.bar(x, y, align='center', width=0.2, color='g', zorder=4)
    barlst[3].set_color('r')
    barlst[4].set_color('r')
    for xi, yi, l in zip(*[x, y, list(map(str, y))]):
        ax.text(xi - len(l) * .02, yi + up, l,
                bbox=dict(facecolor='w', edgecolor='w', alpha=.5))
    ax.set_xticks(x)
    ax.set_xticklabels(
        ['Observed', 'Queried', 'Match', 'No match (obs)',
         'No match (qrd)'])
    ax.tick_params(axis='x', which='major', labelsize=12)

    #
    fig.tight_layout()
    # Generate output plot file.
    plt.savefig('output/' + clust_name + '.png', dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close("all")
