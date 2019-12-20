
import numpy as np
import matplotlib.pyplot as plt


def main(x, y, mag_lim=1.5):
    """
    Predicts the values for a best fit between arrays x and y, defines
    outlier bands around this fit line.
    """

    # To prevent stars in the lower magnitude limit from biasing the fit,
    # we bin the magnitude range and estimate the fit using the median of
    # each range.
    xymin, xymax = max(x.min(), y.min()), min(x.max(), y.max())
    bins = np.linspace(xymin, xymax, 10)
    a1 = np.array(list(zip(*[bins[0::2], bins[1::2]])))
    a2 = np.array(list(zip(*[bins[1::2], bins[2::2]])))
    bb = [item for t in zip(a1, a2) for item in t]
    xy_mean = []
    for l1, l2 in bb:
        msk = (x > l1) & (x < l2)
        medians = np.median([x[msk], y[msk]], axis=1)
        if not np.isnan(medians).any():
            xy_mean.append(medians)

    xm, ym = np.array(xy_mean).T
    if xm.size < 2:
        msk = ~np.isnan(np.array(x)) & ~np.isnan(np.array(y))
        xm, ym = x[msk], y[msk]

    z = np.polyfit(xm, ym, 1)
    fit = np.poly1d(z)

    # polyfit coefficients in 'z' are in the form: y = z[1] + z[0]*x
    # Re-write equation as: a*x + b*y + c = 0, this is:
    # -z[1]*x + 1*y + (-z[0]) = 0, and store these three coefficients in 'z0'
    z0 = (-z[0], 1, -z[1])
    # Obtain the point-line distance for all (x,y) pairs.
    dist = abs(x * z0[0] + y * z0[1] + z0[2]) / np.sqrt(z0[0]**2 + z0[1]**2)

    # Fill masked values with small mag distance value. These stars are kept,
    # since there is no magnitude information to reject them.
    dist = dist.filled(0.)

    # In-out mask.
    msk_in = dist <= mag_lim
    x_out, y_out = x[~msk_in], y[~msk_in]

    # Transform the perpendicular distance to a distance in y, to define the
    # lower|upper bands.
    d_lim = mag_lim * z[0] * np.sqrt(1 + (1. / z[0]**2))

    zl = (z[0], z[1] - d_lim)
    fit_l = np.poly1d(zl)

    zu = (z[0], z[1] + d_lim)
    fit_u = np.poly1d(zu)

    return [msk_in.data, fit, fit_l, fit_u, x_out, y_out]


if __name__ == '__main__':
    # Define example data
    x = np.linspace(1, 15, 200)
    y = x * 4 + 2.5
    x = x + np.random.random_sample(size=x.shape) * 20
    y = y + np.random.random_sample(size=x.shape) * 20

    # Fit x to y
    msk_in, fit, fit_l, fit_u, x_out, y_out = main(x, y)
    x_in, y_in = x[msk_in], y[msk_in]

    p_x = np.linspace(x.min(), x.max(), 10)
    p_y = fit(p_x)
    lower = fit_l(p_x)
    upper = fit_u(p_x)

    plt.scatter(x_in, y_in, c='g', s=10)
    plt.scatter(x_out, y_out, c='grey', s=5)
    plt.plot(p_x, p_y, 'b-')
    plt.plot(p_x, lower, 'r--')
    plt.plot(p_x, upper, 'r--')
    plt.show()
