"""Numerical stability tests."""
from starry import Map
import matplotlib.pyplot as pl
import numpy as np
from tqdm import tqdm
cmap = pl.get_cmap('plasma')


def color(l, lmax=8):
    """Return the color for the spherical harmonic degree `l`."""
    return cmap(0.1 + 0.8 * l / lmax)


def impact_param(ax, lmax=8):
    """Test the stability as a function of b."""
    npts = 200
    yo = 0
    cutoff = 2.5e-16
    barr = np.logspace(-5, np.log10(2.0), npts)
    xo = np.sqrt(barr ** 2 - yo ** 2)

    # Double precision
    ylm = Map(lmax)
    ylm.taylor = True

    # Quad precision (~exact)
    ylm128 = Map(lmax)
    ylm128.use_mp = True

    # Compute
    for l in range(0, lmax + 1):

        # Set the constant term
        # so we don't divide by zero when
        # computing the fractional error below
        ylm.reset()
        ylm[0, 0] = 1
        ylm128.reset()
        ylm128[0, 0] = 1

        # Set all odd terms
        for m in range(-l, l + 1):
            # Skip even terms (they're fine!)
            if (l + m) % 2 == 0:
                continue
            ylm[l, m] = 1
            ylm128[l, m] = 1

        ro = 0.1
        # Compute
        flux = np.array(ylm.flux(xo=xo, yo=yo, ro=ro))
        flux128 = np.array(ylm128.flux(xo=xo, yo=yo, ro=ro))
        error = np.abs(flux / flux128 - 1)

        # HACK to make it prettier.
        error[error < cutoff] = cutoff
        # HACK Believe it or not, quadruple precision is not good
        # enough when the impact parameter is smaller than about 0.001
        # at high values of l! While our taylor expansions are fine, the
        # high precision comparison flux is wrong, so we get a spurious
        # increase in the error towards smaller b. Let's just trim it
        # for simplicity.
        error[xo < 0.002] = cutoff

        ax.plot(barr, error, ls='-',
                label=r'$l = %d$' % l, color=color(l, lmax))

    ax.legend(loc="upper left", ncol=3)
    ax.axhline(1e-3, color='k', ls='--', alpha=0.75, lw=0.5)
    ax.axhline(1e-6, color='k', ls='--', alpha=0.75, lw=0.5)
    ax.axhline(1e-9, color='k', ls='--', alpha=0.75, lw=0.5)
    ax.annotate("ppt", xy=(1e-5, 1e-3), xycoords="data", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", alpha=0.75)
    ax.annotate("ppm", xy=(1e-5, 1e-6), xycoords="data", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", alpha=0.75)
    ax.annotate("ppb", xy=(1e-5, 1e-9), xycoords="data", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", alpha=0.75)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(5e-17, 20.)
    ax.set_xlim(1e-5, np.log10(2.0))
    ax.set_xlabel("Impact parameter", fontsize=16)
    ax.set_ylabel("Fractional error", fontsize=16)


def occultor_radius(ax, lmax=8):
    """Test the stability as a function of occultor radius."""
    # Knobs
    eps = 1e-6
    npts = 50
    yo = 0
    rarr = np.logspace(-3, 3, npts)

    # Double precision
    ylm = Map(lmax)

    # Quad precision (~exact)
    ylm128 = Map(lmax)
    ylm128.use_mp = True

    # Loop over the degrees
    for l in tqdm(range(lmax + 1)):
        ylm.reset()
        ylm128.reset()
        # Set the coefficients for all orders
        for m in range(-l, l + 1):
            ylm[l, m] = 1
            ylm128[l, m] = 1
        # Occultor radius loop
        error = np.zeros_like(rarr)
        for i, ro in enumerate(rarr):
            xo0 = 0.5 * ((ro + 1) + np.abs(ro - 1))
            xo = np.linspace(xo0 - 25 * eps, xo0 + 25 * eps, 50)
            flux = np.array(ylm.flux(xo=xo, yo=yo, ro=ro))
            flux128 = np.array(ylm128.flux(xo=xo, yo=yo, ro=ro))
            error[i] = np.max(np.abs((flux / flux128 - 1)))
        ax.plot(rarr, error, '-', color=color(l, lmax),
                label=r"$l=%d$" % l)

    ax.legend(loc="upper left", ncol=3)
    ax.axhline(1e-3, color='k', ls='--', alpha=0.75, lw=0.5)
    ax.axhline(1e-6, color='k', ls='--', alpha=0.75, lw=0.5)
    ax.axhline(1e-9, color='k', ls='--', alpha=0.75, lw=0.5)
    ax.annotate("ppt", xy=(1e-3, 1e-3), xycoords="data", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", alpha=0.75)
    ax.annotate("ppm", xy=(1e-3, 1e-6), xycoords="data", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", alpha=0.75)
    ax.annotate("ppb", xy=(1e-3, 1e-9), xycoords="data", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", alpha=0.75)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(5e-17, 20.)
    ax.set_xlim(1e-3, 1e3)
    ax.set_xlabel("Occultor radius", fontsize=16)
    ax.set_ylabel("Fractional error", fontsize=16)



if __name__ == "__main__":

    # Set up
    fig, ax = pl.subplots(2, figsize=(9, 8))
    impact_param(ax[0])
    occultor_radius(ax[1])
    fig.savefig("stability.pdf", bbox_inches='tight')
