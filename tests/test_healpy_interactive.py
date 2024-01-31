#
#   Test healpy functionality.
#   2020 Fabian Jankowski
#

import pytest

from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

# astropy generates members dynamically, pylint therefore fails
# disable the corresponding pylint test for now
# pylint: disable=E1101


@pytest.mark.interactive
def test_healpy_functionality():
    print("{0:<6} {1:>10} {2:>12} {3:>12}".format("Expon", "Nside", "Res", "Npix"))
    print("{0:<6} {1:>10} {2:>12} {3:>12}".format("", "", "(arcsec)", "(Gpix)"))

    for expon in range(6, 20):
        iside = 2**expon

        print(
            "{0:<6} {1:>10} {2:>12.1f} {3:>12.1f}".format(
                expon,
                iside,
                60 * hp.nside2resol(iside, arcmin=True),
                hp.nside2npix(iside) * 1e-9,
            )
        )

    npix = hp.nside2npix(2**13)
    print("Number of pixels: {0:.1f} Gpix".format(npix * 1e-9))

    print("Array size in memory")
    for expon in range(10, 15):
        iside = 2**expon
        npix = hp.nside2npix(iside)
        test = np.zeros(npix, dtype=np.float32)

        print("{0:<6} {1:8.1f} GB".format(expon, test.nbytes / 1024**3))

        del test

    # make test map and populate it
    NSIDE = 2**13
    NPIX = hp.nside2npix(NSIDE)

    print(NSIDE, NPIX)

    m = np.full(NPIX, hp.UNSEEN, dtype=np.float32)
    print(m)
    print("Size in memory: {0:.1f} GB".format(m.nbytes / 1024**3))

    coords = SkyCoord(
        ra=["04:37:15.8961737", "05:34:31.973", "08:35:20.61149", "16:44:49.273"],
        dec=["-47:15:09.110714", "+22:00:52.06", "-45:10:34.8751", "-45:59:09.71"],
        unit=(u.hour, u.deg),
        frame="icrs",
    )

    for item in coords:
        # theta [0, 2*pi], phi [0, pi]
        ra_rad = item.ra.radian
        dec_rad = 0.5 * np.pi - item.dec.radian

        print(ra_rad, dec_rad)

        vec = hp.ang2vec(dec_rad, ra_rad)

        mask = hp.query_disc(
            nside=NSIDE,
            vec=vec,
            radius=np.radians(0.58),
            # radius=np.radians(0.5 * 12 / 60.0)
        )

        print(len(mask))

        m[mask] = 1.0

    hp.cartview(
        m,
        badcolor="lightgray",
        cmap="Reds",
        coord=["C"],
        rot=(0, 0, 0),
        title="",
        unit="min",
        xsize=3200,
    )
    # hp.gnomview(m)
    # hp.cartview(m)
    # hp.orthview(m)
    hp.graticule()

    fig = plt.gcf()
    # ax = fig.gca()

    # ax.set_xlabel('RA (hr)')
    # ax.set_ylabel('Dec (deg)')

    # fig.savefig("skymap.pdf", bbox_inches="tight")
    # np.save("skymap.npy", m)

    plt.show()


if __name__ == "__main__":
    pytest.main([__file__])
