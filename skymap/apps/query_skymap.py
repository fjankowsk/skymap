#
#   2022 Fabian Jankowski
#   Query a Skymap
#

import argparse

import astropy.units as units
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import pandas as pd

from skymap import Skymap


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    parser = argparse.ArgumentParser(
        description="Query Skymap data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "filename", type=str, help="The name of the input sky map file."
    )

    parser.add_argument(
        "ra",
        type=str,
        help="The right ascension of the coordinate to query.",
    )

    parser.add_argument(
        "dec",
        type=str,
        help="The declination of the coordinate to query.",
    )

    parser.add_argument(
        "-p",
        "--plot",
        dest="plot",
        action="store_true",
        default=False,
        help="Plot the coordinate on the sky map.",
    )

    args = parser.parse_args()

    return args


#
# MAIN
#


def main():
    args = parse_args()

    m = Skymap.load_from_file(args.filename)

    coords = SkyCoord(
        ra=[args.ra],
        dec=[args.dec],
        unit=(units.hourangle, units.deg),
        frame="icrs",
    )

    exposure = m.query(coords, [0.1 for _ in range(len(coords))])
    print("{0} {1}".format(exposure, m.unit))

    df_sources = pd.DataFrame(
        {"ra": [args.ra], "dec": [args.dec], "name": [""], "type": [""]}
    )

    if args.plot:
        m.show(coordinates="galactic", sources=df_sources)
        plt.show()


if __name__ == "__main__":
    main()
