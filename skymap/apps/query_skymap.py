#
#   2022 Fabian Jankowski
#   Query a Skymap
#

import astropy.units as units
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

from skymap import Skymap

m = Skymap.load_from_file("skymap_coherent.pkl")

coords = SkyCoord(
    ra=["08:35:20.61149"],
    dec=["-45:10:34.8751"],
    unit=(units.hourangle, units.deg),
    frame="icrs",
)

exposure = m.query(coords, [0.1 for _ in range(len(coords))])
print(exposure)

m.show()

plt.show()
