#
#   2020 Fabian Jankowski
#

import os

from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from numpy.testing import assert_equal, assert_raises

from skymap import Skymap

# astropy generates members dynamically, pylint therefore fails
# disable the corresponding pylint test for now
# pylint: disable=E1101


def test_creation():
    nside = 2**8
    quantity = "time"
    unit = "min"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m)


def test_addition_success():
    nside = 2**8
    quantity = "time"
    unit = "min"

    m1 = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m1)

    m2 = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m2)

    mtot = m1 + m2
    print(mtot)

    for param in ["arrangement", "coordinate", "dtype", "nside", "quantity", "unit"]:
        assert getattr(mtot, param) == getattr(m1, param)
        assert getattr(mtot, param) == getattr(m2, param)

    assert_equal(mtot.data, m1.data + m2.data)


def test_addition_fail():
    nside = 2**8
    quantity = "time"
    unit = "min"

    m1 = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m1)

    # nside
    m2 = Skymap(nside=2**4, quantity=quantity, unit=unit)
    print(m2)

    with assert_raises(RuntimeError):
        m1 + m2

    # quantity
    m3 = Skymap(nside=nside, quantity="bla", unit=unit)
    print(m3)

    with assert_raises(RuntimeError):
        m1 + m3

    # unit
    m4 = Skymap(nside=nside, quantity=quantity, unit="bla")
    print(m4)

    with assert_raises(RuntimeError):
        m1 + m4

    # nside and quantity
    m5 = Skymap(nside=2**4, quantity="bla", unit=unit)
    print(m5)

    with assert_raises(RuntimeError):
        m1 + m5

    # quantity and unit
    m6 = Skymap(nside=nside, quantity="bla", unit="bla")
    print(m6)

    with assert_raises(RuntimeError):
        m1 + m6

    # nside, quantity and unit
    m7 = Skymap(nside=2**4, quantity="bla", unit="bla")
    print(m7)

    with assert_raises(RuntimeError):
        m1 + m7


def test_private_access():
    nside = 2**8
    quantity = "time"
    unit = "min"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)

    with assert_raises(AttributeError):
        m.__data

    with assert_raises(AttributeError):
        m.data = None


def test_save():
    nside = 2**8
    quantity = "time"
    unit = "min"
    filename = "skymap_save_test.pkl.bz2"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m)

    m.save_to_file(filename)

    if not os.path.isfile(filename):
        raise RuntimeError("File does not exist.")
    else:
        # remove test file
        os.remove(filename)


def test_load():
    nside = 2**8
    quantity = "time"
    unit = "min"
    filename = "skymap_save_test.pkl.bz2"

    m1 = Skymap(nside=nside, quantity=quantity, unit=unit)

    m1.save_to_file(filename)

    # from instance
    m2 = m1.load_from_file(filename)

    for param in ["arrangement", "coordinate", "dtype", "nside", "quantity", "unit"]:
        assert getattr(m2, param) == getattr(m1, param)

    assert_equal(m2.data, m1.data)

    # using the class method
    m3 = Skymap.load_from_file(filename)

    for param in ["arrangement", "coordinate", "dtype", "nside", "quantity", "unit"]:
        assert getattr(m3, param) == getattr(m1, param)

    assert_equal(m3.data, m1.data)


def test_save_to_fits():
    nside = 2**8
    quantity = "time"
    unit = "min"
    filename = "skymap_save_test.fits"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m)

    m.save_to_fits(filename)

    if not os.path.isfile(filename):
        raise RuntimeError("File does not exist.")
    else:
        # remove test file
        os.remove(filename)


def test_size():
    nside = 2**10
    quantity = "time"
    unit = "min"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m.size)


def test_add_exposure():
    nside = 2**10
    quantity = "time"
    unit = "min"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)

    coords = SkyCoord(
        ra=["04:37:15.8961737", "05:34:31.973", "08:35:20.61149", "16:44:49.273"],
        dec=["-47:15:09.110714", "+22:00:52.06", "-45:10:34.8751", "-45:59:09.71"],
        unit=(u.hour, u.deg),
        frame="icrs",
    )

    radius = np.full(len(coords), 0.58)
    length = np.full(len(coords), 10.0)

    m.add_exposure(coords, radius, length)
    print(m)


def test_comments():
    nside = 2**8
    quantity = "time"
    unit = "min"
    comments = ["Hello, test!", "Bla", "blub"]

    m = Skymap(nside=nside, quantity=quantity, unit=unit)

    for item in comments:
        m.add_comment(item)

    assert len(m.comments) == len(comments)
    assert m.comments == comments


def test_save_comments():
    nside = 2**8
    quantity = "time"
    unit = "min"
    filename = "skymap_save_test.fits"
    comments = ["Hello", "Bla", "Test"]

    m = Skymap(nside=nside, quantity=quantity, unit=unit)
    print(m)

    for item in comments:
        m.add_comment(item)

    assert len(m.comments) == len(comments)
    assert m.comments == comments

    m.save_to_fits(filename)

    if not os.path.isfile(filename):
        raise RuntimeError("File does not exist.")
    else:
        # remove test file
        os.remove(filename)


def test_low_map_resolution():
    nside = 2**5
    quantity = "time"
    unit = "min"

    m = Skymap(nside=nside, quantity=quantity, unit=unit)

    coords = SkyCoord(
        ra=["04:37:15.8961737", "05:34:31.973", "08:35:20.61149", "16:44:49.273"],
        dec=["-47:15:09.110714", "+22:00:52.06", "-45:10:34.8751", "-45:59:09.71"],
        unit=(u.hour, u.deg),
        frame="icrs",
    )

    radius = np.full(len(coords), 0.0058)
    length = np.full(len(coords), 10.0)

    with assert_raises(RuntimeError):
        m.add_exposure(coords, radius, length)


if __name__ == "__main__":
    import nose2

    nose2.main()
