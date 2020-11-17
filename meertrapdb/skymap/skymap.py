# -*- coding: utf-8 -*-
#
#   2020 Fabian Jankowski
#   A sky exposure map.
#

import copy
import logging
import os.path

import healpy as hp
import numpy as np


class Skymap(object):
    """
    A sky exposure map.
    """

    def __init__(self, nside, unit):
        """
        A sky exposure map.

        Parameters
        ----------
        nside: int
            The HEALPIX `nside` parameter for the map.
        unit: str
            The unit of the map data.
        """

        self._arrangement = 'ring'
        self._coordinate = 'icrs'
        self._dtype = np.float32
        self.log = logging.getLogger('meertrapdb.skymap.skymap')
        self._nside = nside
        self._npix = hp.nside2npix(self._nside)
        self._unit = unit

        # create empty map
        #self._data = np.full(self._npix, hp.UNSEEN, dtype=self._dtype)
        self._data = np.zeros(self._npix, dtype=self._dtype)

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {
            'coordinate': self._coordinate,
            'dtype': self._dtype,
            'nside': self._nside,
            'unit': self._unit,
            'min': self.min,
            'mean': self.mean,
            'max': self.max,
            'sum': self.sum
        }

        info_str = '{0}'.format(info_dict)

        return info_str

    def __str__(self):
        """
        String representation of the object.
        """

        info_str = 'Skymap, {0}'.format(repr(self))

        return info_str

    def load_from_file(self, filename):
        """
        Load sky exposure data from file.

        Parameters
        ----------
        filename: str
            The name of the file.

        Raises
        ------
        RuntimeError
            If the file does not exist.
        """

        if not os.path.isfile(filename):
            raise RuntimeError('The file does not exist: {0}'.format(filename))

        # XXX: check that meta parameters match
        self._data = np.load(filename)

    def save_to_file(self, filename):
        """
        Save sky exposure data to file.

        Parameters
        ----------
        filename: str
            The name of the file.
        """

        np.save(filename, self._data)

    def gen_from_database(self):
        """
        Generate the skymap data from the database.
        """
        pass

    def __add__(self, other):
        """
        Add two Skymap objects together, adding the exposure data.

        Parameters
        ----------
        other: ~Skymap
            The `other` Skymap object to add.

        Returns
        -------
        total: ~Skymap
            The updated total exposure skymap.

        Raises
        ------
        RuntimeError
            In case of a mismatch between the Skymap objects.
        """

        if self._arrangement == other._arrangement \
        and self._coordinate == other._coordinate \
        and self._dtype == other._dtype \
        and self._nside == other._nside \
        and self._unit == other._unit:
            pass

        else:
            raise RuntimeError('The Skymap objects are incompatible.')

        total = copy.deepcopy(self)
        total._data = self._data + other._data

        return total

    def __radd__(self, other):
        """
        Reverse addition of Skymap objects.
        """

        return self.__add__(other)

    @property
    def size(self):
        """
        The size of the sky map in GB.
        """

        return self._data.nbytes / 1024**3

    @property
    def min(self):
        """
        The minimum value of the sky map.
        """

        return np.min(self._data)

    @property
    def mean(self):
        """
        The mean value of the sky map.
        """

        return np.mean(self._data)

    @property
    def median(self):
        """
        The median value of the sky map.
        """

        return np.median(self._data)

    @property
    def max(self):
        """
        The maximum value of the sky map.
        """

        return np.max(self._data)

    @property
    def sum(self):
        """
        The sum of the sky map.
        """

        return np.sum(self._data, dtype=np.float128)

    def add_exposure(self, coords, radii, lengths):
        """
        Add exposure to the sky map.

        Parameters
        ----------
        coords: ~astropy.SkyCoord
            The centre of the beam.
        radii: ~np.array of float
            The radius of the beam in degrees.
        lengths: ~np.array of float
            The exposure length in `unit`.
        """

        for item, radius, length in zip(coords, radii, lengths):
            # theta [0, 2*pi], phi [0, pi]
            ra_rad = item.ra.radian
            dec_rad = 0.5 * np.pi - item.dec.radian

            self.log.info('RA, Dec: {0}, {1}'.format(ra_rad, dec_rad))

            vec = hp.ang2vec(dec_rad, ra_rad)

            mask = hp.query_disc(
                nside=self._nside,
                vec=vec,
                radius=np.radians(radius)
            )

            self.log.info('Number of HEAL pixels: {0}'.format(len(mask)))

            self._data[mask] = self._data[mask] + length

    def show(self):
        """
        """

        import matplotlib.pyplot as plt

        hp.mollzoom(
            self._data,
            #badcolor='lightgray',
            cmap='Reds',
            coord=['C'],
            rot=(0, 0, 0),
            title='',
            unit=self._unit,
            xsize=1600
        )

        #hp.gnomview(m)
        #hp.cartview(m)
        #hp.orthview(m)
        hp.graticule()

        plt.show()
