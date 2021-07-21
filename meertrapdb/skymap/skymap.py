# -*- coding: utf-8 -*-
#
#   2020 Fabian Jankowski
#   A sky exposure map.
#

import bz2
import copy
import logging
import pickle
import os.path

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.time import Time
import healpy as hp
from healpy.visufunc import (projscatter, projtext)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

# astropy.units generates members dynamically, pylint therefore fails
# disable the corresponding pylint test for now
# pylint: disable=E1101


class Skymap(object):
    """
    A sky exposure map.
    """

    name = 'Skymap'

    def __init__(self, nside, quantity, unit):
        """
        A sky exposure map.

        Parameters
        ----------
        nside: int
            The HEALPIX `nside` parameter for the map.
        quantity: str
            The name of the quantity that is stored in the map (e.g. 'time').
        unit: str
            The unit of the map data (e.g. 'min').
        """

        self.__arrangement = 'ring'
        self.__comments = []
        self.__coordinate = 'icrs'
        self.__dtype = np.float32
        self.__exposures = 0
        self.__log = logging.getLogger('meertrapdb.skymap.skymap')
        self.__nside = nside
        self.__npix = hp.nside2npix(self.__nside)
        self.__quantity = quantity
        self.__unit = unit

        # create empty map
        #self.__data = np.full(self.__npix, hp.UNSEEN, dtype=self.__dtype)
        self.__data = np.zeros(self.__npix, dtype=self.__dtype)

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {
            'coordinate': self.coordinate,
            'dtype': self.dtype,
            'nside': self.nside,
            'quantity': self.quantity,
            'unit': self.unit,
            'exposures': self.exposures,
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

        info_str = '{0}: {1}'.format(self.name, repr(self))

        return info_str

    @classmethod
    def load_from_file(cls, filename):
        """
        Load Skymap from file.

        Parameters
        ----------
        filename: str
            The name of the file.

        Returns
        -------
        loaded: ~Skymap
            The loaded Skymap object.

        Raises
        ------
        RuntimeError
            If the file does not exist.
        """

        if not os.path.isfile(filename):
            raise RuntimeError('The file does not exist: {0}'.format(filename))

        with bz2.open(filename, 'rb') as fd:
            loaded = pickle.load(fd)

        return loaded

    def save_to_file(self, filename):
        """
        Save Skymap to file.

        Parameters
        ----------
        filename: str
            The name of the file.
        """

        with bz2.open(filename, 'wb') as fd:
            pickle.dump(self, fd, protocol=pickle.DEFAULT_PROTOCOL)

        self.__log.info('Saved skymap to file: {0}'.format(filename))

    def save_to_fits(self, filename):
        """
        Save Skymap to a FITS file.

        Parameters
        ----------
        filename: str
            The name of the file.

        Raises
        ------
        NotImplementedError if the coordinate system is unknown.
        """

        if self.arrangement == 'nest':
            nest = True
        else:
            nest = False

        if self.coordinate == 'icrs':
            coord = 'C'
        elif self.coordinate == 'galactic':
            coord = 'G'
        else:
            raise NotImplementedError('Coordinate system unknown: {0}'.format(self.coordinate))

        extra_header = [
            ('COMMENT', 'Generated using Skymap'),
            ('COMMENT', 'Created on {0}'.format(Time.now().iso))
        ]

        for item in self.comments:
            extra_header.append(
                ('COMMENT', item)
            )

        hp.write_map(
            filename=filename,
            m=self.data,
            nest=nest,
            coord=coord,
            column_names=[self.quantity],
            column_units=[self.unit],
            extra_header=extra_header,
            dtype=self.dtype
        )

        self.__log.info('Saved skymap to FITS file: {0}'.format(filename))

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

        if self.arrangement == other.arrangement \
        and self.coordinate == other.coordinate \
        and self.dtype == other.dtype \
        and self.nside == other.nside \
        and self.quantity == other.quantity \
        and self.unit == other.unit:
            pass

        else:
            raise RuntimeError('The Skymap objects are incompatible.')

        total = copy.deepcopy(self)
        total.__data = self.data + other.data

        return total

    def __radd__(self, other):
        """
        Reverse addition of Skymap objects.
        """

        return self.__add__(other)

    @property
    def arrangement(self):
        """
        The HEALPIX pixel arrangement of the sky map.
        """

        return self.__arrangement

    @property
    def comments(self):
        """
        The comments on the sky map.
        """

        return self.__comments

    def add_comment(self, comment):
        """
        Add a comment to the sky map.

        Parameters
        ----------
        comment: str
            The comment to add to the sky map.
        """

        self.__comments.append(comment)

    @property
    def coordinate(self):
        """
        The coordinate frame of the sky map.
        """

        return self.__coordinate

    @property
    def data(self):
        """
        The sky map data array.
        """

        return self.__data

    @property
    def dtype(self):
        """
        The data type of the sky map.
        """

        return self.__dtype

    @property
    def exposures(self):
        """
        The number of exposures added to the sky map.
        """

        return self.__exposures

    @property
    def nside(self):
        """
        The HEALPIX `nside` parameter.
        """

        return self.__nside

    @property
    def quantity(self):
        """
        The quantity stored in the sky map.
        """

        return self.__quantity

    @property
    def unit(self):
        """
        The unit of the sky map.
        """

        return self.__unit

    @property
    def size(self):
        """
        The size of the sky map in GB.
        """

        return self.data.nbytes / 1024**3

    @property
    def min(self):
        """
        The minimum value of the sky map.
        """

        return np.min(self.data)

    @property
    def mean(self):
        """
        The mean value of the sky map.
        """

        return np.mean(self.data)

    @property
    def median(self):
        """
        The median value of the sky map.
        """

        return np.median(self.data)

    @property
    def max(self):
        """
        The maximum value of the sky map.
        """

        return np.max(self.data)

    @property
    def sum(self):
        """
        The sum of the sky map.
        """

        return np.sum(self.data, dtype=np.float128)

    def add_exposure(self, coords, radii, lengths):
        """
        Add exposure to the sky map.

        Parameters
        ----------
        coords: ~astropy.SkyCoord
            The centres of the beams.
        radii: ~np.array of float
            The radii of the beams in degrees.
        lengths: ~np.array of float
            The exposure lengths in `unit`.
        """

        for item, radius, length in zip(coords, radii, lengths):
            # theta [0, 2*pi], phi [0, pi]
            ra_rad = item.ra.radian
            dec_rad = 0.5 * np.pi - item.dec.radian

            self.__log.debug('RA, Dec: {0}, {1}'.format(ra_rad, dec_rad))

            vec = hp.ang2vec(dec_rad, ra_rad)

            mask = hp.query_disc(
                nside=self.nside,
                vec=vec,
                radius=np.radians(radius)
            )

            self.__log.debug('Number of HEAL pixels: {0}'.format(len(mask)))

            self.__data[mask] = self.data[mask] + length
            self.__exposures = self.exposures + 1

    def query(self, coords, radii):
        """
        Query the exposure for given coordinates.

        Parameters
        ----------
        coords: ~astropy.SkyCoord
            The centres of the beams.
        radii: ~np.array of float
            The radii of the beams in degrees.

        Returns
        -------
        exposures: list of float
            The exposures in units of the sky map.
        """

        exposures = []

        for item, radius in zip(coords, radii):
            # theta [0, 2*pi], phi [0, pi]
            ra_rad = item.ra.radian
            dec_rad = 0.5 * np.pi - item.dec.radian

            self.__log.debug('RA, Dec: {0}, {1}'.format(ra_rad, dec_rad))

            vec = hp.ang2vec(dec_rad, ra_rad)

            mask = hp.query_disc(
                nside=self.nside,
                vec=vec,
                radius=np.radians(radius)
            )

            self.__log.debug('Number of HEAL pixels: {0}'.format(len(mask)))

            sel = np.copy(self.data[mask])
            sel = sel[sel > 0]

            if len(sel) > 0:
                exposure = np.median(sel)
            else:
                exposure = np.nan

            exposures.append(exposure)

        return exposures

    def show(self, coordinates='equatorial', sources=None, shownames=False):
        """
        Visualise the Skymap exposure data and sources.

        Parameters
        ----------
        coordinates: str (default: equatorial)
            The coordinate system to use (equatorial, galactic).
        sources: pandas.DataFrame
            Source positions to overplot (default: None).
        shownames: bool (default: False)
            Whether to show the names of the sources on the sky map.

        Raises
        ------
        NotImplementedError
            In case the `coordinates` coordinate system is not implemented.
        """

        if coordinates == 'equatorial':
            coord = ['C']
            rot = (180.0, 0, 0)
        elif coordinates == 'galactic':
            coord = ['C', 'G']
            rot = (0, 0, 0)
        else:
            raise NotImplementedError('Coordinate system is not available.')

        # mask all empty areas
        masked = np.copy(self.data)
        masked[masked < 0.01] = np.nan

        cmap = copy.copy(plt.get_cmap('Reds'))
        cmap.set_under('white')

        fig = plt.figure()

        hp.mollview(
            masked,
            badcolor='white',
            bgcolor='white',
            cmap=cmap,
            coord=coord,
            fig=fig.number,
            norm=LogNorm(vmin=np.nanmin(masked), vmax=np.nanmax(masked)),
            rot=rot,
            title='',
            unit=self.unit,
            xsize=6400
        )

        hp.graticule()

        # add markers
        if coordinates == 'equatorial':
            ras = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22]

            coords = SkyCoord(
                ra=ras,
                dec=[5 for _ in range(len(ras))],
                unit=(units.hourangle, units.deg),
                frame='icrs'
            )

            for i in range(len(coords)):
                projtext(
                    coords[i].ra.deg,
                    coords[i].dec.deg,
                    s='{0} h'.format(ras[i]),
                    lonlat=True,
                    coord='C',
                    clip_on=True,
                    color='black',
                    fontfamily='serif',
                    fontsize='medium',
                    horizontalalignment='center',
                    verticalalignment='center',
                    snap=True,
                    zorder=3
                )

            decs = [-60, -30, 30, 60]

            coords = SkyCoord(
                ra=[5.0 for _ in range(len(decs))],
                dec=decs,
                unit=(units.hourangle, units.deg),
                frame='icrs'
            )

            for i in range(len(coords)):
                projtext(
                    coords[i].ra.deg,
                    coords[i].dec.deg,
                    s=r'${0:+}\degree$'.format(decs[i]),
                    lonlat=True,
                    coord='C',
                    clip_on=True,
                    color='black',
                    fontfamily='serif',
                    fontsize='medium',
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    snap=True,
                    zorder=3
                )
        elif coordinates == 'galactic':
            gl = [30, 60, 90, 120, 150, 210, 240, 270, 300, 330]

            coords = SkyCoord(
                l=gl,
                b=[15.0 for _ in range(len(gl))],
                unit=(units.deg, units.deg),
                frame='galactic'
            )

            for i in range(len(coords)):
                projtext(
                    coords[i].l.deg,
                    coords[i].b.deg,
                    s=r'${0}\degree$'.format(gl[i]),
                    lonlat=True,
                    coord='G',
                    clip_on=True,
                    color='black',
                    fontfamily='serif',
                    fontsize='medium',
                    horizontalalignment='center',
                    verticalalignment='center',
                    snap=True,
                    zorder=3
                )

            gb = [-60, -30, 30, 60]

            coords = SkyCoord(
                l=[105.0 for _ in range(len(gb))],
                b=gb,
                unit=(units.deg, units.deg),
                frame='galactic'
            )

            for i in range(len(coords)):
                projtext(
                    coords[i].l.deg,
                    coords[i].b.deg,
                    s=r'${0:+}\degree$'.format(gb[i]),
                    lonlat=True,
                    coord='G',
                    clip_on=True,
                    color='black',
                    fontfamily='serif',
                    fontsize='medium',
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    snap=True,
                    zorder=3
                )

        # highlight sources
        if sources is not None:
            types = np.unique(sources['type'])
            colors = ['tab:olive', 'tab:blue', 'tab:purple', 'tab:green']

            for i, item in enumerate(types):
                mask = (sources['type'] == item)
                sel = sources[mask]

                coords = SkyCoord(
                    ra=sel['ra'],
                    dec=sel['dec'],
                    unit=(units.hourangle, units.deg),
                    frame='icrs'
                )

                color = colors[i % len(colors)]

                # no need to convert the coordinates manually
                # healpy does that for us automatically
                projscatter(
                    coords.ra.deg,
                    coords.dec.deg,
                    lonlat=True,
                    coord='C',
                    marker='*',
                    facecolor=color,
                    edgecolor='black',
                    lw=0.5,
                    s=50,
                    zorder=5
                )

                # show the source names
                if shownames:
                    for i in range(len(sel)):
                        projtext(
                            coords[i].ra.deg,
                            coords[i].dec.deg,
                            s=sel['name'].iloc[i],
                            lonlat=True,
                            coord='C',
                            clip_on=True,
                            color='black',
                            fontfamily='sans-serif',
                            fontsize='xx-small',
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            snap=True,
                            zorder=6
                        )

        fig.savefig(
            'skymap_{0}.png'.format(coordinates),
            bbox_inches='tight',
            dpi=300
        )

        fig.savefig(
            'skymap_{0}.pdf'.format(coordinates),
            bbox_inches='tight',
            dpi=300
        )

        plt.draw()

    def show_interactive(self):
        """
        Interactively visualise the Skymap exposure data.
        """

        hp.mollzoom(
            self.data,
            cmap='Reds',
            coord=['C'],
            norm=LogNorm(),
            rot=(0, 0, 0),
            title='',
            unit=self.unit,
            xsize=6400
        )

        hp.graticule()

        plt.show()
