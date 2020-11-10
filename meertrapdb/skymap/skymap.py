# -*- coding: utf-8 -*-
#
#   2020 Fabian Jankowski
#   A sky exposure map.
#

import logging


class Skymap(object):
    """
    A sky exposure map.
    """

    def __init__(self):
        """
        A sky exposure map.
        """

        self.log = logging.getLogger('meertrapdb.skymap.skymap')

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {
            'skymap': 'XXX'
        }

        info_str = '{0}'.format(info_dict)

        return info_str

    def __str__(self):
        """
        String representation of the object.
        """

        info_str = 'Skymap'

        return info_str

    def load_from_file(self, filename):
        """
        Load sky exposure data from file.

        Parameters
        ----------
        filename: str
            The name of the file.
        """
        pass

    def save_to_file(self, filename):
        """
        Save sky exposure data to file.

        Parameters
        ----------
        filename: str
            The name of the file.
        """
        pass

    def gen_from_database(self):
        """
        Generate the skymap data from the database.
        """
        pass

    def merge(self, other):
        """
        Merge sky exposure data, updating the current data.
        """
        pass

    def __add__(self, other):
        """
        Add two Skymap objects together.

        Returns
        -------
        total: ~Skymap
            The updated total exposure skymap.
        """
        pass
