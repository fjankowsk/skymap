# -*- coding: utf-8 -*-
#
#   2020 Fabian Jankowski
#   A multi-beam tiling.
#

import logging


class Tiling(object):
    """
    A multi-beam tiling.
    """

    def __init__(self):
        """
        A multi-beam tiling.
        """

        self.log = logging.getLogger('meertrapdb.skymap.tiling')

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {
            'tiling': 'XXX'
        }

        info_str = '{0}'.format(info_dict)

        return info_str

    def __str__(self):
        """
        String representation of the object.
        """

        info_str = 'Tiling'

        return info_str
