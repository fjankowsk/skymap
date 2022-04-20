#
#   2020 Fabian Jankowski
#   A multi-beam tiling.
#

import logging


class Tiling(object):
    """
    A multi-beam tiling.
    """

    name = "Tiling"

    def __init__(self):
        """
        A multi-beam tiling.
        """

        self.__log = logging.getLogger("skymap.tiling")
        self.__beams = []

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {"nbeam": self.nbeam}

        info_str = "{0}".format(info_dict)

        return info_str

    def __str__(self):
        """
        String representation of the object.
        """

        info_str = "{0}: {1}".format(self.name, repr(self))

        return info_str

    @property
    def nbeam(self):
        """
        The number of beams in the tiling.
        """

        return len(self.__beams)

    def show(self):
        """
        Display the tiling visually.
        """

        pass
