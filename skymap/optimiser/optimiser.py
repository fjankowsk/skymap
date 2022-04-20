#
#   2020 Fabian Jankowski
#   Optimise sky exposure using various methods.
#

import copy
import logging


class Optimiser(object):
    """
    A sky exposure optimiser.
    """

    name = "Optimiser"

    # available methods
    __optimisation_methods = ["watershed"]

    def __init__(self):
        """
        A sky exposure optimiser.
        """

        self.__skymap = None
        self.__log = logging.getLogger("skymap.optimiser.optimiser")

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {"optimiser": "XXX"}

        info_str = "{0}".format(info_dict)

        return info_str

    def __str__(self):
        """
        String representation of the object.
        """

        info_str = "{0}: {1}".format(self.name, repr(self))

        return info_str

    def load_skymap(self, skymap):
        """
        Load skymap exposure data.
        """

        self.__skymap = copy.deepcopy(skymap)

    def optimise(self, pointing, params, method="watershed"):
        """
        Run the sky exposure optimisation.

        Parameters
        ----------
        pointing: ~astropy.SkyCoord
            The pointing position in equatorial coordinates to optimise for.
        params: dict
            The sub-array, requested beam and LST parameters.
        method: str
            The optimisation method to use.

        Raises
        ------
        NotImplementedError
            If the given optimisation `method` is not implemented.
        """

        if method not in self.__optimisation_methods:
            raise NotImplementedError(
                "Optimisation method not implemented: {0}".format(method)
            )

    @property
    def beam_placement(self):
        """
        Retrieve the optimal beam placement after optimisation was run.

        Returns
        -------
        tiling: ~Tiling
            The optimum tiling for the given pointing.
        """
        pass
