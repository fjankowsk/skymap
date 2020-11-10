# -*- coding: utf-8 -*-
#
#   2020 Fabian Jankowski
#   Optimise sky exposure using various methods.
#

import logging


class Optimiser(object):
    """
    A sky exposure optimiser.
    """

    optimisation_methods = ['watershed']

    def __init__(self):
        """
        A sky exposure optimiser.
        """

        self._skymap = None
        self.log = logging.getLogger('meertrapdb.optimiser.optimiser')

    def __repr__(self):
        """
        Representation of the object.
        """

        info_dict = {
            'optimiser': 'XXX'
        }

        info_str = '{0}'.format(info_dict)

        return info_str

    def __str__(self):
        """
        String representation of the object.
        """

        info_str = 'Optimiser'

        return info_str

    def load_skymap(self, skymap):
        """
        Load skymap exposure data.
        """

        self._skymap = skymap
    
    def optimise(self, pointing, params, method='watershed'):
        """
        Run the sky exposure optimisation.

        Parameters
        ----------
        pointing: ~astropy.SkyCoord
            The pointing position in equatorial coordinates to optimise for.
        params: dict
            The sub-array configuration and LST parameters.
        method: str
            The optimisation method to use.

        Raises
        ------
        NotImplementedError
            If the given optimisation `method` is not implemented.
        """

        if method not in self.optimisation_methods:
            raise NotImplementedError('Optimisation method not implemented: {0}'.format(method))

    def get_beam_placement(self):
        """
        Retrieve the optimal beam placement after optimisation was run.

        Returns
        -------
        tiling: ~Tiling
            The optimum tiling for the given pointing.
        """
        pass
