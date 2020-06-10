# Copyright 2017-2020 Andrew Myers, David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

class Particles(Bucket):
    """
    This keeps a dictionary of the species added.
    """
    def __init__(self):
        Bucket.__init__(self, 'particles', _species_dict={})

    def new_species(self, name, **defaults):
        if name is None:
            name = 'species{}'.format(len(self._species_dict))
        try:
            result = self._species_dict[name]
        except KeyError:
            result = Bucket(name, **defaults)
            self._species_dict[name] = result
            self.nspecies = len(self._species_dict)
            self.species_names = self._species_dict.keys()
        return result


particles = Particles()

