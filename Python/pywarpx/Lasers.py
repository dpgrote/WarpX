# Copyright 2019-2020 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

class Lasers(Bucket):
    """
    This keeps a dictionary of the lasers added.
    """
    def __init__(self):
        Bucket.__init__(self, 'lasers', _lasers_dict={})

    def new_laser(self, name):
        if name is None:
            name = 'laser{}'.format(len(self._lasers_dict))
        try:
            result = self._lasers_dict[name]
        except KeyError:
            result = Bucket(name)
            self._lasers_dict[name] = result
            self.nlasers = len(self._lasers_dict)
        return result

lasers = Lasers()

