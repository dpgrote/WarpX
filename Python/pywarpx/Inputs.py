# Copyright 2019-2020 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket
from .WarpX import WarpX
from .Particles import Particles
from .Lasers import Lasers
from .Diagnostics import Diagnostics
from .Constants import Constants

class Inputs(object):

    def __init__(self):
        self.warpx = WarpX()
        self.amr = Bucket('amr')
        self.geometry = Bucket('geometry')
        self.algo = Bucket('algo')
        self.interpolation = Bucket('interpolation')
        self.particles = Particles()
        self.psatd = Bucket('psatd')
        self.lasers = Lasers()
        self.diagnostics = Diagnostics()
        self.my_constants = Constants()

    def create_argv_list(self):
        argv = []
        argv += self.warpx.attrlist()
        argv += self.my_constants.attrlist()
        argv += self.amr.attrlist()
        argv += self.geometry.attrlist()
        argv += self.algo.attrlist()
        argv += self.interpolation.attrlist()
        argv += self.psatd.attrlist()

        argv += self.particles.attrlist()
        for specie in self.particles._species_dict.values():
            argv += specie.attrlist()

        argv += self.lasers.attrlist()
        for laser in self.lasers._lasers_dict.values():
            argv += laser.attrlist()

        argv += self.diagnostics.attrlist()
        for diagnostic in self.diagnostics._diagnostics_dict.values():
            diagnostic.species = diagnostic._species_dict.keys()
            argv += diagnostic.attrlist()
            for species_diagnostic in diagnostic._species_dict.values():
                argv += species_diagnostic.attrlist()

        return argv

    def write_inputs(self, filename='inputs', **kw):
        argv = self.create_argv_list()
        with open(filename, 'w') as ff:

            for k, v in kw.items():
                if v is not None:
                    ff.write('{0} = {1}\n'.format(k, v))

            for arg in argv:
                ff.write('{0}\n'.format(arg))

