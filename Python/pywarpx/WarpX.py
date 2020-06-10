# Copyright 2016-2020 Andrew Myers, David Grote, Maxence Thevenet
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket
from .Geometry import geometry

class WarpX(Bucket):
    """
    A Python wrapper for the WarpX C++ class
    """
    def __init__(self):
        Bucket.__init__(self, 'warpx')

    def init(self, inputs):
        argv = ['warpx'] + inputs.create_argv_list()
        # --- These values are needed by the import of wx
        geometry.prob_lo = inputs.geometry.prob_lo
        geometry.coord_sys = inputs.geometry.coord_sys
        from . import wx
        wx.initialize(argv)

    def evolve(self, nsteps=-1):
        from . import wx
        wx.evolve(nsteps)

    def finalize(self, finalize_mpi=1):
        from . import wx
        wx.finalize(finalize_mpi)

    def getProbLo(self, direction):
        from . import wx
        return wx.libwarpx.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        from . import wx
        return wx.libwarpx.warpx_getProbHi(direction)

warpx = WarpX()

