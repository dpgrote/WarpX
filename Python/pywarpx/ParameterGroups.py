# Copyright 2025 David Grote
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Collisions import collisions  # noqa
from .Constants import my_constants  # noqa
from .Diagnostics import diagnostics, reduced_diagnostics  # noqa
from .Lasers import lasers  # noqa
from .Particles import particles  # noqa
from .WarpX import warpx

algo = warpx.get_parametergroup("algo")
amr = warpx.get_parametergroup("amr")
amrex = warpx.get_parametergroup("amrex")
boundary = warpx.get_parametergroup("boundary")
eb2 = warpx.get_parametergroup("eb2")
external_vector_potential = warpx.get_parametergroup("external_vector_potential")
geometry = warpx.get_parametergroup("geometry")
hybridpicmodel = warpx.get_parametergroup("hybridpicmodel")
interpolation = warpx.get_parametergroup("interpolation")
psatd = warpx.get_parametergroup("psatd")
