# Copyright 2019-2025 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .WarpX import warpx

lasers = warpx.get_parametergroup("lasers", names=[])


def new_laser(name):
    result = warpx.get_parametergroup(name)
    if name not in lasers.names:
        lasers.names.append(name)
    return result
