# Copyright 2025 Modern Electron, David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .WarpX import warpx

collisions = warpx.get_parametergroup("collisions", collision_names=[])


def new_collision(name):
    result = warpx.get_parametergroup(name)
    if name not in collisions.collision_names:
        collisions.collision_names.append(name)
    return result
