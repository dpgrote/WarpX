# Copyright 2017-2025 Andrew Myers, David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .ParameterGroup import ParameterGroup
from .WarpX import warpx

particles = warpx.get_parametergroup(
    "particles", species_names=[], rigid_injected_species=[]
)


def new_species(name):
    result = warpx.get_parametergroup(name, Species)
    if name not in particles.species_names:
        particles.species_names.append(name)
    return result


def valid_species(name):
    return name in particles.species_names


class Species(ParameterGroup):
    def valid_subgroup_name(self, name):
        return (name in self.injection_sources) or name in ["attribute"]
