# Copyright 2017-2025 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .ParameterGroup import ParameterGroup
from .Particles import valid_species
from .WarpX import warpx

diagnostics = warpx.get_parametergroup("diagnostics", diags_names=[])
reduced_diagnostics = warpx.get_parametergroup("warpx", diags_names=[])


def new_diagnostic(name):
    diag = warpx.get_parametergroup(name, Diagnostic)
    if name not in diagnostics.diags_names:
        diagnostics.diags_names.append(name)
    return diag


class Diagnostic(ParameterGroup):
    """
    This is the same as a ParameterGroup, but checks that any attributes are always given the same value.
    """

    def add_new_attr_with_check(self, name, value):
        if name.startswith("_"):
            self._localsetattr(name, value)
        else:
            if name in self.argvattrs:
                assert value == self.argvattrs[name], Exception(
                    f"Diagnostic attributes not consistent for "
                    f'"{self.instancename}": '
                    f'"{value}" != "{self.argvattrs[name]}"'
                )
            self.argvattrs[name] = value

    def valid_subgroup_name(self, name):
        return valid_species(name)

    def __setattr__(self, name, value):
        self.add_new_attr_with_check(name, value)

    def set_or_replace_attr(self, name, value):
        """
        Explicitly set or replace an existing attribute
        (since __setattr__ cannot be used for replacing
        as it would raise an Exception)
        """
        assert not name.startswith("_")
        self.argvattrs[name] = value
