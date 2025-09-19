# Copyright 2016-2025 Andrew Myers, David Grote, Maxence Thevenet
# Remi Lehe, Lorenzo Giacomel
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import re
import sys

from ._libwarpx import libwarpx
from .ParameterGroup import ParameterGroup


class WarpX(ParameterGroup):
    """
    A Python wrapper for the WarpX C++ class
    """

    def create_argv_list(self, **kw):
        argv = []

        for k, v in kw.items():
            if v is not None:
                argv.append(f"{k} = {v}")

        argv += self.attrlist()

        for parametergroup in self._parametergroup_dict.values():
            argv += parametergroup.attrlist()

        return argv

    def get_parametergroup(self, parametergroup_name, group_class=ParameterGroup, **kw):
        try:
            return self._parametergroup_dict[parametergroup_name]
        except KeyError:
            parametergroup = group_class(parametergroup_name, **kw)
            self._parametergroup_dict[parametergroup_name] = parametergroup
            return parametergroup

    def init(self, mpi_comm=None, **kw):
        # note: argv[0] needs to be an absolute path so it works with AMReX backtraces
        # https://github.com/AMReX-Codes/amrex/issues/3435
        argv = [sys.executable] + self.create_argv_list(**kw)
        libwarpx.initialize(argv, mpi_comm=mpi_comm)

    def evolve(self, nsteps=-1):
        libwarpx.warpx.evolve(nsteps)

    def finalize(self, finalize_mpi=1):
        libwarpx.finalize(finalize_mpi)

    def getProbLo(self, direction):
        return libwarpx.libwarpx_so.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        return libwarpx.libwarpx_so.warpx_getProbHi(direction)

    def write_inputs(self, filename="inputs", **kw):
        argv = self.create_argv_list(**kw)

        # Sort the argv list to make it more human readable
        argv.sort()

        with open(filename, "w") as ff:
            prefix_old = ""
            for arg in argv:
                # This prints the name of the input group (prefix) as a header
                # before each group to make the input file more human readable
                prefix_new = re.split(r" |\.", arg)[0]
                if prefix_new != prefix_old:
                    if prefix_old != "":
                        ff.write("\n")
                    ff.write(f"# {prefix_new}\n")
                    prefix_old = prefix_new

                ff.write(f"{arg}\n")


warpx = WarpX("warpx", _parametergroup_dict={})
