# Copyright 2017-2020 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

class Diagnostic(Bucket):
    """
    This is the same as a Bucket, but checks that any attributes are always given the same value.
    """
    def __init__(self, name):
        Bucket.__init__(self, name, _species_dict={})

    def add_new_attr_with_check(self, name, value):
        if name.startswith('_'):
            self._localsetattr(name, value)
        else:
            if name in self.argvattrs:
                assert value == self.argvattrs[name], Exception(f'Diagnostic attributes not consistent for {self.instancename}')
            self.argvattrs[name] = value

    def __setattr__(self, name, value):
        self.add_new_attr_with_check(name, value)


class Diagnostics(Bucket):
    """
    This keeps a dictionary of diagnostics
    """
    def __init__(self):
        Bucket.__init__(self, 'diagnostics', _diagnostics_dict={})

    def new_diagnostic(self, name):
        if name is None:
            name = 'diag{}'.format(len(self._diagnostics_dict))
        try:
            result = self._diagnostics_dict[name]
        except KeyError:
            result = Diagnostic(name)
            self._diagnostics_dict[name] = result
            self.diags_names = self._diagnostics_dict.keys()
        return result


diagnostics = Diagnostics()

