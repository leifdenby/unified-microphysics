# coding: utf-8
import unified_microphysics
import pyclouds.cloud_microphysics
from pyclouds.common import Var
import numpy as np
import random
from test_common import um_constants

um_fortran = unified_microphysics.fortran
um_register = um_fortran.microphysics_register
um_common = um_fortran.microphysics_common
pylib = um_fortran.microphysics_pylib

def test_init():
    pylib.init('dummy')

    assert um_register.n_gases == 2
    assert um_register.n_solids == 4

    pylib.init('dummy')

    assert um_register.n_gases == 2
    assert um_register.n_solids == 4
