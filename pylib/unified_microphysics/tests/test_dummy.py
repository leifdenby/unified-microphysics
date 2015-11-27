# coding: utf-8
import unified_microphysics
import pyclouds.cloud_microphysics
from pyclouds.common import Var
import numpy as np
import random

um_fortran = unified_microphysics.fortran
um_register = um_fortran.microphysics_register
um_common = um_fortran.microphysics_common
um_constants = unified_microphysics.constants
pylib = um_fortran.microphysics_pylib

def test_init():
    pylib.init('dummy', 'isobaric')

    assert um_register.n_compressible_species == 1
    assert um_register.n_incompressible_species == 4

    pylib.init('dummy', 'isobaric')

    assert um_register.n_compressible_species == 1
    assert um_register.n_incompressible_species == 4

def test_state_mapping():
    pylib.init('dummy', 'isobaric')
    state_mapping = unified_microphysics.utils.PyCloudsUnifiedMicrophysicsStateMapping()
    F = np.zeros((Var.NUM))

    if um_register.idx_water_vapour != 0:
        F[Var.q_v] = random.random()
    if um_register.idx_cwater != 0:
        F[Var.q_l] = random.random()
    if um_register.idx_rain != 0:
        F[Var.q_r] = random.random()
    if um_register.idx_cice != 0:
        F[Var.q_i] = random.random()
    F[Var.T] = random.random()

    p = random.random()

    y = state_mapping.pycloud_um(F=F, p=p)
    F2, p2 = state_mapping.um_pycloud(y=y)

    assert np.all(F == F2)
    assert p == p2

def test_state_mapping2():
    pylib.init('dummy', 'isobaric')
    state_mapping = unified_microphysics.utils.PyCloudsUnifiedMicrophysicsStateMapping()

    y = np.zeros((um_register.n_variables))

    # fortran indexing starts at 1
    if um_register.idx_water_vapour != 0:
        y[um_register.idx_water_vapour-1] = random.random()
    if um_register.idx_cwater != 0:
        y[um_register.idx_cwater-1] = random.random()
    if um_register.idx_rain != 0:
        y[um_register.idx_rain-1] = random.random()
    if um_register.idx_cice != 0:
        y[um_register.idx_cice-1] = random.random()

    y[um_register.idx_temp-1] = random.random()
    y[um_register.idx_pressure-1] = random.random()

    F, p2 = state_mapping.um_pycloud(y=y)
    y2 = state_mapping.pycloud_um(F=F, p=p2)

    assert np.all(y == y2)
