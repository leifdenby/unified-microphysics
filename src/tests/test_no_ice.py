# coding: utf-8
import unified_microphysics
import pyclouds.cloud_microphysics
from pyclouds.common import Var
import numpy as np
import random

pylib = unified_microphysics.microphysics_pylib
um_register = unified_microphysics.microphysics_register
um_common = unified_microphysics.microphysics_common

def test_init():
    pylib.init('no_ice')

    assert um_register.n_gases == 2
    assert um_register.n_solids == 2

def test_state_mapping():
    pylib.init('no_ice')
    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()
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

    q_g, q_tr, T = state_mapping.pycloud_um(F=F)
    F2 = state_mapping.um_pycloud(q_g=q_g, q_tr=q_tr, T=T)

    assert np.all(F == F2)

def test_cond_evap():
    from test_common import um_constants

    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    mu_model = unified_microphysics.mphys_no_ice

    T = 300.
    p = 101325.
    rho = 1.125
    rho_g = rho*0.9  # TODO: generally gas has a lower density that whole mixture I think
    qv = 0.017
    ql = 0.01

    dqldt_1 = pyclouds_model.dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)
    dqldt_2 = mu_model.dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, t=T, p=p)

    assert abs(dqldt_1 - dqldt_2) < 1.0e-16


def _test_full():
    from test_common import um_constants

    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()
    F = np.zeros((Var.NUM))
    F[Var.q_v] = 0.017
    F[Var.T] = 285.0
    p = 101325.0

    q_g, q_tr, T = state_mapping.pycloud_um(F=F)

    dqdt_g, dqdt_tr, dTdt, _ = pylib.dqdt(q_g=q_g, q_tr=q_tr, temp=T, pressure=p)
    dFdz1 = state_mapping.um_pycloud(q_g=dqdt_g, q_tr=dqdt_tr, T=dTdt)

    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    dFdz2 = pyclouds_model.dFdt(F=F, t=None, p=p)

    Var.print_formatted(dFdz1)
    Var.print_formatted(dFdz2)

    assert np.all(dFdz1 == dFdz2)


if __name__ == "__main__":
    test_cond_evap()
