# coding: utf-8
import unified_microphysics
import pyclouds.cloud_microphysics
from pyclouds.common import Var
import numpy as np
import random
from test_common import um_constants

um_fortran = unified_microphysics.fortran

pylib = um_fortran.microphysics_pylib
um_register = um_fortran.microphysics_register
um_common = um_fortran.microphysics_common

def test_init():
    pylib.init('no_ice')

    assert um_register.n_gases == 1
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

def test_state_mapping2():
    pylib.init('no_ice')
    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()

    q_g = np.zeros((um_register.n_gases))
    q_tr = np.zeros((um_register.n_solids, um_register.n_moments__max))

    # fortran indexing starts at 1
    if um_register.idx_water_vapour != 0:
        q_g[um_register.idx_water_vapour-1] = random.random()
    if um_register.idx_cwater != 0:
        q_tr[um_register.idx_cwater-1,0] = random.random()
    if um_register.idx_rain != 0:
        q_tr[um_register.idx_rain-1,0] = random.random()
    if um_register.idx_cice != 0:
        q_tr[um_register.idx_cice-1,0] = random.random()

    T = random.random()

    F = state_mapping.um_pycloud(q_g=q_g, q_tr=q_tr, T=T)
    q_g2, q_tr2, T2 = state_mapping.pycloud_um(F=F)

    assert np.all(q_g == q_g2)
    assert np.all(q_tr == q_tr2)
    assert T == T2


def test_cond_evap():
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    mu_model = um_fortran.mphys_no_ice

    T = 300.
    p = 101325.
    qv = 0.017
    ql = 0.01
    qr = 0.0
    qd = 1.0 - ql - qv
    rho = mu_model.rho_f(qd=qd, ql=ql, qv=qv, qr=qr, p=p, temp=T)
    rho_g = mu_model.rho_f(qd=qd, ql=0.0, qv=qv, qr=0.0, p=p, temp=T)

    dqldt_1 = pyclouds_model.dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)
    dqldt_2 = mu_model.dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, t=T, p=p)

    assert abs(dqldt_1 - dqldt_2) < 1.0e-16

    T = 300.
    p = 101325.
    qv = 0.017
    ql = 0.00
    qr = 0.0
    qd = 1.0 - ql - qv
    rho = mu_model.rho_f(qd=qd, ql=ql, qv=qv, qr=qr, p=p, temp=T)
    rho_g = mu_model.rho_f(qd=qd, ql=0.0, qv=qv, qr=0.0, p=p, temp=T)

    dqldt_1 = pyclouds_model.dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)
    dqldt_2 = mu_model.dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, t=T, p=p)

    assert abs(dqldt_1 - dqldt_2) < 1.0e-16

    T = 300.
    p = 101325.
    qv = 0.013
    ql = 0.003
    qr = 0.0
    qd = 1.0 - ql - qv
    rho = mu_model.rho_f(qd=qd, ql=ql, qv=qv, qr=qr, p=p, temp=T)
    rho_g = mu_model.rho_f(qd=qd, ql=0.0, qv=qv, qr=0.0, p=p, temp=T)

    dqldt_1 = pyclouds_model.dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)
    dqldt_2 = mu_model.dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, t=T, p=p)

    assert abs(dqldt_1 - dqldt_2) < 1.0e-16


def _test_full():
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


def test_equation_of_state():
    pylib.init('no_ice')
    assert um_register.n_gases == 1
    assert um_register.n_solids == 2

    mu_model = um_fortran.mphys_no_ice
    mu_common = um_fortran.microphysics_common

    # setup a state
    ql = 0.0010
    qv = 0.01
    qr = 0.0011
    qi = 0.000  # the `no_ice` model doesn't represent ice
    qd = 1.0 - ql - qv - qr - qi
    T = 285.0
    p = 101325.0

    # init pyclouds model
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    rho_1 = mu_model.rho_f(qd=qd, ql=ql, qv=qv, qr=qr, p=p, temp=T)
    rho_2 = pyclouds_model.calc_mixture_density(qd=qd, qv=qv, ql=ql, qi=qi, qr=qr, p=p, T=T)

    assert abs(rho_1 - rho_2) < 1.0e-16

    # setup a state
    ql = 0.0010
    qv = 0.00
    qr = 0.000
    qi = 0.000  # the `no_ice` model doesn't represent ice
    qd = 1.0 - ql - qv - qr - qi
    T = 285.0
    p = 101325.0

    # init pyclouds model
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    rho_1 = mu_model.rho_f(qd=qd, ql=ql, qv=qv, qr=qr, p=p, temp=T)
    rho_2 = pyclouds_model.calc_mixture_density(qd=qd, qv=qv, ql=ql, qi=qi, qr=qr, p=p, T=T)

    assert abs(rho_1 - rho_2) < 1.0e-16


if __name__ == "__main__":
    test_state_mapping2()
