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
mphys_no_ice = um_fortran.mphys_no_ice
# um_integration = um_fortran.microphysics_integration

def test_init():
    pylib.init('no_ice')

    assert um_register.n_compressible_species == 1
    assert um_register.n_incompressible_species == 2

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

    p = random.random()

    y = state_mapping.pycloud_um(F=F, p=p)
    F2, p2 = state_mapping.um_pycloud(y=y)

    assert np.all(F == F2)
    assert p == p2

def test_state_mapping2():
    pylib.init('no_ice')
    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()

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

    T = 288.
    p = 88676.
    qv = 1.2e-2
    ql = 2.0e-7
    qr = 0.0
    qd = 1.0 - ql - qv
    rho = mu_model.rho_f(qd=qd, ql=ql, qv=qv, qr=qr, p=p, temp=T)
    rho_g = mu_model.rho_f(qd=qd, ql=0.0, qv=qv, qr=0.0, p=p, temp=T)

    dqldt_1 = pyclouds_model.dql_dt__cond_evap(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=T, p=p)
    dqldt_2 = mu_model.dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, t=T, p=p)

    assert abs(dqldt_1 - dqldt_2) < 1.0e-16


def test_full1():
    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)


    F = np.zeros((Var.NUM))
    F[Var.q_v] = 0.017
    F[Var.T] = 285.0
    p = 101325.0

    y = state_mapping.pycloud_um(F=F, p=p)

    dydt = mphys_no_ice.dydt(y=y, t=0.0)
    dFdz1, _ = state_mapping.um_pycloud(y=dydt)

    dFdz2 = pyclouds_model.dFdt(F=F, t=None, p=p)

    assert np.all(dFdz1 == dFdz2)

def test_full2():
    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    # sub-saturation state
    F = np.zeros((Var.NUM))
    T = 288.
    p = 88676.
    qv = 1.1e-2
    ql = 2.0e-4
    F[Var.q_v] = qv
    F[Var.q_l] = ql
    F[Var.T] = T

    y = state_mapping.pycloud_um(F=F, p=p)

    dydt = mphys_no_ice.dydt(y=y, t=0.0)
    dFdz1, _ = state_mapping.um_pycloud(y=dydt)

    dFdz2 = pyclouds_model.dFdt(F=F, t=None, p=p)

    Sw = qv/pyclouds_model.qv_sat(T=T, p=p)

    assert Sw < 1.0

    Var.print_formatted(dFdz1, '%.10g')
    Var.print_formatted(dFdz2, '%.10g')
    Var.print_formatted(dFdz2-dFdz1, '%.10g')

    assert np.all(dFdz1 == dFdz2)
    assert dFdz1[Var.q_l] < 0.0


def test_equation_of_state():
    pylib.init('no_ice')
    assert um_register.n_compressible_species == 1
    assert um_register.n_incompressible_species == 2

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


def test_heat_capacity():
    state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()

    # init pyclouds model
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    # sub-saturation state
    F = np.zeros((Var.NUM))
    T = 288.
    p = 88676.
    qv = 1.2e-2
    ql = 2.0e-7
    qr = 0.0
    qd = 1.0 - ql - qv
    F[Var.q_v] = qv
    F[Var.T] = T
    F[Var.q_l] = ql

    y = state_mapping.pycloud_um(F=F, p=p)

    cp_m__1 = pyclouds_model.cp_m(F=F)
    cp_m__2 = mphys_no_ice.cp_m(y)

    assert abs(cp_m__1 - cp_m__2) < 1.0e-16


    # sub-saturation state
    F = np.zeros((Var.NUM))
    T = 288.
    p = 88676.
    qv = 1.2e-2
    ql = 5.8e-8
    qr = 0.0
    qd = 1.0 - ql - qv
    F[Var.q_v] = qv
    F[Var.T] = T
    F[Var.q_l] = ql

    y = state_mapping.pycloud_um(F=F, p=p)

    cp_m__1 = pyclouds_model.cp_m(F=F)
    cp_m__2 = mphys_no_ice.cp_m(y)

    assert abs(cp_m__1 - cp_m__2) < 1.0e-16

# def test_integration():
    # pylib.init('no_ice')
    # assert um_register.n_compressible_species == 1
    # assert um_register.n_incompressible_species == 2

    # mu_model = um_fortran.mphys_no_ice
    # mu_common = um_fortran.microphysics_common

    # state_mapping = pyclouds.cloud_microphysics.PyCloudsUnifiedMicrophysicsStateMapping()

    # F = np.zeros((Var.NUM))
    # F[Var.q_v] = 0.01
    # F[Var.T] = 300.
    # p = 101325.0

    # y = state_mapping.pycloud_um(F=F, p=p)

    # dy, _ = um_integration.calc_dy_with_message(y, dt=1.0)


if __name__ == "__main__":
    test_state_mapping2()
