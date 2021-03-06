# coding: utf-8
import unified_microphysics
import pyclouds.cloud_microphysics
from pyclouds.common import Var
import numpy as np
import random

um_fortran = unified_microphysics.fortran

pylib = um_fortran.microphysics_pylib
um_register = um_fortran.microphysics_register
um_common = um_fortran.microphysics_common
mphys_no_ice = um_fortran.mphys_no_ice

um_constants = unified_microphysics.constants

def test_init():
    pylib.init('no_ice', 'isobaric')

    assert um_register.n_compressible_species == 1
    assert um_register.n_incompressible_species == 2


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
    model_constraint = 'isobaric'
    pylib.init('no_ice', model_constraint)

    state_mapping = unified_microphysics.utils.PyCloudsUnifiedMicrophysicsStateMapping()
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants, model_constraint=model_constraint)


    F = np.zeros((Var.NUM))
    F[Var.q_v] = 0.017
    F[Var.T] = 285.0
    F[Var.p] = 101325.0

    y = state_mapping.pycloud_um(F=F)

    # calc heat capacity of mixture
    c_m = pylib.mixture_heat_capacity(y)

    dydt = mphys_no_ice.dydt(y=y, t=0.0, c_m=c_m)
    dFdz1 = state_mapping.um_pycloud(y=dydt)

    dFdz2 = pyclouds_model.dFdt(F=F, t=None)

    assert np.all(dFdz1 == dFdz2)

def test_full2():
    model_constraint = 'isobaric'
    pylib.init('no_ice', model_constraint)

    state_mapping = unified_microphysics.utils.PyCloudsUnifiedMicrophysicsStateMapping()
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants, model_constraint=model_constraint)

    # sub-saturation state
    F = np.zeros((Var.NUM))
    T = 288.
    p = 88676.
    qv = 1.1e-2
    ql = 2.0e-4
    F[Var.q_v] = qv
    F[Var.q_l] = ql
    F[Var.T] = T
    F[Var.p] = p

    y = state_mapping.pycloud_um(F=F)

    # calc heat capacity of mixture
    c_m = pylib.mixture_heat_capacity(y)

    dydt = mphys_no_ice.dydt(y=y, t=0.0, c_m=c_m)
    dFdz1 = state_mapping.um_pycloud(y=dydt)

    dFdz2 = pyclouds_model.dFdt(F=F, t=None)

    Sw = qv/pyclouds_model.qv_sat(T=T, p=p)

    assert Sw < 1.0

    Var.print_formatted(dFdz1, '%.10g')
    Var.print_formatted(dFdz2, '%.10g')
    Var.print_formatted(dFdz2-dFdz1, '%.10g')

    assert np.all(dFdz1 == dFdz2)
    assert dFdz1[Var.q_l] < 0.0


def test_equation_of_state():
    pylib.init('no_ice', 'isobaric')
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


def test_heat_capacity_constant_pressure():
    pylib.init('no_ice', 'isobaric')
    state_mapping = unified_microphysics.utils.PyCloudsUnifiedMicrophysicsStateMapping()

    # init pyclouds model
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    # sub-saturation state
    F = Var.make_state(T=288., q_v=1.2e-2, q_l=2.0e-7)
    p = 88676.

    y = state_mapping.pycloud_um(F=F)

    cp_m__1 = pyclouds_model.cp_m(F=F)
    cp_m__2 = pylib.mixture_heat_capacity(y)

    assert abs(cp_m__1 - cp_m__2) < 1.0e-16


    # sub-saturation state
    F = Var.make_state(T=288., q_v=1.2e-2, q_l=5.8e-8, q_r=0.0)
    p = 88676.
    y = state_mapping.pycloud_um(F=F)

    cp_m__1 = pyclouds_model.cp_m(F=F)
    cp_m__2 = pylib.mixture_heat_capacity(y)

    assert abs(cp_m__1 - cp_m__2) < 1.0e-16


def test_heat_capacity_constant_volume():
    pylib.init('no_ice', 'isometric')
    state_mapping = unified_microphysics.utils.PyCloudsUnifiedMicrophysicsStateMapping()

    # init pyclouds model
    pyclouds_model = pyclouds.cloud_microphysics.FiniteCondensationTimeMicrophysics(constants=um_constants)

    F = Var.make_state(T=288., q_v=1.2e-2, q_l=2.0e-7)
    p = 88676.

    y = state_mapping.pycloud_um(F=F)

    cv_m__1 = pyclouds_model.cv_m(F=F)
    cv_m__2 = pylib.mixture_heat_capacity(y)

    assert abs(cv_m__1 - cv_m__2) < 1.0e-16


def test_long_integration_step():
    """ Make sure that integrator converges to the same solution independently of
    whether the host model requests integration of many small steps or few
    large steps
    """
    pylib.init('no_ice', 'isobaric')

    # get the tolerance values defined in fortran out
    abs_tol = um_constants.integration_abs_tol
    rel_tol = um_constants.integration_rel_tol

    # super-saturated state
    F = np.zeros((Var.NUM))
    F[Var.q_v] =  1.2E-002
    F[Var.q_l] = 3.0E-004
    F[Var.T] = 287.5
    p0 = 87360.0  # [Pa]
    F[Var.p] = p0

    t_max = 300.

    n_1 = 100
    t = np.linspace(0.0, t_max, n_1)
    F1, _ = unified_microphysics.utils.multistep_integration(F0=F, t=t)

    n_2 = n_1/50
    t = np.linspace(0.0, t_max, n_2)
    F2, _ = unified_microphysics.utils.multistep_integration(F0=F, t=t)

    n_3 = n_1*50
    t = np.linspace(0.0, t_max, n_3)
    F3, _ = unified_microphysics.utils.multistep_integration(F0=F, t=t)

    # reference solution
    f1 = F1[-1]

    # make sure we don't have any nans in the solution
    assert not np.any(np.isnan(f1))

    # compare with large steps
    df = F1[-1]-F2[-1]
    rel_err = np.ma.masked_array(df/f1, f1 == 0.0).filled(0.0)
    assert np.linalg.norm(rel_err) < rel_tol

    # compare with small steps
    df = F1[-1]-F2[-1]
    rel_err = np.ma.masked_array(df/f1, f1 == 0.0).filled(0.0)
    assert np.linalg.norm(rel_err) < rel_tol

    # TODO: this is a hack, I'm not sure we can compare errors like this
    abs_err = np.linalg.norm(df)*float(n_2)/float(n_1)*0.01
    print abs_err
    assert abs_err < abs_tol


    # check for changes larger than initial value since this is an indicator of instability
    vars = [Var.q_v, Var.q_l]
    q_max = sum([F[n] for n in vars])

    dF1 = F1[1:] - F1[:-1]
    for n in vars:
        assert np.all(np.abs(dF1[:,n]) < q_max)

    dF2 = F2[1:] - F2[:-1]
    for n in vars:
        assert np.all(np.abs(dF2[:,n]) < q_max)

    dF3 = F3[1:] - F3[:-1]
    for n in vars:
        assert np.all(np.abs(dF3[:,n]) < q_max)

    print "done"

def test_long_integration_very_small():
    """
    Make sure that the integrator can handle very small concentrations
    """
    pylib.init('no_ice', 'isobaric')

    # get the tolerance values defined in fortran out
    abs_tol = um_constants.integration_abs_tol
    rel_tol = um_constants.integration_rel_tol

    # super-saturated state
    F = np.zeros((Var.NUM))
    F[Var.q_v] =  1.3E-002
    F[Var.q_l] = 1.7E-184
    F[Var.T] = 297.019
    p0 = 99835.578483693796 # [Pa]

    # t = [0., 10., 20.]
    # F1, _ = unified_microphysics.utils.multistep_integration(F0=F, t=t)

    # low concentration state
    F = np.zeros((Var.NUM))
    F[Var.q_v] =  1.3E-200
    F[Var.q_l] = 1.7E-184
    F[Var.T] = 297.019
    p0 = 99835.578483693796 # [Pa]

    t = [0., 10., 20.]
    # F1, _ = unified_microphysics.utils.multistep_integration(F0=F, t=t)

    states_dt = []
    # states_dt.append("27861.407803556413        235.09957042102803        1.5134365769422074E-004   2.1718649357875098E-228")
    # states_dt.append("99814.018498220757        297.52241095914763        1.3299859120580425E-002   1.4262881033833411E-024   0.0000000000000000")
    # states_dt.append("95518.458265905574        293.74472467007939        1.1847279113138988E-002   6.5939621073124898E-024   0.0000000000000000")
    states_dt.append((21.261768343566246, "87356.354491443504        287.29050207276362        1.1620359396091596E-002   2.3304008211767397E-004   0.0000000000000000"))
    states_dt.append((0.0, "99835.008685662935        297.02016991797183        1.3598797808482679E-002   0.0000000000000000        0.0000000000000000"))
    states_dt.append((8.517770963152503E-008, "1.361900412522894E-002   291.183518903837        91095.1638538913     7.523163845262640E-037  0.000000000000000E+000"))
    states_dt.append((13.8219455653082, "1.332737186612671E-002   291.285612560467        91094.7037112492     3.673419846319648E-040  0.000000000000000E+000"))

    for dt, state_str in states_dt:
        state = [float(s) for s in state_str.split()]
        F = np.zeros((Var.NUM))
        F[Var.q_v] = state[2]
        F[Var.q_l] = state[3]
        F[Var.T] = state[1]
        F[Var.p] = state[0]

        t = [0., dt, 2.*dt]
        F1, _ = unified_microphysics.utils.multistep_integration(F0=F, t=t)

if __name__ == '__main__':
    test_long_integration_step()
