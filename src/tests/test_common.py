# coding: utf-8

import pyclouds.parameterisations
import pyclouds.common
import unified_microphysics

# dictionary-based mapping of constants in fortran code to easy testing with
# equations written in python
um_constants_fortran = unified_microphysics.microphysics_constants
um_constants = {
    'cp_d': um_constants_fortran.cp_d,
    'cv_d': um_constants_fortran.cv_d,
    'cp_v': um_constants_fortran.cp_v,
    'cv_v': um_constants_fortran.cv_v,
    'pv_sat': {
        'p0vs': um_constants_fortran.p0vs,
        'a0_lq': um_constants_fortran.a0_lq,
        'a1_lq': um_constants_fortran.a1_lq,
        'a0_ice': um_constants_fortran.a0_ice,
        'a1_ice': um_constants_fortran.a1_ice,
    },
    'Ka': {
        'a_K': um_constants_fortran.a_k,
        'b_K': um_constants_fortran.b_k,
    },
    'Dv': {
        'a_D': um_constants_fortran.a_d,
        'b_D': um_constants_fortran.b_d,
    },
    'L_v': um_constants_fortran.l_cond,
    'rho_l': um_constants_fortran.rho_w,
    'rho_i': um_constants_fortran.rho_i,
}

assert_close = lambda v1, v2: abs(v2 - v1)/v1 < 0.05

def test_pvsat():
    f1 = unified_microphysics.microphysics_common.saturation_vapour_pressure

    pyclouds_parameterisations = pyclouds.parameterisations.ParametersationsWithSpecificConstants(constants=um_constants)
    f2 = pyclouds_parameterisations.pv_sat

    assert f1(300.) == f2(300.)
    assert f1(273.) == f2(273.)
    assert f1(250.) == f2(250.)

    # use constants defined in unified microphysics code
    f3 = pyclouds.parameterisations.pv_sat
    assert_close(f1(300.), f3(300.))
    assert_close(f1(273.), f3(273.))
    assert_close(f1(250.), f3(250.))

def test_thermal_conductivity_parameterisation():
    f1 = unified_microphysics.microphysics_common.thermal_conductivity

    pyclouds_parameterisations = pyclouds.parameterisations.ParametersationsWithSpecificConstants(constants=um_constants)
    f2 = pyclouds_parameterisations.Ka

    assert f1(300.) == f2(300.)
    assert f1(273.) == f2(273.)
    assert f1(250.) == f2(250.)

    # use constants defined in unified microphysics code
    f3 = pyclouds.parameterisations.Ka
    assert_close(f1(300.), f3(300.))
    assert_close(f1(273.), f3(273.))
    assert_close(f1(250.), f3(250.))
