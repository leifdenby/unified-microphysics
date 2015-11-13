# coding: utf-8

import pyclouds.parameterisations
import pyclouds.common

import unified_microphysics

um_fortran = unified_microphysics.fortran
um_constants = unified_microphysics.constants

assert_close = lambda v1, v2: abs(v2 - v1)/v1 < 0.05

def test_pvsat():
    f1 = um_fortran.microphysics_common.saturation_vapour_pressure

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
    f1 = um_fortran.microphysics_common.thermal_conductivity 
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

def test_water_vapour_diffusion_parameterisation():
    p0 = 101325.0
    f1 = um_fortran.microphysics_common.water_vapour_diffusivity
    pyclouds_parameterisations = pyclouds.parameterisations.ParametersationsWithSpecificConstants(constants=um_constants)
    f2 = pyclouds_parameterisations.Dv

    assert f1(300., p0) == f2(300., p0)
    assert f1(273., p0) == f2(273., p0)
    assert f1(250., p0) == f2(250., p0)

    # use constants defined in unified microphysics code
    f3 = pyclouds.parameterisations.Dv
    assert_close(f1(300., p0), f3(300., p0))
    assert_close(f1(273., p0), f3(273., p0))
    assert_close(f1(250., p0), f3(250., p0))
