# coding: utf-8

import pyclouds.parameterisations
import pyclouds.common
import unified_microphysics


def test_pvsat():
    f1 = unified_microphysics.microphysics_common.saturation_vapour_pressure

    um_constants = {
        'cp_d': unified_microphysics.microphysics_constants.cp_d,
        'cv_d': unified_microphysics.microphysics_constants.cv_d,
        'cp_v': unified_microphysics.microphysics_constants.cp_v,
        'cv_v': unified_microphysics.microphysics_constants.cv_v,
        'pv_sat': {
            'p0vs': unified_microphysics.microphysics_constants.p0vs,
            'a0_lq': unified_microphysics.microphysics_constants.a0_lq,
            'a1_lq': unified_microphysics.microphysics_constants.a1_lq,
            'a0_ice': unified_microphysics.microphysics_constants.a0_ice,
            'a1_ice': unified_microphysics.microphysics_constants.a1_ice,
        }
    }

    pyclouds_parameterisations = pyclouds.parameterisations.ParametersationsWithSpecificConstants(constants=um_constants)
    f2 = pyclouds_parameterisations.pv_sat

    def test(T):
        assert f1(T) == f2(T)

    test(300.)
    test(273.)
    test(250.)
