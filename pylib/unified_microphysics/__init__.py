import utils
import unified_microphysics_fortran as fortran

# dictionary-based mapping of constants in fortran code to ease testing with
# equations written in python
um_constants_fortran = fortran.microphysics_constants
constants = {
    'cp_d': um_constants_fortran.cp_d,
    'cv_d': um_constants_fortran.cv_d,
    'cp_v': um_constants_fortran.cp_v,
    'cv_v': um_constants_fortran.cv_v,
    'cp_l': um_constants_fortran.cp_l,
    'cv_l': um_constants_fortran.cv_l,
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

__all__ = ['utils', 'fortran', 'constants', ]
