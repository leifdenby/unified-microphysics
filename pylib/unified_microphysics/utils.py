import warnings
import numpy as np

from pyclouds.common import Var
import unified_microphysics_fortran as um_fortran

class PyCloudsUnifiedMicrophysicsStateMapping():
    """
    Utility for providing mapping between pycloud state representation and
    unified-microphysics model representation.

    TODO: Make pyclouds use indexing from `unified microphysics`
    
    XXX: this is kinda hacky, I only reassign state variables that are actually
    used in the initialized model in `unified microphysics`
    """
    def __init__(self):
        register = um_fortran.microphysics_register
        # Fortran indexing starts at 1
        self.idx_water_vapour = register.idx_water_vapour-1
        self.idx_cwater = register.idx_cwater-1
        self.idx_rain = register.idx_rain-1
        self.idx_cice = register.idx_cice-1
        self.idx_temp = register.idx_temp-1
        self.idx_pressure = register.idx_pressure-1

        self.n_vars = register.n_variables

        if register.idx_temp == 0:
            raise Exception("Init hasn't been called on the fortran backend")

    def um_pycloud(self, y):
        F = np.zeros((Var.NUM))
        if self.idx_water_vapour != -1:
            F[Var.q_v] = y[self.idx_water_vapour]
        if self.idx_cwater != -1:
            F[Var.q_l] = y[self.idx_cwater]
        if self.idx_rain != -1:
            F[Var.q_r] = y[self.idx_rain]
        if self.idx_cice != -1:
            F[Var.q_i] = y[self.idx_cice]
        F[Var.T] = y[self.idx_temp]
        F[Var.p] = y[self.idx_pressure]
        return F

    def pycloud_um(self, F):
        y = np.zeros((self.n_vars,))

        if self.idx_water_vapour != -1:
            y[self.idx_water_vapour] = F[Var.q_v]
        if self.idx_cwater != -1:
            y[self.idx_cwater] = F[Var.q_l]
        if self.idx_rain != -1:
            y[self.idx_rain] = F[Var.q_r]
        if self.idx_cice != -1:
            y[self.idx_cice] = F[Var.q_i]
        y[self.idx_temp] = F[Var.T]
        y[self.idx_pressure] = F[Var.p]
        return y

def multistep_integration(F0, t, with_debug=False):
    warnings.warn("Creating a integration wrapper that is hardcoded for `no_ice`")
    state_mapping = PyCloudsUnifiedMicrophysicsStateMapping()

    y = state_mapping.pycloud_um(F0)
    F = [state_mapping.um_pycloud(y=y),]
    t_ = [t[0],]

    n_steps = 0
    for tn in range(len(t)-1):
        # modifies `y` in-place
        um_fortran.microphysics_pylib.integrate_microphysics(y=y, t=t[tn], t_end=t[tn+1])
        # TODO: re-implement getting out the total number of steps
        # n_steps += m_total
        F.append(state_mapping.um_pycloud(y=y))
        t_.append(t[tn+1])

    if with_debug:
        return np.array(F), np.array(t_), n_steps
    else:
        return np.array(F), np.array(t_)
