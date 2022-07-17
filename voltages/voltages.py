from compounds.Li.Li import *
from compounds.Li2O2.Li2O2 import *
from compounds.Li2O.Li2O import *
from compounds.LiO2.LiO2 import *
from convergence import get_ecut_convergence, get_kpts_convergence


def get_eq_voltage(xc):
    # calculate the equilibrium potential for the case of a Li2O/Li metal battery from the intercalation energy of Li in
    # Li2O. For simplicity, use that assumption that all vibrational energies and entropic terms cancel each other.

    # get converged calculations
    Li_ecut = get_ecut_convergence(get_Li, xc)
    Li_kpts = get_kpts_convergence(get_Li, xc)
    Li2O2_ecut = get_ecut_convergence(get_Li2O2, xc)
    Li2O2_kpts = get_kpts_convergence(get_Li2O2, xc)
    Li2O_ecut = get_ecut_convergence(get_Li2O, xc)
    Li2O_kpts = get_kpts_convergence(get_Li2O, xc)
    LiO2_ecut = get_ecut_convergence(get_LiO2, xc)
    LiO2_kpts = get_kpts_convergence(get_LiO2, xc)

    Li = get_Li(xc, Ecut=Li_ecut, kpts=(Li_kpts, Li_kpts, Li_kpts))
    Li2O2 = get_Li2O2(xc, Ecut=Li2O2_ecut, kpts=(Li2O2_kpts, Li2O2_kpts, Li2O2_kpts))
    Li2O = get_Li2O(xc, Ecut=Li2O_ecut, kpts=(Li2O_kpts, Li2O_kpts, Li2O_kpts))
    LiO2 = get_LiO2(xc, Ecut=LiO2_ecut, kpts=(LiO2_kpts, LiO2_kpts, LiO2_kpts))

    epot_Li_cell = Li.get_potential_energy()
    epot_Li2O2_cell = Li2O2.get_potential_energy()
    epot_Li2O_cell = Li2O.get_potential_energy()
    epot_LiO2_cell = LiO2.get_potential_energy()


    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li = epot_Li_cell / 1  # len is 2, we want to get Li out of Li2, so we divide by 1. 1/2=0.5
    epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_Li2O = epot_Li2O_cell / 4  # len is 12, we want to get Li2O out of Li8O4, so we divide by 4. 3/12=0.25
    epot_LiO2 = epot_LiO2_cell / 4  # len is 12, we want to get LiO2 out of Li4O2, so we divide by 4. 3/12=0.25

    # calculate equilibrium voltage for both reaction ∆U/e
    V_eq_w_Li2O2 = (2 * epot_Li2O - epot_Li2O2 - 2 * epot_Li) / 2
    V_eq_w_LiO2 = (2 * epot_Li2O - epot_LiO2 - 3 * epot_Li) / 3

    V_eq_avg = (V_eq_w_LiO2 + V_eq_w_Li2O2) / 2
    print('V_eq_w_Li2O2 =', V_eq_w_Li2O2)
    print('V_eq_w_LiO2 =', V_eq_w_LiO2)
    print('avg(V_eq) =', V_eq_avg)

    return V_eq_avg


get_eq_voltage('PBE')

