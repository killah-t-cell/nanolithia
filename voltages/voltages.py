
def get_eq_voltage(A, B, n: int):
    # calculate the equilibrium potential for the case of a Li2O/Li metal battery from the intercalation energy of Li in
    # Li2O. For simplicity, use that assumption that all vibrational energies and entropic terms cancel each other.

    # Load Li
    epot_Li_cell = 1

    # Load potential energies (from db)
    epot_A_cell = A.get_potential_energy()
    epot_B_cell = B.get_potential_energy()

    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li = epot_Li_cell / 1  # len is 2, we want to get Li out of Li2, so we divide by 1. 1/2=0.5
    epot_A = epot_A_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_B = epot_B_cell / 4  # len is 12, we want to get Li2O out of Li8O4, so we divide by 4. 3/12=0.25

    # calculate equilibrium voltage for both reaction ∆U/e
    V_eq = (2 * epot_A - epot_B - n * epot_Li) / n
    print('V_eq =', V_eq)

    return V_eq
