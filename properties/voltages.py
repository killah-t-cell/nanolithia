def get_eq_voltage(epot_A, epot_B, n: int):
    # calculate the equilibrium potential for the case of a Li2O/Li metal battery from the intercalation energy of Li in
    # Li2O. For simplicity, use that assumption that all vibrational energies and entropic terms cancel each other.

    # Li Todo this is a bit hacky, it's better to get the data from the DB
    epot_Li = -1.90

    # calculate equilibrium voltage for both reaction ∆U/e
    V_eq = - ((epot_A - epot_B - n * epot_Li) / n)
    print('V_eq =', V_eq)

    return V_eq
