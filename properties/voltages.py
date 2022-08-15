def get_eq_voltage(gibbs_A, gibbs_B, n):
    # calculate the equilibrium potential for the case of a Li2O/Li metal battery.
    # calculate equilibrium voltage for both reaction ∆G/n where n is the number of electrons transferred
    V_eq = - ((gibbs_A - gibbs_B) / n)
    print('V_eq =', V_eq)

    return V_eq
