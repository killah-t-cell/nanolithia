from ase.build import bulk
from gpaw import GPAW, PW
from matplotlib import pyplot as plt
from compounds.Li2O.Li2O import get_Li2O


# This function takes an atom call function and tries to converge it wrt to ecut or kpts
# it outputs the converged k and ecut based on a given tol
# and prints/saves the plt as a side effect
def get_ecut_convergence(atoms, xc, tol=1e-2, k=8, ecut_max=800):
    ecut_convergence = {}
    for ecut in range(100, ecut_max, 50):
        atom = atoms(xc, Ecut=ecut, kpts=(k, k, k))
        e = atom.get_potential_energy()
        ecut_convergence[ecut] = e

    lists = sorted(ecut_convergence.items())  # sorted by key, return a list of tuples
    ecut_list, e_list = zip(*lists)  # unpack a list of pairs into two tuples
    plt.plot(ecut_list, e_list)
    plt.show()

    return ecut_list[next((i for i in range(1, len(lists)) if lists[i][1] - lists[i + 1][1] < tol), lists[-1])]


def get_kpts_convergence(atoms, xc, tol=1e-2, ecut=500):
    k_convergence = {}
    for k in [4, 8, 12, 16, 20, 24]:
        atom = atoms(xc, Ecut=ecut, kpts=(k, k, k))
        e = atom.get_potential_energy()
        k_convergence[k] = e

    lists = sorted(k_convergence.items())  # sorted by key, return a list of tuples
    k_list, e_list = zip(*lists)  # unpack a list of pairs into two tuples
    plt.plot(k_list, e_list)
    plt.show()

    return k_list[next(i for i in range(1, len(lists)) if lists[i][1] - lists[i + 1][1] < tol)]


Li2O = get_Li2O('PBE', Ecut=800, kpts=(12, 12, 12))
Li2O.get_potential_energy() / 12
