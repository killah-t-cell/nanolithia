# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import ase.db

from compounds.Li.Li import get_Li
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from compounds.LiO2.LiO2 import get_LiO2

from matplotlib import pyplot as plt

def converge(db, xc, tol=1e-4, k_range=range(4, 18), ecut_range=range(350, 801, 50), *args):
    """
    :param k_range:
    :param ecut_range:
    :param args: AtomsRow
    :return:
    """
    for arg in args:
        # if not converged
        # converge
        id = db.reserve(name=arg.name, converged=False)
        if id is not None:  # skip calculation if already done
            ecut_convergence = {}
            nkpts_convergence = {}
            for k in k_range:
                atom = arg(db, xc, nkpts=k)
                e = atom.get_potential_energy()
                ecut_convergence[ecut] = e
                # store converged tol
            for ecut in ecut_range:
                arg(db, xc, ecut=ecut)

            # plot results
            ecut_lists = sorted(ecut_convergence.items())  # sorted by key, return a list of tuples
            nkpts_lists = sorted(nkpts_convergence.items())  # sorted by key, return a list of tuples
            ecut_list, ecut_e_list = zip(*ecut_lists)  # unpack a list of pairs into two tuples
            nkpts_list, nkpts_e_list = zip(*nkpts_lists)  # unpack a list of pairs into two tuples

            plt.plot(ecut_list, ecut_e_list)
            plt.plot(nkpts_list, nkpts_e_list)

            # save results
            # update nkpts for arg in db
            # update ecut for arg in db
            # update tol for arg in db
            # if tol is not achieved, warn and save in db as not converged to tol, but update tol

            # save in plots folder with proper names
            plt.show()
        else:
            print(arg, "is already converged")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    db = ase.db.connect('compounds.db')
    xc = 'PBE'
    Li = get_Li(db, xc)
    Li2O = get_Li2O(db, xc)
    Li2O2 = get_Li2O2(db, xc)
    LiO2 = get_LiO2(db, xc)

    # converge(get_Li, get_Li2O, get_Li2O2, get_LiO2)

    # load voltages
    # if not there
    # run voltage experiments
