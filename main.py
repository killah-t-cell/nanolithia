# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import ase.db

from compounds.Li.Li import get_Li
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from compounds.LiO2.LiO2 import get_LiO2


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
            for k in k_range:
                arg(db, xc, nkpts=k)
                # store converged tol
                # update nkpts for arg in db
            for ecut in ecut_range:
                arg(db, xc, ecut=ecut)
                # update ecut for arg in db


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

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
