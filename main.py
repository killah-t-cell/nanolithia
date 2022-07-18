# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import ase.db

from compounds.Li.Li import get_Li
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from compounds.LiO2.LiO2 import get_LiO2

from matplotlib import pyplot as plt
import matplotlib as mpl

from convergence import converge

if __name__ == '__main__':
    db = ase.db.connect('compounds.db')
    xc = 'PBE'
    Li = get_Li(db, xc)
    Li2O = get_Li2O(db, xc)
    Li2O2 = get_Li2O2(db, xc)
    LiO2 = get_LiO2(db, xc)

    # compounds_to_converge = (get_Li, get_Li2O, get_Li2O2, get_LiO2)
    # converge(db, xc, *compounds_to_converge)

    # load voltages
    # if not there
    # run voltage experiments
    # del db[db.get(name=f'Li-{xc}-8x8x8-500', converged=True).id]
