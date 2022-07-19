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
from voltages.voltages import get_eq_voltage

if __name__ == '__main__':
    # get compounds
    db = ase.db.connect('compounds.db')
    xc = 'PBE'
    Li = get_Li(db, xc)
    Li2O = get_Li2O(db, xc)
    Li2O2 = get_Li2O2(db, xc)
    LiO2 = get_LiO2(db, xc)

    # get potential energies
    epot_Li2O_cell = Li2O.get_potential_energy()
    epot_Li2O2_cell = Li2O2.get_potential_energy()
    epot_LiO2_cell = LiO2.get_potential_energy()

    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_Li2O = epot_Li2O_cell / 4  # len is 12, we want to get Li2O out of Li8O4, so we divide by 4. 3/12=0.25
    epot_LiO2 = epot_LiO2_cell / 2
    print(epot_Li2O)
    print(epot_Li2O2)
    print(epot_LiO2)


    get_eq_voltage(2*epot_Li2O, epot_Li2O2, 2)
    get_eq_voltage(2*epot_Li2O, epot_LiO2, 3)
    # compounds_to_converge = (get_Li2O2, get_LiO2, get_Li, get_Li2O)
    # converge(db, xc, *compounds_to_converge)
