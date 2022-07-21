# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import ase.db
from ase.calculators.dftd3 import DFTD3
from gpaw import PW, GPAW

from compounds.Li.Li import get_Li
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from compounds.LiO2.LiO2 import get_LiO2

from matplotlib import pyplot as plt
import matplotlib as mpl

from compounds.LiO3.LiO3 import get_LiO3
from compounds.LiO8.LiO8 import get_LiO8
from compounds.O2.O2 import get_O2
from convergence import converge
from voltages.voltages import get_eq_voltage

# get all energies
# get formation energies
# construct hull
# calculate 0K voltage profile

def get_formation_energy(db, xc, compound_epot, x, y):
    Li = get_Li(db, xc)
    epot_Li_cell = Li.get_potential_energy()
    epot_Li = epot_Li_cell / 2

    O = get_O2(db, xc)
    epot_O_cell = O.get_potential_energy()
    epot_O = epot_O_cell / 4

    return compound_epot - (x * epot_Li - y * epot_O)

if __name__ == '__main__':
    # get compounds
    db = ase.db.connect('compounds.db')
    xc = 'PBE'

    Li = get_Li(db, xc)
    Li2O = get_Li2O(db, xc)
    Li2O2 = get_Li2O2(db, xc)
    LiO2 = get_LiO2(db, xc)
    LiO3 = get_LiO3(db, xc)
    LiO8 = get_LiO8(db, xc)
    O2 = get_O2(db, xc)

    # get potential energies
    epot_Li2O_cell = Li2O.get_potential_energy()
    epot_Li2O2_cell = Li2O2.get_potential_energy()
    epot_LiO2_cell = LiO2.get_potential_energy()
    epot_LiO3_cell = LiO3.get_potential_energy()
    epot_LiO8_cell = LiO8.get_potential_energy()
    epot_Li_cell = Li.get_potential_energy()

    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_Li2O = epot_Li2O_cell / 4  # len is 12, we want to get Li2O out of Li8O4, so we divide by 4. 3/12=0.25
    epot_LiO2 = epot_LiO2_cell / 2
    epot_LiO3 = epot_LiO3_cell / 2
    epot_LiO8 = epot_LiO8_cell / 2
    epot_Li = epot_Li_cell / 2
    print(epot_Li2O)
    print(epot_Li2O2)
    print(epot_LiO2)

    # get voltages
    get_eq_voltage(2 * epot_Li2O, epot_Li2O2, 2)
    get_eq_voltage(2 * epot_Li2O, epot_LiO2, 3)
    # compounds_to_converge = (get_Li2O2, get_LiO2, get_Li, get_Li2O)
    # converge(db, xc, *compounds_to_converge)

    # get convex hull TODO this is super hacky and brittle, should make it sturdier, but let's get a convex hull first
    efLiO2 = get_formation_energy(db, xc, epot_LiO2, 1, 2)
    efLi2O2 = get_formation_energy(db, xc, epot_Li2O2, 2, 2)
    efLi2O = get_formation_energy(db, xc, epot_Li2O, 2, 1)
    efLiO3 = get_formation_energy(db, xc, epot_LiO3, 1, 3)
    efLiO8 = get_formation_energy(db, xc, epot_LiO8, 1, 8)

    db.update(id=4, formation_energy=efLi2O)
    db.update(id=6, formation_energy=efLi2O2)
    db.update(id=8, formation_energy=efLiO2)
    db.update(id=12, formation_energy=efLiO3)
    db.update(id=14, formation_energy=efLiO8)




