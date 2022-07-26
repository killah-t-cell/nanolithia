# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import ase.db
from ase.calculators.dftd3 import DFTD3
from ase.dft import DOS
from ase.phasediagram import PhaseDiagram, Pourbaix
from gpaw import PW, GPAW

import numpy as np

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
from properties.voltages import get_eq_voltage

# TODO
# get correct formation energies
# validate 0K voltage profile

def get_formation_energy(db, xc, compound_epot, x, y):
    Li = get_Li(db, xc).toatoms()
    epot_Li_cell = Li.get_potential_energy()
    epot_Li = epot_Li_cell / 2

    O = get_O2(db, xc, ecut=750).toatoms()
    epot_O_cell = O.get_potential_energy()
    epot_O = epot_O_cell / 4

    return compound_epot - (x * epot_Li + y * epot_O)


if __name__ == '__main__':
    # set high quality plotting
    mpl.rcParams['figure.dpi'] = 300  # increase plot dpi

    # get compounds
    db = ase.db.connect('compounds.db')
    xc = 'PBE'

    Li2O = get_Li2O(db, xc).toatoms()
    Li2O2 = get_Li2O2(db, xc).toatoms()
    LiO2 = get_LiO2(db, xc).toatoms()
    LiO3 = get_LiO3(db, xc).toatoms()
    LiO8 = get_LiO8(db, xc).toatoms()
    Li = get_Li(db, xc).toatoms()
    O2 = get_O2(db, xc).toatoms()

    # get potential energies
    epot_Li2O_cell = Li2O.get_potential_energy()
    epot_Li2O2_cell = Li2O2.get_potential_energy()
    epot_LiO2_cell = LiO2.get_potential_energy()
    epot_LiO3_cell = LiO3.get_potential_energy()
    epot_LiO8_cell = LiO8.get_potential_energy()
    epot_O2_cell = O2.get_potential_energy()


    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_Li2O = epot_Li2O_cell / 4  # len is 12, we want to get Li2O out of Li8O4, so we divide by 4. 3/12=0.25
    epot_LiO2 = epot_LiO2_cell / 2
    epot_LiO3 = epot_LiO3_cell / 2
    epot_LiO8 = epot_LiO8_cell / 2

    # get properties
    get_eq_voltage(2 * epot_Li2O, epot_Li2O2, 2)
    get_eq_voltage(2 * epot_Li2O, epot_LiO2, 3)

    # converge
    print('starting convergence study')
    compounds_to_converge = (get_Li2O2, get_LiO2, get_O2)
    converge(db, 'LDA', *compounds_to_converge)

    # get convex hull
    efLiO2 = get_formation_energy(db, xc, epot_LiO2, 1, 2)
    efLi2O2 = get_formation_energy(db, xc, epot_Li2O2, 2, 2)
    efLi2O = get_formation_energy(db, xc, epot_Li2O, 2, 1)
    efLiO3 = get_formation_energy(db, xc, epot_LiO3, 1, 3)
    efLiO8 = get_formation_energy(db, xc, epot_LiO8, 1, 8)

    # print formation energies
    print('epot_Li2O=', epot_Li2O / 3)
    print('efLi2O=', efLi2O / 3)
    print('efLi2O2=', efLi2O2 / 4)
    print('efLiO2=', efLiO2 / 3)
    print('efLiO8=', efLiO8 / 9)

    # save new formation energies to db
    db.update(id=get_Li2O(db, xc).id, formation_energy=efLi2O)
    db.update(id=get_Li2O2(db, xc).id, formation_energy=efLi2O2)
    db.update(id=get_LiO2(db, xc).id, formation_energy=efLiO2)
    db.update(id=get_LiO3(db, xc).id, formation_energy=efLiO3)
    db.update(id=get_LiO8(db, xc).id, formation_energy=efLiO8)

    refs = [('Li2O', get_Li2O(db, xc).formation_energy),
            ('Li2O2', get_Li2O2(db, xc).formation_energy),
            ('LiO2', get_LiO2(db, xc).formation_energy),
            ('LiO3', get_LiO3(db, xc).formation_energy),
            ('LiO8', get_LiO8(db, xc).formation_energy),
            ('Li', 0.0),
            ('O', 0.0),
            ]

    # plot convex hull
    pd = PhaseDiagram(refs)
    pd.plot()
    plt.savefig(f'plots/convex-hull-0K.png')
    pd.plot(show=True)

    # get voltage profile
    v_points = [get_eq_voltage(2 * epot_Li2O, epot_Li2O2, 2),
                get_eq_voltage(epot_Li2O2, epot_LiO2, 1),
                get_eq_voltage(4 * epot_LiO2, epot_LiO8, 3)]
    concentrations = (1 / 3, 1 / 2, 2 / 3)

    plt.plot(concentrations,v_points)
    plt.ylabel('Voltage (V)')
    plt.xlabel('O')
    plt.savefig(f'plots/voltage-profile-0K.png')
    plt.show()
