# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import ase.db
from matplotlib import pyplot as plt

from compounds.compounds import Compound

import matplotlib as mpl

from global_vars import MOL_MASS_LI2O, MOL_VOL_LI2O
from properties.densities import get_mass_energy_density, get_volumetric_energy_density, get_specific_capacity
from properties.voltages import get_eq_voltage

if __name__ == '__main__':
    # set high quality plotting
    mpl.rcParams['figure.dpi'] = 300  # increase plot dpi

    # get compounds
    db = ase.db.connect('db/compounds.db')
    xc = 'PBE'

    Li2O = Compound('Li2O', 'mp-1960', db, xc, nbands=20, converged=True)
    Li2O2 = Compound('Li2O2', 'mp-841', db, xc, nkpts=6, ecut=900, U_correction={'O': ':p,0.96,0'}, converged=True)
    LiO2 = Compound('LiO2', 'mp-1018789', db, xc, nkpts=4, ecut=500, U_correction={'O': ':p,0.33,0'}, converged=True)
    Li = Compound('Li', '1', db, xc, nbands=-10, converged=True)
    O2 = Compound('O2', 'mp-12957', db, xc, nkpts=4, ecut=900, U_correction={'O': ':p,0.75,0'}, converged=True)

    # get potential energies
    epot_Li2O_cell = Li2O.set_energy()
    epot_Li2O2_cell = Li2O2.set_energy()
    epot_LiO2_cell = LiO2.set_energy()
    epot_O2_cell = O2.set_energy()
    epot_Li_cell = Li.set_energy()

    print(epot_Li2O_cell)
    print(epot_Li2O2_cell)
    print(epot_LiO2_cell)
    print(epot_O2_cell)
    print(epot_Li_cell)

    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_Li2O = epot_Li2O_cell
    epot_LiO2 = epot_LiO2_cell / 2

    # get properties
    vol1 = get_eq_voltage(2 * epot_Li2O, epot_Li2O2, 2)
    vol2 = get_eq_voltage(2 * epot_Li2O, epot_LiO2, 3)
    vol = (vol1+vol2)/2
    x = 1.5

    get_mass_energy_density(vol, MOL_MASS_LI2O, x, watt_hour=True)
    get_volumetric_energy_density(vol, MOL_VOL_LI2O, x, watt_hour=True)
    get_specific_capacity(MOL_MASS_LI2O, x)

    # get electronic structure
    for compound in [Li2O, Li2O2, LiO2]:
        compound.get_bandgap()
        compound.set_pdos()
        compound.set_band_structure(emax=13)
        compound.set_dos()
        compound.set_ldos()

    # ####
    # Li2O = Compound('Li2O', 'mp-1960', db, xc, nkpts=8, ecut=575, nbands=20, converged=True, extension='cif')
    # Li2O2 = Compound('Li2O2', 'mp-841', db, xc, nkpts=6, ecut=903, U_correction={'O': ':p,0.55,0'}, converged=True)
    # LiO2 = Compound('LiO2', 'mp-1018789', db, xc, nkpts=4, ecut=500, U_correction={'O': ':p,0.36,0'}, converged=True)
    # Li = Compound('Li', 'mp-1', db, xc, nbands=-10, converged=True)
    # O2 = Compound('O2', 'mp-12957', db, xc, nkpts=4, ecut=900, U_correction={'O': ':p,0.75,0'}, converged=True)
    #
    # epot_Li2O_cell = Li2O.set_energy()
    # epot_Li2O2_cell = Li2O2.set_energy()
    # epot_LiO2_cell = LiO2.set_energy()
    # epot_O2_cell = O2.set_energy()
    # epot_Li_cell = Li.set_energy()
    #
    # epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    # epot_Li2O = epot_Li2O_cell
    # epot_LiO2 = epot_LiO2_cell / 2
    #
    #
    # def get_formation_energy(compound_epot, x, y):
    #     epot_Li = epot_Li_cell / 4
    #
    #     epot_O_cell = O2.set_energy()
    #     epot_O = epot_O_cell / 4
    #
    #     return compound_epot - (x * epot_Li + y * epot_O)
    #
    #
    # ef_li2O = get_formation_energy(epot_Li2O, 2, 1)
    # print('ef li2o:', ef_li2O / 3)  # goal: -2.067 eV
    #
    # ef_li2O2 = get_formation_energy(epot_Li2O2, 2, 2)
    # print('ef li2o2:', ef_li2O2 / 4)  # goal: -1.651 eV
    #
    # ef_liO2 = get_formation_energy(epot_LiO2, 1, 2)
    # print('ef lio2:', ef_liO2 / 3)  # goal: -1.033 eV
    #
    # ####