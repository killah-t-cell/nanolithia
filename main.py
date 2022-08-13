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
from properties.profiles import gibbs_energy, LI2O_ENTROPY, LI2O2_ENTROPY, LIO2_ENTROPY

if __name__ == '__main__':
    # set high quality plotting
    mpl.rcParams['figure.dpi'] = 300  # increase plot dpi

    # get compounds
    xc = 'PBE'
    #  TODO get right band gap with HSE06 values of the
    #  Corrections are slightly adapted from A Facile Mechanism for Recharging
    #  Li2O2 in Li−O2 Batteries paper. The adaptation is to be expected because DFT+U is implementation specific
    Li2O = Compound('Li2O', 'mp-1960', xc, magmoms=[0.6, 0.6, 0.6], converged=True, setups={'O': ':p, -1.04,0'}, extension='cif')
    Li2O2 = Compound('Li2O2', 'mp-841', xc, magmoms=[0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6], setups={'O': ':p, -0.72,0'}, converged=True)
    LiO2 = Compound('LiO2', 'mp-1018789', xc, magmoms=[0.6, 0.6, 0.6, 0.6, 0.6, 0.6],setups={'O': ':p, -0.24,0'}, converged=True)
    Li = Compound('Li', 'mp-1', xc, ecut=530, converged=True)
    O2 = Compound('O2', 'mp-12957', xc, magmoms=[0.6, 0.6, 0.6, 0.6], ecut=530, converged=True)

    # get potential energies
    epot_Li2O_cell = Li2O.get_energy()
    epot_Li2O2_cell = Li2O2.get_energy()
    epot_LiO2_cell = LiO2.get_energy()
    epot_O2_cell = O2.get_energy()

    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li2O2 = epot_Li2O2_cell / 2
    epot_Li2O = epot_Li2O_cell
    epot_LiO2 = epot_LiO2_cell / 2

    print('epot_Li2O', epot_Li2O)
    print('epot_Li2O2', epot_Li2O2_cell)
    print('epot_LiO2', epot_LiO2_cell)

    # get properties
    Li2O_gibbs = gibbs_energy(epot_Li2O, LI2O_ENTROPY, 2, 1)
    Li2O2_gibbs = gibbs_energy(epot_Li2O2, LI2O2_ENTROPY, 2, 2)
    LiO2_gibbs = gibbs_energy(epot_LiO2, LIO2_ENTROPY, 1, 2)
    vol1 = get_eq_voltage(2 * Li2O_gibbs, Li2O2_gibbs, 2)
    vol2 = get_eq_voltage(2 * Li2O_gibbs, LiO2_gibbs, 3)
    vol = (vol1+vol2)/2
    x = 1.5

    get_mass_energy_density(vol, MOL_MASS_LI2O, x, watt_hour=True)
    get_volumetric_energy_density(vol, MOL_VOL_LI2O, x, watt_hour=True)
    get_specific_capacity(MOL_MASS_LI2O, x)

    # get electronic structure
    for compound in [Li2O, Li2O2, LiO2]:
        compound.get_bandgap()
        compound.set_pdos()
        compound.set_ldos()
        compound.set_band_structure(emax=13)
        compound.set_dos()

