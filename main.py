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

    # del db[db.get(formula='Li2O4', ecut=500, xc=xc, nkpts=4).id]

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

    print(epot_Li2O_cell)
    print(epot_Li2O2_cell)
    print(epot_LiO2_cell)

    # The calculated energies are for the full cells. Convert them to the energy per formula unit.
    epot_Li2O2 = epot_Li2O2_cell / 2  # len is 8, we want to get Li2O2 out of Li4O4, so we divide by 2. 4/8=0.5
    epot_Li2O = epot_Li2O_cell / 4  # len is 12, we want to get Li2O out of Li8O4, so we divide by 4. 3/12=0.25
    epot_LiO2 = epot_LiO2_cell / 2

    print(epot_Li2O)
    print(epot_Li2O2)
    print(epot_LiO2)

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
        compound.set_pdos()
        compound.set_band_structure(emax=13)
        compound.set_dos()
        compound.set_ldos()

