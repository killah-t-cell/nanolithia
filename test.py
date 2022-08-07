# So methods like DFT+vdW correct (or try to) for the lack of dispersion forces, and methods like DFT+U try to
# correct for the lack of correlation
import os

# Li2O is linear
# Li2O2 is nonlinear
# LiO2 is nonlinear
from ase import Atoms
from ase.io import read
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
from gpaw import GPAW, PW
nkpts = 4
ecut = 550
xc = 'PBE'
parameters = dict(mode=PW(ecut), kpts=(nkpts, nkpts, nkpts), xc=xc, symmetry='off')

positions = [[0., 0., 3.851761], [0., 0., 0.], [1.590181, -0.91809317, 5.7776415], [1.590181, 0.91809317, 1.9258805],
             [1.590181, -0.91809317, 1.15021287], [1.590181, 0.91809317, 5.00197387],
             [1.590181, 0.91809317, 6.55330913], [1.590181, -0.91809317, 2.70154813]]
Li2O2 = Atoms(symbols='Li4O4', pbc=True,
              cell=[[1.590181, -2.754274, 0.0], [1.590181, 2.754274, 0.0], [0.0, 0.0, 7.703522]], positions=positions)

calc = GPAW(**parameters)
Li2O2.calc = calc
potentialenergy = Li2O2.get_potential_energy()

vib = Vibrations(Li2O2)
vib.run()
vib_energies = vib.get_energies(method='Frederiksen') # TODO Frederiksen vs standard?
thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=Li2O2,
                        geometry='nonlinear') # TODO add geometry data to compound
H = thermo.get_enthalpy(temperature=298.15)


## THIS WORKS
# nkpts = 4
# ecut = 100
# xc = 'PBE'
# parameters = dict(mode=PW(ecut), kpts=(nkpts, nkpts, nkpts), xc=xc, symmetry='off')
#
# positions = [[0., 0., 3.851761], [0., 0., 0.], [1.590181, -0.91809317, 5.7776415], [1.590181, 0.91809317, 1.9258805],
#              [1.590181, -0.91809317, 1.15021287], [1.590181, 0.91809317, 5.00197387],
#              [1.590181, 0.91809317, 6.55330913], [1.590181, -0.91809317, 2.70154813]]
# Li2O2 = Atoms(symbols='Li4O4', pbc=True,
#               cell=[[1.590181, -2.754274, 0.0], [1.590181, 2.754274, 0.0], [0.0, 0.0, 7.703522]], positions=positions)
#
# calc = GPAW(**parameters)
# Li2O2.calc = calc
# potentialenergy = Li2O2.get_potential_energy()
#
# vib = Vibrations(Li2O2)
# vib.run()
# vib_energies = vib.get_energies(method='Frederiksen')
# thermo = IdealGasThermo(vib_energies=vib_energies,
#                         potentialenergy=potentialenergy,
#                         atoms=Li2O2,
#                         geometry='linear',
#                         symmetrynumber=2, spin=0)
# H = thermo.get_enthalpy(temperature=298.15)

# if not os.path.exists('vib_test_Li2O2.gpw'):
#     nkpts = 4
#     ecut = 100
#     xc = 'PBE'
#     parameters = dict(mode=PW(ecut), kpts=(nkpts, nkpts, nkpts), xc=xc)
#     Li2O2 = read('structures/Li2O2_mp-841.poscar')
#     calc = GPAW(**parameters)
#     Li2O2.calc = calc
#
#     potentialenergy = Li2O2.get_potential_energy()
#     calc.write('vib_test_Li2O2.gpw')
# else:
#     calc = GPAW('vib_test_Li2O2.gpw')
#     Li2O2 = calc.get_atoms()
#     potentialenergy = Li2O2.get_potential_energy()
#
#
# vib = Vibrations(Li2O2)
# vib.run()
# vib_energies = vib.get_energies()
# thermo = IdealGasThermo(vib_energies=vib_energies,
#                         potentialenergy=potentialenergy,
#                         atoms=Li2O2,
#                         geometry='linear',
#                         symmetrynumber=2, spin=0)
# H = thermo.get_enthalpy(temperature=298.15)

####
# nkpts = 6
# ecut = 500
# xc = 'PBE'
# parameters = dict(mode=PW(ecut), kpts=(nkpts, nkpts, nkpts), xc=xc)
#
# positions = [[0., 0., 3.851761],[0., 0., 0.],[1.590181, -0.91809317, 5.7776415],[1.590181, 0.91809317, 1.9258805],[1.590181, -0.91809317, 1.15021287],[1.590181, 0.91809317, 5.00197387],[1.590181, 0.91809317, 6.55330913],[1.590181, -0.91809317, 2.70154813]]
# Li2O2 = Atoms(symbols='Li4O4', pbc=True,
#               cell=[[1.590181, -2.754274, 0.0], [1.590181, 2.754274, 0.0], [0.0, 0.0, 7.703522]], positions=positions)
#
# calc = GPAW(**parameters)
# Li2O2.calc = calc
# potentialenergy = Li2O2.get_potential_energy()
# calc.write('vib_test_Li2O2.gpw')
#
# calc = GPAW('vib_test_Li2O2.gpw')
# Li2O2 = calc.get_atoms()
# vib = Vibrations(Li2O2)
# vib.run()
# vib_energies = vib.get_energies()
#
# thermo = IdealGasThermo(vib_energies=vib_energies,
#                         potentialenergy=potentialenergy,
#                         atoms=Li2O2,
#                         geometry='linear',
#                         symmetrynumber=2, spin=0)
# G = thermo.get_enthalpy(temperature=298.15)

# ### EXAMPLE
# atoms = molecule('N2')
# atoms.calc = EMT()
# dyn = QuasiNewton(atoms)
# dyn.run(fmax=0.01)
# potentialenergy = atoms.get_potential_energy()
#
# vib = Vibrations(atoms)
# vib.run()
# vib_energies = vib.get_energies(read_cache=False)
#
# thermo = IdealGasThermo(vib_energies=vib_energies,
#                         potentialenergy=potentialenergy,
#                         atoms=atoms,
#                         geometry='linear',
#                         symmetrynumber=2, spin=0)
# G = thermo.get_enthalpy(temperature=298.15)
