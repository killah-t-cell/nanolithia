import os

import numpy as np
from ase.calculators.emt import EMT
from ase.dft import DOS
from ase.dft.bandgap import bandgap
from ase.io import read
from ase.phonons import Phonons
from gpaw import GPAW, PW
from matplotlib import pyplot as plt

from global_vars import ROOT_DIR


class Compound:
    def __init__(self, formula, structure, xc, kpts=None, ecut=540, nbands=None, converged=False, tol='null'
                 , setups=None, magmoms=None, extension='.poscar'):
        if kpts is None:
            kpts = {'size': (5, 5, 7), 'gamma': True}
        self.ef = None
        self.bs = None
        self.bandgap = None
        self.dos_weights = None
        self.pot_energy = None
        self.dos_energies = None
        self.magmoms = magmoms
        self.formula = formula
        self.kpts = kpts
        self.ecut = ecut
        self.nbands = nbands
        self.xc = xc
        self.tol = tol
        self.converged = converged
        self.structure = structure
        self.setups = setups
        self.name = f'{formula}-{structure}-{xc}-{ecut:.0f}-{kpts}'
        # load crystal
        self.atoms = read(os.path.join(ROOT_DIR, f'structures/{formula}_{structure}{extension}'))

        if self.magmoms is not None:
            self.atoms.set_initial_magnetic_moments(self.magmoms)

        # attach calculator
        self.calc = GPAW(mode=PW(ecut), kpts=kpts, xc=xc, txt=os.path.join(ROOT_DIR, f'logs/{self.name}.txt'))

        if self.setups is not None:
            self.calc.set(setups=self.setups)

        if self.nbands is not None:
            self.calc.set(nbands=nbands)

        self.volume = self.atoms.get_volume()
        self.atoms.calc = self.calc

    ## --------------- FUNDAMENTAL CALCULATIONS --------------- ##
    def get_energy(self):
        if os.path.exists(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw')):
            self.calc = GPAW(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw'))
            self.ef = self.calc.get_fermi_level()
            self.atoms = self.calc.get_atoms()
            self.pot_energy = self.atoms.get_potential_energy()
        else:
            # get potential energy
            self.pot_energy = self.atoms.get_potential_energy()
            self.ef = self.calc.get_fermi_level()
            self.calc.write(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw'))

        return self.pot_energy

    def converge(self):
        return

    ## --------------- ELECTRONIC STRUCTURES --------------- ##
    def get_bandgap(self, exact_exchange=False, direct=False):
        self.get_energy()
        self.bandgap = bandgap(self.calc, direct=direct)
        #
        # energies = non_self_consistent_energy('foo.gpaw', xcname = 'HSE06')
        # bandgap = bandgap('foo.gpaw', direct=False)

        return self.bandgap


    def set_dos(self):
        self.get_energy()
        dos = DOS(self.calc, npts=800, width=0.2)
        self.dos_energies = dos.get_energies()
        self.dos_weights = dos.get_dos()

        plt.plot(self.dos_energies - self.ef, self.dos_weights)
        plt.xlabel('Energy (eV)')
        plt.ylabel(f'{self.formula}_DOS')
        plt.savefig(os.path.join(ROOT_DIR, f'plots/{self.name}-DOS.png'))
        plt.show()

        return self.dos_energies, self.dos_weights

    def set_pdos(self, npts=201):
        self.get_energy()

        for atom in [0, -1]:
            energy, pdos = self.calc.get_orbital_ldos(a=atom, spin=0, angular='spd', npts=npts, width=None)
            I = np.trapz(pdos, energy)
            center = np.trapz(pdos * energy, energy) / I
            width = np.sqrt(np.trapz(pdos * (energy - center) ** 2, energy) / I)
            plt.plot(energy - self.ef, pdos, label=f'{self.atoms[atom].symbol}', lw=2, alpha=0.7)

        plt.legend(loc='best')
        plt.xlabel('Energy (eV)')
        plt.ylabel(f'projected {self.formula}-DOS on atoms')
        plt.savefig(os.path.join(ROOT_DIR, f'plots/{self.name}-PDOS.png'))
        plt.show()

        return center, width

    def set_ldos(self, npts=201):
        self.get_energy()

        for atom in [0, -1]:  # 0 is Li, -1 is O
            for c in 'spd':
                energies, ldos = self.calc.get_orbital_ldos(a=0, spin=0, angular=c, width=0.4)
                plt.plot(energies - self.ef, ldos, label=f'{self.atoms[atom].symbol}-{c}' + '-up')

                energies, ldos = self.calc.get_orbital_ldos(a=0, spin=1, angular=c, width=0.4)
                plt.plot(energies - self.ef, ldos, label=f'{self.atoms[atom].symbol}-{c}' + '-down')

        plt.legend(loc='best')
        plt.xlabel('Energy (eV)')
        plt.ylabel(f'Localized {self.formula}-DOS on atoms')
        plt.savefig(os.path.join(ROOT_DIR, f'plots/{self.name}-LDOS.png'))
        plt.show()

    def set_band_structure(self, symmetry='off', emax=10.0):
        self.get_energy()

        # set prerequisites
        # lat = self.atoms.cell.get_bravais_lattice()

        path = self.atoms.cell.bandpath(density=7)
        self.calc.set(kpts=path, fixdensity=True,
                      symmetry=symmetry)
        self.atoms.get_potential_energy()

        # get and save bs
        self.bs = self.calc.band_structure()
        self.bs.write(os.path.join(ROOT_DIR, f'bs/{self.name}_bs.json'))
        self.bs.plot(filename=os.path.join(ROOT_DIR, f'plots/{self.name}-band-structure.png'), show=True, emax=emax)

    ## --------------- PHONON CALCULATIONS --------------- ##
    def get_phonons(self, N=2):
        # Phonon calculator
        ph = Phonons(self.atoms, self.calc, supercell=(N, N, N), delta=0.05)
        ph.run()

        # Read forces and assemble the dynamical matrix
        ph.read(acoustic=True)
        ph.clean()

        path = self.atoms.cell.bandpath(density=7)
        bs = ph.get_band_structure(path)

        dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

        # Plot the band structure and DOS:
        fig = plt.figure(1, figsize=(7, 4))
        ax = fig.add_axes([.12, .07, .67, .85])

        emax = 0.035
        bs.plot(ax=ax, emin=0.0, emax=emax)

        dosax = fig.add_axes([.8, .07, .17, .85])
        dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color='grey',
                           edgecolor='k', lw=1)

        dosax.set_ylim(0, emax)
        dosax.set_yticks([])
        dosax.set_xticks([])
        dosax.set_xlabel("DOS", fontsize=18)

        fig.savefig(os.path.join(ROOT_DIR, f'plots/{self.name}-phonon.png'))
