import os

import numpy as np
from ase.calculators.dftd3 import DFTD3
from ase.dft import DOS
from ase.dft.bandgap import bandgap
from ase.io import read
from gpaw import GPAW, PW
from gpaw.hybrids.energy import non_self_consistent_energy
from matplotlib import pyplot as plt
from ase.spacegroup import get_spacegroup

from global_vars import ROOT_DIR


class Compound:
    def __init__(self, formula, structure, db, xc, nkpts=8, ecut=500, nbands=None, converged=False, tol='null',
                 spinpol=False, U_correction=None, extension='poscar', primitive=True):
        self.primitive = primitive
        self.bandgap = None
        self.ef = None
        self.id = None
        self.bs = None
        self.dos_weights = None
        self.pot_energy = None
        self.dos_energies = None
        self.formula = formula
        self.nkpts = nkpts
        self.ecut = ecut
        self.nbands = nbands
        self.xc = xc
        self.db = db
        self.tol = tol
        self.converged = converged
        self.structure = structure
        self.U_correction = U_correction

        self.name = f'{formula}-{structure}-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'
        self.parameters = dict(kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                               spinpol=spinpol,
                               convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
                                            'energy': 2.0e-4,  # eV / electron
                                            'density': 1.0e-3, },
                               xc=xc)
        # load crystal
        self.atoms = read(os.path.join(ROOT_DIR, f'structures/{formula}_{structure}.{extension}'))

        # attach calculator
        self.calc = GPAW(mode=PW(ecut), txt=os.path.join(ROOT_DIR, f'logs/{self.name}.txt'),**self.parameters)

        if self.U_correction is not None:
            self.calc.set(setups=self.U_correction)

        if self.nbands is not None:
            self.calc.set(nbands=nbands)

        self.atoms.calc = self.calc

    def set_energy(self):
        self.id = self.db.reserve(name=self.name, xc=self.xc, nkpts=self.nkpts, ecut=self.ecut)
        if self.id is not None:  # skip calculation if already done
            # get potential energy
            self.pot_energy = self.atoms.get_potential_energy()
            self.ef = self.calc.get_fermi_level()
            self.calc.write(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw'))
            # save
            del self.db[self.id]
            self.db.write(self.atoms,
                          name=self.name,
                          xc=self.xc,
                          nkpts=self.nkpts,
                          ecut=self.ecut,
                          relaxed=True,
                          calc_parameters=str(self.parameters),
                          converged=self.converged,
                          structure=self.structure,
                          primitive=self.primitive,
                          tol=self.tol)
        else:
            self.calc = GPAW(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw'))
            self.ef = self.calc.get_fermi_level()
            self.atoms = self.calc.get_atoms()
            self.pot_energy = self.atoms.get_potential_energy()

        return self.pot_energy

    def get_bandgap(self, direct=False):
        self.set_energy()
        self.bandgap = bandgap(self.calc, direct=direct)
        return self.bandgap

    def set_dos(self):
        self.set_energy()
        dos = DOS(self.calc, npts=800, width=0)
        self.dos_energies = dos.get_energies()
        self.dos_weights = dos.get_dos()

        plt.plot(self.dos_energies, self.dos_weights)
        plt.xlabel('Energy (eV)')
        plt.ylabel(f'{self.formula}_DOS')
        plt.savefig(os.path.join(ROOT_DIR, f'plots/{self.formula}-{self.xc}-{self.nkpts}-{self.ecut}-DOS.png'))
        plt.show()

        return self.dos_energies, self.dos_weights

    def set_pdos(self, npts=201):
        self.set_energy()

        for atom in [0, -1]:
            energy, pdos = self.calc.get_orbital_ldos(a=atom, spin=0, angular='spd', npts=npts, width=None)
            I = np.trapz(pdos, energy)
            center = np.trapz(pdos * energy, energy) / I
            width = np.sqrt(np.trapz(pdos * (energy - center) ** 2, energy) / I)
            plt.plot(energy, pdos, label=f'{self.atoms[atom].symbol}', lw=2, alpha=0.7)

        plt.legend(loc='best')
        plt.xlabel('Energy (eV)')
        plt.ylabel(f'projected {self.formula}-DOS on atoms')
        plt.savefig(os.path.join(ROOT_DIR, f'plots/{self.formula}-{self.xc}-{self.nkpts}-{self.ecut}-PDOS.png'))
        plt.show()

        return center, width

    def set_ldos(self, npts=201):
        self.set_energy()

        for atom in [0, -1]:  # 0 is Li, -1 is O
            for orbit in ['s', 'p']:
                energy, pdos = self.calc.get_orbital_ldos(a=atom, spin=0, angular=orbit, npts=npts, width=None)
                plt.plot(energy, pdos, label=f'{self.atoms[atom].symbol}-{orbit}', lw=2, alpha=0.7)

        plt.legend(loc='best')
        plt.xlabel('Energy (eV)')
        plt.ylabel(f'Localized {self.formula}-DOS on atoms')
        plt.savefig(os.path.join(ROOT_DIR, f'plots/{self.formula}-{self.xc}-{self.nkpts}-{self.ecut}-LDOS.png'))
        plt.show()

    def set_band_structure(self, symmetry='off', emax=10.0):
        self.set_energy()

        # set band structure calc
        path = self.atoms.cell.bandpath(density=7)
        self.calc.set(kpts=path, fixdensity=True,
                      symmetry=symmetry)
        self.atoms.get_potential_energy()

        # get and save bs
        self.bs = self.calc.band_structure()
        self.bs.write(os.path.join(ROOT_DIR, f'bs/{self.name}_bs.json'))
        self.bs.plot(filename=os.path.join(ROOT_DIR, f'plots/{self.name}-band-structure.png'), show=True, emax=emax)
