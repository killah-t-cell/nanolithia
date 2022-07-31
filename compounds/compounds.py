import os
import pathlib
import sys

import ase.db
from ase.dft import DOS
from ase.io import read, write
from gpaw import GPAW, FermiDirac, PW, restart
from ase.calculators.dftd3 import DFTD3

from global_vars import ROOT_DIR


class Compound:
    def __init__(self, formula, structure, db, xc, nkpts=8, ecut=500, nbands=None, converged=False, tol='null',
                 spinpol=False, U_correction=None):
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
        self.atoms = read(os.path.join(ROOT_DIR, f'structures/{formula}_{structure}.poscar'))
        os.path.join(ROOT_DIR, f'logs/{self.name}.txt'),

        # attach calculator
        if xc == 'DFTD3':
            dft = GPAW(mode=PW(ecut),
                       kpts=nkpts,
                       nbands=nbands,
                       txt=os.path.join(ROOT_DIR, f'logs/{self.name}.txt'),
                       xc='PBE')
            self.calc = DFTD3(dft=dft, xc='PBE')
        else:
            self.calc = GPAW(mode=PW(ecut), txt=os.path.join(ROOT_DIR, f'logs/{self.name}.txt'),
                             **self.parameters)

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
                          tol=self.tol)
        else:
            self.calc = GPAW(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw'))
            self.atoms = self.calc.get_atoms()
            self.pot_energy = self.atoms.get_potential_energy()

        return self.pot_energy

    def set_dos(self):
        self.set_energy()
        dos = DOS(self.calc, npts=800, width=0)
        self.dos_energies = dos.get_energies()
        self.dos_weights = dos.get_dos()

    def set_pdos(self):
        return

    def set_ldos(self):
        return

    def set_band_structure(self, symmetry='off', path='GXWKL', npoints=60, emax=10.0):
        self.set_energy()
        GPAW(os.path.join(ROOT_DIR, f'calculators/{self.name}.gpw')).fixed_density(
            nbands=self.nbands,
            symmetry=symmetry,
            kpts={'path': path, 'npoints': npoints},
            convergence={'bands': self.nbands/2})

        self.bs = self.calc.band_structure()
        self.bs.write(os.path.join(ROOT_DIR, f'bs/{self.name}_bs.json'))

        self.bs.plot(filename=os.path.join(ROOT_DIR, f'plots/{self.name}-band-structure.png'), show=True, emax=emax)
