import pathlib

from gpaw import GPAW, FermiDirac, PW
from ase.optimize import QuasiNewton
from ase.build import bulk
from ase.calculators.dftd3 import DFTD3
from ase.constraints import StrainFilter
from ase.optimize.bfgs import BFGS
import ase.db
from gpaw.xc.bee import BEEFEnsemble


def get_Li(db, xc, nkpts=8, ecut=500, nbands=-10, converged=False, tol='null', structure='bcc'):
    """Define a lithium crystal and save it to the database, if it hasn't already been saved

        db: Database
            Database for collecting results.
        xc: Exchange correlation functional
            Functional to be used for calculations
        kpts: int
            Use a (kpts * kpts * kpts) Monkhorst-Pack grid.
        ecut: float
            Cutoff energy for plane waves.
        nbands: number of bands
            Number of bands to be added to calculator

    Returns Li crystal.
    """
    name = f'Li-{structure}-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'

    parameters = dict(mode=PW(ecut),
                      nbands=nbands,
                      kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                      spinpol=True,
                      convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
                                  'energy': 2.0e-4,  # eV / electron
                                  'density': 1.0e-3, },
                      xc=xc)

    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        Li = bulk('Li', crystalstructure=structure, a=3.51, cubic=True)

        # attach calculator
        if xc == 'DFTD3':
            dft = GPAW(mode=PW(ecut),
                       kpts=nkpts,
                       nbands=nbands,
                       txt=name + '.txt',
                       xc='PBE')
            calc = DFTD3(dft=dft, xc='PBE')
        else:
            calc = GPAW(txt=name + '.txt',
                        **parameters)
        Li.calc = calc

        # relax Li
        sf = StrainFilter(Li, mask=[1, 1, 1, 0, 0, 0])
        opt = BFGS(sf)
        opt.run(fmax=0.01)

        # get potential energy
        Li.get_potential_energy()

        # save
        del db[id]
        db.write(Li,
                 name=name,
                 xc=xc,
                 nkpts=nkpts,
                 ecut=ecut,
                 relaxed=True,
                 converged=converged,
                 structure=structure,
                 tol=tol)

    return db.get(name=name, xc=xc, nkpts=nkpts, ecut=ecut, structure=structure)



