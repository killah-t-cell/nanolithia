import pathlib
from ase.io import read, write
from gpaw import GPAW, FermiDirac, PW
from ase.calculators.dftd3 import DFTD3


def get_Li2O(db, xc, nkpts=8, ecut=500, nbands=20, converged=False, tol='null', structure='mp-1960', spinpol=False):
    """Define a Li2O crystal and save it to the database, if it hasn't already been saved

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

    Returns Li2O crystal.
    """
    name = f'Li2O-{structure}-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'

    parameters = dict(nbands=nbands,
                      kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                      spinpol=spinpol,
                      convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
                                   'energy': 2.0e-4,  # eV / electron
                                   'density': 1.0e-3, },
                      xc=xc)

    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        # load crystal
        Li2O = read(pathlib.Path(__file__).parent / f'Li2O_{structure}.poscar')  # mp-1960

        # attach calculator
        if xc == 'DFTD3':
            dft = GPAW(mode=PW(ecut),
                       kpts=nkpts,
                       nbands=nbands,
                       txt=name + '.txt',
                       xc='PBE')
            calc = DFTD3(dft=dft, xc='PBE')
        else:
            calc = GPAW(mode=PW(ecut), txt=name + '.txt',
                        **parameters)

        Li2O.calc = calc

        # get potential energy
        Li2O.get_potential_energy()

        # save
        del db[id]
        db.write(Li2O,
                 name=name,
                 xc=xc,
                 nkpts=nkpts,
                 ecut=ecut,
                 relaxed=True,
                 calc_parameters=str(parameters),
                 converged=converged,
                 structure=structure,
                 tol=tol)

    return db.get(name=name, xc=xc, nkpts=nkpts, ecut=ecut, structure=structure)
