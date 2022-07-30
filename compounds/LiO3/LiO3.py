import pathlib
from os.path import exists
from ase.io import read, write
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3


def get_LiO3(db, xc, nkpts=8, ecut=500, converged=False, tol='null', structure='mp-1001790'):
    """Define a LiO3 crystal and save it to the database, if it hasn't already been saved

        db: Database
            Database for collecting results.
        xc: Exchange correlation functional
            Functional to be used for calculations
        kpts: int
            Use a (kpts * kpts * kpts) Monkhorst-Pack grid.
        ecut: float
            Cutoff energy for plane waves.

    Returns LiO3 crystal
    """
    name = f'LiO3-{structure}-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'
    U_correction = {'O': ':p,0.33,0'} # got this from the correction {'O': ':p,0.33,0'} for superoxides.

    parameters = dict(kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                      spinpol=True,
                      convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
                                  'energy': 2.0e-4,  # eV / electron
                                  'density': 1.0e-3, },
                      setups=U_correction,
                      xc=xc)

    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        # load crystal
        LiO3 = read(pathlib.Path(__file__).parent / f'LiO3_{structure}.poscar')

        # attach calculator
        if xc == 'DFTD3':
            dft = GPAW(mode=PW(ecut),
                       kpts=nkpts,
                       txt=name + '.txt',
                       setups=U_correction,
                       xc='PBE')
            calc = DFTD3(dft=dft, xc='PBE')
        else:
            calc = GPAW(mode=PW(ecut),txt=name + '.txt',
                        **parameters)

        LiO3.calc = calc

        # get potential energy
        LiO3.get_potential_energy()

        # save
        del db[id]
        db.write(LiO3,
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