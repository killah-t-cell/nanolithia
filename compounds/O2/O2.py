import pathlib
from ase.io import read
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3

def get_O2(db, xc, nkpts=8, ecut=500, converged=False, tol='null', structure='mp-12957'):
    """Define a O2 crystal and save it to the database, if it hasn't already been saved

        db: Database
            Database for collecting results.
        xc: Exchange correlation functional
            Functional to be used for calculations
        kpts: int
            Use a (kpts * kpts * kpts) Monkhorst-Pack grid.
        ecut: float
            Cutoff energy for plane waves.

    Returns O2 crystal
    """
    name = f'O2-{structure}-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'
    U_correction = {'O': ':p,1.05,0'}  # because it's an oxide

    parameters = dict(mode=PW(ecut),
                      kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                      setups=U_correction,
                      xc=xc)

    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        # load crystal
        O2 = read(pathlib.Path(__file__).parent / f'O2_{structure}.poscar')

        # attach calculator
        if xc == 'DFTD3':
            dft = GPAW(mode=PW(ecut),
                       kpts=nkpts,
                       txt=name + '.txt',
                       setups=U_correction,
                       xc='PBE')
            calc = DFTD3(dft=dft, xc='PBE')
        else:
            calc = GPAW(txt=name + '.txt',
                        **parameters)

        O2.calc = calc

        # get potential energy
        O2.get_potential_energy()

        # save
        del db[id]
        db.write(O2,
                 name=name,
                 xc=xc,
                 nkpts=nkpts,
                 ecut=ecut,
                 relaxed=True,
                 converged=converged,
                 structure=structure,
                 tol=tol)

    return db.get(name=name, xc=xc, nkpts=nkpts, ecut=ecut, structure=structure)