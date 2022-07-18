import pathlib
from ase.io import read, write
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3


def get_Li2O2(db, xc, nkpts=8, ecut=500, converged=False, tol='null'):
    """Define a Li2O2 crystal and save it to the database, if it hasn't already been saved

        db: Database
            Database for collecting results.
        xc: Exchange correlation functional
            Functional to be used for calculations
        kpts: int
            Use a (kpts * kpts * kpts) Monkhorst-Pack grid.
        ecut: float
            Cutoff energy for plane waves.

    Returns Li202 crystal.
    """
    name = f'Li2O2-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'
    U_correction = {'O': ':d,0.76'}

    parameters = dict(mode=PW(ecut),
                      kpts={'size': (nkpts, nkpts, nkpts)},
                      setups=U_correction,
                      xc='PBE')

    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        # load crystal
        Li2O2 = read(pathlib.Path(__file__).parent / 'Li2O2.poscar')

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

        Li2O2.calc = calc

        # get potential energy
        Li2O2.get_potential_energy()

        # save
        del db[id]
        db.write(Li2O2,
                 name=name,
                 xc=xc,
                 nkpts=nkpts,
                 ecut=ecut,
                 relaxed=True,
                 converged=converged,
                 tol=tol)
        return Li2O2
    else:
        return db.get_atoms(name=name, xc=xc, nkpts=nkpts, ecut=ecut)

