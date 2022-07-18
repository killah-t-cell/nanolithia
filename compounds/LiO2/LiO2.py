import pathlib
from os.path import exists
from ase.io import read, write
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3


def get_LiO2(db, xc, nkpts=8, ecut=500, converged=False, tol='null'):
    """Define a LiO2 crystal and save it to the database, if it hasn't already been saved

        db: Database
            Database for collecting results.
        xc: Exchange correlation functional
            Functional to be used for calculations
        kpts: int
            Use a (kpts * kpts * kpts) Monkhorst-Pack grid.
        ecut: float
            Cutoff energy for plane waves.

    Returns LiO2 crystal (which is mostly an unstable crystal).
    """
    name = f'LiO2-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'
    U_correction = {'O': ':d,0.33'}

    parameters = dict(mode=PW(ecut),
                      kpts={'size': (nkpts, nkpts, nkpts)},
                      setups=U_correction,
                      xc='PBE')

    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        # load crystal
        LiO2 = read(pathlib.Path(__file__).parent / 'LiO2.poscar')

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

        LiO2.calc = calc

        # get potential energy
        LiO2.get_potential_energy()

        # save
        del db[id]
        db.write(LiO2,
                 name=name,
                 xc=xc,
                 nkpts=nkpts,
                 ecut=ecut,
                 relaxed=True,
                 converged=converged,
                 tol=tol)
        return LiO2
    else:
        return db.get_atoms(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
