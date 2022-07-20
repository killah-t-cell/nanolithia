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
    U_correction = {'O': ':p,0.33,0'} # the original paper has this at {'O': ':p,0.33,0'} for superoxides.
                                      # 0.76 with ecut=575 works too

    parameters = dict(mode=PW(ecut),
                      kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                      spinpol=True,
                      convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
                                  'energy': 2.0e-4,  # eV / electron
                                  'density': 1.0e-3, },
                      setups=U_correction,
                      xc=xc)

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


# nkpts=8
# ecut=500
# U_correction = {'O': ':p,0.33,0'}
# parameters = dict(mode=PW(ecut),kpts={'size': (nkpts, nkpts, nkpts)}, setups=U_correction, xc='PBE')
# LiO2 = read(pathlib.Path(__file__).parent / 'LiO2.poscar')
# calc = GPAW(**parameters)
# LiO2.calc = calc
# # get potential energy
# e = LiO2.get_potential_energy()
# e


# goal -54.744
# 575, 8, PBE -> -50.25838214237331

# currently running without setups, ecut=900, k=8 -> -53.043477

# goal -13.686
# new structure, with setups 0.33 -> -14.36
# new structure, without setups -> -14.8
# new structure, with setups 0.33,ecut=500 -> -13.67 !!!!!! THIS WORKS !!!!!!
# new structure, with setups 0.76,ecut=500 -> -13.07
# new structure, with setups 0.76,ecut=575 -> -13.67 !!!!!! THIS WORKS !!!!!!
# new structure, with setups 0.76,ecut=600 -> -13.76
# new structure, with setups 0.76,ecut=700 -> -13.9156691899