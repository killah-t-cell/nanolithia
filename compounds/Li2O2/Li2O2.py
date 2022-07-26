import pathlib
from ase.io import read, write
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3


def get_Li2O2(db, xc, nkpts=8, ecut=575, converged=False, tol='null', structure='mp-841', spinpol=False):
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
    name = f'Li2O2-{structure}-{xc}-{nkpts}x{nkpts}x{nkpts}-{ecut:.0f}'
    U_correction = {'O': ':p,0.76,0'}

    parameters = dict(mode=PW(ecut),
                      kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
                      spinpol=spinpol,
                      convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
                                  'energy': 2.0e-4,  # eV / electron
                                  'density': 1.0e-3, },
                      setups=U_correction,
                      xc=xc)


    id = db.reserve(name=name, xc=xc, nkpts=nkpts, ecut=ecut)
    if id is not None:  # skip calculation if already done
        # load crystal
        Li2O2 = read(pathlib.Path(__file__).parent / f'Li2O2_{structure}.poscar')

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
                 structure=structure,
                 tol=tol)

    return db.get(name=name, xc=xc, nkpts=nkpts, ecut=ecut, structure=structure)

# goal -36.972
#550, 8, unnormalized, -36.7
#575, 8, unnormalized, -37.0
#600, 8, unnormalized, -37.18
#800, 8, unnormalized, -37.4


# LDA, with U correction
# -41.700487

# LDA without U correction
# -41.70544355635925

# LDA without U correction, unnormalized
# -41.70544355635925

# PBE with U correction {'O': ':p,0.76'} k=8, ecut=500
# -36.4

# PBE with U correction {'O': ':p,0.76'} k=8, ecut=600

# PBE without U correction, k=8, ecut=600