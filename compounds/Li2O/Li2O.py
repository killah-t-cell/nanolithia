import pathlib
from ase.io import read, write
from gpaw import GPAW, FermiDirac, PW
from ase.calculators.dftd3 import DFTD3


def get_Li2O(xc, kpts=(8, 8, 8), Ecut=500, nbands=20):
    Li2O = read(pathlib.Path(__file__).parent / 'Li2O.poscar')

    if xc == 'DFTD3':
        dft = GPAW(mode=PW(Ecut),
                   kpts=kpts,
                   nbands=nbands,
                   txt=pathlib.Path(__file__).parent / f'Li2O-{xc}.log',
                   xc='PBE')
        calc = DFTD3(dft=dft, xc='PBE')
    else:
        calc = GPAW(mode=PW(Ecut),
                    nbands=nbands,
                    kpts=kpts,
                    txt=pathlib.Path(__file__).parent / f'Li2O-{xc}.log',
                    xc=xc)

    Li2O.calc = calc
    write(pathlib.Path(__file__).parent / f'Li2O-{xc}.traj', Li2O)
    return Li2O
