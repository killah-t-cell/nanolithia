import pathlib
from os.path import exists
from ase.io import read, write
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3

def get_LiO2(xc, kpts=(8,8,8), Ecut=500):
    LiO2 = read(pathlib.Path(__file__).parent / 'LiO2.poscar')

    if xc == 'DFTD3':
        dft = GPAW(mode=PW(Ecut),
                   kpts=kpts,
                   txt=pathlib.Path(__file__).parent / f'LiO2-{xc}.log',
                   xc='PBE')
        calc = DFTD3(dft=dft, xc='PBE')
    else:
        calc = GPAW(mode=PW(Ecut),
                    kpts=kpts,
                    txt=pathlib.Path(__file__).parent / f'LiO2-{xc}.log',
                    xc=xc)

    LiO2.calc = calc
    write(pathlib.Path(__file__).parent / f'LiO2-{xc}.traj', LiO2)
    return LiO2
