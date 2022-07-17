import os
import pathlib
from os.path import exists

from ase.calculators.dftd3 import DFTD3
from ase.io import read, write
from gpaw import GPAW, PW


def set_atoms(xc, file, elementName):
    if exists(pathlib.Path(__file__).parent / f'{elementName}-{xc}.traj'):
        atoms = read(pathlib.Path(__file__).parent / f'{elementName}-{xc}.traj')
    else:
        atoms = read(pathlib.Path(__file__).parent / file)
        if xc == 'DFTD3':
            dft = GPAW(mode=PW(500),
                       kpts=(8, 8, 8),
                       nbands=-10,
                       txt=pathlib.Path(__file__).parent / f'{elementName}-{xc}.log',
                       xc='PBE')
            calc = DFTD3(dft=dft, xc='PBE')
        else:
            calc = GPAW(mode=PW(500),
                        kpts=(8, 8, 8),
                        nbands=-10,
                        txt=pathlib.Path(__file__).parent / f'{elementName}-{xc}.log',
                        xc=xc)

        atoms.calc = calc
        write(pathlib.Path(__file__).parent / f'{elementName}-{xc}.traj', atoms)
    return atoms