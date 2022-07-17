import pathlib
from ase.io import read, write
from gpaw import GPAW, PW
from ase.calculators.dftd3 import DFTD3


def get_Li2O2(xc, kpts=(8, 8, 8), Ecut=500):
    Li2O2 = read(pathlib.Path(__file__).parent / 'Li2O2.poscar')

    if xc == 'DFTD3':
        dft = GPAW(mode=PW(Ecut),
                   kpts=kpts,
                   txt=pathlib.Path(__file__).parent / f'Li2O2-{xc}.log',
                   xc='PBE')
        calc = DFTD3(dft=dft, xc='PBE')
    else:
        calc = GPAW(mode=PW(Ecut),
                    kpts=kpts,
                    txt=pathlib.Path(__file__).parent / f'Li2O2-{xc}.log',
                    xc=xc)

    Li2O2.calc = calc
    write(pathlib.Path(__file__).parent / f'Li2O2-{xc}.traj', Li2O2)
    return Li2O2

# xc = 'PBE'
# Li2O2 = read(pathlib.Path(__file__).parent / 'Li2O2.poscar')
# if True:
#     calc = GPAW(mode=PW(500),
#                 kpts=(8, 8, 8),
#                 nbands=-10,
#                 txt=pathlib.Path(__file__).parent / f'Li2O2-{xc}.log',
#                 xc=xc)
# Li2O2.calc = calc
# print(Li2O2.get_potential_energy())
