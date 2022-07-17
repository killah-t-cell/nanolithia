import pathlib
from os.path import exists
from ase import Atoms
from ase.io import read, write
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, PW
from ase.optimize import QuasiNewton
from ase.build import bulk
from ase.calculators.dftd3 import DFTD3
from ase.constraints import StrainFilter
from ase.optimize.bfgs import BFGS
from gpaw.xc.bee import BEEFEnsemble

def get_Li(xc, kpts=(8,8,8), Ecut=500, nbands=-10):
    if exists(pathlib.Path(__file__).parent / f'Li-{xc}.traj'):
        Li = read(pathlib.Path(__file__).parent / f'Li-{xc}.traj')
    else:
        Li = bulk('Li', crystalstructure='bcc', a=3.51, cubic=True)

    if xc == 'DFTD3':
        dft = GPAW(mode=PW(Ecut),
                   kpts=kpts,
                   nbands=nbands,
                   txt=pathlib.Path(__file__).parent / f'Li-{xc}.log',
                   xc='PBE')
        calc = DFTD3(dft=dft, xc='PBE')
    else:
        calc = GPAW(mode=PW(Ecut),
                    kpts=kpts,
                    nbands=nbands,
                    txt=pathlib.Path(__file__).parent / f'Li-{xc}.log',
                    xc=xc)

    Li.calc = calc
    sf = StrainFilter(Li, mask=[1, 1, 1, 0, 0, 0])
    opt = BFGS(sf)
    opt.run(fmax=0.01)
    write(pathlib.Path(__file__).parent / f'Li-{xc}.traj', Li)

    return Li


def get_emsemble(Li):
    ens = BEEFEnsemble(Li.calc)
    Li_ens_cell = ens.get_ensemble_energies(2000)
    with paropen(pathlib.Path(__file__).parent / f'ensemble_Li.dat', 'a') as result:
        for e in Li_ens_cell:
            print(e, file=result)


