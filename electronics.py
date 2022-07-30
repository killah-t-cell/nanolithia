import ase.db
import matplotlib.pyplot as plt
from gpaw import GPAW
from ase.dft.dos import DOS

from compounds.Li2O.Li2O import get_Li2O
from compounds.LiO2.LiO2 import get_LiO2

db = ase.db.connect('compounds.db')
xc = 'PBE'
calc = get_Li2O(db, xc)

dos = DOS(calc, npts=800, width=0)
energies = dos.get_energies()
weights = dos.get_dos()
#
# ax = plt.gca()
# ax.plot(energies, weights)
# ax.set_xlabel('Energy [eV]')
# ax.set_ylabel('DOS [1/eV]')
# plt.savefig('dos.png')
# plt.show(