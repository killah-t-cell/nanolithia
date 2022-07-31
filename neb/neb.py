from ase.constraints import FixAtoms
import ase.db
from ase.neb import NEB
from ase.optimize import BFGS
from ase.visualize import view
from gpaw import GPAW, PW

from compounds.compounds import get_Li2O

db = ase.db.connect('compounds.db')
xc = 'PBE'

initial = get_Li2O(db, xc).toatoms()
final = initial.copy()

# change position of one lithium
cell = final.get_cell()
pos = final.get_positions()
pos[6] = pos[6] + cell[1] / 3 + cell[0] / 3
# pos[7] = pos[7] + cell[1] / 3 - cell[0] / 3
final.set_positions(pos)

view(final)

# make 7 bands
images = [initial]
images += [initial.copy() for i in range(5)]  # These will become the minimum energy path images.
images += [final]

# make an initial guess for images
neb = NEB(images)
neb.interpolate()

# view images
view(images)

# calculate all energy for all images
for image in images[0:7]:
    calc = GPAW(mode=PW(500), kpts=(5, 5, 6), xc='LDA', symmetry={'point_group': False})
    image.calc = calc
    image.set_constraint(FixAtoms(mask=[atom.symbol == 'O' for atom in image]))

# get potential energies of first and last image
images[0].get_potential_energy()
images[0].get_forces()
images[6].get_potential_energy()
images[6].get_forces()

optimizer = BFGS(neb, trajectory='Li2O_neb.traj', logfile='Li2O_neb.log')
optimizer.run(fmax=0.10)

# get energy barrier
IS_image = images[0]
TS_image = images[3]

epot_IS = IS_image.get_potential_energy()
epot_TS = TS_image.get_potential_energy()

barrier = epot_TS - epot_IS
print('Energy barrier:', barrier) # Energy barrier: 0.8683292062632475


