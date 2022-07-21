
# Li, O
# what are the ways to combine those fuckers
# cheating a bit by  stealing from docs, seeing if it works
from pathlib import Path
import random

from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.data import atomic_numbers, reference_states
from ase.ga.data import PrepareDB
from ase.ga import set_raw_score
from gpaw import GPAW, PW

parameters = dict(mode=PW(500),
                      kpts={'size': (8, 8, 8), 'gamma': True},
                      xc='LDA')

def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)

cathode_elements = ['Li', 'O']
# Use experimental lattice constants
lattice_constants = dict((m, reference_states[atomic_numbers[m]]['a' if m == 'Li' else 'd'])
                         for m in cathode_elements)

# Create the references (pure slabs) manually
pure_slabs = []
refs = {}
print('Reference energies:')
for m in cathode_elements:
    slab = fcc111(m, size=(2, 4, 3), a=lattice_constants[m],
                  vacuum=5, orthogonal=True, periodic=True)
    slab.calc = GPAW(**parameters)

    # We save the reference energy as E_A / N
    e = slab.get_potential_energy()
    e_per_atom = e / len(slab)
    refs[m] = e_per_atom
    print('{0} = {1:.3f} eV/atom'.format(m, e_per_atom))

    # The mixing energy for the pure slab is 0 by definition
    set_raw_score(slab, 0.0)
    pure_slabs.append(slab)
    
# The population size should be at least the number of different compositions
pop_size = 2 * len(slab)

# We prepare the db and write a few constants that we are going to use later
target = Path('genetic_hull.db')
if target.exists():
    target.unlink()
db = PrepareDB(target, population_size=pop_size,
               reference_energies=refs, metals=cathode_elements,
               lattice_constants=lattice_constants)

# We add the pure slabs to the database as relaxed because we have already
# set the raw_score
for slab in pure_slabs:
    db.add_relaxed_candidate(slab, atoms_string=''.join(slab.get_chemical_symbols()))

# Now we create the rest of the candidates for the initial population
for i in range(pop_size - 2):
    # How many of each metal is picked at random, making sure that
    # we do not pick pure slabs
    nA = random.randint(0, len(slab) - 2)
    nB = len(slab) - 2 - nA
    symbols = [cathode_elements[0]] * nA + [cathode_elements[1]] * nB + cathode_elements

    # Making a generic slab with the correct lattice constant
    slab = fcc111('X', size=(2, 4, 3),
                  a=get_avg_lattice_constant(symbols),
                  vacuum=5, orthogonal=True, periodic=True)

    # Setting the symbols and randomizing the order
    slab.set_chemical_symbols(symbols)
    random.shuffle(slab.numbers)

    # Add these candidates as unrelaxed, we will relax them later
    atoms_string = ''.join(slab.get_chemical_symbols())
    db.add_unrelaxed_candidate(slab, atoms_string=atoms_string)
    
import numpy as np
from ase.ga.population import RankFitnessPopulation
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.slab_operators import (CutSpliceSlabCrossover,
                                   RandomSlabPermutation,
                                   RandomCompositionMutation)
from ase.ga import set_raw_score

from ase.calculators.emt import EMT

# Connect to the database containing all candidates
db = DataConnection('genetic_hull.db')

# Retrieve saved parameters
pop_size = db.get_param('population_size')
refs = db.get_param('reference_energies')
metals = db.get_param('metals')
lattice_constants = db.get_param('lattice_constants')

def get_mixing_energy(atoms):
    # Set the correct cell size from the lattice constant
    new_a = get_avg_lattice_constant(atoms.get_chemical_symbols())
    # Use the orthogonal fcc cell to find the current lattice constant
    current_a = atoms.cell[0][0] / np.sqrt(2)
    atoms.set_cell(atoms.cell * new_a / current_a, scale_atoms=True)

    # Calculate the energy
    atoms.calc = GPAW(**parameters)
    e = atoms.get_potential_energy()

    # Subtract contributions from the pure element references
    # to get the mixing energy
    syms = atoms.get_chemical_symbols()
    for m in set(syms):
        e -= syms.count(m) * refs[m]
    return e


def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)


def get_comp(atoms):
    return atoms.get_chemical_formula()

# Specify the number of generations this script will run
num_gens = 10

# Specify the procreation operators for the algorithm and
# how often each is picked on average
# The probability for an operator is the prepended integer divided by the sum
# of integers
oclist = [(3, CutSpliceSlabCrossover()),
          (1, RandomSlabPermutation()),
          (1, RandomCompositionMutation())
          ]
operation_selector = OperationSelector(*zip(*oclist))

# Pass parameters to the population instance
# A variable_function is required to divide candidates into groups here we use
# the chemical composition
pop = RankFitnessPopulation(data_connection=db,
                            population_size=pop_size,
                            variable_function=get_comp)

# Evaluate the starting population
# The only requirement of the evaluation is to set the raw_score
# Negative mixing energy means more stable than the pure slabs
# The optimization always progress towards larger raw score,
# so we take the negative mixing energy as the raw score
print('Evaluating initial candidates')
while db.get_number_of_unrelaxed_candidates() > 0:
    a = db.get_an_unrelaxed_candidate()
    set_raw_score(a, -get_mixing_energy(a))
    db.add_relaxed_step(a)
pop.update()

# Below is the iterative part of the algorithm
gen_num = db.get_generation_number()
for i in range(num_gens):
    print('Creating and evaluating generation {0}'.format(gen_num + i))
    new_generation = []
    for _ in range(pop_size):
        # Select parents for a new candidate
        parents = pop.get_two_candidates()

        # Select an operator and use it
        op = operation_selector.get_operator()
        offspring, desc = op.get_new_individual(parents)
        # An operator could return None if an offspring cannot be formed
        # by the chosen parents
        if offspring is None:
            continue

        set_raw_score(offspring, -get_mixing_energy(offspring))
        new_generation.append(offspring)

    # We add a full relaxed generation at once, this is faster than adding
    # one at a time
    db.add_more_relaxed_candidates(new_generation)

    # update the population to allow new candidates to enter
    pop.update()

import numpy as np
import matplotlib.pyplot as plt
from ase.phasediagram import PhaseDiagram
from ase.db import connect
from ase.io import write

db = connect('genetic_hull.db')

# Select the evaluated candidates and retrieve the chemical formula and mixing
# energy for the phase diagram
refs = []
dcts = list(db.select('relaxed=1'))
for dct in dcts:
    refs.append((dct.formula, -dct.raw_score))

pd = PhaseDiagram(refs)
ax = pd.plot(show=not True,  # set to True to show plot
             only_label_simplices=True)
plt.savefig('hull.png', dpi=300)

# View the simplices of the convex hull
simplices = []
toview = sorted(np.array(dcts)[pd.hull], key=lambda x: x.mass)
for dct in toview:
    simplices.append(dct.toatoms())

write('hull.traj', simplices)
