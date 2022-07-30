from ase.build import molecule
from gpaw import GPAW

a = 8.0
h = 0.2

energies = {}
with open(f'results-{h:.2f}.txt', 'w') as resultfile:
    for name in ['O2', 'O']:
        system = molecule(name)
        system.set_cell((a, a, a))
        system.center()

        calc = GPAW(h=h,
                    txt=f'gpaw-{name}-{h:.2f}.txt')
        if name == 'O':
            calc.set(hund=True)

        system.calc = calc

        energy = system.get_potential_energy()
        energies[name] = energy
        print(name, energy, file=resultfile)

    e_atomization = energies['O2'] - 2 * energies['O']
    print(e_atomization, file=resultfile)



    # parameters = dict(mode=PW(ecut),
    #                   kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
    #                   setups=U_correction,
    #                   xc=xc)



# get_O2() right now gives = -8.22145087379
# which has a binding energy of -8.22145087379 - (2*-1.8302341313767456) =-4.56098261104
# the binding energy we want is -4.95
# so we are underestimating O2. We need O2 to give energy -8.61