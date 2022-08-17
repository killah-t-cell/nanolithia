import numpy as np
from ase.build import stack
from ase.phasediagram import PhaseDiagram
from matplotlib import pyplot as plt

from compounds.compounds import Compound
from properties.profiles import LIO2_ENTROPY, LI2O2_ENTROPY, LI2O_ENTROPY, gibbs_energy, formation_energy
from properties.voltages import get_eq_voltage


def get_profiles(xc):
    # create supercell
    Li2O_conventional_standard = Compound('Li2O', 'mp-1960_conventional_standard', xc,
                                          magmoms=[0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6],
                                          converged=True, setups={'O': ':p, -1.05,0'},
                                          extension='.cif')
    Li16O8 = stack(Li2O_conventional_standard.atoms, Li2O_conventional_standard.atoms, reorder=True)

    # set up corrections
    li_states = [4, 8, 16]
    U_corrections = [-0.33, -0.75, -1.05]
    entropies = [LIO2_ENTROPY, LI2O2_ENTROPY, LI2O_ENTROPY]
    intermediates = [{'name': 'Li16O8', 'num_li': 16, 'gibbs': -6.156413617038097, 'formation': -6.189181117038097,
                      'compound': Li2O_conventional_standard}]

    # get intermediates, their energies, and their formation and gibbs energies.
    for num_li in range(15, 3, -1):
        # linearly interpolate U and entropy corrections
        inter_U_correction = np.interp(num_li, li_states, U_corrections)
        inter_entropy = np.interp(num_li, li_states, entropies)

        # set up intermediate cell
        magmoms = [0.6] * (num_li + 8)
        name = f'Li{num_li}O8'
        del Li16O8[0]
        Li16O8.write(f'structures/{name}_intermediate.cif')
        intermediate = Compound(name, 'intermediate', xc, magmoms=magmoms,
                                converged=True, setups={'O': f':p, {inter_U_correction},0'},
                                extension='.cif')

        # get intermediate energies
        epot_intermediate_cell = intermediate.get_energy()  # TODO why is it not converging?
        intermediate_gibbs = gibbs_energy(epot_intermediate_cell / 8, inter_entropy, num_li / 8, 1)
        intermediate_ef = formation_energy(epot_intermediate_cell / 8, num_li / 8, 1)
        intermediates.append({'name': name, 'num_li': num_li, 'gibbs': intermediate_gibbs,
                              'formation': intermediate_ef, 'compound': intermediate})

    refs = [('Li', 0.0), ('O', 0.0)]
    voltages = []
    concentrations = [inter['num_li'] for inter in intermediates]
    print(concentrations)
    volumes = [inter['compound'].volume for inter in intermediates]

    for inter, next_inter in zip(intermediates, intermediates[1:] + [intermediates[0]]):
        refs.append((inter['name'], inter['formation']))
        voltages.append(get_eq_voltage(inter['gibbs'], next_inter['gibbs'], 0.125))  # not sure if 0.125 or 1

    # plot convex hull
    pd = PhaseDiagram(refs)
    pd.plot()
    plt.savefig(f'plots/convex-hull-0K.png')
    pd.plot(show=True)

    # plot voltage profile
    plt.plot(concentrations, voltages)
    plt.ylabel('Voltage (V)')
    plt.xlabel('Li')
    plt.savefig(f'plots/voltage-profile.png')
    plt.show()

    # plot volume expansion
    plt.plot(concentrations, volumes)
    plt.ylabel('Volume')
    plt.xlabel('Li')
    plt.savefig(f'plots/volume-expansion-profile.png')
    plt.show()
