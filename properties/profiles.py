from compounds.compounds import Compound

# Entropies
O2_ENTROPY = 2.13
LI_ENTROPY = 0.30
LI2O2_ENTROPY = 0.59
LI2O_ENTROPY = 0.39
LIO2_ENTROPY = 1.21
O2_ENERGY = Compound('O2', 'mp-12957', 'PBE', magmoms=[0.6, 0.6, 0.6, 0.6], ecut=530).get_energy() / 2  # or -9.896
LI_ENERGY = Compound('Li', 'mp-1', 'PBE', ecut=530).get_energy() / 4  # or 1.90


def formation_energy(unit_formula_energy, a, b):
    return unit_formula_energy - a * LI_ENERGY - b * O2_ENERGY / 2


def gibbs_energy(unit_formula_energy, atoms_entropy, a, b, T=0.0257):  # T is 298.15 K in eV
    return unit_formula_energy - T * atoms_entropy - a * (LI_ENERGY - T * LI_ENTROPY) - b / 2 * (
            O2_ENERGY - T * O2_ENTROPY)


def hybrid_energy():
    return
