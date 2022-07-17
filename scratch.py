from compounds.Li.Li import *
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from convergence import get_ecut_convergence

epot_Li = -1.9169006433926967
epot_Li2O2 = -18.937235272212817
epot_Li2O = -14.01310969105073
epot_LiO2 = -12.450224352472798

v1 = - (2*epot_Li2O - epot_Li2O2 - 2*epot_Li)/2

Li = get_Li('PBE')
print(Li.get_potential_energy()/2)

Li2O = get_Li2O('PBE', ecut=600, nkpts=(16, 16, 16))
print(Li2O.get_potential_energy()/12)


Li2O2_ecut = get_ecut_convergence(get_Li2O2, 'PBE')
Li2O2_ecut
print(Li2O2_ecut)


# V_eq_w_Li2O2 = (2 * epot_Li2O - epot_Li2O2 - 2 * epot_Li) / 2

'''
if I know V_eq_w_Li2O2 = 2.87 V, then I can get the (approximate) solution for Li2O2 and LiO2
V_eq_w_Li2O2 = 2*epot_Li2O/2 - epot_Li2O2/2 - 2* epot_Li/2

V_eq_w_Li2O2 + epot_Li2O2/2 = 2*epot_Li2O/2 - 2* epot_Li/2
epot_Li2O2/2 = 2*epot_Li2O/2 - 2* epot_Li/2 - V_eq_w_Li2O2
epot_Li2O2 = 2* epot_Li2O - 2* epot_Li - 2*V_eq_w_Li2O2
'''

# -61.153762361981805/12

# # This script will optimize lattice constant of metallic lithium
# for xc in ['LDA', 'PBE']:
#     Li_metal = bulk('Li', crystalstructure='bcc', a=3.51)
#
#     calc = GPAW(mode=PW(900),
#                 kpts=(16, 16, 16),
#                 txt=f'Li-metal-{xc}.log',
#                 xc=xc)
#
#     Li_metal.calc = calc
#
#     epot = Li_metal.get_potential_energy()
#     print(epot)
#
#
#
# for xc in ['LDA', 'PBE']:
#     Li_metal = bulk('Li', crystalstructure='bcc', a=3.51)
#
#     if xc == 'DFTD3':
#         dft = GPAW(mode=PW(500),
#                    kpts=(8, 8, 8),
#                    nbands=-10,
#                    txt=f'Li-metal-{xc}.log',
#                    xc='PBE')
#         calc = DFTD3(dft=dft, xc='PBE')
#     else:
#         calc = GPAW(mode=PW(500),
#                     kpts=(8, 8, 8),
#                     nbands=-10,
#                     txt=f'Li-metal-{xc}.log',
#                     xc=xc)
#
#     Li_metal.calc = calc
#
#     print(len(Li_metal))
#     e = Li_metal.get_potential_energy()
#     print(e)
