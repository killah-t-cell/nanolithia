from compounds.Li.Li import *
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from compounds.LiO2.LiO2 import get_LiO2
from convergence import get_ecut_convergence, get_kpts_convergence

xc = 'PBE'
tol = 1e-4

# Li (-1.9 eV / atom)
Li_ecut = get_ecut_convergence(get_Li, xc, tol=tol)
Li_kpts = get_kpts_convergence(get_Li, xc, tol=tol)
Li = get_Li(xc, Ecut=Li_ecut, kpts=(Li_kpts, Li_kpts, Li_kpts))
print(Li.get_potential_energy() / len(Li)) # -1.9166736223308312

# Li2O (-4.771 eV / atom)
Li2O_ecut = get_ecut_convergence(get_Li2O, xc, tol=tol)
Li2O_kpts = get_kpts_convergence(get_Li2O, xc, tol=tol)
Li2O = get_Li2O('PBE', Ecut=Li2O_ecut, kpts=(Li2O_kpts, Li2O_kpts, Li2O_kpts))
print(Li2O.get_potential_energy() / len(Li2O))

# These correction terms thus obtained are 1.05 eV/O2, 0.76 eV/O22−, and 0.33 eV/O2− for oxides, peroxides,
# and superoxides, respectively, in HSE. With these new corrections and experimentally reported entropies,
# our calculations predict the voltages of Li2O and Li2O2 to be 2.93 and 2.97 V, respectively, which are in excellent
# agreement with experimental voltages of 2.91 and 2.96 V for Li2O and Li2O2, respectively.48
# 1. get approximate value based on known voltage
# 2. turn getLi2O2 into GGA+U and add the right U value
# 3. Compare both results
# Li2O2 (needs U correction)
Li2O2_ecut = get_ecut_convergence(get_Li2O2, xc, tol=tol)
Li2O2_kpts = get_kpts_convergence(get_Li2O2, xc, tol=tol)
Li2O2 = get_Li2O2(xc, Ecut=Li2O2_ecut, kpts=(Li2O2_kpts, Li2O2_kpts, Li2O2_kpts))
print(Li2O2.get_potential_energy() / len(Li2O2))

# LiO2 (needs U correction)
LiO2_ecut = get_ecut_convergence(get_LiO2, xc, tol=tol)
LiO2_kpts = get_kpts_convergence(get_LiO2, xc, tol=tol)
LiO2 = get_LiO2(xc, Ecut=LiO2_ecut, kpts=(LiO2_kpts, LiO2_kpts, LiO2_kpts))
print(LiO2.get_potential_energy() / len(LiO2))
