from ase.build import bulk

from matplotlib import pyplot as plt
import matplotlib as mpl

# TODO move this to Compound() at some point
def converge(db, xc, *atoms, tol=1e-4, k_range=range(4, 12), ecut_range=range(350, 901, 50)):
    """
    :param db:
    :param xc:
    :param tol:
    :param k_range:
    :param ecut_range:
    :param atoms: list of get_compound functions
    :return:
    """
    for atom in atoms:
        formula = atom.__name__.replace('get_', '')
        ecut_convergence = {}
        nkpts_convergence = {}
        for ecut in ecut_range:
            at = atom(db, xc, ecut=ecut)
            e = at.toatoms().get_potential_energy()
            ecut_convergence[ecut] = e
        for k in k_range:
            at = atom(db, xc, nkpts=k)
            e = at.toatoms().get_potential_energy()
            nkpts_convergence[k] = e

        # process logs
        ecut_lists = sorted(ecut_convergence.items())  # sorted by key, return a list of tuples
        nkpts_lists = sorted(nkpts_convergence.items())  # sorted by key, return a list of tuples
        ecut_list, ecut_e_list = zip(*ecut_lists)  # unpack a list of pairs into two tuples
        nkpts_list, nkpts_e_list = zip(*nkpts_lists)  # unpack a list of pairs into two tuples
        ecut_converged = next(
            (i for i, j, _i, _j in zip(ecut_list, ecut_list[1:], ecut_e_list, ecut_e_list[1:]) if _i - _j <= abs(tol)),
            ecut_list[-1])
        nkpts_converged = next(
            (i for i, j, _i, _j in zip(nkpts_list, nkpts_list[1:], nkpts_e_list, nkpts_e_list[1:]) if
             _i - _j <= abs(tol)),
            nkpts_list[-1])

        # TODO if tol is not achieved, warn and save in db as not converged to tol, but update tol

        # write converged result into db
        atom(db, xc, nkpts=nkpts_converged, ecut=ecut_converged, converged=True, tol=tol)

        # plot logs
        mpl.rcParams['figure.dpi'] = 300  # increase plot dpi
        fig, (ax_ecut, ax_nkpts) = plt.subplots(2)
        fig.suptitle(f'{formula} {xc} convergence study')

        ax_ecut.ticklabel_format(useOffset=False)
        ax_ecut.plot(ecut_list, ecut_e_list)
        ax_ecut.plot(ecut_converged, ecut_convergence[ecut_converged], 'o')
        ax_ecut.set(xlabel='Ecut (eV)', ylabel='Potential energy (eV)')

        ax_nkpts.plot(nkpts_list, nkpts_e_list)
        ax_nkpts.plot(nkpts_converged, nkpts_convergence[nkpts_converged], 'o')
        ax_nkpts.set(xlabel='k-points', ylabel='Potential energy (eV)')

        plt.tight_layout()
        plt.savefig(f'plots/{formula}-{xc}-convergence.png')
        plt.show()

        print(formula, 'ecut conv=',ecut_converged, 'nkpts conv=',nkpts_converged)
