# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import ase.db

from compounds.Li.Li import get_Li
from compounds.Li2O.Li2O import get_Li2O
from compounds.Li2O2.Li2O2 import get_Li2O2
from compounds.LiO2.LiO2 import get_LiO2

from matplotlib import pyplot as plt
import matplotlib as mpl


def converge(db, xc, *atom_rows, tol=1e-4, k_range=range(4, 18), ecut_range=range(350, 801, 50)):
    """
    :param db:
    :param xc:
    :param tol:
    :param k_range:
    :param ecut_range:
    :param atom_rows:
    :return:
    """
    for row, callback in atom_rows:
        id = db.reserve(name=row.name, converged=True)
        if id is not None:  # skip convergence if already done
            ecut_convergence = {}
            nkpts_convergence = {}
            for k in k_range:
                atom = callback(db, xc, nkpts=k)
                e = atom.get_potential_energy()
                nkpts_convergence[k] = e
                # store converged tol
            for ecut in ecut_range:
                atom = callback(db, xc, ecut=ecut)
                e = atom.get_potential_energy()
                ecut_convergence[ecut] = e

            # process results
            ecut_lists = sorted(ecut_convergence.items())  # sorted by key, return a list of tuples
            nkpts_lists = sorted(nkpts_convergence.items())  # sorted by key, return a list of tuples
            ecut_list, ecut_e_list = zip(*ecut_lists)  # unpack a list of pairs into two tuples
            nkpts_list, nkpts_e_list = zip(*nkpts_lists)  # unpack a list of pairs into two tuples

            # save results
            ecut_converged = next(
                (i for i in range(1, len(ecut_lists)) if ecut_lists[i][1] - ecut_lists[i + 1][1] <= tol),
                ecut_lists[-1])
            nkpts_converged = next(
                (i for i in range(1, len(nkpts_lists)) if nkpts_lists[i][1] - nkpts_lists[i + 1][1] <= tol),
                nkpts_lists[-1])

            # TODO if tol is not achieved, warn and save in db as not converged to tol, but update tol

            # update db
            db.update(row.id,
                      nkpts=nkpts_converged,
                      ecut=ecut_converged,
                      converged=True,
                      convergence_tol=tol,
                      atoms=callback(db, xc, nkpts=nkpts_converged, ecut=ecut_converged))

            # plot results
            mpl.rcParams['figure.dpi'] = 300  # increase plot dpi
            fig, (ax_ecut, ax_nkpts) = plt.subplots(2, 1)
            ax_ecut.plot(ecut_list, ecut_e_list, 'r')
            ax_ecut.set(xlabel='Ecut (eV)', ylabel='Potential energy (eV)')
            ax_ecut.set_title(f'{row.formula}-{xc}-ecut')

            ax_nkpts.plot(nkpts_list, nkpts_e_list, 'r')
            ax_nkpts.set(xlabel='k-points', ylabel='Potential energy (eV)')
            ax_nkpts.set_title(f'{row.formula}-{xc}-nkpts')

            plt.show()
            ax_ecut.savefig(f'plots/{row.formula}-{row.xc}-ecut-convergence.png')
            ax_nkpts.savefig(f'plots/{row.formula}-{row.xc}-kpts-convergence.png')
        else:
            print(row, "is already converged")


if __name__ == '__main__':
    db = ase.db.connect('compounds.db')
    xc = 'PBE'
    Li = get_Li(db, xc)
    Li2O = get_Li2O(db, xc)
    Li2O2 = get_Li2O2(db, xc)
    LiO2 = get_LiO2(db, xc)

    # del db[db.get(name=f'Li-{xc}-8x8x8-500', converged=True).id]
    compounds_to_converge = ((db.get(name=f'Li-{xc}-8x8x8-500'), get_Li),
                             (db.get(name=f'Li2O-{xc}-8x8x8-500'), get_Li2O),
                             (db.get(name=f'LiO2-{xc}-8x8x8-500'), get_LiO2),
                             (db.get(name=f'Li2O2-{xc}-8x8x8-500'), get_Li2O2))

    converge(db, xc, *compounds_to_converge)

    # load voltages
    # if not there
    # run voltage experiments
