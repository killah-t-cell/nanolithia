# def construct_convex(db, xc, ecut=500, nkpts=8):
#     # create convex db
#     convex_db = ase.db.connect('hull.db')
#
#     # get Li8O4 cell
#     Li2O = get_Li2O(db, xc)
#
#     # set calculator parameters
#     parameters = dict(mode=PW(ecut),
#                       kpts={'size': (nkpts, nkpts, nkpts), 'gamma': True},
#                       spinpol=True,
#                       convergence={'eigenstates': 1.0e-4,  # eV^2 / electron
#                                   'energy': 2.0e-4,  # eV / electron
#                                   'density': 1.0e-3, },
#                       xc=xc)
#
#     # get all structures
#     for i in range(1, 9):
#         name = str('Li' + str(8 - i) + 'O4' if (8 - i > 0) else 'O4') + f'-{xc}-{ecut}-{nkpts}'
#
#         id = convex_db.reserve(name=name, xc=xc)
#         if id is not None:  # skip calculation if already done
#             # copy structure of Li2O and pop a Li
#             structure = Li2O.copy()
#             [structure.pop(0) for _ in range(0, i)]
#
#             # attach calculator
#             if xc == 'DFTD3':
#                 dft = GPAW(mode=PW(ecut),
#                            kpts=nkpts,
#                            txt=name + '.txt',
#                            xc='PBE')
#                 calc = DFTD3(dft=dft, xc='PBE')
#             else:
#                 calc = GPAW(txt=name + '.txt',
#                             **parameters)
#
#             structure.calc = calc
#
#             # get potential energy
#             e = structure.get_potential_energy()
#             print(e)
#
#             # save structure with name and potential energy
#             del convex_db[id]
#             convex_db.write(structure,
#                      name=name,
#                      num_Li=i,
#                      xc=xc)
#
#     # calculate formation energy for each structure
#     LiO4 = convex_db.get_atoms(name=f'Li1O4-{xc}-{ecut}-{nkpts}')
#     O4 = convex_db.get_atoms(name=f'O4-{xc}-{ecut}-{nkpts}')
#     e_LiO4 = LiO4.get_potential_energy()
#     e_O4 = O4.get_potential_energy()
#
#     for structure_row in convex_db.select():
#         structure = convex_db.get_atoms(name=structure_row.name, xc=xc)
#         se = structure.get_potential_energy()
#         ef = se - (structure_row.num_Li*e_LiO4) - (1 - structure_row.num_Li)*e_O4
#         convex_db.update(structure_row.id, formation_energy=ef)
#         print(ef)
#
#     # construct convex hull
#
#     # plot convex hull
