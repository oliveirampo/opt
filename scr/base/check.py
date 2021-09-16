import sys


from scr.base import optimize


def prm_ratio(molecules, atomTypes):
    """Compute parameter ratio.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
    """

    prm_ratio_EE(molecules, atomTypes)
    prm_ratio_LJ(molecules, atomTypes)


def prm_ratio_LJ(molecules, atomTypes):
    """Compute parameter ratio for LJ types.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
    """

    types_map = {}
    unique_iac_count = {}
    prmsToOptimize = optimize.getParametersToBeOptimized(atomTypes)
    for prm_list in prmsToOptimize:
        if prm_list[0].typ == 'eps':
            unique_iac_count[prm_list[0].iac] = []
            for prm in prm_list:
                types_map[prm.iac] = prm_list[0].iac

    iac_count = count(types_map, molecules)
    for iac in unique_iac_count:
        unique_iac_count[iac] = iac_count[iac]

    print('#LJ types')
    print_prm_ratio(unique_iac_count, atomTypes)


def prm_ratio_EE(molecules, atomTypes):
    """Compute parameter ratio for EE types.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
    """

    types_map = {iac: iac for iac in atomTypes}
    iac_count = count(types_map, molecules)

    print('#EE types')
    print_prm_ratio(iac_count, atomTypes)


def count(types_map, molecules):
    """Count number of molecules, dns and hvp values for each type.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
    """

    iac_count = {key: [0, 0, 0] for key in types_map}
    used_mol = {}

    for cod in molecules:
        mol = molecules[cod]
        run = mol.run

        if not run:
            continue

        atoms = mol.atoms
        properties = mol.properties

        used_iac = []
        for idx in atoms:
            iac = atoms[idx].iac

            if iac not in types_map:
                continue
            iac = types_map[iac]

            if iac in used_iac:
                continue
            used_iac.append(iac)

            if iac not in used_mol:
                used_mol[iac] = []

            if cod[:-1] not in used_mol[iac]:
                used_mol[iac].append(cod[:-1])
                iac_count[iac][0] = iac_count[iac][0] + 1

            if len(properties) != 2:
                sys.exit('TODO')

            dns_wei = properties[0].wei
            hvp_wei = properties[1].wei

            if dns_wei:
                iac_count[iac][1] = iac_count[iac][1] + 1
            if hvp_wei:
                iac_count[iac][2] = iac_count[iac][2] + 1

    return iac_count


def print_prm_ratio(iac_mol, atomTypes):
    for iac in iac_mol:
        # print(iac, atomTypes[iac].nam)
        nam = atomTypes[iac].nam
        nMol = iac_mol[iac][0]
        nDns = iac_mol[iac][1]
        nHvp = iac_mol[iac][2]
        print('\t{:4} {:10} {:3} {:3} {:3}'.format(iac, nam, nMol, nDns, nHvp))
