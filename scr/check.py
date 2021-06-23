import sys


def prm_ratio(molecules, atomTypes):
    """Compute parameter ratio.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
    """

    iac_mol = {}
    used_mol = {}

    for iac in atomTypes:
        iac_mol[iac] = [0, 0, 0]

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
            if iac in used_iac:
                continue
            used_iac.append(iac)

            if iac not in used_mol:
                used_mol[iac] = []

            if cod[:-1] not in used_mol[iac]:
                used_mol[iac].append(cod[:-1])
                iac_mol[iac][0] = iac_mol[iac][0] + 1

            if len(properties) != 2:
                sys.exit('TODO')

            dns_wei = properties[0].wei
            hvp_wei = properties[1].wei

            if dns_wei:
                iac_mol[iac][1] = iac_mol[iac][1] + 1
            if hvp_wei:
                iac_mol[iac][2] = iac_mol[iac][2] + 1

    for iac in iac_mol:
        # print(iac, atomTypes[iac].nam)
        nam = atomTypes[iac].nam
        nMol = iac_mol[iac][0]
        nDns = iac_mol[iac][1]
        nHvp = iac_mol[iac][2]
        print('\t{:4} {:10} {:3} {:3} {:3}'.format(iac, nam, nMol, nDns, nHvp))

    # get pamrater ratio for LJ parameters
    iac_LJ = {13: 13, 14: 14, 15: 15, 16: 16, 17:17,
              45: 13, 39: 14, 69: 15, 68: 16, 32: 32, 33: 33, 34: 34, 67: 67,
              127: 127, 103: 13, 104: 14, 105: 15, 106: 16,
              128: 128, 123: 123, 101: 101,
              125: 123, 100: 101,
              129: 128, 124: 123, 126: 127, 102: 101, 115: 13, 116: 14, 117: 15, 118: 16,
              21: 21, 3: 3, 111: 13, 112: 14, 113: 15, 114: 16,
              20: 128, 1: 123, 12: 101, 97: 21, 98: 3,
              145: 21, 7: 7, 144: 13, 142: 14, 140: 15, 138: 16,
              4: 123, 18: 101, 150: 21, 6: 6, 146: 13, 147: 14, 148: 15, 149: 16}

    # TODO
