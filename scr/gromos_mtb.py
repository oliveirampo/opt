from collections import OrderedDict
import shutil
import sys
import os

import createNewIFP


def updateMTB(newIacDict, atomTypes, molecules, mtbDir, mtbGROMOSDir):
    """reates GROMOS (mtb) topology files.

    :param newIacDict: (dict) Dictionary that maps IAC idx to first idx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param mtbDir: (str) Directory with template mtb files.
    :param mtbGROMOSDir: (str) Directory where new mtb will be saved.
    :return:
    """

    copyMTB(mtbDir, mtbGROMOSDir)
    updateParameters(newIacDict, atomTypes, molecules, mtbGROMOSDir)


def copyMTB(mtbDir, mtbGROMOSDir):
    """Copies template of MTB file.

    :param mtbDir: (str) Directory with template MTB file.
    :param mtbGROMOSDir: (str) Directory where new mtb will be saved.
    """

    if not os.path.exists(mtbGROMOSDir):
        os.makedirs(mtbGROMOSDir)

    mtbSrc = '{}/mol.mtb'.format(mtbDir)
    mtbTrg = '{}/mol.mtb'.format(mtbGROMOSDir)

    if not os.path.exists(mtbSrc):
        sys.exit('\n\tNo such file: {}\n'.format(mtbSrc))

    shutil.copyfile(mtbSrc, mtbTrg)


def updateParameters(newIacDict, atomTypes, molecules, mtbGROMOSDir):
    """Updates values of atomic partial charges.

    :param newIacDict: (dict) Dictionary that maps IAC idx to first idx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    :param mtbGROMOSDir: (str) Directory where new mtb will be saved.
    """

    charges = OrderedDict()

    for cod in molecules:
        mol = molecules[cod]

        if not cod[-1].isdigit():
            cod = cod[:-1]

        if cod in charges:
            continue

        atoms = mol.atoms
        molCharges = []

        for idx in atoms:
            chg = atoms[idx].charge.cur
            molCharges.append(chg)

        if len(molCharges) != len(atoms):
            sys.exit('\nERROR: wrong number of atoms.')

        charges[cod] = molCharges

    mtbFile = '{}/mol.mtb'.format(mtbGROMOSDir, cod)
    rows, mtb = readMTBCHG(mtbFile)

    # update charges
    updateCharges(rows, mtb, charges)

    # update atom type indexes
    updateIAC(rows, mtb, newIacDict)

    writeMTB(mtb, mtbFile)


def readMTBCHG(fileName):
    """Reads 'MTBUILDBLSOLUTE' blocks from mtb file.

    :param fileName: (str) Name of mtb file.
    :return:
        rows: (dict) Dictionary with rows from 'MTBUILDBLSOLUTE' blocks for each molecule.
        metb: (list) Rows of mtb file.
    """

    with open(fileName, 'r') as f:
        mtb = f.readlines()

    rows = OrderedDict()
    for i in range(len(mtb)):
        if mtb[i].startswith('MTBUILDBLSOLUTE'):
            cod = mtb[i + 2].strip()
            rows[cod] = {}

            for j in range(i + 1, len(mtb)):
                if mtb[j].startswith('END'):
                    break

                if not mtb[j].startswith('#'):
                    if len(mtb[j].strip().split()) > 6:
                        # To check, try to convert 1st (atom names) to float.
                        try:
                            float(mtb[j].strip().split()[1])
                        except ValueError:
                            rows[cod][j] = mtb[j]

    return rows, mtb


def formatStrCHG(r):
    """Returns formatted string with charge.

    :param r: (list) Row.
    :return:
        s: (str) Formatted string.
    """

    charge = float(r[4])
    r[4] = '{:.3f}'.format(charge)

    nb = '     '.join(r[7:])
    s = '{:>5} {:>2} {:>6} {:>4} {:>7} {:>2} {:>3}     {}\n'.format(r[0], r[1], r[2], r[3], r[4], r[5], r[6], nb)
    return s


def updateCharges(rows, mtb, charges):
    """Helper function to update charges in mtb file.

    :param rows: (dict) Dictionary with rows from 'MTBUILDBLSOLUTE' block.
    :param mtb: (list) Rows of topology file.
    :param charges: (list) List with charges.
    """

    for cod in charges:
        selectedRows = rows[cod]
        selectedCharges = charges[cod]

        if len(selectedRows) != len(selectedCharges):
            print(selectedRows)
            print(selectedCharges)
            sys.exit('\nERROR: wrong number of atoms.')

        count = 0
        for r in selectedRows:
            line = selectedRows[r].strip().split()
            idx = line[0]
            idx = int(idx)

            count += 1
            if idx != count:
                sys.exit('\nERROR: idx')

            pos = idx - 1
            line[4] = selectedCharges[pos]
            s = formatStrCHG(line)

            selectedRows[r] = s
            mtb[r] = s


def updateIAC(rows, mtb, newIacDict):
    """Helper function to update atom type indexes in mtb file.

    :param rows: (dict) Dictionary with rows from 'MTBUILDBLSOLUTE' block.
    :param mtb: (list) Rows of topology file.
    :param newIacDict: (dict) Dictionary that maps IAC idx to first idx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    """

    for cod in rows:
        selectedRows = rows[cod]

        count = 0
        for r in selectedRows:
            line = selectedRows[r].strip().split()
            idx = line[0]
            idx = int(idx)

            count += 1
            if idx != count:
                sys.exit('\nERROR: idx')

            iac = line[2]
            iac = int(iac)

            if iac not in newIacDict:
                print(cod)
                print(line)
                sys.exit('ERROR: No such atom type: {}'.format(iac))

            newIac = newIacDict[iac]
            line[2] = newIac
            s = formatStrCHG(line)

            selectedRows[r] = s
            mtb[r] = s


def writeMTB(mtb, fileName):
    """Writes rows from list to topology file.

    :param mtb: (list) Rows of mtb file.
    :param fileName: (str) Name of topology file.
    """

    with open(fileName, 'w') as out:
        for row in mtb:
            out.write(row)
