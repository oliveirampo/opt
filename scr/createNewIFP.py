"""Creates new IFP file with new dummy atom types.

Classes:
    AtomType

Methods:

"""

import numpy as np
import math
import sys

import effectiveParameter


class AtomType(object):
    """Class to hold information about atom types in IFP file."""

    def __init__(self, atomNum, atomName, sqrtC06, sqrtC12_1, sqrtC12_2, sqrtC12_3, sqrtC06NB, sqrtC12NB, matrix):
        """

        :param atomNum: (str) Atom number.
        :param atomName: (str) Atom name.
        :param sqrtC06: (str) Squared root of C06.
        :param sqrtC12_1: (str) Squared root of C12 (type 1).
        :param sqrtC12_2: (str) Squared root of C12 (type 2).
        :param sqrtC12_3: (str) Squared root of C12 (type 3).
        :param sqrtC06NB: (str) Squared root of C06 (1-4 LJ interaction).
        :param sqrtC12NB: (str) Squared root of C06 (1-4 LJ interaction).
        :param matrix: (list) Selection matrix for the repulsive Lennard-Jones interaction parameters.
        """

        self.atomNum = atomNum
        self.atomName = atomName
        self.sqrtC06 = sqrtC06
        self.sqrtC12_1 = sqrtC12_1
        self.sqrtC12_2 = sqrtC12_2
        self.sqrtC12_3 = sqrtC12_3
        self.sqrtC06NB = sqrtC06NB
        self.sqrtC12NB = sqrtC12NB
        self.matrix = matrix

    def __str__(self):
        return '[' + self.atomNum + '  ' + self.atomName + '  ' + self.sqrtC06 + '  ' + self.sqrtC12_1 + ']'

    def addToMatrix(self, n):
        """Adds element to matrix.

        :param n: (str) Matrix element.
        """

        self.matrix.append(n)


def create(conf, parameters):
    """Creates new IFP file with new dummy atom types.

    :param conf: (configuration.Conf) Configuration object
    :param parameters: (collections.OrderedDict) Ordered dictionary of atom types.
    :return
        newIacDict: (dict) Dictionary that maps IAC idx to first idx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    """

    fileName = 'prm/2016H66_upd.ifp'
    outName = 'prm/2016H66_upd_NEW.ifp'

    row1, row2, block = readFile(fileName)

    atomTypes = readBlock(block)

    minimalSetIAC, newIacDict = getLJAtomTypes(parameters)

    N = len(minimalSetIAC)
    addAtomType(conf, atomTypes, N, minimalSetIAC, newIacDict, parameters)

    check(atomTypes)

    writeOut(row1, row2, atomTypes, fileName, outName)

    # add other IACs.
    addOtherIAC(newIacDict)

    return newIacDict


def addOtherIAC(newIacDict):
    """Adds indexes of other IAC to dictionary.

    :param newIacDict:
    :return:
        newIacDict: (dict) Dictionary that maps IAC idx to first idx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    """

    # maxIac = 0
    # for key in newIacDict:
    #     iac = newIacDict[key]
    #     if iac > maxIac:
    #         maxIac = iac
    #
    # nIac = maxIac - len(minimalSetIAC)
    #
    # for iac in range(1, nIac + 1):
    #     if iac not in newIacDict:
    #         newIacDict[iac] = iac

    # add united-atom
    for iac in [13, 14, 15, 16]:
        if iac not in newIacDict:
            newIacDict[iac] = iac

    # specific for O+N family
    for pair in [[103, 13], [104, 14], [105, 15], [106, 16],
                 [115, 13], [116, 14], [117, 15], [118, 16], [119, 13], [120, 14], [121, 15], [122, 16],
                 [21, 21], [111, 13], [112, 14], [113, 15], [114, 16],
                 [97, 21],
                 [145, 21], [144, 13], [142, 14], [140, 15], [138, 16],
                 [150, 21], [146, 13], [147, 14], [148, 15], [149, 16]]:

        oldIAC = pair[0]
        newIAC = pair[1]
        newIacDict[oldIAC] = newIAC


def getLJAtomTypes(parameters):
    """Returns unique set of IAC parameters.
    Ignore those whose rng_sig == 0.
    Only use first parameter if they share the same sig and eps values.

    :param parameters: (collections.OrderedDict) Ordered dictionary of atom types.
    :return:
        minimalSetIAC: (list) List with indexes of parameters that should be added to new IFP file
                        and have non-similar sig values.
        newIacDict: (dict) Dictionary that maps IAC idx to first idx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    """

    newIacDict = {}
    minimalSetIAC = []
    symmetrySigTypes = {}

    # check for symmetry
    for iacIdx in parameters:
        atmTyp = parameters[iacIdx]
        if not atmTyp.sig.rng:
            continue

        if atmTyp.sig.hasSymmetricIAC():
            symm = atmTyp.sig.symmetry

            if symm in symmetrySigTypes:
                symmetrySigTypes[symm].append(iacIdx)
            else:
                symmetrySigTypes[symm] = [iacIdx]

    for iacIdx in parameters:
        atmTyp = parameters[iacIdx]
        if not atmTyp.sig.rng:
            continue

        symm = atmTyp.sig.symmetry
        if symm:
            newIacDict[iacIdx] = symmetrySigTypes[symm][0]
            minimalSetIAC.append(symmetrySigTypes[symm][0])
        else:
            newIacDict[iacIdx] = iacIdx
            minimalSetIAC.append(iacIdx)

    minimalSetIAC = list(set(minimalSetIAC))

    return minimalSetIAC, newIacDict


def readFile(fileName):
    """Reads IFP file.

    :param fileName: (str) Name of IFP file.
    :return:
        row1: (int) Number of row that starts 'SINGLEATOMLJPAIR' block.
        row2: (int) Number of row that ends 'SINGLEATOMLJPAIR' block.
        block: (list) List with rows from 'SINGLEATOMLJPAIR' block.
    """

    with open(fileName, 'r') as f:
        lines = f.readlines()

    lines = [row.strip().split() for row in lines]

    for i in range(len(lines)):
        if lines[i][0] == 'SINGLEATOMLJPAIR':
            row1 = i
            row2, block = getBlock(lines, i + 1)

            return row1, row2, block


def getBlock(lines, row):
    """Helper function to read 'SINGLEATOMLJPAIR' block.

    :param lines: (list) Rows of IFP file.
    :param row: (int) Number of row that starts 'SINGLEATOMLJPAIR' block.
    :return:
        i: (int) Number of row that ends 'SINGLEATOMLJPAIR' block.
        block: (list) List with rows from 'SINGLEATOMLJPAIR' block.
    """

    block = []

    for i in range(row, len(lines)):
        if lines[i][0] == 'END':
            return i, block
        else:
            if lines[i][0] != '#':
                block.append(lines[i])


def readBlock(block):
    """Extrack data from 'SINGLEATOMLJPAIR' block in list format.

    :param block: (list) List with rows from 'SINGLEATOMLJPAIR' block.
    :return: (list) List of AtomType objects.
    """

    atomTypes = []

    for i in range(len(block)):
        # if block[i][0] == '#number':
        #     nAtoms = block[i + 1][0]

        if block[i][0] == '#CS6':
            atm = getAtomType(block, i - 1)
            atomTypes.append(atm)

    return atomTypes


def getAtomType(block, row):
    """Extracts atom type information from 'SINGLEATOMLJPAIR' block in list format.

    :param block: (list) List with rows from 'SINGLEATOMLJPAIR' block.
    :param row: (int) Row number in block list.
    :return:
        atm: (AtomType)
    """

    atomBlock = []

    for i in range(row, len(block)):
        if block[i][0] == '#---':
            break
        else:
            if block[i][0] != '#CS6':
                atomBlock.append(block[i])

    atomNum = atomBlock[0][0]
    atomName = atomBlock[0][1]
    sqrtC06 = atomBlock[0][2]
    sqrtC12_1 = atomBlock[0][3]
    sqrtC12_2 = atomBlock[0][4]
    sqrtC12_3 = atomBlock[0][5]

    sqrtC06NB = atomBlock[1][0]
    sqrtC12NB = atomBlock[1][1]

    matrix = []
    for i in range(2, len(atomBlock)):
        for j in range(len(atomBlock[i])):
            matrix.append(atomBlock[i][j])

    atm = AtomType(atomNum, atomName, sqrtC06, sqrtC12_1, sqrtC12_2, sqrtC12_3, sqrtC06NB, sqrtC12NB, matrix)

    return atm


def addDummyAtomType(atomTypes, N):
    """Adds AtomType object to

    :param atomTypes: (list) List of AtomType objects.
    :param N: (int) Number of new dummy atom types.
    """

    for i in range(1, N + 1):
        # extend matrix of all existing atom type by 1
        for a in atomTypes:
            a.addToMatrix(1)

        # create dummy atom and add it to list of atomTypes
        numAtoms = len(atomTypes)
        num = str(numAtoms + 1)

        name = 'MA' + num
        val = '0.000000e-00'

        matrix = [1 for _ in range(numAtoms)]
        matrix.append(1)
        atomTypes.append(AtomType(num, name, val, val, val, val, val, val, matrix))


def addAtomType(conf, atomTypes, N, minimalSetIAC, newIacDict, parameters):
    """Adds AtomType object to

    :param conf: (configuration.Conf) Configuration object
    :param atomTypes: (list) List of AtomType objects.
    :param N: (int) Number of new dummy atom types.
    :param minimalSetIAC: (list) List with indexes of parameters that should be added to new IFP file
                        and have non-similar sig values.
    :param newIacDict: (dict) Dictionary that maps IAC idx to first dx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    :param parameters: (collections.OrderedDict) Ordered dictionary of atom types.
    """

    cr = conf.cr
    matrix = conf.matrix
    scl_sig_NEI = conf.scl_sig_NEI
    scl_eps_NEI = conf.scl_eps_NEI
    nrmParameter = effectiveParameter.NRM()
    neiParameter = effectiveParameter.NEI()

    for i in range(0, N):
        # extend matrix of all existing atom type by 1
        for a in atomTypes:
            a.addToMatrix(1)

        # create dummy atom and add it to list of atomTypes
        numAtoms = len(atomTypes)
        newIac = str(numAtoms + 1)

        iacIdx = minimalSetIAC[i]
        symmetricIacIdx = newIacDict[iacIdx]
        iac = parameters[symmetricIacIdx]

        # update IAC to its new value.
        oldIac = newIac
        for key in newIacDict:
            value = newIacDict[key]
            if value == iacIdx:
                newIacDict[key] = int(newIac)

        name = iac.typ

        sigi = iac.sig.cur
        epsi = iac.eps.cur
        sigij = cr.getSigma(sigi, sigi)
        epsij = cr.getEpsilon(epsi, epsi, sigi, sigi)

        c12_1 = 4.0 * epsij * math.exp(12.0 * math.log(sigij))
        c6, c12_2 = nrmParameter.computeCR(cr, scl_sig_NEI, scl_eps_NEI, iac, iac, matrix)
        c6_nei, c12_nei = neiParameter.computeCR(cr, scl_sig_NEI, scl_eps_NEI, iac, iac, matrix)
        c12_3 = 0.0

        c6 = toScientific(np.sqrt(c6))
        c12_1 = toScientific(np.sqrt(c12_1))
        c12_2 = toScientific(np.sqrt(c12_2))
        c12_3 = toScientific(np.sqrt(c12_3))
        c6_nei = toScientific(np.sqrt(c6_nei))
        c12_nei = toScientific(np.sqrt(c12_nei))

        currentMatrix = [1 for _ in range(numAtoms)]
        currentMatrix.append(1)

        atomTypes.append(AtomType(newIac, name, c6, c12_1, c12_2, c12_3, c6_nei, c12_nei, currentMatrix))

    fixMatrixEntries(matrix, atomTypes, minimalSetIAC, newIacDict)


def fixMatrixEntries(matrix, atomTypes, minimalSetIAC, newIacDict):
    """Updates entries of matrices for each atom type according to matrix.dat file.

    :param matrix: (Matrix) Matrix with usage of C12(II) parameter.
    :param atomTypes: (list) List of AtomType objects.
    :param minimalSetIAC: (list) List with indexes of parameters that should be added to new IFP file
                        and have non-similar sig values.
    :param newIacDict: (dict) Dictionary that maps IAC idx to first dx in list of symmetric IAC indexes.
                        If there are none, the idx is mapped to itself.
    :return:
    """

    for oldIAC1 in newIacDict:
        newIac1 = newIacDict[oldIAC1]

        if oldIAC1 in matrix:
            pos1 = newIac1 - 1
            atomTyp1 = atomTypes[pos1]
            currentMatrix = atomTyp1.matrix

            if newIac1 != int(atomTyp1.atomNum):
                sys.exit('ERROR with index of atoms: {} != {}'.format(newIac1, atomTyp1.atomNum))

            listIAC = matrix[oldIAC1]
            for oldIAC2 in listIAC:
                newIAC2 = newIacDict[oldIAC2]
                pos2 = newIAC2 - 1

                currentMatrix[pos2] = 2


def toScientific(val):
    return '{:.6e}'.format(val)


def check(atomTypes):
    """Checks if the matrices of AtomType objects have correct size.

    :param atomTypes: (list) List of AtomType objects.
    """

    N = len(atomTypes)

    for a in atomTypes:
        if len(a.matrix) != N:
            print(a)
            sys.exit('Matrix with wrong size\n')


def writeOut(row1, row2, atomTypes, fileName, outName):
    """Writes new IFP file with new dummy atom types.

    :param row1 (int) Number of row that starts 'SINGLEATOMLJPAIR' block.
    :param row2 (int) Number of row that ends 'SINGLEATOMLJPAIR' block.
    :param atomTypes: (list) List of AtomType objects.
    :param fileName: (str) Name of IFP file.
    :param outName: (str) Name of new IFP file.
    """

    N = len(atomTypes)

    with open(fileName, 'r') as f:
        rawLines = f.readlines()

    lines = [raw.strip().split() for raw in rawLines]

    out = open(outName, 'w')
    for i in range(0, row2):
        if i <= row1:
            writeLine(rawLines[i], out)
        else:
            if lines[i][0] == '#' or lines[i][0] == '#number':
                if len(lines[i]) == 2:
                    if lines[i][1] == 'NRATT':
                        writeLine(rawLines[i], out)
                        out.write('{0:>10}\n'.format(N))
                else:
                    writeLine(rawLines[i], out)

    writeAtomTypes(atomTypes, out)

    for i in range(row2, len(lines)):
        writeLine(rawLines[i], out)

    out.close()


def writeLine(line, file):
    """Writes lines to file."""

    for item in line:
        file.write('%s' % item)


def writeAtomTypes(atomTypes, out):
    """Writes AtomType objects to file.

    :param atomTypes: (list) List of AtomType objects.
    :param out: (_io.TextIOWrapper) Output object.
    :return:
    """

    for a in atomTypes:
        out.write("{0:>9} {1:>7} {2:>14} {3:>14} {4:>13} {5:>13}\n"
                  .format(a.atomNum, a.atomName, a.sqrtC06, a.sqrtC12_1, a.sqrtC12_2, a.sqrtC12_3))

        out.write('#CS6 CS12 parameters LJ14PAIR\n')
        out.write("   {0:<12} {1:>14}\n".format(a.sqrtC06NB.strip(), a.sqrtC12NB))

        matrix = a.matrix

        for i in range(1, len(matrix) + 1):
            if (i % 20 == 0) and (i != 0):
                out.write("   %s\n" % matrix[i - 1])
            else:
                out.write("   %s" % matrix[i - 1])

        out.write('\n#---\n')
