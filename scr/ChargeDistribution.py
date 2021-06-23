"""Methods for charge distribution
and assignment of charge groups.

Classes:
--------
    ChargeDistributionMethod
    BondChargeDistributionMethod
    NONE
    EEM
    QEqAtomic
    QEqBond
    QEqBond
    AACT
    SQE

    Charge_group_type
    Atomic
    Halo
    AA_Alk
    O_N
"""

from abc import ABC, abstractmethod
from collections import Counter
import numpy as np
import math
import sys

import myExceptions


class ChargeDistributionMethod(ABC):
    """Base class for charge distribution methods.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
        solve(electronegativity, JMatrix, n_atoms, mol):
            Solves system of equations.
        coulomb(cg, atomTypes, mol):
            Computes coulomb interaction terms between any two atoms.
        get_object(n):
            Returns object of class that implements base class ChargeDistributionMethod.
    """

    @abstractmethod
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order: (int) Max order.
        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :param kappa: (float) Kappa parameter (not used here).
        :param lam: (float) Lambda parameter (not used here).
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        pass

    @abstractmethod
    def solve(self, electronegativity, JMatrix, n_atoms, mol):
        """Solves system of equations.

        :param electronegativity: (numpy.ndarray) Array of electronegativities.
        :param JMatrix: (numpy.ndarray)
        :param n_atoms: (int) Number of atoms.
        :param mol: (Molecule) Molecule.
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        pass

    @staticmethod
    def get_object(n):
        """Returns object of class that implements base class ChargeDistributionMethod.

        :param n: (str) Code for class that implements base class ChargeDistributionMethod.
        :return: One of the classes that implements ChargeDistributionMethod.
        """

        classes = {'NONE': NONE, 'EEM': EEM, 'QEqAtomic': QEqAtomic, 'QEqBond': QEqBond, 'AACT': AACT, 'SQE': SQE}

        if n not in classes:
            raise myExceptions.ClassNotImplemented(n, 'ChargeDistribution.ChargeDistributionMethod')
        cr = classes[n]()
        return cr

    def checkDim(self, obj, N):
        """Generic material for dimensionality checks.

        :param obj: (arr) 1-D or 2-D array.
        :param N: (int) Array dimension.
        """

        if len(obj.shape) == 1:
            assert len(obj) == N, "Error: got vector of length" + str(len(obj)) + ", expected " + str(N)
        else:  # vec is a matrix
            assert obj.shape == (N, N), "Error: got matrix of shape" + str(obj.shape) + ", expected (" + str(
                N) + "," + str(N) + ")"

    def getParameters(self, mol, cg, atom_types):
        """Returns diameter, hardness and electronegativity of atoms in the given charge group.

        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :return:
            selected_connectivity: (numpy.ndarray)
            selected_distance_matrix: (numpy.ndarray)
            diameters: (numpy.ndarray)
            hardness: (numpy.ndarray)
            electronegativity (numpy.ndarray)
        """

        atoms = mol.atoms
        nAtoms = cg.shape[0]
        hardness = np.zeros(nAtoms)
        diameters = np.zeros(nAtoms)
        electronegativity = np.zeros(nAtoms)
        selected_distance_matrix = np.zeros((nAtoms, nAtoms))
        selected_connectivity = np.zeros((nAtoms, nAtoms), dtype=bool)

        distance_matrix = mol.distance_matrix
        connectivity = mol.connectivity

        for i in range(cg.shape[0]):
            idx1 = cg[i]
            atom1 = atoms[idx1]
            iac1 = atom1.iac

            atom_type1 = atom_types[iac1]
            hrd1 = atom_type1.hrd.cur
            eln1 = atom_type1.eln.cur
            sig1 = atom_type1.sig.cur
            vdw = atom_type1.vdw

            hardness[i] = hrd1
            diameters[i] = vdw
            electronegativity[i] = eln1

            for j in range(i, cg.shape[0]):
                idx2 = cg[j]

                pos1 = idx1 - 1
                pos2 = idx2 - 1

                selected_distance_matrix[i, j] = distance_matrix[pos1, pos2]
                selected_distance_matrix[j, i] = distance_matrix[pos1, pos2]

                selected_connectivity[i, j] = connectivity[pos1, pos2]
                selected_connectivity[j, i] = connectivity[pos1, pos2]

        return selected_connectivity, selected_distance_matrix, diameters, hardness, electronegativity

    def assignCharges(self, cg, mol, charges):
        """Assigns atomic partial charges to atoms in the given charge group.

        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param mol: (Molecule) Molecule.
        :param charges: (numpy.ndarray) Array of charges.
        :return:
        """

        atoms = mol.atoms
        n_atoms = cg.shape[0]

        for i in range(len(cg)):
            idx = cg[i]
            atom = atoms[idx]
            # atom.charge.cur = Decimal(charges[i])
            atom.charge.cur = charges[i]

        # use only 3 digits for charge
        # qtot = Decimal(0.0)
        qtot = 0.0
        for idx in cg:
            qtot += atoms[idx].charge.cur

        for idx in cg:
            atoms[idx].charge.cur -= qtot / n_atoms

        for idx in cg:
            # atoms[idx].charge.cur = int(atoms[idx].charge.cur * Decimal(1000.0)) / Decimal(1000.0)
            atoms[idx].charge.cur = int(atoms[idx].charge.cur * 1000.0) / 1000.0

        # qtot = Decimal(0.0)
        qtot = 0.0
        for idx in cg:
            qtot += atoms[idx].charge.cur

        # atoms[cg[0]].charge.cur -= qtot
        # add extra (-) charge to first distinct charge.
        charges = []
        for idx in cg:
            charges.append(atoms[idx].charge.cur)
        count = Counter(charges)
        for idx in cg:
            chg = atoms[idx].charge.cur
            n = count[chg]
            if n == 1:
                atoms[idx].charge.cur -= qtot
                break

        charges = []
        for idx in cg:
            chg = atoms[idx].charge.cur
            charges.append(chg)
        qtot = sum(charges)
        if abs(qtot) > 1.00E-10:
            sys.exit('Q = {} for {}'.format(qtot, mol.cod))


    def coulombIntegrals(self, maxOrder, N, connectivity, distance_matrix, diameters):
        """Computes coulomb integrals.

        :param maxOrder: (float) Max order.
        :param N: (int) number of atoms.
        :param connectivity: (numpy.ndarray) Connectivity matrix.
        :param distance_matrix (numpy.ndarray) Distance matrix.
        :param diameters: (numpy.ndarray) Diameter of atoms.
        :return: (numpy.ndarray) Coulomb integrals.
        """

        coulomb = np.zeros((N, N))

        # warnings.filterwarnings('error')
        # try:
        for i in range(0, N):
            for j in range(i+1, N):
                if connectivity[i,j] <= maxOrder:
                    distance = distance_matrix[i, j]
                    diameter_sum = diameters[i] + diameters[j]

                    if distance == 0.0:
                        coulomb[i, j] = 0.0

                    elif diameter_sum == 0.0:
                        coulomb[i, j] = 0.0

                    else:
                        coulomb[i,j] = 1.44 * 1. / distance_matrix[i, j] * math.erf(distance_matrix[i,j] \
                            / np.sqrt(diameters[i]**2 + diameters[j]**2))

                    coulomb[j, i] = coulomb[i, j]

        # except Warning:
        #     print(mol.cod)
        #     sys.exit(123)

        return coulomb


class BondChargeDistributionMethod(ChargeDistributionMethod):
    """Base class for charge distribution methods with bond contribution.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
        solve(bondElneg, bondJMatrix, B, N):
            Solves system of equations.
        bondVars(charge_transfer_topology):
        toAtomicCharges(self, bondCharges, bVars, netCharge, N):
        bondElectronegativity(self, electronegativity, bVars, B):

    """

    @abstractmethod
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order:
        :param mol:
        :param cg:
        :param atom_types:
        :param kappa:
        :param lam:
        :return:
        """

        pass

    # B x B system of equations
    # returns bond charges
    def solve(self, bondElneg, bondJMatrix, B, N):
        """Solves system of equations.

        :param bondElneg: (numpy.ndarray)
        :param bondJMatrix: (bondJMatrix)
        :param B: (int)
        :param N: (int)
        :return:
            bondCharges: (numpy.ndarray) Bond charges.
        """

        # system only has an unique solution for N-1 bond variables
        # otherwise we use SVD to get a particular solution instead
        if B <= N - 1:
            bondCharges = np.linalg.solve(bondJMatrix, -bondElneg)
        else:
            # https://stackoverflow.com/questions/59292279/solving-linear-systems-of-equations-with-svd-decomposition
            U, s, Vh = np.linalg.svd(bondJMatrix)
            c = np.dot(U.T, -bondElneg)
            w = np.dot(np.diag(1 / s), c)
            bondCharges = np.dot(Vh.conj().T, w)

        return bondCharges

    @staticmethod
    def bondVars(charge_transfer_topology):
        """Gets bond variable definitions as pairs of indices.

        :param charge_transfer_topology: (numpy.ndarray)
        :return:
            bVars: (numpy.ndarray)
            B: (int)
        """

        bVars = np.argwhere(charge_transfer_topology)
        upperTriangle = np.where(bVars[:, 0] < bVars[:, 1])
        bVars = bVars[upperTriangle]
        B = len(bVars)
        return bVars, B

    def toAtomicCharges(self, bondCharges, bVars, netCharge, N):
        """Maps bond charges to atoms.

        :param bondCharges: (numpy.ndarray)
        :param bVars: (numpy.ndarray)
        :param netCharge: (int) Net charge.
        :param N: (int) Number of atoms.
        :return:
            charges (numpy.ndarray)
        """
        charges = np.repeat(netCharge / N, N)
        for b, [i, j] in enumerate(bVars):
            charges[i] = charges[i] - bondCharges[b]
            charges[j] = charges[j] + bondCharges[b]

        return charges

    def bondElectronegativity(self, electronegativity, bVars, B):
        """Electronegativity vector in bond variables.

        :param electronegativity: (numpy.ndarray) Array of electronegativities.
        :param bVars: (numpy.ndarray)
        :param B: (int)
        :return:
            bondElneg: (numpy.ndarray)
        """

        bondElneg = np.zeros(B)
        for b, [i, j] in enumerate(bVars):
            bondElneg[b] = electronegativity[j] - electronegativity[i]
        return bondElneg

    def calcBondJMatrix(self, JMatrix, bVars, B):
        """Transforms J Matrix to bond space
        hardness: 1D vector (N)
        res: B x B (B: #entries in upper triangle of CTT matrix)

        :param JMatrix: (numpy.ndarray)
        :param bVars: (numpy.ndarray)
        :param B: (int)
        :return:
            bondJMatrix: (numpy.ndarray)
        """

        # initialize B x B hardness matrix
        bondJMatrix = np.zeros((B, B))

        # construct hardness matrix
        for b1, [i, j] in enumerate(bVars):
            for b2, [k, l] in enumerate(bVars):
                bondJMatrix[b1, b2] = JMatrix[i, k] - JMatrix[i, l] - JMatrix[j, k] + JMatrix[j, l]

        return bondJMatrix


class NONE(ChargeDistributionMethod):
    """Do not perform charge distribution.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Ignores charge distribution.
        solve(electronegativity, JMatrix, n_atoms, mol):
            Ignores charge distribution.
    """

    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Ignores charge distribution."""
        pass

    def solve(self, electronegativity, JMatrix, n_atoms, mol):
        """Ignores charge distribution."""
        pass


class EEM(ChargeDistributionMethod):
    """Compute atomic partial charges via the Electrogenative Equalization (EE) method.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
        solve(electronegativity, JMatrix, n_atoms, mol):
            Solves system of equations.
    """

    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order: (int) Max order.
        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :param kappa: (float) Kappa parameter (not used here).
        :param lam: (float) Lambda parameter (not used here).
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        n_atoms = cg.shape[0]
        selected_connectivity, selected_distance_matrix, diameters, hardness, electronegativity = self.getParameters(
            mol, cg, atom_types)

        # atomic J Matrix and Coulomb integrals
        coulomb = self.coulombIntegrals(max_order, n_atoms, selected_connectivity, selected_distance_matrix, diameters)
        JMatrix = np.diag(hardness) + coulomb
        charges = self.solve(electronegativity, JMatrix, n_atoms, mol)
        self.assignCharges(cg, mol, charges)
        return charges

    def solve(self, electronegativity, JMatrix, n_atoms, mol):
        """Solves system of equations.

        :param electronegativity: (numpy.ndarray) Array of electronegativities.
        :param JMatrix: (numpy.ndarray)
        :param n_atoms: (int) Number of atoms.
        :param mol: (Molecule) Molecule.
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        # prepare augmented hardness matrix
        X = np.zeros((n_atoms + 1, n_atoms + 1,))
        X[:,:-1][:-1] = JMatrix
        X[-1] = 1
        X[:, -1] = -1
        X[-1, -1] = 0

        # prepare vector with electronegativities
        Y = np.zeros(n_atoms + 1)
        Y[:-1] = -electronegativity
        Y[-1] = mol.net_charge

        res = np.linalg.solve(X, Y)
        charges = res[:-1]
        # electronegativityEq = res[-1]

        return charges


class QEqAtomic(ChargeDistributionMethod):
    """Compute atomic partial charges via the QEqAtomic method.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
        solve(electronegativity, JMatrix, n_atoms, mol):
            Solves system of equations.
    """

    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order: (int) Max order.
        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :param kappa: (float) Kappa parameter (not used here).
        :param lam: (float) Lambda parameter (not used here).
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        # Same as for EEM.
        n_atoms = cg.shape[0]
        selected_connectivity, selected_distance_matrix, diameters, hardness, electronegativity = self.getParameters(
            mol, cg, atom_types)

        # compute Coulomb integrals
        coulomb = self.coulombIntegrals(max_order, n_atoms, selected_connectivity, selected_distance_matrix, diameters)
        JMatrix = np.diag(hardness) + coulomb
        charges = self.solve(electronegativity, JMatrix, n_atoms, mol)

        self.assignCharges(cg, mol, charges)

        return charges

    def solve(self, electronegativity, JMatrix, n_atoms, mol):
        """Solves system of equations.

        :param electronegativity: (numpy.ndarray) Array of electronegativities.
        :param JMatrix: (numpy.ndarray)
        :param n_atoms: (int) Number of atoms.
        :param mol: (Molecule) Molecule.
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        # prepare hardness matrix
        X = JMatrix.copy()
        X = X - X[0]
        X[0] = 1

        # prepare vector with electronegativity differences
        Y = -electronegativity.copy()
        Y = Y - Y[0]
        Y[0] = mol.net_charge

        # solve system of equations
        charges = np.linalg.solve(X, Y)
        return charges


class QEqBond(BondChargeDistributionMethod):
    """Compute atomic partial charges via the QEqBond method.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
    """

    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order: (int) Max order.
        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :param kappa: (float) Kappa parameter.
        :param lam: (float) Lambda parameter.
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        n_atoms = cg.shape[0]
        net_charge = mol.net_charge
        charge_transfer_topology = mol.charge_transfer_topology
        selected_connectivity, selected_distance_matrix, diameters, hardness, electronegativity = self.getParameters(
            mol, cg, atom_types)

        # atomic J Matrix
        coulomb = self.coulombIntegrals(max_order, n_atoms, selected_connectivity, selected_distance_matrix, diameters)
        JMatrix = np.diag(hardness) + coulomb

        # transform to bond variables
        bVars, B = self.bondVars(charge_transfer_topology)
        bondElneg = self.bondElectronegativity(electronegativity, bVars, B)
        bondJMatrix = self.calcBondJMatrix(JMatrix, bVars, B)

        # solve system
        bondCharges = self.solve(bondElneg, bondJMatrix, B, n_atoms)
        charges = self.toAtomicCharges(bondCharges, bVars, net_charge, n_atoms)

        self.assignCharges(cg, mol, charges)
        return charges


class AACT(BondChargeDistributionMethod):
    """Compute atomic partial charges via the AACT method.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
    """

    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order: (int) Max order.
        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :param kappa: (float) Kappa parameter.
        :param lam: (float) Lambda parameter.
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        n_atoms = cg.shape[0]
        net_charge = mol.net_charge
        charge_transfer_topology = mol.charge_transfer_topology
        selected_connectivity, selected_distance_matrix, diameters, hardness, electronegativity = self.getParameters(
            mol, cg, atom_types)
        bondHardness = mol.get_bond_hardness(hardness)

        # here, the atomic J matrix has 0 on the diagonal
        JMatrix = self.coulombIntegrals(max_order, n_atoms, selected_connectivity, selected_distance_matrix, diameters)

        # transform to bond variables
        bVars, B = self.bondVars(charge_transfer_topology)
        # self.checkDim(self.bondHardness, self.B)
        bondElneg = self.bondElectronegativity(electronegativity, bVars, B)
        bondJMatrix = self.calcBondJMatrix(JMatrix, bVars, B)

        # add bond hardness on the diagonal
        bondJMatrix = bondJMatrix + 2 * np.diag(bondHardness)

        # solve system
        bondCharges = self.solve(bondElneg, bondJMatrix, B, n_atoms)
        charges = self.toAtomicCharges(bondCharges, bVars, net_charge, n_atoms)

        self.assignCharges(cg, mol, charges)
        return charges


class SQE (BondChargeDistributionMethod):
    """Compute atomic partial charges via the SQE method.

    Methods:
    --------
        compute(max_order, mol, cg, atom_types, kappa, lam):
            Computes charges via implemented charge distribution method.
    """

    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        """Computes charges via implemented charge distribution method.

        :param max_order: (int) Max order.
        :param mol: (Molecule) Molecule.
        :param cg: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
        :param atom_types: (collections.OrderedDict) Ordered dictionary of atom types.
        :param kappa: (float) Kappa parameter.
        :param lam: (float) Lambda parameter.
        :return:
            charges: (numpy.ndarray) Computed charges.
        """

        n_atoms = cg.shape[0]
        net_charge = mol.net_charge
        charge_transfer_topology = mol.charge_transfer_topology
        selected_connectivity, selected_distance_matrix, diameters, hardness, electronegativity = self.getParameters(
            mol, cg, atom_types)
        bondHardness = mol.get_bond_hardness(hardness)

        # atomic J matrix with diagonal scaled by lam^2
        coulomb = self.coulombIntegrals(max_order, n_atoms, selected_connectivity, selected_distance_matrix, diameters)
        scalingFactor1 = lam * lam
        JMatrix = scalingFactor1 * np.diag(hardness) + coulomb

        # transform to bond variables
        bVars, B = self.bondVars(charge_transfer_topology)
        # self.checkDim(self.bondHardness, self.B)
        bondElneg = self.bondElectronegativity(electronegativity, bVars, B)
        bondJMatrix = self.calcBondJMatrix(JMatrix, bVars, B)

        # add bond hardness on the diagonal, scaled by kappa^2
        scalingFactor2 = kappa * kappa
        bondJMatrix = bondJMatrix + scalingFactor2 * 2 * np.diag(bondHardness)

        # solve system
        bondCharges = self.solve(bondElneg, bondJMatrix, B, n_atoms)
        charges = self.toAtomicCharges(bondCharges, bVars, n_atoms, n_atoms)

        self.assignCharges(cg, mol, charges)
        return charges


class Charge_group_type(ABC):
    """Assigns atoms to the charge group.

    Methods:
    --------
        setCG(mol):
            Defines charge groups of molecule.
    """

    @abstractmethod
    def setCG(self, mol):
        """Defines charge groups of molecule.

        :param mol: (Molecule) Molecule.
        """
        pass

    @staticmethod
    def get_object(n):
        """Returns object of class that implements base class Charge_group_type.

        :param n: (str) Code for class that implements base class Charge_group_type.
        :return: One of the classes that implements Charge_group_type.
        """

        classes = {'ATOMIC': Atomic, 'HALO': Halo, 'AA-ALK': AA_Alk, 'O_N': O_N, 'MIX': Mix}

        if n not in classes:
            raise myExceptions.ClassNotImplemented(n, 'ChargeDistribution.Charge_group_type')
        cr = classes[n]()
        return cr


class Atomic(Charge_group_type):
    """Defines one charge group for the whole molecule.

    Methods:
    --------
        setCG(mol):
            Defines charge groups of molecule.
    """

    def setCG(self, mol):
        """Defines charge groups of molecule.

        :param mol: (Molecule) Molecule.
        """

        CGs = []
        atoms = mol.atoms

        # add all atoms
        indexes = []
        for idx1 in atoms:
            indexes.append(idx1)

        CGs.append(indexes)
        mol.CGs = CGs


class Halo(Charge_group_type):
    """Defines charge groups for HAL family.

    Methods:
    --------
        setCG(mol):
            Defines charge groups of molecule.
    """

    def setCG(self, mol):
        """Defines charge groups of molecule.

        :param mol: (Molecule) Molecule.
        """

        mol.CGs = []


class AA_Alk(Charge_group_type):
    """Defines charge groups for ALK family.

    Methods:
    --------
        setCG(mol):
            Defines charge groups of molecule.
    """

    def setCG(self, mol):
        """Defines charge groups of molecule.

        :param mol: (Molecule) Molecule.
        """

        CGs = []
        atoms = mol.atoms

        # add carbon atoms
        for idx1 in atoms:
            atom = atoms[idx1]
            iac = atom.iac
            if iac != 20:
                CGs.append([idx1])

            # add hydrogen atom
            else:
                for i in range(len(CGs)):
                    for j in range(len(CGs[i])):
                        idx2 = CGs[i][j]
                        if atom.isNrmNB(idx2):
                            CGs[i].append(idx1)

        mol.CGs = CGs


class O_N(Charge_group_type):
    """Defines charge groups for O+N family.

    Methods:
    --------
        setCG(mol):
            Defines charge groups of molecule.
    """

    RO = [100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118]
    CC_EST = [119, 120, 121, 122]
    C_EST = [102]

    CARBON_OH = [111, 112, 113, 114]
    CARBON_NH = [138, 139, 140, 141, 142, 143, 144]
    CARBON_NOH = [146, 147, 148, 149]
    N_amid = 6
    N_amin = 7
    O_double_C_acid = 1  # O=C
    C_double_O_acid = 12  # C=O
    OH_aci = 98
    O_double_C_amid = 4  # O=C
    C_double_O_amid = 18  # C=O

    def setCG(self, mol):
        """Defines charge groups of molecule.

        :param mol: (Molecule) Molecule.
        """

        atoms = mol.atoms

        indexes = [[]]
        for idx1 in atoms:
            a = atoms[idx1]
            a.used = 0

        # add first non fixed atom to first array
        for idx1 in atoms:
            a = atoms[idx1]
            if not a.ignore:
                indexes[0].append(idx1)
                a.used = 1
                break

        for idx1 in atoms:
            atm1 = atoms[idx1]

            if atm1.ignore is False:

                iac1 = atm1.iac
                pos = -1
                # loop over list of lists
                for i in range(len(indexes)):
                    # loop over current list
                    for j in range(len(indexes[i])):

                        idx2 = indexes[i][j]
                        bond_1_2 = mol.are_bonded(idx1, idx2)

                        if idx1 == idx2:
                            continue

                        elif iac1 == self.C_double_O_acid:
                            break

                        elif iac1 == self.OH_aci:
                            break

                        elif iac1 == self.O_double_C_amid:
                            break

                        elif iac1 == self.N_amid:
                            break

                        elif iac1 == self.N_amin or iac1 in self.CARBON_NH:
                            # check if there is a nb already in indexes
                            for x in range(len(indexes)):
                                for y in range(len(indexes[x])):
                                    idx3 = indexes[x][y]
                                    bond_1_3 = mol.are_bonded(idx1, idx3)

                                    if bond_1_3:
                                        iac3 = atoms[idx3].iac
                                        if iac1 in self.CARBON_NH and iac3 in self.CARBON_NH:
                                            continue
                                        atm1.used = 1
                                        pos = x
                                        break

                            if pos == -1:
                                break

                        elif bond_1_2:
                            atm2 = atoms[idx2]
                            iac2 = atm2.iac

                            if iac1 in self.CARBON_OH and iac2 in self.CARBON_OH:
                                continue

                            elif iac1 == self.C_double_O_amid and iac2 == self.C_double_O_amid:
                                continue

                            elif (iac1 != self.O_double_C_acid) and (iac2 == self.C_double_O_acid):
                                continue

                            elif iac1 in self.RO and iac2 in self.RO:
                                continue

                            elif iac1 in self.CC_EST:
                                continue

                            atm1.used = 1
                            pos = i
                            break

                        # 2 esters = C=O should be in the same CG
                        # check is nbs of C=O are already connected to any other ester group
                        elif iac1 in self.C_EST:
                            neighbors_indexes = mol.get_neighbors_of_atom(idx1)

                            for idx_nb in neighbors_indexes:
                                atom_nb = atoms[idx_nb]
                                iac_nb = atom_nb.iac

                                for x in range(len(indexes)):
                                    for y in range(len(indexes[x])):
                                        other_idx_tmp = indexes[x][y]
                                        are_bonded = mol.are_bonded(idx_nb, other_idx_tmp)

                                        if idx_nb == other_idx_tmp:
                                            continue

                                        elif atom_nb.ignore:
                                            break

                                        elif iac_nb in self.C_EST:
                                            break

                                        elif iac_nb in self.CC_EST:
                                            break

                                        elif are_bonded:
                                            atm1.used = 1
                                            pos = x
                                            break

                    # else:
                    # nothing to do

                    if pos != -1:
                        indexes[pos].append(idx1)
                        break

                if atm1.used == 0:
                    indexes.append([idx1])

        mol.CGs = indexes


class Mix(Charge_group_type):
    """Defines charge groups for MIX family (ALk + HAL + O+N + MIX).

    Methods:
    --------
        setCG(mol):
            Defines charge groups of molecule.
    """

    X = [32, 33, 34, 67]
    CX = [45, 39, 69, 68]

    RO = [100, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118]
    CC_EST = [119, 120, 121, 122]
    C_EST = [102]

    CARBON_OH = [111, 112, 113, 114]
    CARBON_NH = [138, 139, 140, 141, 142, 143, 144]
    CARBON_NOH = [146, 147, 148, 149]
    N_amid = 6
    N_amin = 7
    O_double_C_acid = 1  # O=C
    C_double_O_acid = 12  # C=O
    OH_aci = 98
    O_double_C_amid = 4  # O=C
    C_double_O_amid = 18  # C=O

    def setCG(self, mol):
        """Defines charge groups of molecule.

        :param mol: (Molecule) Molecule.
        """

        atoms = mol.atoms

        indexes = []
        for idx1 in atoms:
            a = atoms[idx1]
            a.used = 0

        # Try first to add all CX atoms
        for idx1 in atoms:
            a = atoms[idx1]
            if (not a.ignore) and (a.iac in self.CX):
                indexes.append([idx1])
                a.used = 1

        # If list is still empty, add first non fixed atom to first array
        if len(indexes) == 0:
            for idx1 in atoms:
                a = atoms[idx1]
                if not a.ignore:
                    indexes.append([idx1])
                    a.used = 1
                    break

        for idx1 in atoms:
            atm1 = atoms[idx1]

            if atm1.ignore is False:

                iac1 = atm1.iac
                pos = -1

                if atm1.used == 1:
                    continue

                # loop over list of lists
                for i in range(len(indexes)):

                    # loop over current list
                    for j in range(len(indexes[i])):

                        idx2 = indexes[i][j]
                        bond_1_2 = mol.are_bonded(idx1, idx2)

                        if idx1 == idx2:
                            continue

                        elif iac1 == self.C_double_O_acid:
                            break

                        elif iac1 == self.OH_aci:
                            break

                        elif iac1 == self.O_double_C_amid:
                            break

                        elif iac1 == self.N_amid:
                            break

                        elif iac1 == self.N_amin or iac1 in self.CARBON_NH:
                            # check if there is a nb already in indexes
                            for x in range(len(indexes)):
                                for y in range(len(indexes[x])):
                                    idx3 = indexes[x][y]
                                    bond_1_3 = mol.are_bonded(idx1, idx3)

                                    if bond_1_3:
                                        iac3 = atoms[idx3].iac
                                        if iac1 in self.CARBON_NH and iac3 in self.CARBON_NH:
                                            continue
                                        atm1.used = 1
                                        pos = x
                                        break

                            if pos == -1:
                                break

                        elif bond_1_2:
                            atm2 = atoms[idx2]
                            iac2 = atm2.iac

                            if iac1 in self.CARBON_OH and iac2 in self.CARBON_OH:
                                continue

                            elif iac1 == self.C_double_O_amid and iac2 == self.C_double_O_amid:
                                continue

                            elif (iac1 != self.O_double_C_acid) and (iac2 == self.C_double_O_acid):
                                continue

                            elif iac1 in self.RO and iac2 in self.RO:
                                continue

                            elif (iac1 in self.CARBON_OH and iac2 in self.C_EST) or \
                                (iac2 in self.CARBON_OH and iac1 in self.C_EST):
                                continue

                            elif (iac1 in self.RO and iac2 in self.CARBON_OH) or \
                                    (iac2 in self.RO and iac1 in self.CARBON_OH):
                                continue

                            elif (iac1 in self.CARBON_OH and iac2 in self.CARBON_NH) or \
                                    (iac2 in self.CARBON_OH and iac1 in self.CARBON_NH):
                                continue

                            elif iac1 in self.CC_EST:
                                continue

                            atm1.used = 1
                            pos = i
                            break

                        # 2 esters = C=O should be in the same CG
                        # check is nbs of C=O are already connected to any other ester group
                        elif iac1 in self.C_EST:
                            neighbors_indexes = mol.get_neighbors_of_atom(idx1)

                            for idx_nb in neighbors_indexes:
                                atom_nb = atoms[idx_nb]
                                iac_nb = atom_nb.iac

                                for x in range(len(indexes)):
                                    for y in range(len(indexes[x])):
                                        other_idx_tmp = indexes[x][y]
                                        are_bonded = mol.are_bonded(idx_nb, other_idx_tmp)

                                        if idx_nb == other_idx_tmp:
                                            continue

                                        elif atom_nb.ignore:
                                            break

                                        elif iac_nb in self.C_EST:
                                            break

                                        elif iac_nb in self.CC_EST:
                                            break

                                        elif are_bonded:
                                            atm1.used = 1
                                            pos = x
                                            break

                    # else:
                    # nothing to do

                    if pos != -1:
                        indexes[pos].append(idx1)
                        break

                if atm1.used == 0:
                    indexes.append([idx1])

        mol.CGs = indexes
        # print(indexes, '\n')
