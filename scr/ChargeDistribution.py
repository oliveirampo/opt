"""Compute atomic partial charges via the Electrogenative Equalization (EE) method.

Classes:

    EEM
    NONE
    Halo
    AA_Alk
    O_N
"""

from abc import ABC, abstractmethod
import numpy as np
import math
import sys


class ChargeDistributionMethod(ABC):
    """ Base class for the EE method.

    Methods:
        setCG(mol):
            Assign charge groups of molecule.
        computeEEM(cg, atomTypes, mol):
            Compute charges via EE method.
        coulomb(cg, atomTypes, mol):
            Compute coulomb interaction between two atoms.
    """

    @abstractmethod
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        pass

    @abstractmethod
    def solve(self, electronegativity, JMatrix, n_atoms, mol):
        pass

    def checkDim(self, obj, N):
        '''Generic material for dimensionality checks.
        :param obj:
        :param N:
        :return:
        '''
        if len(obj.shape) == 1:
            assert len(obj) == N, "Error: got vector of length" + str(len(obj)) + ", expected " + str(N)
        else:  # vec is a matrix
            assert obj.shape == (N, N), "Error: got matrix of shape" + str(obj.shape) + ", expected (" + str(
                N) + "," + str(N) + ")"

    def getParameters(self, mol, cg, atom_types):
        atoms = mol.atoms
        nAtoms = cg.shape[0]
        hardness = np.zeros(nAtoms)
        diameters = np.zeros(nAtoms)
        electronegativity = np.zeros(nAtoms)

        for i in range(cg.shape[0]):
            idx1 = cg[i]
            atom1 = atoms[idx1]
            iac1 = atom1.iac

            atom_type1 = atom_types[iac1]
            hrd1 = atom_type1.hrd.cur
            eln1 = atom_type1.eln.cur
            sig1 = atom_type1.sig.cur

            hardness[i] = hrd1
            diameters[i] = sig1
            electronegativity[i] = eln1

        return diameters, hardness, electronegativity

    def assignCharges(self, cg, mol, charges):
        atoms = mol.atoms
        n_atoms = cg.shape[0]

        for i in range(len(cg)):
            idx = cg[i]
            atom = atoms[idx]
            atom.charge.cur = charges[i]

        # use only 3 digits for charge
        qtot = 0.0
        for idx in cg:
            qtot += atoms[idx].charge.cur

        for idx in cg:
            atoms[idx].charge.cur -= qtot / n_atoms

        for idx in cg:
            atoms[idx].charge.cur = int(atoms[idx].charge.cur * 1000.0) / 1000.0

        qtot = 0.0
        for idx in cg:
            qtot += atoms[idx].charge.cur

        atoms[cg[0]].charge.cur -= qtot

    def coulombIntegrals(self, maxOrder, mol, N, diameters):
        connectivity = mol.connectivity
        distance_matrix = mol.distance_matrix

        coulomb = np.zeros((N, N))
        for i in range(0, N):
            for j in range(i+1, N):
                if connectivity[i,j] <= maxOrder:
                    distance = distance_matrix[i, j]
                    if distance == 0.0:
                        coulomb[i, j] = 0.0
                    else:
                        coulomb[i,j] = 1. / distance_matrix[i,j] * math.erf(distance_matrix[i,j] \
                            / np.sqrt(diameters[i]**2 + diameters[j]**2 ))
                    coulomb[j,i] = coulomb[i,j]

        return coulomb


class BondChargeDistributionMethod(ChargeDistributionMethod):
    @abstractmethod
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        pass

    # B x B system of equations
    # returns bond charges
    def solve(self, bondElneg, bondJMatrix, B, N):

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

    # get bond variable definitions as pairs of indices
    @staticmethod
    def bondVars(charge_transfer_topology):
        bVars = np.argwhere(charge_transfer_topology)
        upperTriangle = np.where(bVars[:, 0] < bVars[:, 1])
        bVars = bVars[upperTriangle]
        B = len(bVars)
        return bVars, B

    # map bond charges to atoms
    def toAtomicCharges(self, bondCharges, bVars, netCharge, N):
        charges = np.repeat(netCharge / N, N)
        for b, [i, j] in enumerate(bVars):
            charges[i] = charges[i] - bondCharges[b]
            charges[j] = charges[j] + bondCharges[b]
        return charges

    # electronegativity vector in bond variables
    def bondElectronegativity(self, electronegativity, bVars, B):
        bondElneg = np.zeros(B)
        for b, [i, j] in enumerate(bVars):
            bondElneg[b] = electronegativity[j] - electronegativity[i]
        return bondElneg

    # transform J Matrix to bond space
    # hardness: 1D vector (N)
    # res: B x B (B: #entries in upper triangle of CTT matrix)
    def calcBondJMatrix(self, JMatrix, bVars, B):

        # initialize B x B hardness matrix
        bondJMatrix = np.zeros((B, B))

        # construct hardness matrix
        for b1, [i, j] in enumerate(bVars):
            for b2, [k, l] in enumerate(bVars):
                bondJMatrix[b1, b2] = JMatrix[i, k] - JMatrix[i, l] - JMatrix[j, k] + JMatrix[j, l]

        return bondJMatrix


class NONE(ChargeDistributionMethod):
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        pass

    def solve(self, electronegativity, JMatrix, n_atoms, mol):
        pass


class EEM(ChargeDistributionMethod):
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        n_atoms = cg.shape[0]
        diameters, hardness, electronegativity = self.getParameters(mol, cg, atom_types)

        # atomic J Matrix and Coulomb integrals
        coulomb = self.coulombIntegrals(max_order, mol, n_atoms, diameters)
        JMatrix = np.diag(hardness) + coulomb
        charges, electronegativityEq = self.solve(electronegativity, JMatrix, n_atoms, mol)
        self.assignCharges(cg, mol, charges)
        return charges

    def solve(self, electronegativity, JMatrix, n_atoms, mol):
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
        electronegativityEq = res[-1]

        return charges, electronegativityEq


class QEqAtomic(ChargeDistributionMethod):
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        # Same as for EEM.
        n_atoms = cg.shape[0]
        diameters, hardness, electronegativity = self.getParameters(mol, cg, atom_types)

        # compute Coulomb integrals
        coulomb = self.coulombIntegrals(max_order, mol, n_atoms, diameters)
        JMatrix = np.diag(hardness) + coulomb
        charges = self.solve(electronegativity, JMatrix, n_atoms, mol)

        self.assignCharges(cg, mol, charges)

        return charges

    def solve(self, electronegativity, JMatrix, n_atoms, mol):
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
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        n_atoms = cg.shape[0]
        net_charge = mol.net_charge
        charge_transfer_topology = mol.charge_transfer_topology
        diameters, hardness, electronegativity = self.getParameters(mol, cg, atom_types)

        # atomic J Matrix
        coulomb = self.coulombIntegrals(max_order, mol, n_atoms, diameters)
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
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        n_atoms = cg.shape[0]
        net_charge = mol.net_charge
        charge_transfer_topology = mol.charge_transfer_topology
        diameters, hardness, electronegativity = self.getParameters(mol, cg, atom_types)
        bondHardness = mol.get_bond_hardness(hardness)

        # here, the atomic J matrix has 0 on the diagonal
        JMatrix = self.coulombIntegrals(max_order, mol, n_atoms, diameters)

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
    def compute(self, max_order, mol, cg, atom_types, kappa, lam):
        n_atoms = cg.shape[0]
        net_charge = mol.net_charge
        charge_transfer_topology = mol.charge_transfer_topology
        diameters, hardness, electronegativity = self.getParameters(mol, cg, atom_types)
        bondHardness = mol.get_bond_hardness(hardness)

        # atomic J matrix with diagonal scaled by lam^2
        coulomb = self.coulombIntegrals(max_order, mol, n_atoms, diameters)
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
    @abstractmethod
    def setCG(self, mol):
        pass


class Atomic(Charge_group_type):
    def setCG(self, mol):
        CGs = []
        atoms = mol.atoms

        # add all atoms
        indexes = []
        for idx1 in atoms:
            indexes.append(idx1)

        CGs.append(indexes)
        mol.CGs = CGs


class Halo(Charge_group_type):
    def setCG(self, mol):
        mol.CGs = []


class AA_Alk(Charge_group_type):
    def setCG(self, mol):
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

                                    if atm1.isNrmNB(idx3):
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
