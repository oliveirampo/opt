from abc import ABC, abstractmethod
import numpy as np
import sys
import math


class EEM(ABC):
	@abstractmethod
	def setCG(self, mol):
		pass


	def computeEEM(self, cg, atomTypes, mol):
		atoms = mol.atoms
		nAtoms = len(cg)
		X = np.empty((nAtoms + 1, nAtoms + 1,))
		X[:] = np.zeros(nAtoms + 1)
		Y = np.empty((nAtoms + 1, 1,))

		for i in range(len(cg)):
			idx1 = cg[i]
			atom1 = atoms[idx1]
			iac1 = atom1.iac

			atomType1 = atomTypes[iac1]
			hrd1 = atomType1.hrd.cur
			eln1 = atomType1.eln.cur
			sig1 = atomType1.sig.cur

			X[i, i] = hrd1
			Y[i, 0] = -eln1

			for j in range(i + 1, len(cg)):
				idx2 = cg[j]
				atom2 = atoms[idx2]
				iac2 = atom2.iac

				atomType2 = atomTypes[iac2]
				sig2 = atomType2.sig.cur

				coul = self.coulomb(atom1, atom2, sig1, sig2)
				X[i, j] = coul
				X[j, i] = coul

		for i in range(len(cg)):
			X[nAtoms, i] = 1
			X[i, nAtoms] = -1

		res = np.linalg.solve(X, Y)
		EEM_ene = res[nAtoms, 0]  # TODO: save somewhere

		# print('{}\n{}\n\n'.format(X, Y))
		# print('res = \n{}'.format(res[:-1]))
		# sys.exit(1)

		for i in range(len(cg)):
			idx = cg[i]
			atom = atoms[idx]
			atom.charge.cur = res[i, 0]

		# use only 3 digits for charge
		qtot = 0.0
		for idx in cg:
			qtot += atoms[idx].charge.cur

		for idx in cg:
			atoms[idx].charge.cur -= qtot / nAtoms

		for idx in cg:
			atoms[idx].charge.cur = int(atoms[idx].charge.cur * 1000.0) / 1000.0

		qtot = 0.0
		for idx in cg:
			qtot += atoms[idx].charge.cur

		atoms[cg[0]].charge.cur -= qtot


	def coulomb(self, atom1, atom2, sig1, sig2):
		dist = 0.0

		# check if iac2 is first neighbor of iac1
		if atom1.isNrmNB(atom2.idx):
			dist = atom1.getNrmBndDist(atom2.idx)

		# otherwise use distance as second neighbor
		else:
			if atom1.isNeiNB(atom2.idx):
				dist = atom1.getNeiBndDist(atom2.idx)
			else:
				return 0.0

		sig_sum = sig1 + sig2
		if sig_sum == 0:
			return 0.0
		return math.erf(abs(dist) / (math.sqrt(sig1 * sig1 + sig2 * sig2))) * 1.0 / abs(dist)


class NONE(EEM):
	def setCG(self, mol):
		mol.CGs = []


class Halo(EEM):
	def setCG(self, mol):
		mol.CGs = []


class AA_Alk(EEM):
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


class O_N(EEM):
	def setCG(self, mol):
		CGs = []
		atoms = mol.atoms

		# add all atoms
		indexes = []
		for idx1 in atoms:
			indexes.append(idx1)

		CGs.append(indexes)
		mol.CGs = CGs

		atoms = mol.atoms
		for idx in atoms:
			atom = atoms[idx]
			if not atom.ignore:
				charge = atom.charge
				mol.addParameter(charge)
