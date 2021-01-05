from abc import ABC, abstractmethod
from decimal import Decimal as dec
import numpy as np
import math

from scr import effectiveParameter


class EEM(ABC):
	def setCG(self, mol):
		mol.CGs = []


	@abstractmethod
	def setCG(self, mol):
		pass


	def computeEEM(self, atomTypes, cg, atoms):
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
			atom.curChg = dec(res[i, 0])

		# use only 3 digits for charge
		qtot = dec(0.0)
		for idx in cg:
			qtot += atoms[idx].curChg

		for idx in cg:
			atoms[idx].curChg -= qtot / nAtoms

		for idx in cg:
			atoms[idx].curChg = dec(int(atoms[idx].curChg * dec(1000.0)) / dec(1000.0))

		qtot = dec(0.0)
		for idx in cg:
			qtot += atoms[idx].curChg

		atoms[cg[0]].curChg -= qtot

		# for idx in cg:
		# 	print(idx, atoms[idx].curChg)


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
		for indexes in CGs:
			# cg = chgCG(indexes)
			cg = effectiveParameter.createEffectiveParameterFactory('CHG', indexes, '', '', '', '', '', '', '', '')
			mol.addParameter(cg)








