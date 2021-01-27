from abc import ABC, abstractmethod
import numpy as np
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

		# print(X)
		# print(Y)
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

			if not atm1.ignore:

				iac1 = atm1.iac
				pos = -1
				# loop over list of lists
				for i in range(len(indexes)):
					# loop over current list
					for j in range(len(indexes[i])):

						idx2 = indexes[i][j]
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

						elif atm1.isNrmNB(idx2):
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
							nbs = atm1.bnd_nrm
							for idx_nb in nbs:
								atom_nb = atoms[idx_nb]
								iac_nb = atom_nb.iac

								for x in range(len(indexes)):
									for y in range(len(indexes)):
										other_idx_tmp = indexes[x][y]

										if idx_nb == other_idx_tmp:
											continue

										elif atom_nb.ignore:
											break

										elif iac_nb in self.C_EST:
											break

										elif iac_nb in self.CC_EST:
											break

										elif atm1.isNrmNB(other_idx_tmp):
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

		atoms = mol.atoms
		for idx in atoms:
			atom = atoms[idx]
			# print(atom.idx, atom.iac, atom.nam)
			# if not atom.ignore:
			charge = atom.charge
			mol.addParameter(charge)
