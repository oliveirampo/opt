import pandas as pd
import numpy as np
import sys

from effectiveParameter import C6
from effectiveParameter import C12
from effectiveParameter import NEI
from effectiveParameter import NRM
import parameter_utils, effectiveParameter
from sensitivity import Sensitivity


class Molecule:
	def __init__(self, cod, frm, run, pre_sim, tem_sim, properties, mlp_ref, blp_ref, eps_ref):
		self._cod = cod
		self._frm = frm

		self._run = False
		run = float(run)
		if run == 1:
			self._run = True

		self._pre_sim = float(pre_sim)
		self._tem_sim = float(tem_sim)

		self._properties = properties

		self._mlp_ref = float(mlp_ref)
		self._blp_ref = float(blp_ref)
		self._eps_ref = float(eps_ref)

		self._atoms = {}
		self._connectivity = np.empty([0, 0])
		self._distance_matrix = np.empty([0, 0])
		self._CGs = []
		self._parameters = []

		self._sens = Sensitivity(pd.DataFrame())

	@property
	def cod(self):
		return self._cod

	@property
	def frm(self):
		return self._frm

	@property
	def run(self):
		return self._run

	@property
	def pre_sim(self):
		return self._pre_sim

	@property
	def tem_sim(self):
		return self._tem_sim

	@property
	def properties(self):
		return self._properties

	@property
	def mlp_ref(self):
		return self._mlp_ref

	@property
	def blp_ref(self):
		return self._blp_ref

	@property
	def eps_ref(self):
		return self._eps_ref

	@property
	def atoms(self):
		return self._atoms

	@property
	def connectivity(self):
		return self._connectivity

	@property
	def distance_matrix(self):
		return self._distance_matrix

	@property
	def parameters(self):
		return self._parameters

	@property
	def CGs(self):
		return self._CGs

	@property
	def sens(self):
		return self._sens

	@CGs.setter
	def CGs(self, n):
		self._CGs = n

	@run.setter
	def run(self, n):
		n = float(n)
		if n == 1:
			self._run = True

	@connectivity.setter
	def connectivity(self, matrix):
		if matrix.shape[0] != len(self._atoms):
			print('Inconsistency with number of atoms in connectivity matrix.')
			sys.exit(123)
		self._connectivity = matrix[:]

	@distance_matrix.setter
	def distance_matrix(self, matrix):
		if matrix.shape[0] != len(self._atoms):
			print('Inconsistency with number of atoms in distance matrix.')
			sys.exit(123)
		self._distance_matrix = matrix[:]

	@sens.setter
	def sens(self, df):
		self._sens = df

	def addAtom(self, atom, conf):
		try:
			idx = atom.idx
		except KeyError:
			sys.exit('\tERROR: Atom has no idx')

		if conf.ignoreAtom(atom.iac):
			atom.ignore = True

		self._atoms[idx] = atom

	def getAtom(self, idx):
		return self._atoms[idx]

	def checkAtoms(self):
		if len(self._atoms) == 0:
			sys.exit('No atoms were found for {}'.format(self._cod))

	def nAtoms(self):
		return len(self._atoms)

	def are_bonded(self, idx1, idx2):
		pos1 = idx1 - 1
		pos2 = idx2 - 1
		return self._connectivity[pos1, pos2]

	def get_bond_distance(self, idx1, idx2):
		pos1 = idx1 - 1
		pos2 = idx2 - 1
		return self._distance_matrix[pos1, pos2]

	def add_bond_distance(self, idx1, idx2, distance):
		pos1 = idx1 - 1
		pos2 = idx2 - 1
		self._distance_matrix[pos1, pos2] = distance
		self._distance_matrix[pos2, pos1] = distance

	def get_neighbors_of_atom(self, idx1):
		pos1 = idx1 - 1
		neighbors_indexes = np.where(self._connectivity[pos1, :])
		neighbors_indexes = neighbors_indexes[0] + 1
		return neighbors_indexes

	def createLJPairs(self, atomTypes):
		iacList = []
		for idx in self._atoms:
			iac = self._atoms[idx].iac
			if iac not in iacList:
				iacList.append(iac)

		iacList = sorted(iacList)
		pairs = []
		N = len(iacList)
		for i in range(N):
			for j in range(i, N):
				iac1 = iacList[i]
				iac2 = iacList[j]
				pairs.append([iac1, iac2])

		nAtoms = self.nAtoms()
		idx = nAtoms + 1

		for i in range(len(pairs)):
			pair = pairs[i]
			iac1 = pair[0]
			iac2 = pair[1]

			typ1 = parameter_utils.getType(iac1, atomTypes)
			typ2 = parameter_utils.getType(iac2, atomTypes)

			c6_nrm = effectiveParameter.createEffectiveParameterFactory('LJ', idx, C6(), NRM(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c6_nrm)
			idx += 1
			c6_nei = effectiveParameter.createEffectiveParameterFactory('LJ', idx, C6(), NEI(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c6_nei)
			idx += 1
			c12_nrm = effectiveParameter.createEffectiveParameterFactory('LJ', idx, C12(), NRM(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c12_nrm)
			idx += 1
			c12_nei = effectiveParameter.createEffectiveParameterFactory('LJ', idx, C12(), NEI(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c12_nei)
			idx += 1

	def addParameter(self, prm):
		self._parameters.append(prm)

	def writePrmMod(self, f):
		for prm in self._parameters:
			prm.writePrm(f)

	def addTrajectory(self, letter, traj):
		for prop in self._properties:
			if letter == prop.letter:
				prop.trajectory = traj

	def addRunningAverages(self, propCode, avgs):
		for prop in self._properties:
			if propCode == prop.code:
				prop.runningAverages = avgs

	def getNumProps(self):
		return len(self._properties)
