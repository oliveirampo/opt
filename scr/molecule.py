"""Module for molecule object.

Classes:
	Molecule
"""

import pandas as pd
import numpy as np
import sys

from effectiveParameter import C6
from effectiveParameter import C12
from effectiveParameter import NEI
from effectiveParameter import NRM
import parameter_utils
import effectiveParameter
from sensitivity import Sensitivity
from ChargeDistribution import BondChargeDistributionMethod


class Molecule:
	"""Defines a molecule object.

	Attributes:
		_cod: (str) Code.
		_frm: (str) Formula.
		_run: (bool) Flag to consider or not the molecule.
		_pre_sim: (float) Pressure of simulation.
		_tem_sim: (float) Temperature of simulation.
		_properties: (list) List of properties.
		_mlp_ref: (float) Melting point.
		_blp_ref: (float) Boiling point.
		_tem_cri: (float) Critical temperature.
		_eps_ref: (float) Permittivity.

		_atoms = {}
		_connectivity: (numpy.ndarray) Connectivity matrix of booleans.
		_distance_matrix: (numpy.ndarray) Distance matrix.
		_charge_transfer_topology: (numpy.ndarray)
		_bond_hardness: (numpy.ndarray)
		_CGs: (numpy.ndarray) Charge groups - list of list with indexes of atoms in the same charge group.
		_net_charge: (float) Net charge.
		_parameters: (list) List of parameters.

		_sens: (Sensitivity) Sensitivity matrix.

	Methods:
		get_bond_hardness(self, hardness)
		createChargeTransferTopology()
		addAtom(atom, conf)
		getAtom(self, idx)
		checkAtoms()
		nAtoms()
		are_bonded( idx1, idx2)
		get_bond_distance(idx1, idx2)
		add_bond_distance(idx1, idx2, distance)
		get_neighbors_of_atom(idx1)
		createEffectiveAtomicCharges()
		createLJPairs(atomTypes)
		writePrmMod(f)
		addTrajectory(letter, traj)
		addRunningAverages(propCode, avgs)
		getNumProps()
	"""

	def __init__(self, cod, frm, run, pre_sim, tem_sim, properties, mlp_ref, blp_ref, tem_cri, eps_ref):
		"""Constructs all the necessary attributes for the molecule.

		:param cod: (str) Code.
		:param frm: (str) Formula.
		:param run: (bool) Flag to consider or not the molecule.
		:param pre_sim: (float) Pressure of simulation.
		:param tem_sim: (float) Temperature of simulation.
		:param properties: (list) List of properties.
		:param mlp_ref: (float) Melting point.
		:param blp_ref: (float) Boiling point.
		:param tem_cri: (float) Critical temperature.
		:param eps_ref: (float) Permittivity.
		"""

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
		self._tem_cri = float(tem_cri)
		self._eps_ref = float(eps_ref)

		self._atoms = {}
		self._connectivity = np.zeros([0, 0])
		self._distance_matrix = np.zeros([0, 0])
		self._charge_transfer_topology = np.zeros([0, 0])
		self._bond_hardness = np.zeros(0)
		self._CGs = []
		self._net_charge = 0.0
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
	def tem_cri(self):
		return self._tem_cri

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
	def charge_transfer_topology(self):
		return self._charge_transfer_topology

	@property
	def parameters(self):
		return self._parameters

	@property
	def CGs(self):
		return self._CGs

	@property
	def net_charge(self):
		return self._net_charge

	@property
	def sens(self):
		return self._sens

	@CGs.setter
	def CGs(self, list_of_lists):
		self._CGs = np.array([np.array(arr) for arr in list_of_lists], dtype=object)

	@run.setter
	def run(self, n):
		n = float(n)
		if n == 1:
			self._run = True
		else:
			self._run = False

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

	def get_bond_hardness(self, hardness):
		"""Returns bond hardness.
		Assume bond hardnesses = sum of atomic hardnesses

		:param hardness: (ndarray) List of hardness.
		:return:
			bondHardness: (ndarray) bond hardness.
		"""

		charge_transfer_topology = self._charge_transfer_topology

		bVars, B = BondChargeDistributionMethod.bondVars(charge_transfer_topology)
		bondHardness = np.zeros(B)

		for b, [i, j] in enumerate(bVars):
			bondHardness[b] = (hardness[i] + hardness[j]) / 2

		return bondHardness

	def createChargeTransferTopology(self):
		"""Creates charge transfer topology."""

		chargeTransferFilter = lambda x: 1 if x else 0
		chargeTransferTopology = np.vectorize(chargeTransferFilter)(self.connectivity)
		self._charge_transfer_topology = chargeTransferTopology

	def addAtom(self, atom, conf):
		"""Adds atom to dictionary of atoms of molecule.

		:param atom: (Atom)
		:param conf: (configuration.Conf) Configuration object.
		"""

		try:
			idx = atom.idx
		except KeyError:
			sys.exit('\tERROR: Atom has no idx')

		if conf.ignoreAtom(atom.iac):
			atom.ignore = True

		self._atoms[idx] = atom

	def getAtom(self, idx):
		"""Returns atom by index.

		:param idx: (idx) Index of atom in molecule.
		:return:
			atom: (Atom)
		"""

		return self._atoms[idx]

	def checkAtoms(self):
		"""Checks if dictionary of atoms is not empty. Stops pipeline if this is the case."""

		if len(self._atoms) == 0:
			sys.exit('No atoms were found for {}'.format(self._cod))

	def nAtoms(self):
		"""Returns number of atoms in molecule."""

		return len(self._atoms)

	def are_bonded(self, idx1, idx2):
		"""Checks if two atoms are connected given their indexes.

		:param idx1: (int) Index of atom 1.
		:param idx2: (int) Index of atom 2.
		:return:
			:(boolean) True or False.
		"""

		pos1 = idx1 - 1
		pos2 = idx2 - 1
		return self._connectivity[pos1, pos2]

	def get_bond_distance(self, idx1, idx2):
		"""Returns bond distance between two atoms given their indexes.

		:param idx1: (int) Index of atom 1.
		:param idx2: (int) Index of atom 2.
		:return:
			distance: (float) Distance between two atoms.
		"""

		pos1 = idx1 - 1
		pos2 = idx2 - 1
		return self._distance_matrix[pos1, pos2]

	def add_bond_distance(self, idx1, idx2, distance):
		"""Adds bond distance between two atoms to distance matrix.

		:param idx1: (int) Index of atom 1.
		:param idx2: (int) Index of atom 2.
		:param distance: (numpy.ndarray) Distance matrix.
		"""

		pos1 = idx1 - 1
		pos2 = idx2 - 1
		self._distance_matrix[pos1, pos2] = distance
		self._distance_matrix[pos2, pos1] = distance

	def get_neighbors_of_atom(self, idx1):
		"""Returns list with indexes of neighbors.

		:param idx1: (int) Index of atom 1.
		:return:
			neighbors_indexes: (numpy.ndarray) List of indexes of neighbors.
		"""

		pos1 = idx1 - 1
		neighbors_indexes = np.where(self._connectivity[pos1, :])
		neighbors_indexes = neighbors_indexes[0] + 1
		return neighbors_indexes

	def createEffectiveAtomicCharges(self):
		"""Creates (Parameter) charges, and adds them to list of parameters of given molecule."""

		atoms = self._atoms
		for idx in atoms:
			atom = atoms[idx]
			charge = atom.charge
			self._parameters.append(charge)

	def createLJPairs(self, atomTypes):
		"""Creates (Parameter) LJ pairs, and add them to list of parameters of given molecule.

		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		"""

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

	def writePrmMod(self, f):
		"""Writes param.mod file.

		:param f: (str) Output file name.
		:return:
		"""

		for prm in self._parameters:
			prm.writePrm(f)

	def addTrajectory(self, letter, traj):
		"""Extracts instantaneous values of simulation results and adds them to property.

		:param letter: (str) Property letter.
		:param traj: (ndarray) Instantaneous values of simulation results.
		"""

		for prop in self._properties:
			if letter == prop.letter:
				prop.trajectory = traj

	def addRunningAverages(self, propCode, avgs):
		"""Extracts running averages of simulation results and adds them to property.

		:param propCode: (str) Code of property.
		:param avgs: (list) Running averages of simulation results.
		"""

		for prop in self._properties:
			if propCode == prop.code:
				prop.runningAverages = avgs

	def getNumProps(self):
		"""Returns number of properties."""

		return len(self._properties)
