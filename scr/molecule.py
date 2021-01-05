import pandas as pd
import sys

from scr.effectiveParameter import C6
from scr.effectiveParameter import C12
from scr.effectiveParameter import NEI
from scr.effectiveParameter import NRM
from scr import parameter_utils, effectiveParameter
from scr.sensitivity import Sensitivity


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
		idx = int(idx)
		return self._atoms[idx]

	def checkAtoms(self):
		if len(self._atoms) == 0:
			sys.exit('No atoms were found for {}'.format(self._cod))

	def nAtoms(self):
		return len(self._atoms)

	def createEffectivePrms(self, atomTypes, eem, matrix):
		self.createLJPairs(atomTypes, eem, matrix)

	def createLJPairs(self, atomTypes, eem, matrix):
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

			c6_nrm = effectiveParameter.createEffectiveParameterFactory('LJ', [], idx, C6(), NRM(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c6_nrm)
			idx += 1
			c6_nei = effectiveParameter.createEffectiveParameterFactory('LJ', [], idx, C6(), NEI(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c6_nei)
			idx += 1
			c12_nrm = effectiveParameter.createEffectiveParameterFactory('LJ', [], idx, C12(), NRM(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c12_nrm)
			idx += 1
			c12_nei = effectiveParameter.createEffectiveParameterFactory('LJ', [], idx, C12(), NEI(), iac1, iac2, typ1, typ2, 0.0)
			self._parameters.append(c12_nei)
			idx += 1

	def addParameter(self, prm):
		self._parameters.append(prm)

	def writePrmMod(self, f):
		# write charges
		for idx in self._atoms:
			atom = self._atoms[idx]
			# print(atom.curChg)
			f.write('{0:4} {1:>7} {2:>3} {3:>3} {4:>5} {5:3} {6:>15.4f} {7:13.4f}\n'
			.format(idx, 'CHG_ATM', 1, idx, 'MOLEC', atom.nam, 0.0, atom.curChg))

		# write LJ
		for prm in self._parameters:
			prm.writePrm(f)

	def addTrajectory(self, letter, traj):
		for prop in self._properties:
			if letter == prop.letter:
				prop.trajectory = traj

	def addRunningAverages(self, letter, avgs):
		for prop in self._properties:
			if letter == prop.letter:
				prop.runningAverages = avgs

	def getNumProps(self):
		return len(self._properties)
