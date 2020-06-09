from abc import ABC, abstractmethod


class EffectiveParameter(ABC):
	def __init__(self, idx, typ, category, iac1, iac2, typAtm1, typAtm2, val):
		self._idx = int(idx)
		self._typ = typ
		self._category = category
		self._iac1 = int(iac1)
		self._iac2 = int(iac2)
		self._typAtm1 = typAtm1
		self._typAtm2 = typAtm2
		self._oriVal = val
		self._prevVal = val
		self._curVal = val
		self._used = False

	@property
	def idx(self):
		return self._idx
	@property
	def category(self):
		return self._category
	@property
	def typ(self):
		return self._typ
	@property
	def iac1(self):
		return self._iac1
	@property
	def iac2(self):
		return self._iac2
	@property
	def typAtm1(self):
		return self._typAtm1
	@property
	def typAtm2(self):
		return self._typAtm2
	@property
	def oriVal(self):
		return self._oriVal
	@property
	def prevVal(self):
		return self._prevVal
	@property
	def curVal(self):
		return self._curVal
	@property
	def used(self):
		return self._used

	@abstractmethod
	def computeEEM(self, eem, atomTypes, atoms):
		pass

	@abstractmethod
	def writePrm(self, out):
		pass


class chgCG(EffectiveParameter):
	def __init__(self, indexes):
		EffectiveParameter.__init__(self, 1, '-', '', 1, 1, '-', '-', '-')
		self._indexes = indexes

	def computeEEM(self, eem, atomTypes, atoms):
		eem.computeEEM(atomTypes, self._indexes, atoms)

	def writePrm(self, out):
		pass


class LJ(EffectiveParameter):
	def __init__(self, idx, typ, category, iac1, iac2, typAtm1, typAtm2, val):
		EffectiveParameter.__init__(self, idx, typ, category, iac1, iac2, typAtm1, typAtm2, val)

	def computeEEM(self, eem, atomTypes, atoms):
		pass

	def writePrm(self, out):
		idx = self.idx
		typ = self.typ
		category = self.category
		iac1 = self.iac1
		iac2 = self.iac2
		typAtm1 = self.typAtm1
		typAtm2 = self.typAtm2
		val = self.curVal

		prmName = '{}_{}'.format(typ, category)
		out.write('{0:4} {1:>7} {2:>3} {3:>3} {4:<5} {5:3} {6:>15.4f} {7:13.4f}\n'
		.format(idx, prmName, iac1, iac2, typAtm1, typAtm2, 0.0, val))


class C6(LJ):
	def __init__(self, idx, category, iac1, iac2, typAtm1, typAtm2, val):
		LJ.__init__(self, idx, 'C06', category, iac1, iac2, typAtm1, typAtm2, val)


class C12(LJ):
	def __init__(self, idx, category, iac1, iac2, typAtm1, typAtm2, val):
		LJ.__init__(self, idx, 'C12', category, iac1, iac2, typAtm1, typAtm2, val)





