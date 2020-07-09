from abc import ABC, abstractmethod
import math
import sys


class EffectiveParameter(ABC):
	def __init__(self, idx, typ, category, iac1, iac2, typAtm1, typAtm2, val):
		self._idx = int(idx)
		self._typ = typ
		self._category = category
		self._iac1 = int(iac1)
		self._iac2 = int(iac2)
		self._typAtm1 = typAtm1
		self._typAtm2 = typAtm2
		self._val = val
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
	def curVal(self):
		return self._val
	@property
	def used(self):
		return self._used

	@abstractmethod
	def computeEEM(self, eem, atomTypes, atoms):
		pass

	@abstractmethod
	def computeCR(self, cr, atomTypes, matrix):
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

	def computeCR(self, cr, atomTypes, matrix):
		pass

	def writePrm(self, out):
		pass


class LJ(EffectiveParameter):
	def __init__(self, idx, typ, category, iac1, iac2, typAtm1, typAtm2, val):
		EffectiveParameter.__init__(self, idx, typ, category, iac1, iac2, typAtm1, typAtm2, val)
		# C6 or C12 value
		self.LJnrm = 0.0
		self.LJnei = 0.0

	@abstractmethod
	def add_c_value(self, c6, c12):
		pass

	def computeEEM(self, eem, atomTypes, atoms):
		pass

	def computeCR(self, cr, atomTypes, matrix):
		iac1 = atomTypes[self.iac1]
		iac2 = atomTypes[self.iac2]

		sigi = iac1.sig.cur
		sigj = iac2.sig.cur

		epsi = iac1.eps.cur
		epsj = iac2.eps.cur

		if iac1.iac in matrix:
			if iac2.iac in matrix[iac1.iac]:
				sigi = iac1.sig_2.cur
				epsj = iac1.eps_2.cur

		if iac2.iac in matrix:
			if iac1.iac in matrix[iac2.iac]:
				sigj = iac2.sig_2.cur
				epsj = iac2.eps_2.cur

		sigij = cr.getSigma(sigi, sigj)
		epsij = cr.getEpsilon(epsi, epsj, sigi, sigj)

		if sigij == 0.0:
			c6 = 0.0
			c12 = 0.0
		else:
			c6  = 4.0 * epsij * math.exp( 6.0 * math.log(sigij))
			c12 = 4.0 * epsij * math.exp(12.0 * math.log(sigij))

		sigi_nei = sigi
		epsi_nei = epsi
		sigj_nei = sigj
		epsj_nei = epsj
		if iac1.fixed_nei == True:
			sigi_nei = iac1.sig_nei
			epsi_nei = iac1.eps_nei
		if iac2.fixed_nei == True:
			sigj_nei = iac2.sig_nei
			epsj_nei = iac2.eps_nei

		sigij_nei = cr.getSigma(sigi_nei, sigj_nei)
		epsij_nei = cr.getEpsilon(epsi_nei, epsj_nei, sigi_nei, sigj_nei)

		if sigij_nei == 0.0:
			c6_nei = 0.0
			c12_nei = 0.0
		else:
			c6_nei = 4.0 * epsij_nei * math.exp(6.0 * math.log(sigij_nei))
			c12_nei = 4.0 * epsij_nei * math.exp(12.0 * math.log(sigij_nei))

		if self.category == 'NRM':
			self.add_c_value(c6, c12)
		elif self.category == 'NEI':
			self.add_c_value(c6_nei, c12_nei)


	def writePrm(self, out):
		idx = self.idx
		typ = self.typ
		category = self.category
		iac1 = self.iac1
		iac2 = self.iac2
		typAtm1 = self.typAtm1
		typAtm2 = self.typAtm2
		val = self.LJnrm

		prmName = '{}_{}'.format(typ, category)
		out.write('{0:4} {1:>7} {2:>3} {3:>3} {4:<5} {5:3} {6:>15.4f} {7:13.6e}\n'
		.format(idx, prmName, iac1, iac2, typAtm1, typAtm2, 0.0, val))


class C6(LJ):
	def __init__(self, idx, category, iac1, iac2, typAtm1, typAtm2, val):
		LJ.__init__(self, idx, 'C06', category, iac1, iac2, typAtm1, typAtm2, val)

	def add_c_value(self, c6, c12):
		self.LJnrm = c6
		# self.LJnei = c6_nei


class C12(LJ):
	def __init__(self, idx, category, iac1, iac2, typAtm1, typAtm2, val):
		LJ.__init__(self, idx, 'C12', category, iac1, iac2, typAtm1, typAtm2, val)

	def add_c_value(self, c6, c12):
		self.LJnrm = c12
		# self.LJnei = c12_nei





