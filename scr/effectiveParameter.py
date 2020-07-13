from abc import ABC, ABCMeta, abstractmethod
import math


def createEffectiveParameterFactory(type, indexes, idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val):
	class chgCG(EffectiveParameter):
		def __init__(self, indexes):
			self._indexes = indexes

		def computeEEM(self, eem, atomTypes, atoms):
			eem.computeEEM(atomTypes, self._indexes, atoms)

		def computeCR(self, cr, atomTypes, matrix):
			pass

		def writePrm(self, out):
			pass

	class LJ(EffectiveParameter):
		def __init__(self, idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val):
			self._idx = int(idx)
			self._type_c = type_c
			self._type_14 = type_14
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
		def type_14(self):
			return self._type_14

		@property
		def type_c(self):
			return self._type_c

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

		def computeEEM(self, eem, atomTypes, atoms):
			pass

		def computeCR(self, cr, atomTypes, matrix):
			iac1 = atomTypes[self.iac1]
			iac2 = atomTypes[self.iac2]

			sigi_nrm = iac1.sig.cur
			sigj_nrm = iac2.sig.cur

			epsi_nrm = iac1.eps.cur
			epsj_nrm = iac2.eps.cur

			if iac1.iac in matrix:
				if iac2.iac in matrix[iac1.iac]:
					sigi_nrm = iac1.sig_2.cur
					epsj_nrm = iac1.eps_2.cur

			if iac2.iac in matrix:
				if iac1.iac in matrix[iac2.iac]:
					sigj_nrm = iac2.sig_2.cur
					epsj_nrm = iac2.eps_2.cur

			sigi_nei = sigi_nrm
			epsi_nei = epsi_nrm
			sigj_nei = sigj_nrm
			epsj_nei = epsj_nrm
			if iac1.fixed_nei == True:
				sigi_nei = iac1.sig_nei
				epsi_nei = iac1.eps_nei
			if iac2.fixed_nei == True:
				sigj_nei = iac2.sig_nei
				epsj_nei = iac2.eps_nei

			sigi, epsi = self.type_14.getVal(sigi_nrm, epsi_nrm, sigi_nei, epsi_nei)
			sigj, epsj = self.type_14.getVal(sigj_nrm, epsj_nrm, sigj_nei, epsj_nei)

			sigij = cr.getSigma(sigi, sigj)
			epsij = cr.getEpsilon(epsi, epsj, sigi, sigj)

			if sigij == 0.0:
				c6 = 0.0
				c12 = 0.0
			else:
				c6 = 4.0 * epsij * math.exp(6.0 * math.log(sigij))
				c12 = 4.0 * epsij * math.exp(12.0 * math.log(sigij))

			val = self.type_c.getVal(c6, c12)
			self.val = val

		def writePrm(self, out):
			idx = self.idx
			type_c = self.type_c.getType()
			type_14 = self.type_14.getType()
			iac1 = self.iac1
			iac2 = self.iac2
			typAtm1 = self.typAtm1
			typAtm2 = self.typAtm2
			val = self.val

			prmName = '{}_{}'.format(type_c, type_14)
			out.write('{0:4} {1:>7} {2:>3} {3:>3} {4:<5} {5:3} {6:>15.4f} {7:13.6e}\n'
					  .format(idx, prmName, iac1, iac2, typAtm1, typAtm2, 0.0, val))

	if type == 'CHG':
		return chgCG(indexes)
	if type == 'LJ':
		return LJ(idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val)
	assert 0, "Bad shape creation: " + type


class EffectiveParameter(ABC):
	@abstractmethod
	def computeEEM(self, eem, atomTypes, atoms):
		pass

	@abstractmethod
	def computeCR(self, cr, atomTypes, matrix):
		pass

	@abstractmethod
	def writePrm(self, out):
		pass


class LJ_14_Type(metaclass=ABCMeta):
	@classmethod
	def __subclasshook__(cls, subclass):
		return (hasattr(subclass, 'getVal') and
				callable(subclass.getVal) and
				hasattr(subclass, 'getType') and
				callable(subclass.getType) or
				NotImplemented
				)

	@abstractmethod
	def getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		raise NotImplementedError

	@abstractmethod
	def getType(self):
		raise NotImplementedError


class NRM(LJ_14_Type):
	def getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		return sig_nrm, eps_nrm

	def getType(self):
		return 'NRM'


class NEI(LJ_14_Type):
	def getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		return sig_nei, eps_nei

	def getType(self):
		return 'NEI'


class LJ_C_Type(metaclass=ABCMeta):
	@classmethod
	def __subclasshook__(cls, subclass):
		return (hasattr(subclass, 'getVal') and
				callable(subclass.getVal) and
				hasattr(subclass, 'getType') and
				callable(subclass.getType) or
				NotImplemented
				)

	@abstractmethod
	def getVal(self, c6, c12):
		raise NotImplementedError

	@abstractmethod
	def getType(self):
		raise NotImplementedError


class C6(LJ_C_Type):
	def getVal(self, c6, c12):
		return c6

	def getType(self):
		return 'C06'


class C12(LJ_C_Type):
	def getVal(self, c6, c12):
		return c12

	def getType(self):
		return 'C12'






