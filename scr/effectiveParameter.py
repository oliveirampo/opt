"""Effective parameters for non-bonded interactions.

Classes:
--------
	EffectiveParameter
	ParameterType
	NRM
	NEI
	C6
	C12

Methods:
--------
	createEffectiveParameterFactory
"""


from abc import ABC, ABCMeta, abstractmethod
import math

import myExceptions


def createEffectiveParameterFactory(type, idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val):
	"""Creates effective parameter objects given variable 'type'.

	:param type:
	:param idx:
	:param type_c:
	:param type_14:
	:param iac1:
	:param iac2:
	:param typAtm1:
	:param typAtm2:
	:param val:
	:return: Selected EffectiveParameters based on variable 'type' given.
		type = 'LJ': returns LJ.
		type = 'Charge': returns Charge.
	:exception
		myExceptions.ClassNotImplemented
	"""

	class Charge(EffectiveParameter):
		"""LJ effective parameter.

		Methods:
		--------
			computeCR(cr, atomTypes, matrix):
				No combining rule for Charge object.
			writePrm(out):
				Writes parameter values to txt file.
		"""

		def	__init__(self, idx, typAtm, val):
			"""Constructs all the necessary attributes for the given LJ parameter.

			:param idx: (int) Index.
			:param typAtm: (str) Type/name of atom.
			:param val: (float) Value of parameter.
			"""

			self._idx = idx
			self._typAtm = typAtm
			self._ori = float(val)
			self._cur = float(val)
			self._typ = 'CHG_ATM'

			s = '{}_1_{}'.format(self._typ, self._idx)
			self._nam = s

		@property
		def ori(self):
			return self._ori

		@property
		def cur(self):
			return self._cur

		@property
		def nam(self):
			return self._nam

		@ori.setter
		def ori(self, value):
			self._ori = value

		@cur.setter
		def cur(self, value):
			self._cur = value

		def computeCR(self, cr, atomTypes, matrix):
			"""No combining rule for Charge object."""
			pass

		def writePrm(self, out):
			"""Writes parameter values to txt file.

			:param out: (output object) Output file.
			"""

			idx = self._idx
			typ = self._typ
			typAtm = self._typAtm
			val = self._cur

			out.write('{0:>4} {1:>7} {2:>3} {3:>3} {4:<5} {5:3} {6:>15.4f} {7:>13.4f}\n'
					  .format(idx, typ, '1', idx, 'MOLEC', typAtm, 0.0, val))

	class LJ(EffectiveParameter):
		"""LJ effective parameter.

		Methods:
		--------
			computeCR(cr, atomTypes, matrix):
				Computes combining rule.
			writePrm(out):
				Writes parameter values to txt file.
		"""

		def __init__(self, idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val):
			"""Constructs all the necessary attributes for the given LJ parameter.

			:param idx: (int) Index.
			:param type_c: (effectiveParameter.C6 or effectiveParameter.C12) C6 or C12 type.
			:param type_14: (effectiveParameter.NEI or effectiveParameter.NRM) NEI or NRM type.
			:param iac1: (int) Atom type of atom 1.
			:param iac2: (int) Atom type of atom 2.
			:param typAtm1: (str) Type/name of atom 1.
			:param typAtm2: (str) Type/name of atom 2.
			:param val: (float) Value of parameter.
			"""

			self._idx = int(idx)
			self._type_c = type_c
			self._type_14 = type_14
			self._iac1 = int(iac1)
			self._iac2 = int(iac2)
			self._typAtm1 = typAtm1
			self._typAtm2 = typAtm2
			self._ori = float(val)
			self._cur = float(val)
			self._used = False

			s = '{}_{}_{}_{}'.format(self.type_c.getType(), self.type_14.getType(), self.iac1, self._iac2)
			self._nam = s

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
		def ori(self):
			return self._ori

		@property
		def cur(self):
			return self._cur

		@property
		def used(self):
			return self._used

		@ori.setter
		def ori(self, value):
			self._ori = value

		@property
		def nam(self):
			return self._nam

		def computeCR(self, cr, atomTypes, matrix):
			"""Computes combining rule.

			:param cr: (CR) Combining rule.
			:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
			:param matrix: (Matrix) Matrix with usage of C12(II) parameters.
			"""

			iac1 = atomTypes[self.iac1]
			iac2 = atomTypes[self.iac2]

			sigi_nrm = iac1.sig.cur
			sigj_nrm = iac2.sig.cur

			epsi_nrm = iac1.eps.cur
			epsj_nrm = iac2.eps.cur

			if iac1.iac in matrix:
				if iac2.iac in matrix[iac1.iac]:
					sigi_nrm = iac1.sig_2.cur
					epsi_nrm = iac1.eps_2.cur

			if iac2.iac in matrix:
				if iac1.iac in matrix[iac2.iac]:
					sigj_nrm = iac2.sig_2.cur
					epsj_nrm = iac2.eps_2.cur

			sigi_nei = sigi_nrm
			epsi_nei = epsi_nrm
			sigj_nei = sigj_nrm
			epsj_nei = epsj_nrm
			if iac1.fixed_nei:
				sigi_nei = iac1.sig_nei
				epsi_nei = iac1.eps_nei
			if iac2.fixed_nei:
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

			val = self._type_c.getVal(c6, c12)
			self._cur = val

		def writePrm(self, out):
			"""Writes parameter values to txt file.

			:param out: (output object) Output file.
			"""

			idx = self._idx
			type_c = self._type_c.getType()
			type_14 = self._type_14.getType()
			iac1 = self._iac1
			iac2 = self._iac2
			typAtm1 = self._typAtm1
			typAtm2 = self._typAtm2
			val = self._cur

			prmName = '{}_{}'.format(type_c, type_14)
			out.write('{0:4} {1:>7} {2:>3} {3:>3} {4:<5} {5:5} {6:>13.4f} {7:13.6e}\n'
					  .format(idx, prmName, iac1, iac2, typAtm1, typAtm2, 0.0, val))

	if type == 'LJ':
		return LJ(idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val)
	elif type == 'Charge':
		return Charge(idx, typAtm1, val)
	else:
		raise myExceptions.ClassNotImplemented(type, 'effectiveParameter')


class EffectiveParameter(ABC):
	"""Base class for effective (non-bonded) parameters.

	Methods:
	--------
		computeCR(cr, atomTypes, matrix):
			Computes combining rule.
		writePrm(out):
			Writes parameter values to txt file.
	"""

	@abstractmethod
	def computeCR(self, cr, atomTypes, matrix):
		"""Computes combining rule."""
		pass

	@abstractmethod
	def writePrm(self, out):
		"""Writes parameter values to txt file."""
		pass


class ParameterType(metaclass=ABCMeta):
	"""Defines parameter type:
		- NRM or NEI
		- C6 or C12

	Methods:
	--------
		getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
			Returns value of parameter.
		getType(():
			Returns type of parameter.
	"""

	@classmethod
	def __subclasshook__(cls, subclass):
		return (hasattr(subclass, 'getVal') and
				callable(subclass.getVal) and
				hasattr(subclass, 'getType') and
				callable(subclass.getType) or
				NotImplemented
				)


class NRM(ParameterType):
	"""Parameter type = NRM

	Methods:
	--------
	getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		Returns value of parameter.
	getType(():
		Returns type of parameter.
	"""

	def getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		"""Return 'NRM' values of parameter.

		:param sig_nrm: (float) Sigma normal value.
		:param eps_nrm: (float) Epsilon normal value.
		:param sig_nei: (float) Sigma neighbor value.
		:param eps_nei: (float) Epsilon neighbor value.
		:return:
			sig_nrm: (float) Sigma normal value.
			eps_nrm: (float) Epsilon normal value.
		"""

		return sig_nrm, eps_nrm

	def getType(self):
		"""Returns type of parameter."""

		return 'NRM'


class NEI(ParameterType):
	"""Parameter type = NEI

	Methods:
	--------
	getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		Returns value of parameter.
	getType(():
		Returns type of parameter.
	"""

	def getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		"""Return 'NEI' values of parameter.

		:param sig_nrm: (float) Sigma normal value.
		:param eps_nrm: (float) Epsilon normal value.
		:param sig_nei: (float) Sigma neighbor value.
		:param eps_nei: (float) Epsilon neighbor value.
		:return:
			sig_nrm: (float) Sigma neighbor value.
			eps_nrm: (float) Epsilon neighbor value.
		"""

		return sig_nei, eps_nei

	def getType(self):
		"""Returns type of parameter."""

		return 'NEI'


class C6(ParameterType):
	"""Parameter type = C6

	Methods:
	--------
	getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		Returns value of parameter.
	getType(():
		Returns type of parameter.
	"""

	def getVal(self, c6, c12):
		"""Return 'c6' value of parameter.

		:param c6: (float) C6 value.
		:param c12: (float) C12 value.
		:return:
			c6: (float) C6 value.
		"""

		return c6

	def getType(self):
		"""Returns type of parameter."""

		return 'C06'


class C12(ParameterType):
	"""Parameter type = C12

	Methods:
	--------
	getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		Returns value of parameter.
	getType(():
		Returns type of parameter.
	"""

	def getVal(self, c6, c12):
		"""Return 'c12' value of parameter.

		:param c6: (float) C6 value.
		:param c12: (float) C12 value.
		:return:
			c12: (float) C12 value.
		"""

		return c12

	def getType(self):
		"""Returns type of parameter."""

		return 'C12'
