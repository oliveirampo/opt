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
import sys

import myExceptions


def createEffectiveParameterFactory(typ, idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val):
	"""Creates effective parameter objects given variable 'type'.

	:param typ:
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

		def __init__(self, idx_value, typAtm, value):
			"""Constructs all the necessary attributes for the given LJ parameter.

			:param idx_value: (int) Index.
			:param typAtm: (str) Type/name of atom.
			:param value: (float) Value of parameter.
			"""

			self._idx = idx_value
			self._typAtm = typAtm
			self._ori = float(value)
			self._cur = float(value)
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

		def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, atomTypes, matrix):
			"""No combining rule for Charge object."""
			pass

		def writePrm(self, out):
			"""Writes parameter values to txt file.

			:param out: (output object) Output file.
			"""

			idx_value = self._idx
			typ_value = self._typ
			typAtm = self._typAtm
			value = self._cur

			out.write('{0:>4} {1:>7} {2:>3} {3:>3} {4:<5} {5:3} {6:>15.4f} {7:>13.4f}\n'.format(idx_value, typ_value,
																							'1', idx_value, 'MOLEC',
																							typAtm, 0.0, value))

	class LJ(EffectiveParameter):
		"""LJ effective parameter.

		Methods:
		--------
			computeCR(cr, crPrms, atomTypes, matrix):
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

		def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, atomTypes, matrix):
			"""Computes combining rule.

			:param cr: (CR) Combining rule.
			:param crPrms: (dict) Parameters for linear combination of combining rules.
			:param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
			:param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
			:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
			:param matrix: (Matrix) Matrix with usage of C12(II) parameters.
			"""

			iac1 = atomTypes[self.iac1]
			iac2 = atomTypes[self.iac2]

			prm_type = self._type_14  # NRM or NEI
			c6, c12 = prm_type.computeCR(cr, crPrms, scl_sig_NEI, scl_eps_NEI, iac1, iac2, matrix)

			val = self._type_c.getVal(c6, c12, '', '')
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

	if typ == 'LJ':
		return LJ(idx, type_c, type_14, iac1, iac2, typAtm1, typAtm2, val)
	elif typ == 'Charge':
		return Charge(idx, typAtm1, val)
	else:
		raise myExceptions.ClassNotImplemented(typ, 'effectiveParameter')


class EffectiveParameter(ABC):
	"""Base class for effective (non-bonded) parameters.

	Methods:
	--------
		computeCR(cr, crPrms, scl_sig_NEI, scl_eps_NEI, atomTypes, matrix):
			Computes combining rule.
		writePrm(out):
			Writes parameter values to txt file.
	"""

	@abstractmethod
	def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, atomTypes, matrix):
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

	@abstractmethod
	def getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		pass

	@abstractmethod
	def getType(self):
		pass

	@abstractmethod
	def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, iac1, iac2, matrix):
		pass


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

	def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, iac1, iac2, matrix):
		"""Computes combining rule.

		:param cr: (CR) Combining rule.
		:param crPrms: (dict) Parameters for linear combination of combining rules.
		:param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
		:param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
		:param iac1: (IAC) Atom type 1.
		:param iac2: (IAC) Atom type 2.
		:param matrix: (Matrix) Matrix with usage of C12(II) parameters.
		"""

		alpha = crPrms['alpha']['val']

		sigi = iac1.sig.cur
		sigj = iac2.sig.cur
		epsi = iac1.eps.cur
		epsj = iac2.eps.cur

		sigi_2 = sigi
		sigj_2 = sigj
		epsi_2 = epsi
		epsj_2 = epsj

		use_sig_eps_2 = False
		if iac1.iac in matrix:
			if iac2.iac in matrix[iac1.iac]:
				sigi_2 = iac1.sig_2.cur
				epsi_2 = iac1.eps_2.cur
				use_sig_eps_2 = True

		if iac2.iac in matrix:
			if iac1.iac in matrix[iac2.iac]:
				sigj_2 = iac2.sig_2.cur
				epsj_2 = iac2.eps_2.cur
				use_sig_eps_2 = True

		sigi_2_rng = iac1.sig_2.rng
		sigj_2_rng = iac2.sig_2.rng
		epsi_2_rng = iac1.eps_2.rng
		epsj_2_rng = iac2.eps_2.rng

		if use_sig_eps_2:
			sigij = cr.getSigma(sigi, sigj, alpha)
			epsij = cr.getEpsilon(epsi, epsj, sigi, sigj, alpha)

			if (epsi_2_rng != 0) or (epsj_2_rng != 0):
				epsij_2 = cr.getEpsilon(epsi_2, epsj_2, sigi_2, sigj_2, alpha)

				c6 = 4.0 * epsij * math.exp(6.0*math.log(sigij))

				sigij = math.exp((1.0 / 6.0) * math.log(c6 / (4 * epsij_2)))
				epsij = epsij_2

			elif (sigi_2_rng != 0) or (sigj_2_rng != 0):
				sigij_2 = cr.getSigma(sigi_2, sigj_2, alpha)

				c6 = 4.0 * epsij * math.exp(6.0*math.log(sigij))

				sigij = sigij_2
				epsij = c6 / (4 * math.exp(6 * math.log(sigij_2)))

			else:
				sys.exit('computeCR()::TODO')

		else:
			sigij = cr.getSigma(sigi, sigj, alpha)
			epsij = cr.getEpsilon(epsi, epsj, sigi, sigj, alpha)

		if sigij == 0.0:
			c6 = 0.0
			c12 = 0.0
		else:
			c6 = 4.0 * epsij * math.exp(6.0 * math.log(sigij))
			c12 = 4.0 * epsij * math.exp(12.0 * math.log(sigij))

		return c6, c12


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

	def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, iac1, iac2, matrix):
		"""Computes combining rule.

		:param cr: (CR) Combining rule.
		:param crPrms: (dict) Parameters for linear combination of combining rules.
		:param scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
		:param scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
		:param iac1: (IAC) Atom type 1.
		:param iac2: (IAC) Atom type 2.
		:param matrix: (Matrix) Matrix with usage of C12(II) parameters.
		"""

		sigi_nrm = iac1.sig.cur
		sigj_nrm = iac2.sig.cur
		epsi_nrm = iac1.eps.cur
		epsj_nrm = iac2.eps.cur

		sigi_nei = sigi_nrm * scl_sig_NEI
		sigj_nei = sigj_nrm * scl_sig_NEI
		epsi_nei = epsi_nrm * scl_eps_NEI
		epsj_nei = epsj_nrm * scl_eps_NEI

		if iac1.fixed_nei:
			sigi_nei = iac1.sig_nei
			epsi_nei = iac1.eps_nei
		if iac2.fixed_nei:
			sigj_nei = iac2.sig_nei
			epsj_nei = iac2.eps_nei

		alpha = crPrms['alpha']['val']

		sigij = cr.getSigma(sigi_nei, sigj_nei, alpha)
		epsij = cr.getEpsilon(epsi_nei, epsj_nei, sigi_nei, sigj_nei, alpha)

		c6 = 0.0
		c12 = 0.0
		if sigij != 0.0:
			c6 = 4.0 * epsij * math.exp(6.0 * math.log(sigij))
			c12 = 4.0 * epsij * math.exp(12.0 * math.log(sigij))

		return c6, c12


class C6(ParameterType):
	"""Parameter type = C6

	Methods:
	--------
	getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		Returns value of parameter.
	getType(():
		Returns type of parameter.
	"""

	def getVal(self, c6, c12, dummy_1, dummy_2):
		"""Return 'c6' value of parameter.

		:param c6: (float) C6 value.
		:param c12: (float) C12 value.
		:param dummy_1: (float) Dummy variable (not used).
		:param dummy_2: (float) Dummy variable (not used).
		:return:
			c6: (float) C6 value.
		"""

		return c6

	def getType(self):
		"""Returns type of parameter."""

		return 'C06'

	def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, iac1, iac2, matrix):
		raise myExceptions.ClassNotImplemented('C12', 'computeCR()')


class C12(ParameterType):
	"""Parameter type = C12

	Methods:
	--------
	getVal(self, sig_nrm, eps_nrm, sig_nei, eps_nei):
		Returns value of parameter.
	getType(():
		Returns type of parameter.
	"""

	def getVal(self, c6, c12, dummy_1, dummy_2):
		"""Return 'c12' value of parameter.

		:param c6: (float) C6 value.
		:param c12: (float) C12 value.
		:param dummy_1: (float) Dummy variable (not used).
		:param dummy_2: (float) Dummy variable (not used).
		:return:
			c12: (float) C12 value.
		"""

		return c12

	def getType(self):
		"""Returns type of parameter."""

		return 'C12'

	def computeCR(self, cr, crPrms, scl_sig_NEI, scl_eps_NEI, iac1, iac2, matrix):
		raise myExceptions.ClassNotImplemented('C12', 'computeCR()')
