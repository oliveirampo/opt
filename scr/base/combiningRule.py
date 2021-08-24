from abc import ABC, abstractmethod
import math

from scr.base import myExceptions


class CR(ABC):
	"""Base class for combining rules.

	Methods:
	--------
		getSigma(self, si, sj): Returns sigma(i,j).
		getEpsilon(self, ei, ej, si, sj): Returns epsilon(i,j).
		get_object(n): Returns object of class that implements base class CR.
	"""

	@abstractmethod
	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).
		
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		"""
		pass

	@abstractmethod
	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		"""
		pass

	@staticmethod
	def get_object(n):
		"""Returns object of class that implements base class CR.

		:param n: (str) Code for class that implements base class CR.
		:return: One of the classes that implements CR.
		"""
		classes = {"GEOM": GeometricCR, "WH": WaldmanHagler, "LB": LorentzBerthelot, "FH": FenderHalsey,
				   "GEOM_WH": Geometric_WaldmanHagler}

		if n not in classes:
			raise myExceptions.ClassNotImplemented(n, 'combiningRule')

		cr = classes[n]()
		return cr


class GeometricCR(CR):
	"""Implements geometric combining rule."""

	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).

		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			sigma(i,j): (float) Collision diameter of atoms i and j.
		"""

		return math.sqrt(si * sj)

	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			epsilon(i, j): (float) Well depth of atoms i and j.
		"""

		return math.sqrt(ei * ej)


class WaldmanHagler(CR):
	"""Implements Waldman-Hagler combining rule."""

	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).

		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			sigma(i,j): (float) Collision diameter of atoms i and j.
		"""

		return math.pow(0.5 * (math.pow(si, 6) + math.pow(sj, 6)), 1.0/6.0)

	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			epsilon(i, j): (float) Well depth of atoms i and j.
		"""

		if (si + sj) == 0.0:
			return 0.0

		return 2.0 * math.sqrt(ei*ej) * (si**3) * (sj**3) / ((si**6) + (sj**6))


class LorentzBerthelot(CR):
	"""Implements Lorentz-Berthelot combining rule."""

	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).

		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			sigma(i,j): (float) Collision diameter of atoms i and j.
		"""

		return (si + sj) / 2.0

	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			epsilon(i, j): (float) Well depth of atoms i and j.
		"""

		return math.sqrt(ei * ej)


class FenderHalsey(CR):
	"""Implements Fender-Halsey combining rule."""

	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).

		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			sigma(i,j): (float) Collision diameter of atoms i and j.
		"""

		return (si + sj) / 2

	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			epsilon(i, j): (float) Well depth of atoms i and j.
		"""

		if (ei * ej) == 0.0:
			return 0.0

		return (2 * ei * ej) / (ei * ej)


class Geometric_LorentzBerthelot_CR(CR):
	"""Implements linear combination of geometric and LB combining rules."""

	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).

		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			sigma(i,j): (float) Collision diameter of atoms i and j.
		"""

		geom = math.sqrt(si * sj)
		LB = (si + sj) / 2.0
		return alpha * geom + (1 - alpha) * LB

	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			epsilon(i, j): (float) Well depth of atoms i and j.
		"""

		geom = math.sqrt(ei * ej)
		LB = math.sqrt(ei * ej)
		return alpha * geom + (1 - alpha) * LB


class Geometric_WaldmanHagler(CR):
	"""Implements linear combination of geometric and WaldmanHagler combining rules."""

	def getSigma(self, si, sj, alpha):
		"""Returns sigma(i,j).

		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			sigma(i,j): (float) Collision diameter of atoms i and j.
		"""

		geom = math.sqrt(si * sj)
		WH = math.pow(0.5 * (math.pow(si, 6) + math.pow(sj, 6)), 1.0/6.0)
		return alpha * geom + (1 - alpha) * WH

	def getEpsilon(self, ei, ej, si, sj, alpha):
		"""Returns epsilon(i,j).

		:param ei: (float) Well depth of atom i.
		:param ej: (float) Well depth of atom j.
		:param si: (float) Collision diameter of atom i.
		:param sj: (float) Collision diameter of atom j.
		:param alpha: (float) Parameter for linear combination.
		:return:
			epsilon(i, j): (float) Well depth of atoms i and j.
		"""

		WH = 0.0
		if (si + sj) != 0.0:
			WH = 2.0 * math.sqrt(ei*ej) * (si**3) * (sj**3) / ((si**6) + (sj**6))

		geom = math.sqrt(ei * ej)
		return alpha * geom + (1 - alpha) * WH
