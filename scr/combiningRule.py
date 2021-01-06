from abc import ABC, abstractmethod
import math


class CR(ABC):
	@abstractmethod
	def getSigma(self, si, sj):
		pass

	@abstractmethod
	def getEpsilon(self, ei, ej, si, sj):
		pass


class GeometricCR(CR):
	def getSigma(self, si, sj):
		return math.sqrt(si * sj)

	def getEpsilon(self, ei, ej, si, sj):
		return math.sqrt(ei * ej)


class WaldmanHagler(CR):
	def getSigma(self, si, sj):
		return math.pow(0.5 * (math.pow(si, 6) + math.pow(sj, 6)), 1.0/6.0)

	def getEpsilon(self, ei, ej, si, sj):
		return 2.0 * math.sqrt(ei*ej) * (si**3) * (sj**3) / ((si**6) + (sj**6))


class LorentzBerthelot(CR):
	def getSigma(self, si, sj):
		return (si + sj) / 2.0

	def getEpsilon(self, ei, ej, si, sj):
		return math.sqrt(ei * ej)


class FenderHalsey(CR):
	def getSigma(self, si, sj):
		return (si + sj) / 2

	def getEpsilon(self, ei, ej, si, sj):
		return (2 * ei * ej) / (ei * ej)
