import numpy as np
import sys

import IO
import EEM
import combiningRule


class Conf:
	def __init__(self, fileName, it):
		self._it          = int(it)
		self._nJobs       = 4
		self._scl_dns     = 20.0
		self._scl_hvp     = 1.0
		self._opt_nit     = 5000 
		self._rng_scl     = 5.0
		self._eq_stp      = 150000
		self._eq_frq      = 2500
		self._prd_stp     = 150000
		self._prd_frq     = 2500
		self._wall_time   = 12
		self._cr          = combiningRule.GeometricCR()
		self._matrix      = []
		self._eem         = EEM.NONE()
		self._scl_sig_NEI = 1.0
		self._scl_eps_NEI = 1.0
		self._ignoreIAC   = []

		# input files
		self._inpDir       = '00_inp'
		self._molDataFile  = '{}/{}'.format(self._inpDir, "mol.dat")
		self._atomListFile = '{}/{}'.format(self._inpDir, "listAtom.dat")
		self._bondListFile = '{}/{}'.format(self._inpDir, "listBond.dat")
		self._angListFile  = '{}/{}'.format(self._inpDir, "listAng.dat")
		self._matrixFile   = '{}/{}'.format(self._inpDir, "matrix.dat")
		self._symSigFile   = '{}/{}'.format(self._inpDir, "symmetry_sig.dat")
		self._symEpsFile   = '{}/{}'.format(self._inpDir, "symmetry_eps.dat")
		self._prmNeiFile   = '{}/{}'.format(self._inpDir, "prm_NEI.dat")
		self._prmFile      = '{}/{}_{}.dat'.format(self._inpDir, "prm", self._it)
		self._ifpFile = '{}'.format("prm/2016H66_upd.ifp")

		# plot configuration
		#self._iac_name = {}

		self._inpFiles = [self._molDataFile, self._atomListFile, self._bondListFile, self._angListFile, self._matrixFile, self._ifpFile, self._prmNeiFile, self._prmFile]

		IO.readConf(self, fileName)
		# TODO 
		#IO.checkInpDir(self) 
		IO.checkInpFiles(self)

	@property
	def it(self):
		return self._it

	@property
	def nJobs(self):
		return self._nJobs

	@property
	def scl_dns(self):
		return self._scl_dns

	@property
	def scl_hvp(self):
		return self._scl_hvp

	@property
	def opt_nit(self):
		return self._opt_nit

	@property
	def rng_scl(self):
		return self._rng_scl

	@property
	def eq_stp(self):
		return self._eq_stp

	@property
	def eq_frq(self):
		return self._eq_frq

	@property
	def prd_stp(self):
		return self._prd_stp

	@property
	def prd_frq(self):
		return self._prd_frq

	@property
	def wall_time(self):
		return self._wall_time

	@property
	def cr(self):
		return self._cr

	@property
	def matrix(self):
		return self._matrix

	@property
	def eem(self):
		return self._eem

	@property
	def scl_sig_NEI(self):
		return self._scl_sig_NEI

	@property
	def scl_eps_NEI(self):
		return self._scl_eps_NEI

	@property
	def ignoreIAC(self):
		return self._ignoreIAC

	@property
	def inpFiles(self):
		return self._inpFiles

	@property
	def inpDir(self):
		return self._inpDir

	@property
	def molDataFile(self):
		return self._molDataFile

	@property
	def atomListFile(self):
		return self._atomListFile

	@property
	def bondListFile(self):
		return self._bondListFile

	@property
	def angListFile(self):
		return self._angListFile

	@property
	def matrixFile(self):
		return self._matrixFile

	@property
	def symSigFile(self):
		return self._symSigFile

	@property
	def symEpsFile(self):
		return self._symEpsFile

	@property
	def prmNeiFile(self):
		return self._prmNeiFile

	@property
	def prmFile(self):
		return self._prmFile

	@property
	def ifpFile(self):
		return self._ifpFile

	@it.setter
	def it(self, n):
		self._it = n

	@nJobs.setter
	def nJobs(self, n):
		try:
			n = int(n)
			self._nJobs = n

		except ValueError as err:
			print('\nWrong number of jobs')
			print(err)
			sys.exit('\n')

	@scl_dns.setter
	def scl_dns(self, n):
		self._scl_dns = n

	@scl_hvp.setter
	def scl_hvp(self, n):
		self._scl_hvp = n

	@opt_nit.setter
	def opt_nit(self, n):
		self._opt_nit = n

	@rng_scl.setter
	def rng_scl(self, n):
		self._rng_scl = n

	@eq_stp.setter
	def eq_stp(self, n):
		self._eq_stp = n

	@eq_frq.setter
	def eq_frq(self, n):
		self._eq_frq = n

	@prd_stp.setter
	def prd_stp(self, n):
		self._prd_stp = n

	@prd_frq.setter
	def prd_frq(self, n):
		self._prd_frq = n

	@wall_time.setter
	def wall_time(self, n):
		self._wall_time = n

	@cr.setter
	def cr(self, n):
		classes = {"GEOM": combiningRule.GeometricCR, "WH": combiningRule.WaldmanHagler, "LB": combiningRule.LorentzBerthelot, "FH": combiningRule.FenderHalsey}
		cR = classes[n]()
		self._cr = cR

	@matrix.setter
	def matrix(self, n):
		self._matrix = n

	@eem.setter
	def eem(self, n):
		classes = {'HALO': EEM.NONE, 'AA-Alk': EEM.AA_Alk}
		eem = classes[n]()
		self._eem = eem

	@scl_sig_NEI.setter
	def scl_sig_NEI(self, n):
		self._scl_sig_NEI = n

	@scl_eps_NEI.setter
	def scl_eps_NEI(self, n):
		self._scl_eps_NEI = n

	@ignoreIAC.setter
	def ignoreIAC(self, n):
		if isinstance(n, str):
			self._ignoreIAC = [int(n)]
		elif isinstance(n, list):
			self._ignoreIAC = [int(i) for i in n]
		else:
			if np.isnan(n):
				self._ignoreIAC = []
			else:
				sys.exit('ERROR: wrong type for ignoreIAC')

	@inpFiles.setter
	def inpFiles(self, n):
		self._inpFiles = n

	def ignoreAtom(self, iac):
		#print(self._ignoreIAC)
		if iac in self._ignoreIAC:
			return True
		return False

	def __str__(self):
		s = '\t{}:00 {} {}'.format(self._wall_time, self._prd_frq, self._prd_stp)
		return s

