import sys


from parameter import Parameter


class IAC:
	def __init__(self, iac, typ, sig, rng_sig, eps, rng_eps, hrd, rng_hrd, eln, rng_eln, sig_2, rng_sig_2, eps_2, rng_eps_2):
		self._iac = int(iac)
		self._typ = typ
		self._sig = Parameter(iac, "sig", rng_sig, sig, sig, sig)
		self._eps = Parameter(iac, "eps", rng_eps, eps, eps, eps)
		self._hrd = Parameter(iac, "hrd", rng_hrd, hrd, hrd, hrd)
		self._eln = Parameter(iac, "eln", rng_eln, eln, eln, eln)
		self._sig_2 = Parameter(iac, "sig_2", rng_sig_2, sig_2, sig_2, sig_2)
		self._eps_2 = Parameter(iac, "eps_2", rng_eps_2, eps_2, eps_2, eps_2)

		self._sig_nei = Parameter(iac, "sig_nei", rng_sig, sig, sig, sig)
		self._eps_nei = Parameter(iac, "eps_nei", rng_eps, eps, eps, eps)
		self._fixed_nei = False

	@property
	def iac(self):
		return self._iac

	@property
	def typ(self):
		return self._typ

	@property
	def sig(self):
		return self._sig

	@property
	def eps(self):
		return self._eps

	@property
	def hrd(self):
		return self._hrd

	@property
	def eln(self):
		return self._eln

	@property
	def sig_2(self):
		return self._sig_2

	@property
	def eps_2(self):
		return self._eps_2

	@property
	def sig_nei(self):
		return self._sig_nei

	@property
	def eps_nei(self):
		return self._eps_nei

	@property
	def fixed_nei(self):
		return self._fixed_nei

	@sig_nei.setter
	def sig_nei(self, n):
		self._sig_nei = float(n)

	@eps_nei.setter
	def eps_nei(self, n):
		self._eps_nei = float(n)

	@fixed_nei.setter
	def fixed_nei(self, n):
		self._fixed_nei = n

	def __str__(self):
		s = '\t{:3} {:5}'.format(self._iac, self._typ)
		return s

	def updateVal(self, prmTyp, val):
		if prmTyp == 'sig':
			self._sig.cur = val
		else:
			print('No such parameter type: {}'.format(prmTyp))
			sys.exit(1)

	def addSymmetry(self, symTyp, sym):
		if symTyp == 'sig':
			self._sig.symmetry = sym
		elif symTyp == 'eps':
			self._eps.symmetry = sym
		else:
			print('Symmetry type not implemented: {}'.format(symTyp))
			sys.exit(1)
