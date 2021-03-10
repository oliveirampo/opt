"""Atom type object.

Classes:
--------
	IAC
"""


import sys


from parameter import Parameter


class IAC:
	"""Atom type object.

	Methods:
	--------
		addSymmetry(symTyp, sym):
	"""

	def __init__(self, iac, typ, nam, sig, rng_sig, eps, rng_eps, hrd, rng_hrd, eln, rng_eln, sig_2, rng_sig_2, eps_2, rng_eps_2):
		"""Constructs all the necessary attributes for the atom type.

		:param iac: (int) Atom type index.
		:param typ: (str) Type of atom.
		:param nam: (str) Name of atom.
		:param sig: (float) Sigma parameter.
		:param rng_sig: (float) Range of sigma parameter variation.
		:param eps: (float) Epsilon parameter
		:param rng_eps: (float) Range of epsilon parameter variation.
		:param hrd: (float) Hardness.
		:param rng_hrd: (float) Range of hardness parameter variation.
		:param eln: (float) Electronegativity.
		:param rng_eln: (float) Range of electronegativity parameter variation.
		:param sig_2: (float) Sigma 2 - for hydrogen bonding atoms.
		:param rng_sig_2: (float) Range of sigma_2 parameter variation.
		:param eps_2: (float) Epsilon 2 - for hydrogen bonding atoms.
		:param rng_eps_2: (float) Range of epsilon_2 parameter variation.
		"""

		self._iac = int(iac)
		self._typ = typ
		self._nam = nam
		self._sig = Parameter(iac, "sig", rng_sig, sig, sig, sig)
		self._eps = Parameter(iac, "eps", rng_eps, eps, eps, eps)
		self._hrd = Parameter(iac, "hrd", rng_hrd, hrd, hrd, hrd)
		self._eln = Parameter(iac, "eln", rng_eln, eln, eln, eln)
		self._sig_2 = Parameter(iac, "sig_2", rng_sig_2, sig_2, sig_2, sig_2)
		self._eps_2 = Parameter(iac, "eps_2", rng_eps_2, eps_2, eps_2, eps_2)

		self._sig_nei = Parameter(iac, "sig_nei", rng_sig, sig, sig, sig)
		self._eps_nei = Parameter(iac, "eps_nei", rng_eps, eps, eps, eps)
		self._fixed_nei = False

		self._vdw = 0.0

	@property
	def iac(self):
		return self._iac

	@property
	def typ(self):
		return self._typ

	@property
	def nam(self):
		return self._nam

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

	@property
	def vdw(self):
		return self._vdw

	@sig_nei.setter
	def sig_nei(self, n):
		self._sig_nei = float(n)

	@eps_nei.setter
	def eps_nei(self, n):
		self._eps_nei = float(n)

	@fixed_nei.setter
	def fixed_nei(self, n):
		self._fixed_nei = n

	@vdw.setter
	def vdw(self, n):
		self._vdw = n

	def __str__(self):
		s = '\t{:3} {:5}'.format(self._iac, self._typ)
		return s

	def addSymmetry(self, symTyp, sym):
		"""Adds symmetry between atom types.

		:param symTyp: (str) Type of symmetry.
		:param sym: (str) Code that marks atom types as symmetric.

		Symmetric atom types are optimized together
		and share the same value for the given parameter.
		"""

		if symTyp == 'sig':
			self._sig.symmetry = sym
		elif symTyp == 'eps':
			self._eps.symmetry = sym
		else:
			print('Symmetry type not implemented: {}'.format(symTyp))
			sys.exit(1)


def add_vdw_radii(atomTypes, vdw):
	"""Adds van der Waals radii to atom types.

	:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
	:param vdw: (dict) Maps of atom type to van der Waals radius.
	"""

	for iac in atomTypes:
		atomType = atomTypes[iac]
		typ = atomType.typ

		radius = vdw[typ]
		atomType.vdw = radius
